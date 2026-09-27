# get_bal() builds on get_PSM() / get_PSW() for the weights and on cobalt
# (a WeightIt dependency, so always installed) for the balance, so the core
# paths never skip. optmatch ("full") and RegR (save) are Suggests.

bal_data <- function(n = 400L) {
  set.seed(20260927)
  d <- data.frame(x1 = stats::rnorm(n),
                  x2 = stats::rbinom(n, 1, 0.4),
                  x3 = stats::runif(n))
  # controls outnumber treated, so 1:1 and 1:2 matching both find partners
  lp <- -0.9 - 1.1 * d$x1 + 0.8 * d$x2 - 1.2 * d$x3
  d$z <- stats::rbinom(n, 1, stats::plogis(lp))
  d
}

bal_adj <- c("x1", "x2", "x3")

# signed SMD of one covariate under one scheme, from the $balance long table
bal_smd <- function(res, method, variable = "x1") {
  b <- res$balance
  b$smd[b$method == method & b$variable == variable]
}


test_that("get_bal expands shorthand into labelled schemes", {
  d   <- bal_data()
  res <- get_bal(d, treat = "z", adj_var = bal_adj,
                 methods = c("PSM", "ATE", "ATO"))

  labs <- c("PS matching (nearest, ATT)", "IPTW (ATE)",
            "Overlap weighting (ATO)")
  expect_named(res, c("plt", "balance", "data"))
  expect_s3_class(res$plt, "ggplot")
  expect_named(attr(res$plt, "plot_size"), c("width", "height"))
  expect_identical(levels(res$plt$data$Sample), c("Unadjusted", labs))

  expect_named(res$balance, c("variable", "method", "smd"))
  expect_identical(unique(res$balance$method), c("Unadjusted", labs))
  expect_setequal(unique(res$balance$variable), bal_adj)

  expect_identical(nrow(res$data), nrow(d))
  expect_true(all(labs %in% names(res$data)))
  # the "PSM" shorthand is nearest 1:1, ATT, caliper 0.2
  ref <- get_PSM(d, treat = "z", adj_var = bal_adj, method = "nearest",
                 estimand = "ATT", ratio = 1, caliper = 0.2, balance = FALSE)
  expect_identical(res$data[["PS matching (nearest, ATT)"]],
                   ref$data$w_nearest)
})

test_that("the default compares matching with all six weights", {
  res <- get_bal(bal_data(), treat = "z", adj_var = bal_adj)
  expect_identical(
    levels(res$plt$data$Sample),
    c("Unadjusted", "PS matching (nearest, ATT)", "IPTW (ATE)",
      "SMR weighting (ATT)", "SMR weighting (ATC)", "Overlap weighting (ATO)",
      "Matching weighting (ATM)", "Entropy weighting (EW)"))
})

test_that("a custom list takes several matching schemes side by side", {
  d <- bal_data()
  methods <- list(
    `PSM 1:1, caliper 0.2` = list(design = "matching", method = "nearest",
                                  estimand = "ATT", ratio = 1, caliper = 0.2),
    `PSM 1:2, caliper 0.1` = list(design = "matching", method = "nearest",
                                  estimand = "ATT", ratio = 2, caliper = 0.1),
    `Stabilized IPTW (ATE)` = list(design = "weighting", method = "glm",
                                   estimand = "ATE", stabilize = TRUE))
  res <- get_bal(d, treat = "z", adj_var = bal_adj, methods = methods)

  expect_identical(levels(res$plt$data$Sample),
                   c("Unadjusted", names(methods)))
  ref <- get_PSM(d, treat = "z", adj_var = bal_adj, ratio = 2,
                 caliper = 0.1, balance = FALSE)
  expect_identical(res$data[["PSM 1:2, caliper 0.1"]], ref$data$w_nearest)
  expect_false(identical(res$data[["PSM 1:1, caliper 0.2"]],
                         res$data[["PSM 1:2, caliper 0.1"]]))
})

test_that("every scheme is standardised by the unadjusted pooled SD", {
  d   <- bal_data()
  res <- get_bal(d, treat = "z", adj_var = bal_adj,
                 methods = c("PSM", "ATE"))
  x   <- d$x1
  z   <- d$z
  den <- sqrt((stats::var(x[z == 1]) + stats::var(x[z == 0])) / 2)
  hand <- function(w) (stats::weighted.mean(x[z == 1], w[z == 1]) -
                         stats::weighted.mean(x[z == 0], w[z == 0])) / den

  expect_equal(bal_smd(res, "Unadjusted"), hand(rep(1, nrow(d))))
  expect_equal(bal_smd(res, "IPTW (ATE)"), hand(res$data[["IPTW (ATE)"]]))
  expect_equal(bal_smd(res, "PS matching (nearest, ATT)"),
               hand(res$data[["PS matching (nearest, ATT)"]]))
})

test_that("overlap weights balance exactly and stabilizing changes nothing", {
  methods <- list(
    `IPTW`            = list(design = "weighting", estimand = "ATE"),
    `Stabilized IPTW` = list(design = "weighting", estimand = "ATE",
                             stabilize = TRUE),
    `Overlap`         = list(design = "weighting", estimand = "ATO"))
  res <- get_bal(bal_data(), treat = "z", adj_var = bal_adj,
                 methods = methods)

  b <- res$balance
  expect_true(all(abs(b$smd[b$method == "Overlap"]) < 1e-8))
  expect_equal(b$smd[b$method == "Stabilized IPTW"], b$smd[b$method == "IPTW"])
})

test_that("trimmed units take weight 0, not NA", {
  d <- bal_data()
  methods <- list(`Trimmed IPTW` = list(
    design = "weighting", estimand = "ATE",
    trim_args = list(method = "ps", lower = 0.1, upper = 0.9)))
  res <- get_bal(d, treat = "z", adj_var = bal_adj, methods = methods)

  ref <- get_PSW(d, treat = "z", adj_var = bal_adj, estimand = "ATE",
                 trim_args = list(method = "ps", lower = 0.1, upper = 0.9),
                 balance = FALSE)
  w <- res$data[["Trimmed IPTW"]]
  expect_gt(sum(ref$data$.trimmed), 0)
  expect_false(anyNA(w))
  expect_true(all(w[ref$data$.trimmed] == 0))
  expect_identical(w[!ref$data$.trimmed], ref$data$w_ate[!ref$data$.trimmed])
  expect_false(anyNA(res$balance$smd))
})

test_that("literal names work and binary factors keep the variable name", {
  d <- bal_data()
  names(d)[names(d) == "x1"] <- "Age (years)"
  d$sex <- factor(ifelse(d$x2 == 1, "M", "F"))
  d$x2  <- NULL
  adj   <- c("Age (years)", "sex", "x3")

  res <- get_bal(d, treat = "z", adj_var = adj, methods = c("PSM", "ATO"))
  expect_setequal(unique(res$balance$variable), adj)
  expect_true("sex" %in% levels(res$plt$data$var))
})

test_that("incomplete rows are dropped once, for every scheme", {
  d <- bal_data()
  d$x3[c(3, 7)] <- NA
  res <- get_bal(d, treat = "z", adj_var = bal_adj, methods = c("PSM", "ATE"))
  expect_identical(nrow(res$data), nrow(d) - 2L)
  expect_false(anyNA(res$data[["IPTW (ATE)"]]))
})

test_that("method specifications are validated before anything is fitted", {
  d <- bal_data()
  run <- function(methods) get_bal(d, treat = "z", adj_var = bal_adj,
                                   methods = methods)

  expect_error(run(list(`PSM` = list(
    design = "matching", method = "nearest",
    matchit_args = list(ratio = 1, caliper = 0.2)))),
    "`matchit_args`.*ratio = 1, caliper = 0.2")
  expect_error(run(list(`A` = list(design = "weighting", estimand = "ATE",
                                   foo = 1))),
               "`A`.*unknown field.*`foo`")
  expect_error(run(list(`A` = list(estimand = "ATE"))),
               "`A`.*`design`")
  expect_error(run(list(`A` = list(design = "stratify"))),
               "`design` must be \"matching\" or \"weighting\"")
  expect_error(run(list(`A` = list(design = "weighting"))),
               "`A`.*one `estimand`")
  expect_error(run(list(`A` = list(design = "weighting",
                                   estimand = c("ATE", "ATO")))),
               "`A`.*one `estimand`")
  expect_error(run(list(`A` = list(design = "matching",
                                   method = c("nearest", "full")))),
               "`A`.*one matching `method`")
  expect_error(run(list(`A` = list(design = "weighting", estimand = "ATE",
                                   adj_var = "x1"))),
               "`A`.*`adj_var`.*set by get_bal")
  expect_error(run(list(list(design = "weighting", estimand = "ATE"))),
               "named list")
  expect_error(run(stats::setNames(
    rep(list(list(design = "weighting", estimand = "ATE")), 2), c("A", "A"))),
    "duplicated.*`A`")
  expect_error(run(c("PSM", "IPW")), "Unknown shorthand.*\"IPW\"")
  expect_error(run(stats::setNames(
    rep(list(list(design = "weighting", estimand = "ATE")), 15),
    paste("m", 1:15))), "at most 14")

  d$`IPTW (ATE)` <- 1
  expect_error(get_bal(d, treat = "z", adj_var = bal_adj, methods = "ATE"),
               "`IPTW \\(ATE\\)`.*already")
})

test_that("a failing scheme is named in the error", {
  expect_error(
    get_bal(bal_data(), treat = "z", adj_var = bal_adj,
            methods = list(`Bad PSM` = list(design = "matching",
                                             method = "nearest",
                                             estimand = "ATE"))),
    "`Bad PSM`.*not available")
})

test_that("save writes a PDF only when it is a non-empty list", {
  skip_if_not_installed("RegR")
  d   <- bal_data()
  dir <- withr::local_tempdir()

  get_bal(d, "z", bal_adj, methods = "ATO", save = list())
  get_bal(d, "z", bal_adj, methods = "ATO", save = NULL)
  expect_identical(list.files(dir), character(0))

  f   <- file.path(dir, "bal.pdf")
  res <- get_bal(d, "z", bal_adj, methods = "ATO", save = list(filename = f))
  expect_s3_class(res$plt, "ggplot")
  expect_true(file.exists(f))
  expect_gt(file.size(f), 0)

  expect_error(get_bal(d, "z", bal_adj, methods = "ATO", save = "nope.pdf"),
               "`save` must be `NULL` or a list")
})
