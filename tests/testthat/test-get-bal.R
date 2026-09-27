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

# Yang & Dalton (2012) SMD of one factor, written independently of get_bal():
# weighted level proportions, treated minus control, over the mean of the two
# arms' unweighted multinomial covariances, on k - 1 levels with solve()
yd_smd <- function(x, z, w = rep(1, length(x))) {
  x <- factor(x)
  k <- levels(x)[-1L]
  prop <- function(s, wt)
    vapply(k, function(l) sum(wt[s] * (x[s] == l)) / sum(wt[s]), 0)
  covm <- function(s) {
    q <- prop(s, rep(1, length(x)))
    diag(q) - outer(q, q)
  }
  D <- prop(z == 1, w) - prop(z == 0, w)
  sqrt(drop(D %*% solve((covm(z == 1) + covm(z == 0)) / 2, D)))
}

stage_data <- function() {
  d <- bal_data()
  d$sex   <- factor(ifelse(d$x2 == 1, "M", "F"))
  d$stage <- cut(d$x1, c(-Inf, -0.5, 0.5, Inf), labels = c("I", "II", "III"))
  d
}

test_that("cat_smd = \"overall\" gives a factor one Yang & Dalton row", {
  adj <- c("x1", "sex", "stage")
  res <- get_bal(stage_data(), "z", adj, methods = "ATE")

  expect_setequal(unique(res$balance$variable), adj)
  expect_equal(bal_smd(res, "Unadjusted", "stage"),
               yd_smd(res$data$stage, res$data$z))
  expect_equal(bal_smd(res, "IPTW (ATE)", "stage"),
               yd_smd(res$data$stage, res$data$z, res$data[["IPTW (ATE)"]]))
  expect_true("stage" %in% levels(res$plt$data$var))
})

test_that("cat_smd = \"level\" keeps one row per level; two levels are unchanged", {
  d   <- stage_data()
  adj <- c("x1", "sex", "stage")
  lev <- get_bal(d, "z", adj, methods = "ATE", cat_smd = "level")
  ove <- get_bal(d, "z", adj, methods = "ATE")
  lb  <- lev$balance
  ob  <- ove$balance

  expect_setequal(unique(lb$variable),
                  c("x1", "sex", "stage_I", "stage_II", "stage_III"))
  expect_identical(lb$smd[lb$variable %in% c("x1", "sex")],
                   ob$smd[ob$variable %in% c("x1", "sex")])
  # the multivariate SMD is never below the largest level-wise one
  for (m in c("Unadjusted", "IPTW (ATE)"))
    expect_gte(bal_smd(ove, m, "stage"),
               max(abs(lb$smd[lb$method == m & startsWith(lb$variable, "stage_")])))

  expect_error(get_bal(d, "z", adj, methods = "ATE", cat_smd = "total"),
               "should be one of")
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

# the red reference line fmt_ref() adds last, and the point layer's data
bal_ref  <- function(p) {
  v <- Filter(function(l) inherits(l$geom, "GeomVline"), p$layers)
  v <- v[[length(v)]]
  c(as.list(v$data), v$aes_params)     # xintercept lives in the layer data
}
bal_pts <- function(p) ggplot2::layer_data(
  p, which(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), TRUE)))

test_that("love_args defaults keep the RegR-style plot", {
  expect_identical(names(formals(get_bal)),
                   c("data", "treat", "adj_var", "methods", "cat_smd",
                     "tbl", "love_args", "save_plt", "save_tbl"))
  p <- get_bal(bal_data(), "z", bal_adj, methods = c("PSM", "ATE"))$plt

  expect_equal(unname(bal_ref(p)$xintercept), 0.1)
  expect_identical(bal_ref(p)$colour, "red")
  pts <- bal_pts(p)
  expect_setequal(unique(pts$colour), UtilsR::pal_lancet[1:3])
  expect_setequal(unique(pts$shape), c(17, 16, 15))
  expect_true(all(pts$size == 3.5))
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomPath"),
                         TRUE)))
  expect_identical(p$theme$text$size, 16)
  expect_equal(p$theme$legend.position.inside, c(0.99, 0.02))
  expect_equal(p$theme$legend.justification, c(1, 0))
})

test_that("love_args restyles points, lines, theme and legend", {
  p <- get_bal(bal_data(), "z", bal_adj, methods = c("PSM", "ATE"),
               love_args = list(threshold = 0.2, ref_color = "blue",
                                colors = c("black", "grey50", "blue", "red"),
                                shapes = c(1, 2, 5), size = 2, line = FALSE,
                                base_size = 12, legend_position = "bottom"))$plt

  expect_equal(unname(bal_ref(p)$xintercept), 0.2)
  expect_identical(bal_ref(p)$colour, "blue")
  pts <- bal_pts(p)
  expect_setequal(unique(pts$colour), c("black", "grey50", "blue"))
  expect_setequal(unique(pts$shape), c(1, 2, 5))
  expect_true(all(pts$size == 2))
  expect_false(any(vapply(p$layers, function(l) inherits(l$geom, "GeomPath"),
                          TRUE)))
  expect_identical(p$theme$text$size, 12)
  expect_identical(p$theme$legend.position, "bottom")
  # fmt_legend() centres an outside legend on its edge; the long scheme
  # labels are wrapped into two columns rather than one clipped row
  expect_equal(p$theme$legend.justification, c(0.5, 1))
  expect_identical(p$guides$guides$colour$params$ncol, 2)
  expect_identical(p$guides$guides$shape$params$ncol, 2)

  # an explicit justification survives fmt_legend(), which would overwrite it
  p <- get_bal(bal_data(), "z", bal_adj, methods = "ATO",
               love_args = list(legend_justification = c(0, 0)))$plt
  expect_equal(p$theme$legend.justification, c(0, 0))
})

test_that("love_args$var_names relabels the plot, not the table", {
  d <- bal_data()
  d$sex   <- factor(ifelse(d$x2 == 1, "M", "F"))
  d$stage <- factor(rep_len(c("I", "II", "III"), nrow(d)))
  adj <- c("x1", "sex", "stage")
  vn  <- c(x1 = "Age", sex = "Sex", stage = "Stage")
  res <- get_bal(d, "z", adj, methods = "ATO", cat_smd = "level",
                 love_args = list(var_names = vn))

  expect_setequal(levels(res$plt$data$var),
                  c("Age", "Sex", "Stage_I", "Stage_II", "Stage_III"))
  expect_setequal(unique(res$balance$variable),
                  c("x1", "sex", "stage_I", "stage_II", "stage_III"))

  res <- get_bal(d, "z", adj, methods = "ATO", love_args = list(var_names = vn))
  expect_setequal(levels(res$plt$data$var), c("Age", "Sex", "Stage"))
  expect_setequal(unique(res$balance$variable), adj)
})

test_that("love_args is validated before anything is fitted", {
  d <- bal_data()
  run <- function(love_args, methods = c("PSM", "ATE"))
    get_bal(d, "z", bal_adj, methods = methods, love_args = love_args)

  expect_error(get_bal(d, "z", bal_adj, methods = "ATE", threshold = 0.1),
               "unused argument")
  expect_error(run(list(foo = 1)), "`love_args` contains unknown field.*`foo`")
  expect_error(run(list(threshold = "a")),
               "`love_args\\$threshold` must be a single number")
  expect_error(run(list(colors = c("red", "blue"))),
               "`love_args\\$colors` needs at least 3")
  expect_error(run(list(shapes = 1:2)), "`love_args\\$shapes` needs at least 3")
  expect_error(run(list(size = -1)),
               "`love_args\\$size` must be a single positive number")
  expect_error(run(list(line = NA)), "`love_args\\$line` must be TRUE or FALSE")
  expect_error(run(list(var_order = "bogus")),
               "`love_args\\$var_order` must be")
  expect_error(run(list(legend_position = c(1, 2, 3))),
               "`love_args\\$legend_position` must be")
  expect_error(run(list(var_names = c(foo = "Foo"))),
               "`love_args\\$var_names`.*`foo`.*`adj_var`")
  expect_error(run(list(var_names = "Age")),
               "`love_args\\$var_names` must be a named character vector")
})

test_that("save_plt writes a PDF only when it is a non-empty list", {
  skip_if_not_installed("RegR")
  d   <- bal_data()
  dir <- withr::local_tempdir()

  res <- get_bal(d, "z", bal_adj, methods = "ATO", save_plt = list(),
                 save_tbl = list())
  get_bal(d, "z", bal_adj, methods = "ATO", save_plt = NULL, save_tbl = NULL)
  expect_identical(list.files(dir), character(0))
  expect_null(res$tbl)

  f   <- file.path(dir, "bal.pdf")
  res <- get_bal(d, "z", bal_adj, methods = "ATO", save_plt = list(filename = f))
  expect_s3_class(res$plt, "ggplot")
  expect_true(file.exists(f))
  expect_gt(file.size(f), 0)

  expect_error(get_bal(d, "z", bal_adj, methods = "ATO", save_plt = "nope.pdf"),
               "`save_plt` must be `NULL` or a named list")
})

skip_if_no_tbl <- function() {
  for (p in c("gtsummary", "survey", "cardx", "smd")) skip_if_not_installed(p)
}

test_that("tbl = TRUE merges one gtsummary table per scheme; FALSE skips it", {
  skip_if_no_tbl()
  d   <- stage_data()
  adj <- c("x1", "sex", "stage")
  expect_null(get_bal(d, "z", adj, methods = "ATO")$tbl)

  res  <- get_bal(d, "z", adj, methods = c("PSM", "ATE"), tbl = TRUE)
  labs <- c("Unadjusted", "PS matching (nearest, ATT)", "IPTW (ATE)")
  expect_s3_class(res$tbl, "tbl_merge")
  expect_length(res$tbl$tbls, 3L)
  expect_setequal(res$tbl$table_styling$spanning_header$spanning_header,
                  paste0("**", labs, "**"))

  for (i in seq_along(labs)) {
    t <- res$tbl$tbls[[i]]
    b <- t$table_body[t$table_body$row_type == "label", ]
    # gtsummary's own SMD: a factor's Yang & Dalton value equals get_bal's;
    # a continuous one has its variance over n rather than n - 1
    expect_equal(abs(b$estimate[b$variable == "stage"]),
                 bal_smd(res, labs[i], "stage"))
    expect_equal(abs(b$estimate[b$variable == "x1"]),
                 abs(bal_smd(res, labs[i], "x1")), tolerance = 0.01)
    expect_false(anyNA(b$p.value))
    h <- t$table_styling$header
    expect_identical(h$label[h$column == "estimate"], "**SMD**")
    expect_true(h$hide[h$column == "conf.low"])
  }
})

test_that("save_tbl builds and writes the table; each save takes its own fields", {
  skip_if_not_installed("RegR")
  skip_if_not_installed("flextable")
  skip_if_no_tbl()
  d   <- bal_data()
  dir <- withr::local_tempdir()

  expect_error(get_bal(d, "z", bal_adj, methods = "ATO",
                       save_tbl = list(filename = "x.pdf")),
               "`save_tbl` contains unknown field.*`filename`")
  expect_error(get_bal(d, "z", bal_adj, methods = "ATO",
                       save_plt = list(path = dir)),
               "`save_plt` contains unknown field.*`path`")
  expect_error(get_bal(d, "z", bal_adj, methods = "ATO", save_tbl = list("a")),
               "`save_tbl` must be a fully named list")

  # tbl is left FALSE: a non-empty save_tbl builds the table itself
  res <- get_bal(d, "z", bal_adj, methods = "ATO",
                 save_tbl = list(path = dir, title = "bal_tbl"))
  expect_s3_class(res$tbl, "tbl_merge")
  expect_identical(list.files(dir), "bal_tbl.docx")

  f <- file.path(dir, "bal.pdf")
  get_bal(d, "z", bal_adj, methods = "ATO", save_plt = list(filename = f),
          save_tbl = list(path = dir, title = "bal_tbl2"))
  expect_setequal(list.files(dir), c("bal_tbl.docx", "bal_tbl2.docx", "bal.pdf"))
})
