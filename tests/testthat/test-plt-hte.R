# plt_hte_dep() draws from a get_hte() result. grf and sandwich are Suggests;
# every test that needs a forest skips without them.

dep_cache <- new.env()

dep_res <- function() {
  skip_if_not_installed("grf")
  skip_if_not_installed("sandwich")
  if (is.null(dep_cache$res)) {
    set.seed(20260923)
    n <- 500L
    d <- data.frame(age   = round(stats::runif(n, 20, 85)),
                    sex   = factor(sample(c("F", "M"), n, replace = TRUE)),
                    stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
    d$z <- stats::rbinom(n, 1, 0.5)
    d$y <- stats::rbinom(n, 1, stats::plogis(-1 + 0.02 * (d$age - 50) +
                                               d$z * (0.2 + 0.6 * (d$sex == "M"))))
    # sex is added to the covariates with a message
    dep_cache$res <- suppressMessages(get_hte(
      d, cat_var = "z", sub_var = "sex", adj_var = c("age", "stage"),
      surv = "y", grf_args = list(num.trees = 300, seed = 1)))
  }
  dep_cache$res
}

# A patchwork keeps its last panel as the object itself.
panels_of <- function(p)
  if (inherits(p, "patchwork")) c(p$patches$plots, list(p)) else list(p)

strip_of <- function(q) {
  for (l in q$layers)
    if (is.data.frame(l$data) && "panel" %in% names(l$data))
      return(as.character(l$data$panel[1]))
  NA_character_
}

layer_data_of <- function(q, geom, col = NULL) {
  for (l in q$layers)
    if (inherits(l$geom, geom) && (is.null(col) || col %in% names(l$data)))
      return(l$data)
  NULL
}

# Forest CATE averaged over every row with the stage columns set to one level.
manual_pdp <- function(fit, set) {
  X <- fit$X.orig
  for (nm in names(set)) X[, nm] <- set[[nm]]
  mean(stats::predict(fit, X)$predictions)
}


test_that(".hte_dr_var matches grf subsets and a hand-built spline test", {
  res <- dep_res()
  w <- res$fit$W.orig
  z <- stats::qnorm(0.975)

  st <- .hte_dr_var(res$data, "stage", w, z, 3L)
  expect_identical(st$type, "categorical")
  expect_identical(st$df, 2L)
  for (lv in c("I", "II", "III")) {
    ref <- grf::average_treatment_effect(res$fit, subset = res$data$stage == lv)
    row <- st$levels[st$levels$level == lv, ]
    expect_equal(row$estimate, unname(ref[["estimate"]]))
    expect_equal(row$std.error, unname(ref[["std.err"]]))
  }

  age <- .hte_dr_var(res$data, "age", w, z, 3L)
  fit <- stats::lm(.dr_score ~ splines::ns(age, df = 3), data = res$data)
  V   <- sandwich::vcovHC(fit, type = "HC3")
  b   <- stats::coef(fit)[-1]
  expect_identical(age$type, "continuous")
  expect_identical(age$df, 3L)
  expect_equal(age$p_het, stats::pchisq(drop(b %*% solve(V[-1, -1], b)), 3,
                                        lower.tail = FALSE))
})

test_that("get_hte() reports p_het in $importance, equal to the subgroup p_inter", {
  res <- dep_res()
  imp <- res$importance
  expect_named(imp, c("variable", "importance", "n_col", "df", "p_het"))
  expect_equal(imp$p_het[imp$variable == "sex"],
               res$subgroup$p_inter[res$subgroup$measure == "diff"][1])
  expect_identical(imp$df[imp$variable == "age"], 3L)
  expect_identical(imp$df[imp$variable == "stage"], 2L)
})

test_that("plt_hte_dep draws one panel per covariate in importance order", {
  res <- dep_res()
  p <- plt_hte_dep(res)
  expect_s3_class(p, "patchwork")
  strips <- vapply(panels_of(p), strip_of, "")
  expect_identical(sub(" \\(p_het.*", "", strips), res$importance$variable)
  expect_true(all(grepl("(p_het ", strips, fixed = TRUE)))
  expect_identical(names(attr(p, "plot_size")), c("width", "height"))

  fct <- vapply(panels_of(plt_hte_dep(res, x_var = "fct")), strip_of, "")
  expect_identical(sub(" \\(p_het.*", "", fct),
                   setdiff(res$importance$variable, "age"))
  num <- plt_hte_dep(res, x_var = "num")
  expect_false(inherits(num, "patchwork"))
  expect_match(strip_of(num), "^age ")
  expect_identical(strip_of(plt_hte_dep(res, x_var = "stage", display = "cate")),
                   "stage")
})

test_that("the dr layer is the doubly robust summary; pdp averages forest predictions", {
  res <- dep_res()
  p <- plt_hte_dep(res, x_var = "stage", display = c("dr", "pdp"))

  dr <- .hte_dr_var(res$data, "stage", res$fit$W.orig, stats::qnorm(0.975), 3L)
  expect_equal(layer_data_of(p, "GeomPointrange")$estimate, dr$levels$estimate)

  pdp <- layer_data_of(p, "GeomPoint", "estimate")
  expect_equal(pdp$estimate[pdp$x == "II"],
               manual_pdp(res$fit, list(stageI = 0, stageII = 1, stageIII = 0)))

  cont <- plt_hte_dep(res, x_var = "age", display = "pdp",
                      pdp_args = list(grid_n = 5))
  line <- layer_data_of(cont, "GeomLine", "estimate")
  expect_identical(nrow(line), 5L)
  expect_equal(line$estimate[1], manual_pdp(res$fit, list(age = min(res$data$age))))
})

test_that("type = 'heat' tiles the two-way partial dependence", {
  res <- dep_res()
  p <- plt_hte_dep(res, x_var = c("sex", "stage"), type = "heat")
  expect_false(inherits(p, "patchwork"))
  expect_identical(nrow(p$data), 6L)
  cell <- p$data[p$data$sex == "M" & p$data$stage == "III", ]
  expect_equal(cell$estimate,
               manual_pdp(res$fit, list(sexF = 0, sexM = 1, stageI = 0,
                                        stageII = 0, stageIII = 1)))
})

test_that("pdp sets an integer-coded factor through its level code", {
  skip_if_not_installed("grf")
  set.seed(20260923)
  n <- 500L
  d <- data.frame(age   = round(stats::runif(n, 20, 85)),
                  stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
  d$z <- stats::rbinom(n, 1, 0.5)
  d$y <- stats::rbinom(n, 1, stats::plogis(-1 + d$z * (0.2 + 0.4 * (d$stage == "III"))))
  res <- get_hte(d, cat_var = "z", adj_var = c("age", "stage"), surv = "y",
                 factor_encoding = "integer",
                 grf_args = list(num.trees = 300, seed = 1))

  pdp <- layer_data_of(plt_hte_dep(res, x_var = "stage", display = "pdp"),
                       "GeomPoint", "estimate")
  expect_equal(pdp$estimate[pdp$x == "II"], manual_pdp(res$fit, list(stage = 2)))
  heat <- plt_hte_dep(res, x_var = c("age", "stage"), type = "heat",
                      pdp_args = list(grid_n = 3))
  cell <- heat$data[heat$data$stage == "III" & heat$data$age == min(d$age), ]
  expect_equal(cell$estimate,
               manual_pdp(res$fit, list(age = min(d$age), stage = 3)))
})

test_that("a covariate with missing values draws only its observed patients", {
  skip_if_not_installed("grf")
  skip_if_not_installed("sandwich")
  skip_if_not_installed("patchwork")
  set.seed(11)
  n <- 400L
  d <- data.frame(age   = stats::runif(n, 20, 85),
                  stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
  d$age[1:40]    <- NA
  d$stage[41:80] <- NA
  d$z <- stats::rbinom(n, 1, 0.5)
  d$y <- stats::rbinom(n, 1, 0.3 + 0.1 * d$z)
  res <- get_hte(d, cat_var = "z", adj_var = c("age", "stage"), surv = "y",
                 grf_args = list(num.trees = 200, seed = 1))
  expect_identical(nrow(res$data), n)

  p <- plt_hte_dep(res, display = c("cate", "dr", "pdp"))
  for (q in panels_of(p)) expect_no_warning(ggplot2::ggplot_build(q))
  pts <- layer_data_of(plt_hte_dep(res, x_var = "stage", display = "cate"),
                       "GeomPoint")
  expect_identical(nrow(pts), n - 40L)
  expect_false(anyNA(pts$x))
})

test_that("survival: levels and values no arm follows past `time` leave the dr layer", {
  skip_if_not_installed("grf")
  skip_if_not_installed("sandwich")
  skip_if_not_installed("patchwork")
  set.seed(20260924)
  n <- 800L
  d <- data.frame(age   = stats::runif(n, 20, 85),
                  stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
  d$z  <- stats::rbinom(n, 1, 0.5)
  ev   <- stats::rexp(n, 0.02)
  cens <- pmin(stats::rexp(n, 0.01), 120)
  d$time <- pmin(ev, cens)
  d$DSS  <- as.integer(ev <= cens)
  # no treated patient in stage III or over 70 is followed past 60
  cut <- d$z == 1 & (d$stage == "III" | d$age > 70) & d$time > 50
  d$time[cut] <- 50
  d$DSS[cut]  <- 0L
  expect_warning(
    res <- get_hte(d, "z", sub_var = "stage", adj_var = c("age", "stage"),
                   surv = TRUE, time = 60,
                   grf_args = list(num.trees = 300, seed = 1)),
    "stage = III: no patient in one arm")
  imp <- res$importance

  # stage III is NA in $subgroup, so it leaves p_het and the dr layer as well
  expect_equal(imp$p_het[imp$variable == "stage"],
               res$subgroup$p_inter[res$subgroup$measure == "diff"][1])
  expect_identical(imp$df[imp$variable == "stage"], 1L)
  pr <- layer_data_of(plt_hte_dep(res, x_var = "stage", display = "dr"),
                      "GeomPointrange")
  expect_identical(as.character(pr$x), c("I", "II"))

  # the age spline is fit and drawn only where both arms reach past 60
  dd   <- res$data
  past <- dd$time > 60
  lo   <- max(min(dd$age[past & dd$z == 1]), min(dd$age[past & dd$z == 0]))
  hi   <- min(max(dd$age[past & dd$z == 1]), max(dd$age[past & dd$z == 0]))
  line <- layer_data_of(plt_hte_dep(res, x_var = "age", display = "dr"),
                        "GeomLine", "estimate")
  expect_equal(range(line$x), c(lo, hi))
  fit <- stats::lm(.dr_score ~ splines::ns(age, df = 3),
                   data = dd[dd$age >= lo & dd$age <= hi, ])
  V   <- sandwich::vcovHC(fit, type = "HC3")
  b   <- stats::coef(fit)[-1]
  expect_equal(imp$p_het[imp$variable == "age"],
               stats::pchisq(drop(b %*% solve(V[-1, -1], b)), 3,
                             lower.tail = FALSE))
})

test_that("invalid requests stop with a clear message", {
  res <- dep_res()
  expect_error(plt_hte_dep(list()), "hte_res")
  expect_error(plt_hte_dep(res, x_var = "nope"), "nope")
  expect_error(plt_hte_dep(res, x_var = "age", type = "heat"), "exactly two")
  expect_error(plt_hte_dep(res, x_var = c("age", "sex"), type = "heat",
                           display = "dr"), "`display`")
  expect_error(plt_hte_dep(res, dr_args = list(df = 3)), "unknown field")
  expect_error(plt_hte_dep(res, axis_arg = list(share_y = "var")), "share_y")
  expect_error(plt_hte_dep(res, save = "a.pdf"), "`save`")
})

test_that("save writes one PDF and returns the plot unchanged", {
  res <- dep_res()
  expect_s3_class(plt_hte_dep(res, x_var = "stage", save = list()), "ggplot")
  expect_s3_class(plt_hte_dep(res, x_var = "stage", save = NULL), "ggplot")
  skip_if_not_installed("RegR")
  f <- tempfile(fileext = ".pdf")
  p <- plt_hte_dep(res, x_var = "stage", save = list(filename = f))
  expect_true(file.exists(f))
  expect_s3_class(p, "ggplot")
  expect_false(inherits(p, "patchwork"))
})
