plt_test_deps <- function(...) {
  for (pkg in c(...)) skip_if_not_installed(pkg)
}

plt_lm_res <- function() {
  skip_if_not_installed("sensemakr")
  get_sens(sensemakr::darfur,
           treat = "directlyharmed", outcome = "peacefactor",
           adj_var = c("age", "farmer_dar", "herder_dar", "pastvoted",
                       "hhsize_darfur", "female", "village"),
           bench_var = "female", method = "lm",
           bench_args = list(k_treat = 1:3))
}

plt_cox_res <- function(...) {
  skip_if_not_installed("survival")
  skip_if_not_installed("tipr")
  d <- stats::na.omit(survival::lung[, c("time", "status", "sex", "age")])
  d$status <- d$status - 1L
  d$sex <- factor(d$sex, labels = c("male", "female"))
  get_sens(d, treat = "sex", outcome = "status", time = "time",
           adj_var = "age", method = "cox", ...)
}

plt_iv_res <- function(...) {
  skip_if_not_installed("iv.sensemakr")
  e <- new.env()
  utils::data("card", package = "iv.sensemakr", envir = e)
  get_sens(e$card, treat = "educ", outcome = "lwage", instrument = "nearc4",
           adj_var = c("exper", "expersq", "black", "south", "smsa"),
           bench_var = "black", method = "iv", ...)
}

test_that("plt_sens draws a labelled native contour for the linear backend", {
  plt_test_deps("sensemakr")
  res <- plt_lm_res()
  p <- plt_sens(res, type = "contour")

  expect_s3_class(p, "ggplot")
  expect_equal(unname(attr(p, "plot_size")), c(7, 6))
  # Native, not a wrapped grob: the contour levels must survive as text and
  # every benchmark bound must be inside the default window.
  d <- ggplot2::ggplot_build(p)$data
  txt <- unlist(lapply(d, function(l) if ("label" %in% names(l)) l$label))
  expect_true(length(txt) > 3L)
  expect_true(any(grepl("1x female", txt)))
  expect_true(any(grepl("3x female", txt)))
  expect_true(any(grepl("Unadjusted", txt)))
  # the critical contour is the red dashed one
  expect_true(any(vapply(d, function(l)
    any(l$colour == "red" & l$linetype == 2, na.rm = TRUE), logical(1))))
})

test_that("plt_sens reproduces sensemakr adjusted estimates on the contour", {
  plt_test_deps("sensemakr")
  res <- plt_lm_res()
  st <- res$sens$sensitivity_stats
  # The label printed next to each bound must equal sensemakr's own adjusted
  # estimate for that scenario.
  expect_equal(
    as.numeric(sensemakr::adjusted_estimate(st$estimate, st$se, st$dof,
                                            res$bounds$r2_treat,
                                            res$bounds$r2_out)),
    res$bounds$adj_estimate, tolerance = 1e-8)
})

test_that("t-value contours retain the direction of the treatment effect", {
  skip_if_not_installed("sensemakr")
  set.seed(933)
  d <- data.frame(z = rep(0:1, 100), x = rnorm(200))
  d$y <- d$z + 0.2 * d$x + rnorm(200)
  for (direction in c(-1, 1)) {
    dd <- d
    dd$y <- direction * d$y
    res <- get_sens(dd, "z", "y", adj_var = "x", method = "lm")
    p <- plt_sens(res, lim = c(0.6, 0.6),
                  contour_args = list(sensitivity_of = "t-value"))
    layers <- ggplot2::ggplot_build(p)$data
    red <- Filter(function(l) "linetype" %in% names(l) &&
      any(l$colour == "red" & l$linetype == 2, na.rm = TRUE), layers)
    expect_length(red, 1L)
    st <- res$sens$sensitivity_stats
    values <- sensemakr::adjusted_t(st$estimate, st$se, st$dof,
                                     red[[1L]]$x, red[[1L]]$y)
    critical <- direction * abs(stats::qt(0.025, st$dof - 1))
    expect_lt(max(abs(as.numeric(values) - critical)), 0.01)
  }
})

test_that("plt_sens wraps the sensemakr extreme plot as a ggplot", {
  plt_test_deps("sensemakr", "ggplotify")
  p <- plt_sens(plt_lm_res(), type = "extreme", extreme_r2 = c(1, 0.5))
  expect_s3_class(p, "ggplot")
})

test_that("plt_sens defaults to the type that matches the backend", {
  plt_test_deps("sensemakr", "survival", "tipr")
  expect_s3_class(plt_sens(plt_lm_res()), "ggplot")
  p <- plt_sens(plt_cox_res())
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$x, "Confounder-outcome risk ratio")
})

test_that("plt_sens draws a native tipping-point curve for the Cox backend", {
  plt_test_deps("survival", "tipr")
  res <- plt_cox_res()
  p <- plt_sens(res, type = "tip")

  expect_s3_class(p, "ggplot")
  expect_gte(length(p$layers), 5L)
  expect_match(p$labels$title, "Tipping point for sex")
  # the curve must actually cross the null at the reported tipping point
  d <- ggplot2::ggplot_build(p)$data[[2]]
  expect_true(min(d$x) <= res$stats$tip_effect)
  expect_true(max(d$x) >= res$stats$tip_effect)
})

test_that("plt_sens draws the E-value curve through the reported E-value", {
  plt_test_deps("survival", "tipr")
  res <- plt_cox_res()
  p <- plt_sens(res, type = "evalue", title = "custom")

  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "custom")
  expect_equal(unname(attr(p, "plot_size")), c(7.5, 5.5))
  # On the bias-factor curve the point (E, E) must lie on the diagonal.
  b <- 1 / causalR:::.sens_hr_to_rr(res$stats$estimate, FALSE)
  e <- res$stats$evalue_point
  expect_equal(b * (e - 1) / (e - b), e, tolerance = 1e-6)
})

test_that("E-value plots omit the confidence-limit curve when the interval crosses one", {
  plt_test_deps("survival", "tipr")
  set.seed(1)
  d <- data.frame(trt = rep(0:1, each = 80), time = rexp(160),
                  status = rbinom(160, 1, 0.8))
  for (reverse in c(FALSE, TRUE)) {
    dd <- d
    if (reverse) dd$trt <- 1 - dd$trt
    res <- get_sens(dd, "trt", "status", time = "time", method = "cox",
                     evalue_args = list(rare = TRUE))
    expect_lt(res$stats$conf.low, 1)
    expect_gt(res$stats$conf.high, 1)
    expect_equal(res$stats$evalue_ci, 1)
    p <- plt_sens(res, type = "evalue")
    expect_identical(unique(as.character(p$data$which)), "Point estimate")
  }
})

test_that("plt_sens rejects a type the backend cannot draw", {
  plt_test_deps("sensemakr", "survival", "tipr")
  expect_error(plt_sens(plt_lm_res(), type = "tip"),
               "not available for method")
  expect_error(plt_sens(plt_cox_res(), type = "contour"),
               "not available for method")
})

test_that("plt_sens narrows sensitivity_of to what the backend supports", {
  plt_test_deps("iv.sensemakr", "ggplotify")
  res <- plt_iv_res()
  expect_s3_class(plt_sens(res, type = "contour"), "ggplot")
  expect_error(
    plt_sens(res, type = "contour",
             contour_args = list(sensitivity_of = "estimate")),
    "should be one of")
  expect_error(
    plt_sens(res, type = "contour", contour_args = list(nlevels = 4)),
    "unknown field")
})

test_that("IV contours receive the fitted confidence level and manual scenario", {
  plt_test_deps("iv.sensemakr", "ggplotify")
  res <- plt_iv_res(conf_level = 0.90,
                    bench_args = list(bound = c(0.01, 0.02),
                                      bound_label = "Specified scenario"))
  received <- NULL
  local_mocked_bindings(ovb_contour_plot = function(model, ...) {
    received <<- list(...)
    graphics::plot.new()
  }, .package = "iv.sensemakr")
  expect_s3_class(plt_sens(res), "ggplot")
  expect_equal(received$alpha, 0.1)
  expect_equal(received$r2zw.x, 0.01)
  expect_equal(received$r2y0w.zx, 0.02)
  expect_identical(received$bound_label, "Specified scenario")
})

test_that("IV contours use the benchmark names passed to the fitted backend", {
  plt_test_deps("iv.sensemakr", "ggplotify")
  e <- new.env()
  utils::data("card", package = "iv.sensemakr", envir = e)
  d <- e$card
  names(d)[names(d) == "black"] <- "black race"
  res <- get_sens(d, "educ", "lwage", instrument = "nearc4",
                   adj_var = c("exper", "expersq", "black race", "south", "smsa"),
                   bench_var = "black race", method = "iv")
  received <- NULL
  local_mocked_bindings(ovb_contour_plot = function(model, ...) {
    received <<- list(...)
    graphics::plot.new()
  }, .package = "iv.sensemakr")
  expect_s3_class(plt_sens(res), "ggplot")
  expect_identical(attr(res, "analysis")$bench_var, "black race")
  expect_identical(received$benchmark_covariates, "black.race")
})

test_that("plt_sens validates x, lim and save", {
  plt_test_deps("survival", "tipr")
  res <- plt_cox_res()
  expect_error(plt_sens(res$stats), "must be a `sens_res` object")
  expect_error(plt_sens(res, lim = 1), "must be `NULL` or two numbers")
  expect_error(plt_sens(res, save = "out.pdf"), "must be `NULL` or a list")
  expect_s3_class(plt_sens(res, save = list()), "ggplot")
  expect_s3_class(plt_sens(res, save = NULL), "ggplot")
})

test_that("plt_sens writes a PDF through RegR::save_plt when save is non-empty", {
  plt_test_deps("survival", "tipr", "RegR")
  f <- tempfile(fileext = ".pdf")
  on.exit(unlink(f), add = TRUE)
  p <- plt_sens(plt_cox_res(), type = "tip", save = list(filename = f))
  expect_s3_class(p, "ggplot")
  expect_true(file.exists(f))
  expect_gt(file.size(f), 0)
})

test_that("plt_sens honours the binary-confounder parameterisation", {
  plt_test_deps("survival", "tipr")
  res <- plt_cox_res(evalue_args = list(confounder = "binary",
                                        exposed_prev = 0.5,
                                        unexposed_prev = 0.2))
  p <- plt_sens(res, type = "tip")
  expect_match(p$labels$x, "prevalence 0.5 vs 0.2")
})
