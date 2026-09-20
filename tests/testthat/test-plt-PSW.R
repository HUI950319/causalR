# plt_PSW() needs halfmoon for two of its four figures and RegR only to save.

plt_psw_res <- function(balance = TRUE, ...) {
  set.seed(20260921)
  n <- 300L
  d <- data.frame(x1 = stats::rnorm(n),
                  x2 = stats::rbinom(n, 1, 0.4),
                  x3 = stats::runif(n))
  d$z <- stats::rbinom(n, 1, stats::plogis(
    -0.3 + 0.8 * d$x1 - 0.6 * d$x2 + 1.1 * d$x3))
  get_PSW(d, treat = "z", adj_var = c("x1", "x2", "x3"),
          balance = balance, ...)
}


test_that("every type returns a sized ggplot", {
  skip_if_not_installed("halfmoon")
  res <- plt_psw_res()

  for (ty in c("love", "ess", "weight", "ps")) {
    p <- plt_PSW(res, type = ty)
    expect_s3_class(p, "ggplot")
    sz <- attr(p, "plot_size")
    expect_named(sz, c("width", "height"))
    expect_true(all(sz > 0))
    expect_no_error(ggplot2::ggplot_build(p))
  }
})


test_that("the default type follows whether balance was computed", {
  skip_if_not_installed("halfmoon")

  with_bal <- plt_psw_res()
  expect_identical(.psw_plt_spec(with_bal)$type, "love")

  no_bal <- plt_psw_res(balance = FALSE)
  expect_identical(.psw_plt_spec(no_bal)$type, "ess")
  expect_no_error(plt_PSW(no_bal))
  expect_error(plt_PSW(no_bal, type = "love"), "needs the balance table")
})


test_that("the two halfmoon-free types work without it", {
  res <- plt_psw_res(balance = FALSE)
  expect_s3_class(plt_PSW(res, "ess"), "ggplot")
  expect_s3_class(plt_PSW(res, "weight"), "ggplot")
})


test_that("estimand selects which series are drawn", {
  res <- plt_psw_res(balance = FALSE)

  p <- plt_PSW(res, "ess", estimand = c("ATO", "ATE"))
  d <- ggplot2::ggplot_build(p)$plot$data
  expect_setequal(as.character(d$estimand), c("ATE", "ATO"))

  full <- ggplot2::ggplot_build(plt_PSW(res, "ess"))$plot$data
  expect_identical(nrow(full), 6L)

  expect_error(plt_PSW(res, "ess", estimand = "nope"), "should be one of")
})


test_that("the weight figure drops trimmed units", {
  res <- plt_psw_res(balance = FALSE,
                     trim_args = list(method = "ps", lower = 0.3,
                                      upper = 0.7))
  kept <- sum(!res$data$.trimmed)
  expect_lt(kept, nrow(res$data))

  d <- ggplot2::ggplot_build(plt_PSW(res, "weight"))$plot$data
  expect_identical(nrow(d), kept * 6L)
  expect_false(anyNA(d$weight))
})


test_that("save writes a PDF only when it is a non-empty list", {
  skip_if_not_installed("RegR")
  res <- plt_psw_res(balance = FALSE)
  dir <- withr::local_tempdir()

  expect_identical(list.files(dir), character(0))
  plt_PSW(res, "ess", save = list())
  plt_PSW(res, "ess", save = NULL)
  expect_identical(list.files(dir), character(0))

  f <- file.path(dir, "ess.pdf")
  p <- plt_PSW(res, "ess", save = list(filename = f))
  expect_s3_class(p, "ggplot")
  expect_true(file.exists(f))
  expect_gt(file.size(f), 0)

  expect_error(plt_PSW(res, "ess", save = "nope.pdf"),
               "`save` must be `NULL` or a list")
})


test_that("plt_PSW validates its input", {
  res <- plt_psw_res(balance = FALSE)
  expect_error(plt_PSW(list()), "must be a `psw_res` object")
  expect_error(plt_PSW(res, "nope"), "should be one of")
  expect_error(plt_PSW(res, "ess", threshold = c(1, 2)),
               "`threshold` must be a single number")
  expect_error(plt_PSW(res, "ps", bins = 1), "at least 2")
})
