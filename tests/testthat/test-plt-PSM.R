# plt_PSM() needs halfmoon for two of its four figures and RegR only to save.

plt_psm_res <- function(balance = TRUE, ...) {
  set.seed(20260921)
  n <- 300L
  d <- data.frame(x1 = stats::rnorm(n),
                  x2 = stats::rbinom(n, 1, 0.4),
                  x3 = stats::runif(n))
  d$z <- stats::rbinom(n, 1, stats::plogis(
    -0.9 - 1.1 * d$x1 + 0.8 * d$x2 - 1.2 * d$x3))
  get_PSM(d, treat = "z", adj_var = c("x1", "x2", "x3"),
          balance = balance, ...)
}


test_that("every type returns a sized ggplot", {
  skip_if_not_installed("halfmoon")
  res <- plt_psm_res(caliper = 0.2)

  for (ty in c("love", "ess", "weight", "ps")) {
    p <- plt_PSM(res, type = ty)
    expect_s3_class(p, "ggplot")
    sz <- attr(p, "plot_size")
    expect_named(sz, c("width", "height"))
    expect_true(all(sz > 0))
    expect_no_error(ggplot2::ggplot_build(p))
  }
})


test_that("the default type follows whether balance was computed", {
  skip_if_not_installed("halfmoon")

  with_bal <- plt_psm_res()
  expect_identical(.psm_plt_spec(with_bal)$type, "love")

  no_bal <- plt_psm_res(balance = FALSE)
  expect_identical(.psm_plt_spec(no_bal)$type, "ess")
  expect_no_error(plt_PSM(no_bal))
  expect_error(plt_PSM(no_bal, type = "love"), "needs the balance table")
})


test_that("the two halfmoon-free types work without it", {
  res <- plt_psm_res(balance = FALSE)
  expect_s3_class(plt_PSM(res, "ess"), "ggplot")
  expect_s3_class(plt_PSM(res, "weight"), "ggplot")
})


test_that("the weight figure shows only matched units", {
  res  <- plt_psm_res(balance = FALSE, caliper = 0.1)
  kept <- sum(res$data$w_nearest > 0)
  expect_lt(kept, nrow(res$data))

  d <- ggplot2::ggplot_build(plt_PSM(res, "weight"))$plot$data
  expect_identical(nrow(d), kept)
  expect_true(all(d$weight > 0))
})


test_that("method selects which schemes are drawn", {
  skip_if_not_installed("optmatch")
  res <- plt_psm_res(balance = FALSE, method = c("nearest", "full"))

  full_p <- ggplot2::ggplot_build(plt_PSM(res, "ess"))$plot$data
  expect_identical(nrow(full_p), 2L)

  one <- ggplot2::ggplot_build(plt_PSM(res, "ess", method = "full"))$plot$data
  expect_identical(as.character(one$method), "full")

  expect_error(plt_PSM(res, "ess", method = "cem"), "should be one of")
})


test_that("save writes a PDF only when it is a non-empty list", {
  skip_if_not_installed("RegR")
  res <- plt_psm_res(balance = FALSE)
  dir <- withr::local_tempdir()

  expect_identical(list.files(dir), character(0))
  plt_PSM(res, "ess", save = list())
  plt_PSM(res, "ess", save = NULL)
  expect_identical(list.files(dir), character(0))

  f <- file.path(dir, "ess.pdf")
  p <- plt_PSM(res, "ess", save = list(filename = f))
  expect_s3_class(p, "ggplot")
  expect_true(file.exists(f))
  expect_gt(file.size(f), 0)

  expect_error(plt_PSM(res, "ess", save = "nope.pdf"),
               "`save` must be `NULL` or a list")
})


test_that("plt_PSM validates its input", {
  res <- plt_psm_res(balance = FALSE)
  expect_error(plt_PSM(list()), "must be a `psm_res` object")
  expect_error(plt_PSM(res, "nope"), "should be one of")
  expect_error(plt_PSM(res, "ess", threshold = c(1, 2)),
               "`threshold` must be a single number")
  expect_error(plt_PSM(res, "ps", bins = 1), "at least 2")
})
