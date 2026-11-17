# A user-set `axis_args$xlim` narrower than the data clipped those intervals
# without saying so, while `xlim = NULL` takes the range from the data
# (paper-pac pac/TRT: a 183-patient subgroup reached +24.2% on a +-20% axis and
# only a read of the PDF text layer caught it). The plot still clips; it now
# names the rows and the limit that would hold them.

xlim_cache <- new.env()

xlim_res <- function() {
  skip_if_not_installed("grf")
  skip_if_not_installed("sandwich")
  if (is.null(xlim_cache$res)) {
    set.seed(20261117)
    n <- 400L
    d <- data.frame(age = round(stats::runif(n, 20, 85)),
                    sex = factor(sample(c("F", "M"), n, replace = TRUE)))
    d$z <- stats::rbinom(n, 1, 0.5)
    d$y <- stats::rbinom(n, 1, stats::plogis(-1 + 0.02 * (d$age - 50) +
                                               d$z * (0.2 + 0.8 * (d$sex == "M"))))
    xlim_cache$res <- suppressMessages(get_hte(
      d, cat_var = "z", sub_var = "sex", adj_var = "age", surv = "y",
      grf_args = list(num.trees = 200, seed = 1)))
  }
  xlim_cache$res
}

test_that("plt_hte_sub() warns when a given xlim clips a confidence interval", {
  res <- xlim_res()
  expect_warning(
    suppressMessages(plt_hte_sub(res, sub_var = "sex",
                                 axis_args = list(xlim = c(-0.001, 0.001)))),
    "clipped")
})

test_that("a wide enough xlim and the default auto range stay quiet", {
  res <- xlim_res()
  expect_no_warning(
    suppressMessages(plt_hte_sub(res, sub_var = "sex",
                                 axis_args = list(xlim = c(-1, 1)))))
  expect_no_warning(suppressMessages(plt_hte_sub(res, sub_var = "sex")))
})
