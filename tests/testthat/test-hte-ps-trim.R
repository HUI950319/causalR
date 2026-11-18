# get_hte(ps = , trim = ): a fitted propensity for grf, and overlap trimming on it.
# grf grows its own regression forest for e(x) and reaches exactly 0 or 1 wherever
# a covariate region holds one arm only, which leaves the doubly robust scores
# undefined; a logistic score cannot, and `trim` keeps the overlap region.

mk_ps_data <- function(n = 300L, seed = 11L) {
  set.seed(seed)
  d <- data.frame(age = stats::rnorm(n), x2 = stats::rbinom(n, 1L, 0.4))
  d$z    <- stats::rbinom(n, 1L, stats::plogis(-0.4 + 2.5 * d$age))
  d$time <- stats::rexp(n, 0.12 + 0.05 * d$z)
  d$DSS  <- stats::rbinom(n, 1L, 0.85)
  d
}
ps_glm <- function(d) unname(stats::fitted(stats::glm(z ~ age + x2,
                                                      family = stats::binomial(),
                                                      data = d)))

test_that("get_hte() hands grf a fitted propensity and leaves the default alone", {
  skip_if_not_installed("grf")
  skip_if_not_installed("WeightIt")
  d  <- mk_ps_data()
  ga <- list(num.trees = 50, seed = 1)

  res <- get_hte(d, cat_var = "z", adj_var = c("age", "x2"), surv = TRUE,
                 time = 4, ps = "glm", grf_args = ga)
  a <- attr(res, "analysis")
  expect_identical(a$ps, "glm")
  expect_identical(a$n_trimmed, 0L)
  expect_identical(a$n, nrow(d))
  expect_equal(unname(res$fit$W.hat), ps_glm(d), tolerance = 1e-6)

  # Default: grf estimates e(x) with its own regression forest, as before.
  res0 <- get_hte(d, cat_var = "z", adj_var = c("age", "x2"), surv = TRUE,
                  time = 4, grf_args = ga)
  expect_null(attr(res0, "analysis")$ps)
  expect_identical(attr(res0, "analysis")$n_trimmed, 0L)
  expect_false(isTRUE(all.equal(unname(res0$fit$W.hat), ps_glm(d),
                                tolerance = 1e-6)))
})

test_that("trim keeps the overlap region and records how many rows went", {
  skip_if_not_installed("grf")
  skip_if_not_installed("WeightIt")
  d  <- mk_ps_data()
  ga <- list(num.trees = 50, seed = 1)
  p  <- ps_glm(d)
  out <- sum(p < 0.1 | p > 0.9)
  expect_gt(out, 0L)

  res <- get_hte(d, cat_var = "z", adj_var = c("age", "x2"), surv = TRUE,
                 time = 4, ps = "glm", trim = c(0.1, 0.9), grf_args = ga)
  a <- attr(res, "analysis")
  expect_identical(a$n_trimmed, out)
  expect_identical(a$n, nrow(d) - out)
  expect_identical(nrow(res$data), nrow(d) - out)
  expect_true(all(res$fit$W.hat >= 0.1 & res$fit$W.hat <= 0.9))
  expect_match(a$ps, "^glm, trimmed")
})

test_that("a supplied propensity, and the checks on ps / trim", {
  skip_if_not_installed("grf")
  d  <- mk_ps_data()
  ga <- list(num.trees = 50, seed = 1)
  p  <- ps_glm(d)

  # One propensity per row, as a vector and as a column, both reach grf.
  res <- get_hte(d, cat_var = "z", adj_var = c("age", "x2"), surv = TRUE,
                 time = 4, ps = p, grf_args = ga)
  expect_identical(attr(res, "analysis")$ps, "supplied")
  expect_equal(unname(res$fit$W.hat), p, tolerance = 1e-10)

  d2 <- d
  d2$my_ps <- p
  res2 <- get_hte(d2, cat_var = "z", adj_var = c("age", "x2"), surv = TRUE,
                  time = 4, ps = "my_ps", grf_args = ga)
  expect_identical(attr(res2, "analysis")$ps, "column `my_ps`")

  expect_error(
    get_hte(d, cat_var = "z", adj_var = c("age", "x2"), surv = TRUE, time = 4,
            ps = "glm", grf_args = c(ga, list(W.hat = p))),
    "not both")
  expect_error(
    get_hte(d, cat_var = "z", adj_var = c("age", "x2"), surv = TRUE, time = 4,
            trim = c(0.1, 0.9), grf_args = ga),
    "needs `ps`")
  expect_error(
    get_hte(d, cat_var = "z", adj_var = c("age", "x2"), surv = TRUE, time = 4,
            ps = "glm", trim = c(0.9, 0.1), grf_args = ga),
    "lower < upper")
  expect_error(
    get_hte(d, cat_var = "z", adj_var = c("age", "x2"), surv = TRUE, time = 4,
            ps = "nosuch", grf_args = ga),
    "'arg' should be one of|a column of")
})
