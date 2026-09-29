icf_data <- function(n = 1600L, tau = function(d) 2 * d$X1 * d$X3, seed = 11) {
  withr::local_seed(seed)
  d <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rnorm(n),
                  X3 = rbinom(n, 1, 0.5), X4 = rnorm(n),
                  grp = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  d$z <- rbinom(n, 1, plogis(0.4 * d$X1 - 0.3 * d$X2))
  d$y <- d$X2 + d$z * tau(d) + rnorm(n)
  d
}

icf_fast <- list(n_forest = 5L, num_trees = 100L, n_folds = 3L)

icf_call <- function(d, adj_var = c("X1", "X2", "X3", "X4"), ...) {
  get_hte_icf(d, "z", adj_var, surv = "y", rule_args = icf_fast,
              grf_args = list(num.trees = 500L), ...)
}

test_that("get_hte_icf() finds an interaction and estimates it on the other half", {
  skip_if_not_installed("grf")
  d <- icf_data()
  res <- icf_call(d)
  expect_s3_class(res, "hte_icf")
  expect_named(res, c("rules", "cv", "vote", "importance", "est"))
  expect_identical(res$cv$depth[res$cv$selected], 2L)
  expect_true(all(grepl("^(X1|X3) = [01]( & (X1|X3) = [01])?$", res$rules$rule)))
  both <- grepl("X1 = 1", res$rules$rule) & grepl("X3 = 1", res$rules$rule)
  expect_identical(sum(both), 1L)
  expect_gt(res$rules$estimate[both], 1.2)
  expect_true(all(abs(res$rules$estimate[!both]) < 0.6))
  # Honest estimation: the estimation half is disjoint from the discovery half
  expect_s3_class(res$est, "hte_res")
  expect_true(is.factor(res$est$data$.rule))
  expect_identical(sum(res$rules$n_disc) + nrow(res$est$data), nrow(d))
  expect_identical(sum(res$rules$n, na.rm = TRUE), nrow(res$est$data))
  p <- suppressMessages(plt_hte_sub(res$est, sub_var = ".rule", fixed_size = FALSE))
  expect_s3_class(p, "ggplot")
})

test_that("trees splitting X1 and X3 in either order vote for one partition", {
  skip_if_not_installed("grf")
  withr::local_seed(1)
  n <- 1600L
  d <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rnorm(n),
                  X3 = rbinom(n, 1, 0.5), X4 = rnorm(n))
  d$z <- rbinom(n, 1, plogis(0.4 * d$X1 - 0.3 * d$X2))
  d$y <- d$X2 + d$z * (2 * d$X1 + 2 * d$X3) + rnorm(n)
  res <- get_hte_icf(d, "z", c("X1", "X2", "X3", "X4"), surv = "y", depth = 2,
                     rule_args = list(n_forest = 10L, num_trees = 100L,
                                      n_folds = 2L),
                     grf_args = list(num.trees = 500L), seed = 1)
  expect_identical(res$vote$n_leaf, 4L)
  # Before the leaf keys were sorted, 2 of these 10 trees split X1 first and
  # voted apart from the 8 that split X3 first
  expect_identical(res$vote$share, 1)
})

test_that("get_hte_icf() keeps every patient together when the effect is constant", {
  skip_if_not_installed("grf")
  d <- icf_data(tau = function(d) rep(0.5, nrow(d)))
  res <- icf_call(d)
  expect_identical(res$cv$depth[res$cv$selected], 0L)
  expect_identical(res$rules$rule, "All patients")
  expect_true(attr(res, "analysis")$gated)
  # Without the calibration gate the cross-validation alone lets a spurious
  # split through on these data
  open <- get_hte_icf(d, "z", c("X1", "X2", "X3", "X4"), surv = "y",
                      rule_args = c(icf_fast, gate = 1),
                      grf_args = list(num.trees = 500L))
  expect_false(attr(open, "analysis")$gated)
  expect_gt(open$cv$depth[open$cv$selected], 0L)
})

test_that("factor rules read the same under both encodings", {
  skip_if_not_installed("grf")
  d <- icf_data(tau = function(d) 2 * (d$grp == "c"))
  for (enc in c("onehot", "integer")) {
    res <- icf_call(d, adj_var = c("grp", "X2", "X4"), factor_encoding = enc)
    expect_identical(res$cv$depth[res$cv$selected], 1L)
    expect_setequal(res$rules$rule, c("grp = c", "grp != c"))
  }
})

test_that("get_hte_icf() validates its input and restores the RNG", {
  skip_if_not_installed("grf")
  d <- icf_data(n = 200L)
  expect_error(get_hte_icf(d, "z", "X1", surv = "y",
                           rule_args = list(n_tree = 5)), "unknown field")
  expect_error(get_hte_icf(d, "z", "X1", surv = "y",
                           rule_args = list(gate = 0)), "gate")
  expect_error(get_hte_icf(d, "z", "X1", surv = "y",
                           grf_args = list(W.hat = rep(0.5, 200))), "W.hat")
  expect_error(get_hte_icf(d, "z", "X1", surv = "y", depth = 0), "depth")
  expect_error(get_hte_icf(d, "z", "X1", surv = "y", split_frac = 1), "split_frac")
  d$X2[1] <- NA
  expect_error(get_hte_icf(d, "z", c("X1", "X2"), surv = "y"), "complete covariates")
  set.seed(5)
  before <- .Random.seed
  res <- get_hte_icf(icf_data(n = 400L), "z", c("X1", "X3"), surv = "y",
                     depth = 1, rule_args = list(n_forest = 2L, num_trees = 50L,
                                                 n_folds = 2L),
                     grf_args = list(num.trees = 100L))
  expect_identical(.Random.seed, before)
})

test_that("get_hte_icf() runs on a survival outcome", {
  skip_if_not_installed("grf")
  skip_on_cran()
  d <- icf_data(n = 800L)
  withr::local_seed(3)
  ev <- rexp(nrow(d), 0.05 * exp(-d$z * 1.2 * d$X1))
  cens <- rexp(nrow(d), 0.02)
  d$time <- pmin(ev, cens)
  d$DSS <- as.integer(ev <= cens)
  res <- get_hte_icf(d, "z", c("X1", "X2"), time = 12, depth = 1,
                     rule_args = list(n_forest = 2L, num_trees = 50L,
                                      n_folds = 2L),
                     grf_args = list(num.trees = 200L))
  expect_s3_class(res, "hte_icf")
  expect_s3_class(res$est$fit, "causal_survival_forest")
  expect_true(all(c(0L, 1L) %in% res$cv$depth))
})
