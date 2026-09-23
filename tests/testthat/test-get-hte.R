# grf is a Suggests backend; every test that fits a forest skips without it.
# Fits are cached per scenario so each forest is grown once per test run.

hte_args  <- list(num.trees = 300, seed = 1)
hte_cache <- new.env()

hte_bin_data <- function(n = 800L, seed = 20260923) {
  set.seed(seed)
  d <- data.frame(age   = stats::rnorm(n, 60, 10),
                  x2    = stats::rnorm(n),
                  sex   = factor(sample(c("F", "M"), n, replace = TRUE)),
                  stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
  d$z  <- stats::rbinom(n, 1, stats::plogis(0.02 * (d$age - 60) - 0.3 * d$x2))
  d$p0 <- stats::plogis(-1 + 0.5 * d$x2)
  d$p1 <- stats::plogis(-1 + 0.5 * d$x2 + 0.3 + 0.8 * (d$sex == "M"))
  d$y  <- stats::rbinom(n, 1, ifelse(d$z == 1, d$p1, d$p0))
  d
}

hte_surv_data <- function(n = 800L, seed = 20260923) {
  set.seed(seed)
  d <- data.frame(age = stats::rnorm(n, 60, 10),
                  x2  = stats::rnorm(n),
                  sex = factor(sample(c("F", "M"), n, replace = TRUE)))
  d$z  <- stats::rbinom(n, 1, stats::plogis(0.3 * d$x2))
  ev   <- stats::rexp(n, 0.02 * exp(0.3 * d$x2 - 0.5 * d$z))
  cens <- pmin(stats::rexp(n, 0.01), 120)
  d$time <- pmin(ev, cens)
  d$DSS  <- as.integer(ev <= cens)
  d
}

hte_bin <- function() {
  skip_if_not_installed("grf")
  # ratio / OR for ATT, ATC and ATO are skipped with a message
  if (is.null(hte_cache$bin))
    hte_cache$bin <- suppressMessages(get_hte(
      hte_bin_data(), cat_var = "z", adj_var = c("age", "x2", "sex", "stage"),
      surv = "y", estimand = c("ATE", "ATT", "ATC", "ATO"),
      measure = c("diff", "ratio", "OR"), grf_args = hte_args))
  hte_cache$bin
}

hte_surv <- function() {
  skip_if_not_installed("grf")
  if (is.null(hte_cache$surv))
    hte_cache$surv <- get_hte(hte_surv_data(), cat_var = "z",
                              adj_var = c("age", "x2", "sex"), surv = TRUE,
                              time = 60, measure = c("diff", "ratio", "OR"),
                              grf_args = hte_args)
  hte_cache$surv
}

# Arm-specific AIPW means rebuilt from the public fields of a grf forest.
arm_means <- function(fit, idx = TRUE) {
  tau <- as.numeric(predict(fit)$predictions)
  e   <- fit$W.hat
  wc  <- fit$W.orig - e
  r   <- if (inherits(fit, "causal_survival_forest"))
    fit[["_psi"]]$numerator / wc - wc * tau else fit$Y.orig - fit$Y.hat - wc * tau
  g1 <- fit$Y.hat + (1 - e) * tau + fit$W.orig / e * r
  g0 <- fit$Y.hat - e * tau + (1 - fit$W.orig) / (1 - e) * r
  c(mean(g1[idx]), mean(g0[idx]))
}

row_of <- function(tbl, estimand, measure)
  tbl[tbl$estimand == estimand & tbl$measure == measure, , drop = FALSE]


test_that("get_hte returns the documented structure", {
  res <- hte_bin()

  expect_s3_class(res, "hte_res")
  expect_named(res, c("stats", "subgroup", "data", "fit"))
  expect_s3_class(res$fit, "causal_forest")
  expect_null(res$subgroup)
  expect_named(res$stats, c("method", "estimand", "measure", "estimate",
                            "std.error", "conf.low", "conf.high", "p.value",
                            "n", "n_treat"))
  # ratio and OR exist only for the ATE: 4 diff rows + 2 relative rows
  expect_identical(nrow(res$stats), 6L)
  expect_true(all(c(".cate", ".dr_score") %in% names(res$data)))
  expect_identical(attr(res, "analysis")$outcome_type, "binary")
})

test_that("diff reproduces grf for every estimand and $data matches grf", {
  res <- hte_bin()
  map <- c(ATE = "all", ATT = "treated", ATC = "control", ATO = "overlap")

  for (e in names(map)) {
    ref <- grf::average_treatment_effect(res$fit, target.sample = map[[e]])
    row <- row_of(res$stats, e, "diff")
    expect_equal(row$estimate, unname(ref[["estimate"]]))
    expect_equal(row$std.error, unname(ref[["std.err"]]))
  }
  expect_equal(res$data$.cate, as.numeric(predict(res$fit)$predictions))
  expect_equal(res$data$.dr_score, as.numeric(grf::get_scores(res$fit)))
})

test_that("ratio and OR come from arm-specific AIPW scores", {
  res <- hte_bin()
  m   <- arm_means(res$fit)

  rr <- row_of(res$stats, "ATE", "ratio")
  or <- row_of(res$stats, "ATE", "OR")
  expect_equal(rr$estimate, m[1] / m[2])
  expect_equal(or$estimate, (m[1] / (1 - m[1])) / (m[2] / (1 - m[2])))
  # log-scale interval: symmetric around log(estimate)
  expect_equal(log(rr$conf.high) - log(rr$estimate),
               log(rr$estimate) - log(rr$conf.low))
})

test_that("the risk-ratio interval covers the true marginal risk ratio", {
  skip_if_not_installed("grf")
  d   <- hte_bin_data(n = 3000L, seed = 1)
  res <- get_hte(d, cat_var = "z", adj_var = c("age", "x2", "sex"),
                 surv = "y", measure = "ratio",
                 grf_args = list(num.trees = 500, seed = 1))
  truth <- mean(d$p1) / mean(d$p0)
  expect_lt(res$stats$conf.low, truth)
  expect_gt(res$stats$conf.high, truth)
})

test_that("subgroups use the grf subset estimate and a descriptive CATE mean", {
  skip_if_not_installed("grf")
  d <- hte_bin_data()
  expect_message(
    res <- get_hte(d, cat_var = "z", sub_var = c("sex", "stage"),
                   adj_var = c("age", "x2"), surv = "y",
                   measure = c("diff", "ratio"), grf_args = hte_args),
    "sex")
  sub <- res$subgroup

  expect_named(sub, c("sub_var", "level", "estimand", "measure", "n",
                      "n_treat", "estimate", "std.error", "conf.low",
                      "conf.high", "p.value", "cate_mean", "p_inter"))
  expect_identical(nrow(sub), 10L)   # (2 + 3 levels) x (diff, ratio)
  expect_true(all(c("sex", "stage") %in% attr(res, "analysis")$covariates))

  for (lv in c("F", "M")) {
    idx <- d$sex == lv
    ref <- grf::average_treatment_effect(res$fit, subset = idx)
    row <- sub[sub$sub_var == "sex" & sub$level == lv & sub$measure == "diff", ]
    expect_equal(row$estimate, unname(ref[["estimate"]]))
    expect_equal(row$cate_mean, mean(res$data$.cate[idx]))
    m <- arm_means(res$fit, idx)
    rr <- sub[sub$sub_var == "sex" & sub$level == lv & sub$measure == "ratio", ]
    expect_equal(rr$estimate, m[1] / m[2])
    expect_true(is.na(rr$cate_mean))
  }
  expect_true(all(sub$p_inter >= 0 & sub$p_inter <= 1))
  expect_identical(length(unique(sub$p_inter[sub$sub_var == "sex" &
                                               sub$measure == "diff"])), 1L)
})

test_that("survival: S(t) difference matches grf and ratio / OR use the event risk", {
  res <- hte_surv()
  expect_s3_class(res$fit, "causal_survival_forest")
  expect_identical(res$fit$target, "survival.probability")
  expect_identical(res$fit$horizon, 60)

  ref <- grf::average_treatment_effect(res$fit)
  expect_equal(row_of(res$stats, "ATE", "diff")$estimate,
               unname(ref[["estimate"]]))
  expect_equal(res$data$.dr_score, as.numeric(grf::get_scores(res$fit)))

  s <- arm_means(res$fit)
  f <- 1 - s                          # event risk 1 - S(t)
  expect_equal(row_of(res$stats, "ATE", "ratio")$estimate, f[1] / f[2])
  expect_equal(row_of(res$stats, "ATE", "OR")$estimate,
               (f[1] / (1 - f[1])) / (f[2] / (1 - f[2])))
})

test_that("survival RMST is reached through grf_args$target; OR is skipped", {
  skip_if_not_installed("grf")
  expect_message(
    res <- get_hte(hte_surv_data(), cat_var = "z", adj_var = c("age", "x2"),
                   surv = TRUE, time = 60, measure = c("diff", "ratio", "OR"),
                   grf_args = c(hte_args, target = "RMST")),
    "OR")
  expect_identical(res$fit$target, "RMST")
  expect_identical(res$stats$measure, c("diff", "ratio"))
  s <- arm_means(res$fit)
  expect_equal(row_of(res$stats, "ATE", "ratio")$estimate, s[1] / s[2])
})

test_that("continuous outcomes give a ratio of means but no OR", {
  skip_if_not_installed("grf")
  d <- hte_bin_data()
  d$cost <- exp(1 + 0.3 * d$x2 + 0.2 * d$z + stats::rnorm(nrow(d), sd = 0.5))
  res <- get_hte(d, cat_var = "z", adj_var = c("age", "x2"), surv = "cost",
                 measure = c("diff", "ratio"), grf_args = hte_args)
  expect_identical(attr(res, "analysis")$outcome_type, "continuous")
  s <- arm_means(res$fit)
  expect_equal(row_of(res$stats, "ATE", "ratio")$estimate, s[1] / s[2])

  expect_error(get_hte(d, cat_var = "z", adj_var = "age", surv = "cost",
                       measure = "OR", grf_args = hte_args),
               "No requested")
})

test_that("relative measures outside the ATE are skipped with a message", {
  skip_if_not_installed("grf")
  expect_message(
    res <- get_hte(hte_bin_data(), cat_var = "z", adj_var = c("age", "x2"),
                   surv = "y", estimand = c("ATE", "ATT"),
                   measure = c("diff", "ratio"), grf_args = hte_args),
    "ATT")
  expect_identical(paste(res$stats$estimand, res$stats$measure),
                   c("ATE diff", "ATT diff", "ATE ratio"))
})

test_that("invalid requests stop before any forest is grown", {
  skip_if_not_installed("grf")
  d <- hte_bin_data(n = 200L)

  expect_error(get_hte(d, "z", adj_var = "age", surv = FALSE), "competing")
  expect_error(get_hte(d, "z", adj_var = "age", surv = TRUE), "`time` and `DSS`")
  expect_error(get_hte(d, "z", adj_var = "age", surv = "y", time = 60),
               "only applies")
  expect_error(get_hte(d, "z", sub_var = "age", surv = "y"),
               "more than 5 distinct")
  expect_error(get_hte(d, "stage", adj_var = "age", surv = "y"),
               "`cat_var` column `stage`")
  expect_error(get_hte(d, "z", adj_var = "age", surv = "y",
                       grf_args = list(num_trees = 10)),
               "unknown field")
  expect_error(get_hte(d, "z", adj_var = "age", surv = "y", measure = "ratio",
                       grf_args = list(clusters = rep(1:4, 50))),
               "clusters")
  expect_error(get_hte(d, "z", adj_var = "y", surv = "y"), "outcome")
  d$y12 <- d$y + 1
  expect_error(get_hte(d, "z", adj_var = "age", surv = "y12"), "0/1")
})

test_that("print shows the analysis header and the event-risk note", {
  res <- hte_surv()
  out <- capture.output(print(res))
  expect_match(out[1], "<hte_res> causal_survival_forest")
  expect_true(any(grepl("1 - S(t)", out, fixed = TRUE)))
  expect_true(any(grepl(".dr_score", out, fixed = TRUE)))
})
