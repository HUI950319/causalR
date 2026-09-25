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
  expect_named(res, c("stats", "subgroup", "importance", "calibration",
                      "data", "fit"))
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

test_that("every factor enters the forest with one column per level", {
  res <- hte_bin()
  expect_identical(colnames(res$fit$X.orig),
                   c("age", "x2", "sexF", "sexM", "stageI", "stageII", "stageIII"))
  expect_identical(attr(res$fit$X.orig, "assign"), c(1L, 2L, 3L, 3L, 4L, 4L, 4L))
})

test_that("factor_encoding = 'integer' codes each factor as one column of level codes", {
  skip_if_not_installed("grf")
  expect_identical(names(formals(get_hte))[6:7], c("method", "factor_encoding"))
  d <- hte_bin_data()
  d$stage <- factor(d$stage, levels = c("III", "II", "I"))  # codes follow the levels
  d$grp   <- sample(c("b", "a"), nrow(d), replace = TRUE)    # character: alphabetical
  res <- get_hte(d, "z", adj_var = c("age", "sex", "stage", "grp"), surv = "y",
                 factor_encoding = "integer", grf_args = hte_args)
  X <- res$fit$X.orig
  expect_identical(colnames(X), c("age", "sex", "stage", "grp"))
  expect_equal(unname(X[, "stage"]),
               as.numeric(match(as.character(d$stage), c("III", "II", "I"))))
  expect_equal(unname(X[, "grp"]), as.numeric(match(d$grp, c("a", "b"))))
  expect_identical(res$importance$n_col, rep(1L, 4L))
  expect_identical(attr(res, "analysis")$factor_encoding, "integer")

  # subgroups and p_het still work on the original levels
  res <- suppressMessages(get_hte(d, "z", sub_var = "stage", adj_var = "age",
                                  surv = "y", factor_encoding = "integer",
                                  grf_args = hte_args))
  expect_identical(unique(res$subgroup$level), c("III", "II", "I"))
  imp <- res$importance
  expect_equal(imp$p_het[imp$variable == "stage"], res$subgroup$p_inter[1])

  expect_error(get_hte(d, "z", adj_var = "age", surv = "y",
                       factor_encoding = "target"), "should be one of")
})

test_that("$importance sums grf variable importance back to each covariate", {
  res <- hte_bin()
  imp <- res$importance
  vi  <- as.numeric(grf::variable_importance(res$fit))
  src <- c("age", "x2", "sex", "stage")[attr(res$fit$X.orig, "assign")]

  expect_named(imp, c("variable", "importance", "n_col", "df", "p_het"))
  # MLR::plt_bar_per() reads the first two columns: one categorical, one numeric
  expect_type(imp$variable, "character")
  expect_type(imp$importance, "double")
  expect_setequal(imp$variable, c("age", "x2", "sex", "stage"))
  expect_equal(imp$importance, as.numeric(tapply(vi, src, sum)[imp$variable]))
  expect_equal(imp$n_col, unname(c(age = 1L, x2 = 1L, sex = 2L, stage = 3L)[imp$variable]))
  expect_equal(sum(imp$importance), 1)
  expect_false(is.unsorted(rev(imp$importance)))
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

test_that("$calibration reproduces grf::test_calibration() on a causal forest", {
  res <- hte_bin()
  cal <- res$calibration
  ref <- grf::test_calibration(res$fit)
  expect_named(cal, c("term", "estimate", "std.error", "statistic", "p.value"))
  expect_identical(cal$term, rownames(ref))
  expect_equal(unname(as.matrix(cal[-1])), unname(unclass(ref)[, 1:4]))

  # grf's observation weights: sample weights, or equal weight per cluster
  d <- hte_bin_data(n = 400L)
  for (ga in list(list(clusters = rep(1:40, length.out = 400),
                       sample.weights = rep(c(0.5, 1, 1.5), length.out = 400)),
                  list(clusters = rep(1:40, length.out = 400),
                       equalize.cluster.weights = TRUE))) {
    res <- get_hte(d, "z", adj_var = c("age", "x2"), surv = "y",
                   grf_args = c(hte_args, ga))
    expect_equal(unname(as.matrix(res$calibration[-1])),
                 unname(unclass(grf::test_calibration(res$fit))[, 1:4]))
  }
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

test_that("survival diff, ratio and OR intervals cover the true marginal effects", {
  skip_if_not_installed("grf")
  # 100 replicates of this design covered at 0.97-0.99 for every measure
  d   <- hte_surv_data(n = 2000L, seed = 1)
  res <- get_hte(d, cat_var = "z", adj_var = c("age", "x2"), surv = TRUE,
                 time = 60, measure = c("diff", "ratio", "OR"),
                 grf_args = list(num.trees = 500, seed = 1))
  # S_z(60) under hazard 0.02 exp(0.3 x2 - 0.5 z) with x2 ~ N(0, 1)
  S <- function(z) stats::integrate(function(x)
    exp(-0.02 * exp(0.3 * x - 0.5 * z) * 60) * stats::dnorm(x), -Inf, Inf)$value
  f1 <- 1 - S(1)
  f0 <- 1 - S(0)
  truth <- c(diff = S(1) - S(0), ratio = f1 / f0,
             OR = (f1 / (1 - f1)) / (f0 / (1 - f0)))
  st <- res$stats
  expect_identical(st$measure, c("diff", "ratio", "OR"))
  expect_true(all(st$conf.low < truth[st$measure]))
  expect_true(all(st$conf.high > truth[st$measure]))
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

test_that("survival: $calibration detects an effect that differs by sex", {
  skip_if_not_installed("grf")
  # 200 replicates (n = 1000, S(60) difference) with a constant effect
  # rejected at 0.050; grf's own test, on the uncensored 0/1 outcome, at 0.030
  set.seed(1)
  n <- 800L
  d <- data.frame(age = stats::rnorm(n, 60, 10), x2 = stats::rnorm(n),
                  sex = factor(sample(c("F", "M"), n, replace = TRUE)))
  d$z  <- stats::rbinom(n, 1, stats::plogis(0.3 * d$x2))
  ev   <- stats::rexp(n, 0.02 * exp(0.3 * d$x2 -
                                      d$z * (0.1 + 1.2 * (d$sex == "M"))))
  cens <- pmin(stats::rexp(n, 0.01), 120)
  d$time <- pmin(ev, cens)
  d$DSS  <- as.integer(ev <= cens)
  res <- get_hte(d, cat_var = "z", adj_var = c("age", "x2", "sex"),
                 surv = TRUE, time = 60, grf_args = hte_args)
  cal <- res$calibration
  expect_identical(cal$term, c("mean.forest.prediction",
                               "differential.forest.prediction"))
  expect_lt(cal$p.value[2], 0.01)
  expect_true(all(abs(cal$estimate - 1) < stats::qnorm(0.975) * cal$std.error))

  expect_error(.test_calibration(structure(list(), class = "regression_forest")),
               "causal_forest or causal_survival_forest")
})

test_that("survival: a time point past either arm's follow-up stops", {
  skip_if_not_installed("grf")
  d <- hte_surv_data()                          # censoring ends at 120
  expect_error(get_hte(d, "z", adj_var = c("age", "x2"), surv = TRUE,
                       grf_args = hte_args),    # default time = 120
               "followed beyond `time` = 120")

  late <- d$z == 1 & d$time > 50                # treated follow-up ends at 50
  d$time[late] <- 50
  d$DSS[late]  <- 0L
  expect_error(get_hte(d, "z", adj_var = c("age", "x2"), surv = TRUE,
                       time = 60, grf_args = hte_args),
               "No patient with `z` = 1 is followed beyond `time` = 60")
})

test_that("survival: few patients past `time` warn; a subgroup with none is NA", {
  skip_if_not_installed("grf")
  d   <- hte_surv_data()
  cut <- which(d$z == 1 & d$time > 60)[-(1:3)]  # three treated stay past 60
  d$time[cut] <- 59
  d$DSS[cut]  <- 0L
  expect_warning(get_hte(d, "z", adj_var = c("age", "x2"), surv = TRUE,
                         time = 60, grf_args = hte_args),
                 "Few patients are followed beyond `time` = 60: 3 with `z` = 1")

  d   <- hte_surv_data()
  cut <- d$z == 1 & d$sex == "M" & d$time > 50
  d$time[cut] <- 50
  d$DSS[cut]  <- 0L
  expect_warning(
    res <- get_hte(d, "z", sub_var = "sex", adj_var = c("age", "x2"),
                   surv = TRUE, time = 60, grf_args = hte_args),
    "sex = M: no patient in one arm")
  sub <- res$subgroup
  expect_true(is.na(sub$estimate[sub$level == "M"]))
  expect_false(is.na(sub$estimate[sub$level == "F"]))
})

test_that("a propensity of exactly 0 or 1 stops with a positivity message", {
  skip_if_not_installed("grf")
  set.seed(3)
  n <- 1000L
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  d$z <- ifelse(d$x2 > 1, 1L, stats::rbinom(n, 1, 0.5))  # x2 > 1 always treated
  d$y <- stats::rnorm(n) + 0.3 * d$z
  expect_error(get_hte(d, "z", adj_var = c("x1", "x2"), surv = "y",
                       measure = c("diff", "ratio"), grf_args = hte_args),
               "propensity of `z` is exactly 0 or 1")
  # a bounded, known propensity passes
  expect_s3_class(get_hte(d, "z", adj_var = c("x1", "x2"), surv = "y",
                          grf_args = c(hte_args, W.hat = 0.5)), "hte_res")
})

test_that("missing covariates go to grf; only cat_var and the outcome drop rows", {
  skip_if_not_installed("grf")
  d <- hte_bin_data()
  set.seed(7)
  d$grade <- factor(sample(c("G1", "G2", "G3"), nrow(d), replace = TRUE))
  d$grade[sample(nrow(d), 200)] <- NA
  d$age[sample(nrow(d), 80)]    <- NA
  d$y[1:5]  <- NA
  d$z[6:10] <- NA

  res <- suppressMessages(get_hte(d, "z", sub_var = "grade",
                                  adj_var = c("age", "x2", "sex"), surv = "y",
                                  grf_args = hte_args))
  expect_identical(res$stats$n, nrow(d) - 10L)
  expect_identical(nrow(res$data), nrow(d) - 10L)
  expect_true(anyNA(res$fit$X.orig[, "age"]))
  # a missing level leaves every indicator column of that factor missing
  na_grade <- is.na(res$data$grade)
  expect_true(all(is.na(res$fit$X.orig[na_grade,
                                       c("gradeG1", "gradeG2", "gradeG3")])))

  # subgroups and p_het use the patients whose value is observed
  sub <- res$subgroup
  expect_identical(sub$level, c("G1", "G2", "G3"))
  expect_identical(sum(sub$n), sum(!na_grade))
  imp <- res$importance
  expect_equal(imp$p_het[imp$variable == "grade"], sub$p_inter[1])
  expect_false(is.na(imp$p_het[imp$variable == "age"]))
})

test_that("per-row grf_args given for every row of `data` follow the rows kept", {
  skip_if_not_installed("grf")
  d <- hte_bin_data(n = 400L)
  d$y[1:5] <- NA
  res <- get_hte(d, "z", adj_var = c("age", "x2"), surv = "y",
                 grf_args = c(hte_args,
                              list(W.hat    = rep(0.5, nrow(d)),
                                   clusters = rep(1:40, length.out = nrow(d)))))
  expect_identical(nrow(res$data), nrow(d) - 5L)
  expect_equal(res$fit$W.hat, rep(0.5, nrow(d) - 5L))
  expect_length(res$fit$clusters, nrow(d) - 5L)

  # fields already matching the rows analysed pass unchanged
  res <- get_hte(d, "z", adj_var = c("age", "x2"), surv = "y",
                 grf_args = c(hte_args, list(W.hat = rep(0.4, nrow(d) - 5L))))
  expect_equal(res$fit$W.hat, rep(0.4, nrow(d) - 5L))
})

test_that("grf's overlap warnings collapse into one per call", {
  skip_if_not_installed("grf")
  d <- hte_bin_data()
  e <- pmin(pmax(stats::plogis(3 * d$x2), 0.02), 0.98)  # bounded, poor overlap
  set.seed(5)
  d$z <- stats::rbinom(nrow(d), 1, e)
  w <- character()
  res <- withCallingHandlers(
    suppressMessages(get_hte(d, "z", sub_var = c("sex", "stage"),
                             adj_var = c("age", "x2"), surv = "y",
                             estimand = c("ATE", "ATT", "ATC"),
                             grf_args = c(hte_args, list(W.hat = e)))),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  expect_length(w, 1L)
  expect_match(w, "Estimated propensities of `z` range from 0.020 to 0.980")
  expect_equal(attr(res, "analysis")$ps_range, c(0.02, 0.98))
  expect_true(any(grepl("propensity range 0.020 to 0.980",
                        capture.output(print(res)), fixed = TRUE)))
})

test_that("grf_args naming a field get_hte() sets points to the argument", {
  skip_if_not_installed("grf")
  expect_error(get_hte(hte_surv_data(n = 200L), "z", adj_var = "x2",
                       surv = TRUE, time = 60, grf_args = list(horizon = 60)),
               "`grf_args` cannot set `horizon`: get_hte() sets `horizon` from `time`.",
               fixed = TRUE)
  expect_error(get_hte(hte_bin_data(n = 200L), "z", adj_var = "x2", surv = "y",
                       grf_args = list(X = matrix(0), W = 1)),
               "cannot set `X`, `W`: get_hte() sets `X` from `adj_var` / `sub_var` and `W` from `cat_var`.",
               fixed = TRUE)
})

test_that("a date covariate stops instead of becoming one column per date", {
  skip_if_not_installed("grf")
  d <- hte_bin_data(n = 200L)
  d$dx_date <- as.Date("2010-01-01") + seq_len(nrow(d))
  expect_error(get_hte(d, "z", adj_var = c("age", "dx_date"), surv = "y",
                       grf_args = hte_args),
               "`dx_date` must be numeric, factor, character or logical")
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
  expect_true(any(grepl("Calibration (one-sided p):", out, fixed = TRUE)))
  expect_true(any(grepl("differential.forest.prediction", out, fixed = TRUE)))
  expect_true(any(grepl("1 - S(t)", out, fixed = TRUE)))
  expect_true(any(grepl(".dr_score", out, fixed = TRUE)))
})

test_that("unavailable diagnostics preserve valid average effects", {
  skip_if_not_installed("grf")
  d <- hte_bin_data(n = 400L)
  d$tied <- c(rep(0, 394), 1:6)
  expect_warning(r <- get_hte(d, "z", adj_var = c("age", "tied"), surv = "y",
                              grf_args = hte_args), "tied.*unavailable")
  expect_equal(r$stats$estimate,
               unname(grf::average_treatment_effect(r$fit)[["estimate"]]))
  expect_true(is.na(r$importance$p_het[r$importance$variable == "tied"]))
  expect_warning(p <- plt_hte_dep(r, x_var = "tied", display = "dr"),
                 "tied.*unavailable")
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")

  d$y <- 0
  for (v in c("age", "sex")) {
    warnings <- character()
    r <- withCallingHandlers(
      get_hte(d, "z", adj_var = v, surv = "y", grf_args = hte_args),
      warning = function(e) {
        warnings <<- c(warnings, conditionMessage(e))
        invokeRestart("muffleWarning")
      })
    expect_match(warnings, "unavailable", all = TRUE)
    expect_equal(r$stats$estimate, 0)
    expect_equal(r$stats$std.error, 0)
    expect_true(all(is.na(r$importance$p_het)))
    expect_identical(r$calibration$term,
                     c("mean.forest.prediction", "differential.forest.prediction"))
    expect_true(all(is.na(r$calibration$p.value)))
  }
})

test_that("clustered subgroup comparisons include cross-group covariance", {
  skip_if_not_installed("grf")
  set.seed(122)
  cl <- rep(1:40, each = 40)
  d <- data.frame(x = stats::rnorm(1600),
                  g = factor(rep(rep(c("A", "B"), each = 20), 40)),
                  z = stats::rbinom(1600, 1, 0.5))
  d$y <- d$z * (0.3 * (d$g == "B") + stats::rnorm(40, sd = 3)[cl]) +
    stats::rnorm(1600)
  r <- get_hte(d, "z", sub_var = "g", adj_var = c("x", "g"), surv = "y",
               estimand = c("ATE", "ATT", "ATC", "ATO"),
               grf_args = c(hte_args, list(W.hat = 0.5, clusters = cl)))
  a <- r$subgroup[r$subgroup$estimand == "ATE", ]
  u <- sapply(levels(d$g), function(lv) {
    i <- which(d$g == lv)
    v <- numeric(nrow(d))
    v[i] <- (r$data$.dr_score[i] - mean(r$data$.dr_score[i])) / length(i)
    as.numeric(rowsum(v, cl)) * sqrt(40 / 39)
  })
  S <- crossprod(u)
  expect_equal(diag(S), a$std.error^2, ignore_attr = TRUE)
  p <- stats::pchisq(diff(a$estimate)^2 / sum(S * matrix(c(1, -1, -1, 1), 2)),
                     1, lower.tail = FALSE)
  expect_equal(a$p_inter, rep(p, 2))
  expect_lt(p, 0.05)
  expect_equal(r$importance$p_het[r$importance$variable == "g"], p)

  # Independent cluster-level perturbations of each group's residual regression
  # give the cross-covariance for overlap effects.
  u <- sapply(levels(d$g), function(lv) {
    i <- which(d$g == lv)
    m <- stats::lm(I(d$y[i] - r$fit$Y.hat[i]) ~ I(d$z[i] - 0.5))
    scores <- sandwich::estfun(m) %*% sandwich::bread(m) / length(i)
    as.numeric(rowsum(scores[, 2], cl[i])) *
      sqrt(40 / 39 * (length(i) - 1) / (length(i) - 2))
  })
  S <- crossprod(u)
  a <- r$subgroup[r$subgroup$estimand == "ATO", ]
  expect_equal(diag(S), a$std.error^2, ignore_attr = TRUE)
  expect_equal(a$p_inter, rep(stats::pchisq(diff(a$estimate)^2 /
    sum(S * matrix(c(1, -1, -1, 1), 2)), 1, lower.tail = FALSE), 2))
  expect_equal(attr(plt_hte_sub(r, sub_var = "g"), "subgroup")$p_inter,
               r$subgroup$p_inter[r$subgroup$estimand == "ATE"])
})

test_that("heterogeneity inherits weights and clusters from the forest", {
  skip_if_not_installed("grf")
  d <- hte_bin_data(n = 480L)
  d$age[1:8] <- NA
  cl <- rep(seq_len(40), times = rep(c(8, 16), 20))
  for (ga in list(list(sample.weights = seq(0.1, 2, length.out = 480),
                       clusters = cl),
                  list(clusters = cl, equalize.cluster.weights = TRUE))) {
    res <- get_hte(d, "z", adj_var = c("age", "x2", "sex"), surv = "y",
                   grf_args = c(hte_args, ga))
    wt <- if (!is.null(ga$sample.weights)) ga$sample.weights else
      1 / as.numeric(table(cl)[as.character(cl)])
    i <- which(!is.na(d$age))
    m <- stats::lm(.dr_score ~ splines::ns(age, 2), data = res$data[i, ],
                   weights = wt[i])
    V <- sandwich::vcovCL(m, cluster = cl[i], type = "HC3")
    b <- stats::coef(m)[-1L]
    expect_equal(res$importance$p_het[res$importance$variable == "age"],
                 stats::pchisq(drop(b %*% solve(V[-1L, -1L], b)), 2,
                               lower.tail = FALSE))

    mu <- tapply(seq_len(nrow(d)), d$sex,
                 function(i) stats::weighted.mean(res$data$.dr_score[i], wt[i]))
    influence <- sapply(levels(d$sex), function(lv) {
      i <- which(d$sex == lv)
      u <- numeric(nrow(d))
      u[i] <- wt[i] * (res$data$.dr_score[i] - mu[lv]) / sum(wt[i])
      k <- length(unique(cl[i][wt[i] > 0]))
      as.numeric(rowsum(u, cl)) * sqrt(k / (k - 1))
    })
    S <- crossprod(influence)
    expect_equal(res$importance$p_het[res$importance$variable == "sex"],
                 stats::pchisq(diff(mu)^2 / (S[1, 1] + S[2, 2] - 2 * S[1, 2]),
                               1, lower.tail = FALSE), ignore_attr = TRUE)
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p <- plt_hte_dep(res, x_var = "sex", display = "dr")
      layer <- Filter(function(l) inherits(l$geom, "GeomPointrange"), p$layers)[[1]]
      expect_equal(layer$data$estimate, as.numeric(mu))
      expect_equal((layer$data$conf.high - layer$data$conf.low) /
                     (2 * stats::qnorm(0.975)), sqrt(diag(S)), ignore_attr = TRUE)
    }
  }
})
