sens_test_deps <- function(...) {
  for (pkg in c(...)) skip_if_not_installed(pkg)
}

lm_data <- function() {
  skip_if_not_installed("sensemakr")
  sensemakr::darfur
}

lm_adj <- c("age", "farmer_dar", "herder_dar", "pastvoted",
            "hhsize_darfur", "female", "village")

cox_data <- function() {
  skip_if_not_installed("survival")
  set.seed(20260921)
  d <- stats::na.omit(survival::lung[, c("time", "status", "sex", "age", "ph.ecog")])
  d$status <- d$status - 1L
  d$sex <- factor(d$sex, labels = c("male", "female"))
  d
}

dml_data <- function() {
  skip_if_not_installed("dml.sensemakr")
  set.seed(20260921)
  e <- new.env()
  utils::data("pension", package = "dml.sensemakr", envir = e)
  pension <- e$pension
  pension[sample(nrow(pension), 900L), ]
}

iv_data <- function() {
  skip_if_not_installed("iv.sensemakr")
  e <- new.env()
  utils::data("card", package = "iv.sensemakr", envir = e)
  e$card
}

stats_cols <- c("method", "estimand", "term", "scale", "estimate", "std.error",
                "statistic", "conf.low", "conf.high", "dof", "rv_q", "rv_qa",
                "xrv_qa", "r2yd_x", "evalue_point", "evalue_ci",
                "tip_effect", "tip_n")

bounds_cols <- c("method", "estimand", "bound_label", "r2_treat", "r2_out",
                 "adj_estimate", "adj_low", "adj_high", "adj_conf.low",
                 "adj_conf.high", "adj_se", "adj_statistic")

test_that("get_sens returns sensemakr robustness values for the linear backend", {
  sens_test_deps("sensemakr")
  res <- get_sens(lm_data(), treat = "directlyharmed", outcome = "peacefactor",
                  adj_var = lm_adj, bench_var = "female", method = "lm",
                  bench_args = list(k_treat = 1:3))

  expect_s3_class(res, "sens_res")
  expect_named(res, c("stats", "bounds", "fit", "sens"))
  expect_named(res$stats, stats_cols)
  expect_named(res$bounds, bounds_cols)
  expect_identical(attr(res$stats, "method"), "lm")
  expect_equal(nrow(res$stats), 1L)
  expect_equal(nrow(res$bounds), 3L)
  expect_equal(res$stats$estimate, 0.09731582, tolerance = 1e-6)
  expect_equal(res$stats$rv_q,     0.1387764,  tolerance = 1e-6)
  expect_equal(res$stats$rv_qa,    0.07625797, tolerance = 1e-6)
  expect_equal(res$stats$xrv_qa,   0.01705328, tolerance = 1e-6)
  expect_equal(res$stats$r2yd_x,   0.02187309, tolerance = 1e-6)
  expect_equal(res$bounds$bound_label, c("1x female", "2x female", "3x female"))
  expect_equal(res$bounds$adj_estimate[3L], 0.03039602, tolerance = 1e-6)
  expect_true(all(is.na(res$stats[, c("evalue_point", "evalue_ci",
                                      "tip_effect", "tip_n")])))
})

test_that("treatment and benchmark coefficients are matched by model term", {
  skip_if_not_installed("sensemakr")
  set.seed(932)
  d <- data.frame(trt = rep(0:1, 90), trt_age = rnorm(180),
                  age = rnorm(180), age2 = rnorm(180))
  d$y <- d$trt + 0.2 * d$age + 0.4 * d$age2 + rnorm(180)
  for (factor_bench in c(FALSE, TRUE)) {
    dd <- d
    benchmark <- "age"
    if (factor_bench) {
      dd$trt <- factor(dd$trt)
      dd$age <- factor(rep(c("young", "mid", "old"), 60),
                       levels = c("young", "mid", "old"))
      benchmark <- list(age = c("agemid", "ageold"))
    }
    res <- get_sens(dd, "trt", "y", adj_var = c("trt_age", "age", "age2"),
                     bench_var = "age", method = "lm")
    term <- if (factor_bench) "trt1" else "trt"
    ref <- sensemakr::sensemakr(res$fit, treatment = term,
                                benchmark_covariates = benchmark)
    expect_identical(res$stats$term, term)
    expect_equal(res$bounds$r2_treat, ref$bounds$r2dz.x)
    expect_equal(res$bounds$r2_out, ref$bounds$r2yz.dx)
  }
})

test_that("get_sens reproduces the EValue square-root transform on the Cox backend", {
  sens_test_deps("survival", "tipr")
  res <- get_sens(cox_data(), treat = "sex", outcome = "status", time = "time",
                  adj_var = c("age", "ph.ecog"), method = "cox")

  expect_identical(res$stats$scale, "HR")
  expect_identical(res$stats$term, "sexfemale")
  expect_equal(res$stats$estimate,  0.5754446, tolerance = 1e-6)
  expect_equal(res$stats$conf.low,  0.4142130, tolerance = 1e-6)
  expect_equal(res$stats$conf.high, 0.7994351, tolerance = 1e-6)
  # Cross-checked against EValue::evalues.HR(rare = FALSE) in a session where
  # lava is not loaded; see the "Choosing rare" section of get_sens().
  expect_equal(res$stats$evalue_point, 2.2898751, tolerance = 1e-6)
  expect_equal(res$stats$evalue_ci,    1.6103234, tolerance = 1e-6)
  expect_true(res$stats$evalue_ci < res$stats$evalue_point)
  expect_true(is.finite(res$stats$tip_effect))
  expect_true(all(is.na(res$stats[, c("rv_q", "rv_qa", "xrv_qa", "r2yd_x")])))
  expect_null(res$bounds)
})

test_that("get_sens honours rare = TRUE on the Cox backend", {
  sens_test_deps("survival", "tipr")
  res <- get_sens(cox_data(), treat = "sex", outcome = "status", time = "time",
                  adj_var = c("age", "ph.ecog"), method = "cox",
                  evalue_args = list(rare = TRUE))
  expect_equal(res$stats$evalue_point, 2.8700924, tolerance = 1e-6)
  expect_equal(res$stats$evalue_ci,    1.8110848, tolerance = 1e-6)
})

test_that("get_sens returns a bounded confidence region for the DML backend", {
  sens_test_deps("dml.sensemakr")
  res <- get_sens(dml_data(), treat = "e401", outcome = "net_tfa",
                  adj_var = c("age", "inc", "educ", "fsize", "marr",
                              "twoearn", "pira", "hown"),
                  bench_var = "inc", method = "dml",
                  bench_args = list(k_treat = 1:2, bound = c(0.04, 0.03)),
                  dml_args = list(cf_folds = 2L, cf_seed = 20260921L,
                                  dirty_tuning = FALSE))

  expect_named(res$stats, stats_cols)
  expect_identical(res$stats$method, "dml")
  expect_identical(res$stats$estimand, "ate")
  expect_true(is.finite(res$stats$estimate))
  expect_true(res$stats$rv_q >= 0 && res$stats$rv_q < 1)
  expect_true(res$stats$rv_qa <= res$stats$rv_q)
  expect_named(res$bounds, bounds_cols)
  # one manual scenario plus two benchmark multiples
  expect_equal(nrow(res$bounds), 3L)
  expect_true(all(res$bounds$adj_low <= res$bounds$adj_high))
  expect_true(all(is.na(res$stats[, c("evalue_point", "evalue_ci", "r2yd_x")])))
})

test_that("get_sens reports Anderson-Rubin bounds for the IV backend", {
  sens_test_deps("iv.sensemakr")
  res <- get_sens(iv_data(), treat = "educ", outcome = "lwage",
                  instrument = "nearc4",
                  adj_var = c("exper", "expersq", "black", "south", "smsa",
                              "reg661", "reg662", "reg663", "reg664", "reg665",
                              "reg666", "reg667", "reg668", "smsa66"),
                  bench_var = "black", method = "iv",
                  bench_args = list(k_treat = 1:3))

  expect_named(res$stats, stats_cols)
  expect_identical(res$stats$estimand, c("iv", "fs", "rf"))
  expect_equal(res$stats$estimate[1L], 0.1315, tolerance = 1e-3)
  expect_equal(res$stats$rv_qa[1L], 0.00667, tolerance = 1e-3)
  expect_true(all(is.na(res$stats$std.error)))
  expect_true(all(is.na(res$stats$rv_q)))
  expect_equal(nrow(res$bounds), 3L)
  expect_identical(res$bounds$bound_label, c("1x black", "2x black", "3x black"))
})

test_that("get_sens rejects a backend list that does not apply to the method", {
  sens_test_deps("sensemakr")
  expect_error(
    get_sens(lm_data(), treat = "directlyharmed", outcome = "peacefactor",
             adj_var = lm_adj, method = "lm", dml_args = list(cf_folds = 10L)),
    "does not apply to method")
  expect_error(
    get_sens(lm_data(), treat = "directlyharmed", outcome = "peacefactor",
             adj_var = lm_adj, method = "lm", time = "age"),
    "does not apply to method")
})

test_that("get_sens rejects unknown and duplicated named-list fields", {
  sens_test_deps("sensemakr")
  expect_error(
    get_sens(lm_data(), treat = "directlyharmed", outcome = "peacefactor",
             adj_var = lm_adj, method = "lm", bench_args = list(kd = 3)),
    "unknown field")
  expect_error(
    get_sens(lm_data(), treat = "directlyharmed", outcome = "peacefactor",
             adj_var = lm_adj, method = "lm",
             bench_args = list(k_treat = 1, k_treat = 2)),
    "duplicated field")
})

test_that("get_sens validates columns and method-specific requirements", {
  sens_test_deps("sensemakr", "survival", "tipr")
  expect_error(
    get_sens(lm_data(), treat = "nope", outcome = "peacefactor",
             adj_var = lm_adj, method = "lm"),
    "not found in `data`")
  expect_error(
    get_sens(lm_data(), treat = "directlyharmed", outcome = "peacefactor",
             adj_var = "age", bench_var = "female", method = "lm"),
    "must be a subset of `adj_var`")
  expect_error(
    get_sens(cox_data(), treat = "sex", outcome = "status",
             adj_var = "age", method = "cox"),
    "`time` is required")
  expect_error(
    get_sens(cox_data(), treat = "sex", outcome = "status", time = "time",
             adj_var = "age", method = "cox", q = 0.5),
    "`q` must be 1")
})

test_that("get_sens default lists match the internal default constants", {
  expect_identical(eval(formals(get_sens)$bench_args),  causalR:::.SENS_BENCH_DEFAULTS)
  expect_identical(eval(formals(get_sens)$evalue_args), causalR:::.SENS_EVALUE_DEFAULTS)
  expect_identical(eval(formals(get_sens)$dml_args),    causalR:::.SENS_DML_DEFAULTS)
  expect_identical(eval(formals(get_sens)$iv_args),     causalR:::.SENS_IV_DEFAULTS)
})

test_that("print.sens_res echoes a copy-pasteable plt_sens call", {
  sens_test_deps("sensemakr")
  res <- get_sens(lm_data(), treat = "directlyharmed", outcome = "peacefactor",
                  adj_var = lm_adj, bench_var = "female", method = "lm")
  out <- utils::capture.output(print(res))
  expect_true(any(grepl("plt_sens: type", out)))
  expect_true(any(grepl("Not comparable", out)))
  expect_true(any(grepl("sensemakr", out)))
})
