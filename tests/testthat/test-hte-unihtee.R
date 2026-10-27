# get_hte_unihtee(): univariate projections of AIPW pseudo-outcomes. The glm
# route needs only stats + sandwich; forest tests skip without grf.

uni_args <- list(num.trees = 300, seed = 1)

uni_layer <- function(q, geom, col) {
  for (l in q$layers)
    if (inherits(l$geom, geom) && is.data.frame(l$data) && col %in% names(l$data))
      return(l$data)
  NULL
}

# tau = w3 + 1{sex = M}, sd(w3) = 2: the TEM-VIP is 2 per SD of w3, 1 for
# M vs F and 0 for w1, w2 and w4.
uni_cont_data <- function(n = 800L, seed = 510) {
  set.seed(seed)
  d <- data.frame(w1 = stats::rnorm(n), w2 = stats::rnorm(n),
                  w3 = stats::rnorm(n, 0, 2), w4 = stats::rnorm(n),
                  sex = factor(sample(c("F", "M"), n, replace = TRUE)))
  d$a <- stats::rbinom(n, 1, stats::plogis(0.5 * d$w1))
  d$y <- d$w1 + d$w2 + d$a * (d$w3 + (d$sex == "M")) + stats::rnorm(n)
  d
}

uni_bin_data <- function(n = 1500L, seed = 511) {
  set.seed(seed)
  d <- data.frame(w1 = stats::rnorm(n), w2 = stats::rnorm(n),
                  sex = factor(sample(c("F", "M"), n, replace = TRUE)))
  d$a <- stats::rbinom(n, 1, stats::plogis(0.4 * d$w1))
  d$y <- stats::rbinom(n, 1, stats::plogis(-1 + 0.5 * d$w1 +
                                             d$a * (0.2 + 0.6 * d$w2)))
  d
}

test_that("unihtee survival support excludes zero-weight patients", {
  set.seed(321)
  n <- 400L
  d <- data.frame(a = rep(0:1, n / 2L), b = rep(c(0, 0, 1, 1), n / 4L),
                  v = rnorm(n), time = runif(n, 1, 15), DSS = rbinom(n, 1, 0.7))
  no <- d$a == 1 & d$b == 1
  d$time[no] <- runif(sum(no), 1, 4)
  j <- which(no)[1L]
  d$time[j] <- 12
  wt <- rep(1, n); wt[j] <- 0
  for (target in c("survival.probability", "RMST")) {
    expect_warning(out <- get_hte_unihtee(d, "a", tree_var = c("b", "v"), time = 5,
      grf_args = list(num.trees = 100, seed = 7, num.threads = 1, target = target,
                      sample.weights = wt, W.hat = rep(0.5, n))),
      "`b`.*no patient.*followed past")
    expect_true(is.na(out$vip$estimate[out$vip$variable == "b"]))
    expect_true(is.finite(out$vip$estimate[out$vip$variable == "v"]))
    expect_true(is.na(out$hte$importance$p_het[out$hte$importance$variable == "b"]))
    q <- plt_hte_unihtee(out, type = "dep", x_var = "b")
    expect_null(uni_layer(q, "GeomPointrange", "estimate"))
  }
  # A zero-weight extreme must not expand continuous common support.
  x <- seq_len(n); wt[1:4] <- 0
  beyond <- seq_len(n) %in% c(1, 2, 101, 102, 201, 202)
  actual <- .uni_project(x, d$v, d$a, wt, NULL, beyond, TRUE)
  expected <- .uni_project(x, d$v, d$a, wt, NULL, beyond & wt > 0, TRUE)
  expect_equal(actual$i, expected$i)
  expect_equal(actual$range, c(102, 201))
  expect_equal(coef(actual$fit), coef(expected$fit))
  expect_equal(actual$V, expected$V)
})

test_that("unihtee and shared HTE paths retain results after covariance failures", {
  set.seed(1)
  n <- 200L
  cl <- rep(1:4, each = n / 4L)
  d <- data.frame(a = rep(0:1, n / 2L), x = rnorm(n), v = rnorm(n))
  d$x[cl == 1] <- d$x[cl == 1] * 1e8
  d$y <- d$a + d$v + rnorm(n)
  notes <- character()
  out <- withCallingHandlers(get_hte_unihtee(d, "a", tree_var = c("x", "v"),
    surv = "y", grf_args = list(num.trees = 100, seed = 7, num.threads = 1,
      clusters = cl, W.hat = rep(0.5, n), Y.hat = rep(0, n))),
    warning = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("`x`.*(covariance|singular)", notes)))
  expect_true(is.finite(out$vip$std.error[out$vip$variable == "v"]))
  expect_true(is.finite(out$hte$importance$p_het[out$hte$importance$variable == "v"]))
  expect_true(all(is.finite(out$hte$stats$estimate)))
  set.seed(42)
  x <- rnorm(n); x[1L] <- 1e10
  bad <- suppressWarnings(.uni_project(x, d$y, d$a, rep(1, n), NULL, NULL, TRUE))
  expect_match(bad$why, "non-finite HC3 covariance")
})

test_that("shared HTE calibration failure leaves point estimates available", {
  d <- uni_cont_data(n = 200L)
  testthat::local_mocked_bindings(vcovCL = function(...) stop("test covariance failure"),
                                  .package = "sandwich")
  expect_warning(out <- get_hte(d, "a", adj_var = "w3", surv = "y",
    grf_args = list(num.trees = 100, seed = 7, num.threads = 1)),
    "Calibration.*unavailable")
  expect_true(all(is.finite(out$stats$estimate)))
  expect_true(all(is.finite(out$calibration$estimate)))
  expect_true(all(is.na(out$calibration$std.error)))
})

test_that("unihtee volcano keeps finite coordinates when Wald p values underflow", {
  set.seed(1027)
  n <- 400L
  d <- data.frame(a = rep(0:1, n / 2L), x = rnorm(n), v = runif(n))
  set.seed(11)
  d$y <- d$a * (1 + 5 * d$x) + rnorm(n)
  out <- get_hte_unihtee(d, "a", tree_var = c("x", "v"), surv = "y", method = "glm")
  expect_true(any(out$vip$p.value == 0))
  original <- out$vip
  for (measure in c("diff", "ratio")) {
    r <- out
    if (measure == "ratio") {
      for (nm in c("estimate", "conf.low", "conf.high")) r$vip[[nm]] <- exp(r$vip[[nm]])
      r$vip$measure <- measure
      attr(r, "analysis")$measure <- measure
    }
    q <- plt_hte_unihtee(r, type = "volcano")
    pts <- uni_layer(q, "GeomPoint", "neglog10p")
    b <- if (measure == "diff") r$vip$estimate else log(r$vip$estimate)
    expected <- -(log(2) + pnorm(-abs(b / r$vip$std.error), log.p = TRUE)) / log(10)
    expect_equal(pts$neglog10p, expected)
    expect_true(all(is.finite(pts$neglog10p)))
    cut <- uni_layer(q, "GeomHline", "yintercept")$yintercept
    expect_equal(cut, min(expected[r$vip$p.adj < 0.05]))
    expect_true(all(is.finite(ggplot2::ggplot_build(q)$layout$panel_params[[1L]]$y.range)))
  }
  expect_identical(out$vip, original)
})

test_that("unihtee rejects empty candidate selections before drawing", {
  d <- uni_cont_data(n = 200L)
  r <- get_hte_unihtee(d, "a", tree_var = "w3", surv = "y", method = "glm")
  for (type in c("bar", "volcano", "dep")) {
    expect_error(plt_hte_unihtee(r, type = type, x_var = character()),
                  "`x_var` must select at least one candidate")
    expect_s3_class(plt_hte_unihtee(r, type = type), "ggplot")
  }
})

test_that("unihtee detects loss of projection support after deleting a cluster", {
  set.seed(21)
  n <- 200L
  cl <- rep(seq_len(5L), each = n / 5L)
  d <- data.frame(a = rep(0:1, n / 2L), x = as.integer(cl == 1L),
                  v = rnorm(n))
  d$y <- d$v + d$a + rnorm(n)
  expect_warning(out <- get_hte_unihtee(
    d, "a", tree_var = c("x", "v"), surv = "y",
    grf_args = list(num.trees = 100, seed = 1, num.threads = 1,
                    clusters = cl, W.hat = rep(0.5, n), Y.hat = rep(0, n))),
    "`x`.*deleting a cluster")
  expect_true(is.na(out$vip$estimate[out$vip$variable == "x"]))
  expect_true(is.finite(out$vip$std.error[out$vip$variable == "v"]))
  phi <- d$y
  reference <- .uni_project(d$v, phi, d$a, rep(1, n), cl, NULL, TRUE)
  for (j in seq_len(3L)) {
    i <- sample.int(n)
    bad <- .uni_project(d$x[i], phi[i], d$a[i], rep(1, n), cl[i], NULL, FALSE)
    expect_match(bad$why, "deleting a cluster")
    good <- .uni_project(d$v[i], phi[i], d$a[i], rep(1, n), cl[i], NULL, TRUE)
    expect_equal(unname(coef(good$fit)), unname(coef(reference$fit)))
    expect_equal(unname(good$V), unname(reference$V), tolerance = 1e-10)
  }
  # Only one cluster carries variation after zero-weight rows are excluded.
  x <- rep(0, n); x[cl <= 2L] <- seq_len(sum(cl <= 2L))
  expect_match(.uni_project(x, phi, d$a, as.numeric(cl != 2L), cl, NULL, TRUE)$why,
                "deleting a cluster")
})

test_that("unihtee projections are stable under changes of origin and units", {
  d <- uni_cont_data(n = 200L)
  testthat::local_mocked_bindings(
    .uni_glm_fit = function(data, ...) list(e = rep(0.5, nrow(data)),
      mu1 = rep(0, nrow(data)), mu0 = rep(0, nrow(data)), n_clip = 0L),
    .package = "causalR")
  run <- function(d) get_hte_unihtee(d, "a", tree_var = "w3",
                                     surv = "y", method = "glm")
  reference <- run(d)
  reference_line <- uni_layer(plt_hte_unihtee(reference, type = "dep"),
                               "GeomLine", "estimate")
  for (change in list(c(1e8, 1), c(0, 1e-8), c(0, 1e8))) {
    shifted <- d
    shifted$w3 <- change[1L] + change[2L] * d$w3
    expect_warning(out <- run(shifted), NA)
    keep <- c("estimate", "std.error", "conf.low", "conf.high", "p.value")
    expect_equal(out$vip[keep], reference$vip[keep], tolerance = 1e-7)
    line <- uni_layer(plt_hte_unihtee(out, type = "dep"), "GeomLine", "estimate")
    expect_equal(line[c("estimate", "conf.low", "conf.high")],
                  reference_line[c("estimate", "conf.low", "conf.high")], tolerance = 1e-7)
    expect_equal((line$x - change[1L]) / change[2L], reference_line$x,
                  tolerance = 1e-7)
  }
  wt <- rep(c(0, 1, 2, 3), 50L)
  cl <- rep(seq_len(10L), each = 20L)
  i <- wt > 0
  fit <- lm(y ~ w3, data = d[i, ], weights = wt[i])
  expected <- sandwich::vcovCL(fit, cluster = cl[i], type = "HC3")
  for (offset in c(0, 1e8)) {
    pr <- .uni_project(d$w3 + offset, d$y, d$a, wt, cl, NULL, TRUE)
    expect_equal(unname(coef(pr$fit)[2L]), unname(coef(fit)[2L]) * pr$sd,
                  tolerance = 1e-7)
    expect_equal(unname(pr$V[2L, 2L]), unname(expected[2L, 2L]) * pr$sd^2,
                  tolerance = 1e-7)
  }
})

test_that("unihtee rejects finite but non-estimable GLM predictions", {
  d <- uni_cont_data(n = 200L)
  d$arm_x <- d$a * d$w1
  for (binary in c(FALSE, TRUE)) {
    if (binary) d$y <- as.integer(d$y > median(d$y))
    expect_error(get_hte_unihtee(d, "a", tree_var = "w3", adj_var = "arm_x",
                                 surv = "y", method = "glm"),
                  "Fold [0-9]+ control outcome:.*non-estimable")
  }
  d <- uni_cont_data(n = 200L)
  d$duplicate <- d$w3
  expect_warning(out <- get_hte_unihtee(d, "a", tree_var = "w3", adj_var = "duplicate",
                                        surv = "y", method = "glm"), NA)
  ref <- get_hte_unihtee(d, "a", tree_var = "w3", surv = "y", method = "glm")
  expect_equal(out$vip, ref$vip)
  expect_equal(out$data$.score_diff, ref$data$.score_diff)
})

test_that("unihtee plots do not retain unused forest contents", {
  make_plot <- function(payload, type, multiple) {
    set.seed(37)
    n <- 200L
    d <- data.frame(a = rep(0:1, n / 2L), x = rnorm(n), v = rnorm(n))
    d$y <- d$a * d$x + rnorm(n)
    res <- causalR::get_hte_unihtee(d, "a", tree_var = c("x", "v"),
                                   surv = "y", method = "glm")
    res$hte <- list(fit = list(unused = rnorm(payload)))
    p <- causalR::plt_hte_unihtee(res, type = type,
                                 x_var = if (multiple) c("x", "v") else "x")
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(p, path)
    list(bytes = file.info(path)$size,
         layers = ggplot2::ggplot_build(p)$data)
  }
  # A fresh caller avoids retaining testthat fixtures through plot environments.
  environment(make_plot) <- globalenv()
  for (type in c("bar", "volcano", "dep")) {
    for (multiple in c(FALSE, TRUE)) {
      small <- make_plot(0L, type, multiple)
      large <- make_plot(250000L, type, multiple)
      expect_lt(large$bytes - small$bytes, 100000)
      expect_equal(large$layers, small$layers)
    }
  }
})

test_that("unihtee reports unsupported GLM survival before requiring outcome columns", {
  d <- uni_cont_data(n = 200L)
  for (has_survival in c(FALSE, TRUE)) {
    if (has_survival) { d$time <- seq_len(nrow(d)); d$DSS <- d$a }
    expect_error(get_hte_unihtee(d, "a", tree_var = "w3", method = "glm"),
                  'method = "glm" does not handle survival outcomes; use method = "grf"')
  }
  d$time <- d$DSS <- NULL
  expect_error(get_hte_unihtee(d, "a", tree_var = "w3"), "time|DSS")
  expect_warning(out <- get_hte_unihtee(d, "a", tree_var = "w3", surv = "y",
                                        method = "glm"), NA)
  expect_true(is.finite(out$vip$estimate))
})

test_that("unihtee protects input columns from internal names", {
  d <- uni_cont_data(n = 200L)
  for (nm in c(".ps", ".mu1", ".mu0", ".score_diff", ".score_ratio",
               ".score_OR", ".cate", ".dr_score")) {
    z <- d
    names(z)[names(z) == "w1"] <- nm
    for (method in c("glm", "grf"))
      expect_error(get_hte_unihtee(z, "a", tree_var = nm, surv = "y",
                                   method = method), "reserved.*Rename")
  }
  reference <- get_hte_unihtee(d, "a", tree_var = "w1", adj_var = "w3",
                               surv = "y", method = "glm")
  for (nm in c(".w", ".y")) {
    z <- d
    names(z)[names(z) == "w1"] <- nm
    out <- get_hte_unihtee(z, "a", tree_var = nm, adj_var = "w3",
                           surv = "y", method = "glm")
    expect_identical(out$data[[nm]], d$w1)
    expect_equal(out$vip$estimate, reference$vip$estimate)
    expect_equal(out$data$.score_diff, reference$data$.score_diff)
  }
})

test_that("unihtee normalizes labelled missing values and validates storage", {
  d <- uni_cont_data(n = 200L)
  f <- function(z, method = "glm") get_hte_unihtee(
    z, "a", tree_var = c("w1", "w3"), surv = "y", method = method,
    grf_args = if (method == "grf") list(num.trees = 100, seed = 1) else list())
  for (nm in c("w1", "a", "y")) {
    z <- d
    z[[nm]][1:2] <- c(999, 850)
    z[[nm]] <- structure(z[[nm]], class = c("haven_labelled_spss", "haven_labelled"),
                          na_values = 999, na_range = c(800, 900))
    correct <- d
    correct[[nm]][1:2] <- NA
    for (method in c("glm", "grf")) {
      got <- f(z, method)
      ref <- f(correct, method)
      expect_equal(got$vip, ref$vip)
      expect_equal(got$data$.score_diff, ref$data$.score_diff)
      expect_equal(nrow(got$data), nrow(ref$data))
      expect_equal(got$data[[nm]], ref$data[[nm]])
    }
  }
  z <- d
  z$w1 <- structure(z$w1, class = "haven_labelled", labels = c(label = 1))
  expect_equal(f(z)$vip, f(d)$vip)
  z <- d; z$w1 <- factor(replace(d$w1 > 0, 1:2, NA), exclude = NULL)
  correct <- z; correct$w1 <- factor(as.character(z$w1))
  expect_equal(f(z)$vip, f(correct)$vip)
  bad <- list(integer64 = structure(as.double(seq_len(nrow(d))), class = "integer64"),
              matrix = I(cbind(d$w1, d$w3)), infinite = replace(d$w1, 1, Inf),
              list = as.list(d$w1), date = as.Date(seq_len(nrow(d)), origin = "2000-01-01"))
  for (x in bad) {
    z <- d; z$w1 <- x
    expect_error(f(z), "Column `w1`")
  }
  z <- d; names(z)[2] <- "w1"
  expect_error(f(z), "unique")
  expect_error(get_hte_unihtee(d, NULL, tree_var = "w1", surv = "y", method = "glm"),
               "cat_var")
})

test_that("unihtee diagnoses unestimable cross-fitting folds", {
  d <- uni_cont_data(n = 200L)
  d$rare <- factor(c("rare", rep("common", nrow(d) - 1L)))
  expect_error(get_hte_unihtee(d, "a", tree_var = "w3", adj_var = "rare",
                               surv = "y", method = "glm"),
               "Fold [0-9]+.*(propensity|outcome).*`rare`.*level")
  d$y <- 0L
  d$y[c(which(d$a == 0L)[1L], which(d$a == 1L)[1L])] <- 1L
  before <- .Random.seed
  expect_error(get_hte_unihtee(d, "a", tree_var = "w3", surv = "y",
                               method = "glm", measure = c("diff", "ratio", "OR")),
               "Fold [0-9]+.*outcome.*(events|converge)")
  expect_identical(.Random.seed, before)
  d <- uni_cont_data(n = 400L)
  d$a <- as.integer(d$w1 > 0)
  expect_error(get_hte_unihtee(d, "a", tree_var = c("w1", "w3"),
                               surv = "y", method = "glm"),
               "Fold [0-9]+ propensity: GLM did not converge")
})

test_that("unihtee warns about relative-scale instability without clipping", {
  d <- uni_cont_data(n = 200L)
  d$y <- as.integer(d$w1 > 0)
  testthat::local_mocked_bindings(
    .uni_glm_fit = function(data, ...) list(e = rep(0.5, nrow(data)),
      mu1 = rep(1e-12, nrow(data)), mu0 = rep(0.5, nrow(data)), n_clip = 0L),
    .package = "causalR")
  notes <- character()
  r <- withCallingHandlers(get_hte_unihtee(d, "a", tree_var = "w1", surv = "y",
                                           method = "glm", measure = c("ratio", "OR")),
    warning = function(w) {notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning")})
  expect_true(any(grepl("ratio.*near.*boundary", notes)))
  expect_true(any(grepl("OR.*near.*boundary", notes)))
  expect_true(any(grepl("exponentiat", notes)))
  expect_identical(r$data$.mu1, rep(1e-12, nrow(d)))
  g1 <- 1e-12 + d$a / 0.5 * (d$y - 1e-12)
  g0 <- 0.5 + (1 - d$a) / 0.5 * (d$y - 0.5)
  expect_equal(r$data$.score_ratio, log(1e-12 / 0.5) +
                 (g1 - 1e-12) / 1e-12 - (g0 - 0.5) / 0.5)
  expect_equal(r$vip$n, rep(nrow(d), 2L))
})

test_that("unihtee reports non-finite scores before projection", {
  d <- uni_cont_data(n = 200L)
  testthat::local_mocked_bindings(
    .uni_glm_fit = function(data, ...) list(e = c(0, 0, rep(0.5, nrow(data) - 2L)),
      mu1 = rep(0.7, nrow(data)), mu0 = rep(0.5, nrow(data)), n_clip = 0L),
    .package = "causalR")
  expect_warning(r <- get_hte_unihtee(d, "a", tree_var = "w1", surv = "y", method = "glm"),
                  "2 non-finite pseudo-outcomes")
  expect_equal(r$vip$n, 198L)
  expect_equal(sum(!is.finite(r$data$.score_diff)), 2L)
})

test_that("unihtee rejects inactive arguments and records the effective seed", {
  d <- uni_cont_data(n = 200L)
  expect_error(get_hte_unihtee(d, "a", tree_var = "w1", surv = "y", method = "glm",
                               grf_args = list(clusters = rep(1:20, 10))),
               "grf_args.*only.*grf")
  expect_error(get_hte_unihtee(d, "a", tree_var = "w1", surv = "y", method = "grf",
                               glm_args = list(folds = 3)), "glm_args.*only.*glm")
  skip_if_not_installed("grf")
  r <- get_hte_unihtee(d, "a", tree_var = "w1", surv = "y", seed = 123,
                       grf_args = list(num.trees = 100, seed = 7))
  expect_equal(attr(r, "analysis")$seed, r$hte$fit$seed)
  expect_equal(attr(r, "analysis")$seed, 7L)
})

test_that("unihtee keeps levels only for binary candidates", {
  d <- uni_cont_data(n = 200L)
  r <- get_hte_unihtee(d, "a", tree_var = c("w3", "sex"), surv = "y", method = "glm")
  expect_null(attr(r, "analysis")$candidate_levels$w3)
  expect_identical(attr(r, "analysis")$candidate_levels$sex, c("F", "M"))
  old_metadata <- r
  attr(old_metadata, "analysis")$candidate_levels$w3 <- as.character(sort(unique(d$w3)))
  for (v in c("w3", "sex")) {
    actual <- ggplot2::ggplot_build(plt_hte_unihtee(r, type = "dep", x_var = v))$data
    reference <- ggplot2::ggplot_build(plt_hte_unihtee(old_metadata, type = "dep", x_var = v))$data
    expect_equal(actual, reference)
  }
})

test_that("glm route returns the documented table and recovers the truth", {
  skip_if_not_installed("sandwich")
  d   <- uni_cont_data()
  res <- get_hte_unihtee(d, cat_var = "a", adj_var = c("w1", "w2", "w3", "w4", "sex"),
                         surv = "y", method = "glm")
  expect_named(res, c("vip", "data", "hte"))
  expect_null(res$hte)  # no forest on the glm route
  expect_s3_class(res, "hte_unihtee")
  out <- sub(" +$", "", capture.output(print(res)))
  expect_lt(length(out), 18L)
  expect_match(out[1], "^<hte_unihtee> Univariate effect modification \u00b7 glm")
  expect_true(any(grepl("^  w3 +\\+1 SD +", out)))
  expect_true(any(grepl("^# Components: vip, data, and hte", out)))
  expect_named(res$vip, c("variable", "type", "contrast", "n", "sd", "measure",
                          "estimate", "std.error", "conf.low", "conf.high",
                          "p.value", "p.adj"))
  expect_setequal(res$vip$variable, c("w1", "w2", "w3", "w4", "sex"))
  expect_equal(res$vip$p.adj, stats::p.adjust(res$vip$p.value, "BH"))
  expect_false(is.unsorted(res$vip$p.value))
  expect_setequal(res$vip$variable[1:2], c("w3", "sex"))

  v <- res$vip[match(c("w3", "sex"), res$vip$variable), ]
  expect_equal(v$type, c("continuous", "binary"))
  expect_equal(v$contrast, c("+1 SD", "M vs F"))
  expect_lt(abs(v$estimate[1] - 2), 0.3)
  expect_true(v$conf.low[2] < 1 && v$conf.high[2] > 1)

  # the estimate is the least-squares slope of the pseudo-outcome
  dd <- res$data
  expect_equal(v$estimate[1],
               unname(stats::coef(stats::lm(.score_diff ~ I(w3 / stats::sd(w3)),
                                            data = dd))[2]))
  expect_equal(v$estimate[2], mean(dd$.score_diff[dd$sex == "M"]) -
                 mean(dd$.score_diff[dd$sex == "F"]))
  g1 <- dd$.mu1 + dd$a / dd$.ps * (dd$y - dd$.mu1)
  g0 <- dd$.mu0 + (1 - dd$a) / (1 - dd$.ps) * (dd$y - dd$.mu0)
  expect_equal(dd$.score_diff, g1 - g0)
})

test_that("relative measures project log-scale pseudo-outcomes", {
  skip_if_not_installed("sandwich")
  d   <- uni_bin_data()
  res <- get_hte_unihtee(d, cat_var = "a", adj_var = c("w1", "w2", "sex"),
                         surv = "y", method = "glm",
                         measure = c("diff", "ratio", "OR"))
  expect_equal(unique(res$vip$measure), c("diff", "ratio", "OR"))
  for (m in c("diff", "ratio", "OR")) {
    r <- res$vip[res$vip$measure == m, ]
    expect_equal(r$p.adj, stats::p.adjust(r$p.value, "BH"))
    expect_equal(r$variable[1], "w2")
  }
  dd <- res$data
  g1 <- dd$.mu1 + dd$a / dd$.ps * (dd$y - dd$.mu1)
  g0 <- dd$.mu0 + (1 - dd$a) / (1 - dd$.ps) * (dd$y - dd$.mu0)
  expect_equal(dd$.score_ratio, log(dd$.mu1 / dd$.mu0) +
                 (g1 - dd$.mu1) / dd$.mu1 - (g0 - dd$.mu0) / dd$.mu0)
  r <- res$vip[res$vip$measure == "ratio" & res$vip$variable == "w2", ]
  b <- unname(stats::coef(stats::lm(.score_ratio ~ I(w2 / stats::sd(w2)), data = dd))[2])
  expect_equal(r$estimate, exp(b))
  expect_equal(c(r$conf.low, r$conf.high),
               exp(b + c(-1, 1) * stats::qnorm(0.975) * r$std.error))
})

test_that("grf route equals grf::best_linear_projection()", {
  skip_if_not_installed("grf")
  skip_if_not_installed("sandwich")
  d   <- uni_cont_data(n = 600L)
  adj <- c("w1", "w2", "w3", "w4", "sex")
  res <- get_hte_unihtee(d, cat_var = "a", adj_var = adj, surv = "y",
                         grf_args = uni_args)
  hte <- get_hte(d, cat_var = "a", adj_var = adj, surv = "y", grf_args = uni_args)
  expect_equal(res$data$.score_diff, hte$data$.dr_score)
  blp <- grf::best_linear_projection(hte$fit, A = cbind(w3 = d$w3 / stats::sd(d$w3)))
  v <- res$vip[res$vip$variable == "w3", ]
  expect_equal(v$estimate, unname(blp["w3", 1]))
  expect_equal(v$std.error, unname(blp["w3", 2]))

  # The get_hte() result the scores come from is returned with them, and its
  # forest is the one a separate get_hte() call fits.
  expect_s3_class(res$hte, "hte_res")
  expect_equal(res$hte$data$.dr_score, res$data$.score_diff)
  expect_equal(res$hte$data$.cate, hte$data$.cate)
  skip_if_not_installed("forestplot")
  skip_if_not_installed("ggplotify")
  expect_s3_class(plt_hte_sub(res$hte, sub_var = "sex"), "ggplot")
})

test_that("survival outcomes run through the forest on both scales", {
  skip_if_not_installed("grf")
  skip_if_not_installed("sandwich")
  set.seed(20260928)
  n <- 600
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                  grp = sample(c(TRUE, FALSE), n, replace = TRUE))
  d$z  <- stats::rbinom(n, 1, 0.5)
  ev   <- stats::rexp(n, 0.02 * exp(0.3 * d$x1 - d$z * (0.3 + 0.5 * d$x2)))
  cens <- pmin(stats::rexp(n, 0.01), 120)
  d$time <- pmin(ev, cens)
  d$DSS  <- as.integer(ev <= cens)
  res <- get_hte_unihtee(d, cat_var = "z", adj_var = c("x1", "x2", "grp"),
                         time = 60, measure = c("diff", "ratio"),
                         grf_args = uni_args)
  expect_equal(nrow(res$vip), 6L)
  expect_true(all(is.finite(res$vip$estimate)))
  expect_equal(res$vip$contrast[res$vip$variable == "grp"], rep("TRUE vs FALSE", 2))
  expect_true(all(c(".score_diff", ".score_ratio") %in% names(res$data)))
  skip_if_not_installed("patchwork")
  p <- plt_hte_unihtee(res, type = "dep", x_var = c("x2", "grp"), measure = "ratio")
  expect_s3_class(p, "patchwork")
  expect_no_warning(for (q in c(p$patches$plots, list(p))) ggplot2::ggplot_build(q))
})

test_that("tree_var narrows the table and joins the adjustment set", {
  skip_if_not_installed("sandwich")
  d   <- uni_cont_data(n = 400L)
  res <- get_hte_unihtee(d, cat_var = "a", tree_var = "w3",
                         adj_var = c("w1", "w2"), surv = "y", method = "glm")
  expect_equal(res$vip$variable, "w3")
  expect_equal(attr(res, "analysis")$covariates, c("w1", "w2", "w3"))
})

test_that("the glm route leaves the global RNG untouched", {
  skip_if_not_installed("sandwich")
  d <- uni_cont_data(n = 300L)
  set.seed(1)
  before <- .Random.seed
  get_hte_unihtee(d, cat_var = "a", adj_var = c("w1", "w3"), surv = "y",
                  method = "glm")
  expect_identical(.Random.seed, before)
})

test_that("invalid input stops with a pointed message", {
  d <- uni_cont_data(n = 200L)
  d$stage <- factor(sample(c("I", "II", "III"), nrow(d), replace = TRUE))
  d$one   <- 1
  d$time  <- stats::rexp(nrow(d)); d$DSS <- 1L
  f <- function(...) get_hte_unihtee(d, cat_var = "a", surv = "y",
                                     method = "glm", ...)
  expect_error(f(adj_var = c("w1", "stage")), "two levels")
  expect_error(f(adj_var = c("w1", "one")), "constant")
  expect_error(f(adj_var = "w1", glm_args = list(nfold = 5)), "unknown")
  expect_error(f(adj_var = "w1", time = 60), "survival")
  expect_error(f(adj_var = "w1", tree_var = "a"), "cat_var")
  expect_error(get_hte_unihtee(d, cat_var = "a", adj_var = "w1", method = "glm"),
               "method = \"grf\"")
})


# ---- plt_hte_unihtee() -------------------------------------------------------

uni_plot_res <- local({
  cache <- NULL
  function() {
    skip_if_not_installed("sandwich")
    if (is.null(cache))
      cache <<- get_hte_unihtee(uni_cont_data(), cat_var = "a",
                                adj_var = c("w1", "w2", "w3", "w4", "sex"),
                                surv = "y", method = "glm")
    cache
  }
})

test_that("bar draws every signed estimate with its interval, coloured by p.adj", {
  res <- uni_plot_res()
  v   <- res$vip
  p   <- plt_hte_unihtee(res)
  expect_s3_class(p, "ggplot")
  expect_identical(names(attr(p, "plot_size")), c("width", "height"))
  b <- uni_layer(p, "GeomCol", "estimate")
  i <- match(v$variable, b$variable)
  expect_equal(b$estimate[i], v$estimate)
  expect_equal(as.character(b$sig[i]),
               ifelse(v$p.adj < 0.05, ifelse(v$estimate > 0, "Positive", "Negative"), "n.s."))
  e <- uni_layer(p, "GeomErrorbar", "conf.low")
  expect_equal(e$conf.high[match(v$variable, e$variable)], v$conf.high)
  expect_identical(tail(names(formals(plt_hte_unihtee)), 1L), "save")
})

test_that("volcano plots -log10 p and labels only BH-significant candidates", {
  res <- uni_plot_res()
  v   <- res$vip
  p   <- plt_hte_unihtee(res, type = "volcano", sig_level = 0.01)
  pts <- uni_layer(p, "GeomPoint", "neglog10p")
  expect_equal(pts$neglog10p[match(v$variable, pts$variable)], -log10(v$p.value))
  expect_setequal(uni_layer(p, "GeomText", "variable")$variable,
                  v$variable[v$p.adj < 0.01])
})

test_that("dep draws the TEM-VIP line, a spline check and the binary contrast", {
  res <- uni_plot_res()
  v   <- res$vip
  q   <- plt_hte_unihtee(res, type = "dep", x_var = "w3")
  expect_false(inherits(q, "patchwork"))
  ln <- uni_layer(q, "GeomLine", "estimate")
  ln <- ln[order(ln$x), ]
  slope <- (ln$estimate[nrow(ln)] - ln$estimate[1]) / (ln$x[nrow(ln)] - ln$x[1])
  expect_equal(slope * v$sd[v$variable == "w3"], v$estimate[v$variable == "w3"])
  expect_false(is.null(uni_layer(q, "GeomLine", "spline")))
  pr <- uni_layer(plt_hte_unihtee(res, type = "dep", x_var = "sex"),
                  "GeomPointrange", "estimate")
  expect_equal(diff(pr$estimate), v$estimate[v$variable == "sex"])
  skip_if_not_installed("patchwork")
  expect_s3_class(plt_hte_unihtee(res, type = "dep", x_var = c("w3", "sex")),
                  "patchwork")
})

test_that("relative measures draw on the log scale", {
  skip_if_not_installed("sandwich")
  skip_if_not_installed("RegR")
  res <- get_hte_unihtee(uni_bin_data(), cat_var = "a", adj_var = c("w1", "w2", "sex"),
                         surv = "y", method = "glm", measure = c("diff", "ratio"))
  p <- plt_hte_unihtee(res, measure = "ratio")
  b <- uni_layer(p, "GeomCol", "estimate")
  r <- res$vip[res$vip$measure == "ratio", ]
  expect_equal(b$estimate[match(r$variable, b$variable)], r$estimate)
  expect_match(ggplot2::get_labs(p)$x, "ratio")
  q <- plt_hte_unihtee(res, type = "dep", x_var = "w2", measure = "ratio")
  expect_match(ggplot2::get_labs(q)$y, "log")
})

test_that("plt_hte_unihtee rejects invalid input", {
  res <- uni_plot_res()
  expect_error(plt_hte_unihtee(list(a = 1)), "get_hte_unihtee")
  expect_error(plt_hte_unihtee(res, measure = "OR"), "measure")
  expect_error(plt_hte_unihtee(res, x_var = "zz"), "x_var")
  expect_error(plt_hte_unihtee(res, sig_level = 2), "sig_level")
  expect_error(plt_hte_unihtee(res, save = "a.pdf"), "`save`")
  expect_error(plt_hte_unihtee(res, type = "dep", dr_args = list(df = 2)), "unknown")
  expect_error(plt_hte_unihtee(res, dr_args = list(spline_df = 3)), "type = \"dep\"")
})

test_that("save writes one PDF and returns the plot unchanged", {
  res <- uni_plot_res()
  expect_s3_class(plt_hte_unihtee(res, save = list()), "ggplot")
  expect_s3_class(plt_hte_unihtee(res, save = NULL), "ggplot")
  skip_if_not_installed("RegR")
  f <- tempfile(fileext = ".pdf")
  p <- plt_hte_unihtee(res, type = "volcano", save = list(filename = f))
  expect_true(file.exists(f))
  expect_s3_class(p, "ggplot")
})

test_that("dr_args$bins adds quantile-bin means of the pseudo-outcome", {
  res <- uni_plot_res()
  d   <- res$data
  q   <- plt_hte_unihtee(res, type = "dep", x_var = "w3", dr_args = list(bins = 4))
  bn  <- uni_layer(q, "GeomPointrange", "bin")
  expect_equal(bn$estimate, unname(c(tapply(
    d$.score_diff, cut(d$w3, stats::quantile(d$w3, 0:4 / 4), include.lowest = TRUE),
    mean))))
  expect_null(uni_layer(plt_hte_unihtee(res, type = "dep", x_var = "w3"),
                        "GeomPointrange", "bin"))
  expect_error(plt_hte_unihtee(res, type = "dep", dr_args = list(bins = 1)), "bins")
})

test_that("axis_arg$share_y = 'none' gives each dep panel its own range", {
  skip_if_not_installed("patchwork")
  res <- uni_plot_res()
  lim <- function(p) lapply(c(p$patches$plots, list(p)),
                            function(q) q$coordinates$limits$y)
  shared <- lim(plt_hte_unihtee(res, type = "dep", x_var = c("w3", "sex")))
  expect_equal(shared[[1]], shared[[2]])
  own <- lim(plt_hte_unihtee(res, type = "dep", x_var = c("w3", "sex"),
                             axis_arg = list(share_y = "none")))
  expect_false(isTRUE(all.equal(own[[1]], own[[2]])))
  expect_error(plt_hte_unihtee(res, axis_arg = list(share_y = "none")), "type = \"dep\"")
  expect_error(plt_hte_unihtee(res, type = "dep", axis_arg = list(share_y = "x")),
               "share_y")
  expect_identical(tail(names(formals(plt_hte_unihtee)), 4L),
                   c("axis_arg", "var_names", "title", "save"))
})

test_that("var_names relabels every type and rejects clashing labels", {
  res <- uni_plot_res()
  b <- uni_layer(plt_hte_unihtee(res, var_names = c(w3 = "Weight", sex = "Sex")),
                 "GeomCol", "estimate")
  expect_true(all(c("Weight (+1 SD)", "Sex (M vs F)", "w1 (+1 SD)") %in%
                    as.character(b$label)))
  b0 <- uni_layer(plt_hte_unihtee(res, var_names = NULL), "GeomCol", "estimate")
  expect_true("w3 (+1 SD)" %in% as.character(b0$label))
  # a name given twice keeps its first label
  b2 <- uni_layer(plt_hte_unihtee(res, var_names = c(w3 = "Weight", w3 = "Mass")),
                  "GeomCol", "estimate")
  expect_true("Weight (+1 SD)" %in% as.character(b2$label))
  txt <- uni_layer(plt_hte_unihtee(res, type = "volcano", var_names = c(w3 = "Weight")),
                   "GeomText", "display")
  expect_true("Weight" %in% txt$display)
  q <- plt_hte_unihtee(res, type = "dep", x_var = "w3", var_names = c(w3 = "Weight"))
  expect_match(uni_layer(q, "GeomLine", "estimate")$panel[1], "^Weight:")
  expect_identical(ggplot2::get_labs(q)$x, "Weight")
  expect_error(plt_hte_unihtee(res, var_names = c(w3 = "A", w4 = "A")),
               "same display label")
  expect_error(plt_hte_unihtee(res, var_names = list(w3 = 1)), "var_names")
})
