hte_select_data <- function(n = 180L) {
  withr::local_seed(71)
  d <- data.frame(z = rep(0:1, n / 2), pair = rep(seq_len(n / 2), each = 2),
                  x = rnorm(n), group = factor(rep(c("A", "B", "C"), length.out = n)))
  d$ps <- plogis(0.2 * d$x)
  d$y <- (2 * d$z - 1) * (1.8 * d$x + (d$group == "B")) + rnorm(n)
  d$binary <- rbinom(n, 1, plogis((2 * d$z - 1) * 1.5 * d$x))
  event <- rexp(n, exp(-d$z * (0.5 + 1.2 * d$x)))
  censor <- rexp(n, 0.3)
  d$time <- pmin(event, censor)
  d$DSS <- as.integer(event <= censor)
  d
}

hte_select_call <- function(d = hte_select_data(), ...) {
  get_hte_select(d, "z", c("x", "group"), surv = "y", ps_var = "ps",
                 fit_args = list(nfolds = 3L), ...)
}

test_that("validation metrics agree with RATE and effect-scale loss formulas", {
  skip_if_not_installed("grf")
  score <- rep(c(-1, 0, 2, 3), 10)
  evaluation <- list(dr = seq(-3, 4, length.out = 40),
                     y_residual = seq(-2, 2, length.out = 40),
                     w_residual = rep(c(-0.4, 0.6), 20), clusters = NULL)
  value <- causalR:::.hte_select_metrics(score, evaluation, continuous = TRUE)
  for (metric in c("autoc", "qini")) {
    direct <- grf::rank_average_treatment_effect.fit(
      evaluation$dr, score, target = toupper(metric), R = 0)
    expect_equal(value[[metric]], unname(direct$estimate))
  }
  expect_equal(value$r_loss,
               mean((evaluation$y_residual - evaluation$w_residual * 2 * score)^2))
  expect_equal(value$dr_loss, mean((evaluation$dr - 2 * score)^2))
  value <- causalR:::.hte_select_metrics(rep(0, 40), evaluation, continuous = FALSE)
  expect_equal(value$autoc, 0)
  expect_equal(value$qini, 0)
  expect_true(is.na(value$r_loss))
  expect_true(is.na(value$dr_loss))
})

test_that("all outcomes and adjustment routes agree with direct personalized fits", {
  skip_if_not_installed("personalized")
  d <- hte_select_data()
  x <- model.matrix(~x + group, d)[, -1, drop = FALSE]
  for (outcome in c("survival", "continuous", "binary")) {
    y <- switch(outcome, survival = survival::Surv(d$time, d$DSS),
                 continuous = d$y, binary = d$binary)
    selector <- switch(outcome, survival = TRUE, continuous = "y", binary = "binary")
    loss <- switch(outcome, survival = "cox_loss_lasso",
                    continuous = "sq_loss_lasso", binary = "logistic_loss_lasso")
    for (route in c("ps", "match")) {
      matched <- route == "match"
      res <- get_hte_select(d, "z", c("x", "group"), surv = selector,
                             ps_var = if (!matched) "ps" else NULL,
                             match_var = if (matched) "pair" else NULL,
                             fit_args = list(nfolds = 3L), n_select = 1L)
      expect_identical(res$analysis$outcome_type, outcome)
      expect_identical(res$analysis$loss, loss)
      expect_identical(res$selected, res$ranking$variable[1])
      expect_identical(lengths(res$analysis$design_columns), c(x = 1L, group = 2L))
      for (k in seq_len(2L)) {
        variables <- res$forward$variables[[k]]
        expect_identical(variables, res$ranking$variable[seq_len(k)])
        cols <- unlist(res$analysis$design_columns[variables], use.names = FALSE)
        args <- list(x = x[, cols, drop = FALSE], y = y, trt = d$z,
                      loss = loss, method = "weighting", nfolds = 3L,
                      standardize = TRUE)
        if (matched) args$match.id <- factor(d$pair) else {
          args$propensity.func <- function(x, trt) d$ps
          args$foldid <- res$analysis$foldid
        }
        direct <- withr::with_seed(123L, suppressWarnings(
          do.call(personalized::fit.subgroup, args)))
        score <- as.numeric(predict(direct, newx = args$x, type = "benefit.score"))
        expected <- c(sd(score), mean(abs(score)), median(score), mean(score), IQR(score))
        expect_equal(unname(unlist(res$forward[k, c("score_sd", "score_mean_abs",
                                                     "score_median", "score_mean", "score_iqr")])), expected)
        if (k == 1L)
          expect_equal(unname(unlist(res$ranking[1, c("score_sd", "score_mean_abs",
                                                       "score_median", "score_mean", "score_iqr")])), expected)
      }
      expect_identical(res$forward$n, rep(nrow(d), 2))
      expect_s3_class(res$plots$ranking, "ggplot")
      expect_s3_class(res$plots$forward, "ggplot")
      expect_silent(ggplot2::ggplot_build(res$plots$ranking))
      expect_silent(ggplot2::ggplot_build(res$plots$forward))
      expect_identical(names(res), c("ranking", "forward", "selected", "plots", "analysis"))
      expect_identical(res$analysis$fit_args, list(nfolds = 3L, standardize = TRUE))
    }
  }
})

test_that("ranking ties are stable and accumulation never greedily reorders", {
  skip_if_not_installed("personalized")
  d <- hte_select_data()
  d$a <- seq_len(nrow(d))
  d$b <- rev(d$a)
  d$c <- d$a * 2
  seen <- list()
  local_mocked_bindings(.hte_select_fit = function(x, y, trt, ps, match_id,
                                                   foldid, loss, fit_args, seed, context) {
    seen[[length(seen) + 1L]] <<- list(cols = colnames(x), n = nrow(x), ps = ps,
                                     folds = foldid)
    # a and b tie; c would appear best if a greedy algorithm tried it after a.
    value <- if (ncol(x) == 1L) if (colnames(x) == "c") 1 else 2 else 99
    list(statistics = list(score_sd = value, score_mean_abs = 1, score_median = 0,
                           score_mean = 0),
         warnings = "Recorded backend warning")
  })
  res <- get_hte_select(d, "z", c("b", "a", "c"), surv = "y", ps_var = "ps")
  expect_identical(res$ranking$variable, c("b", "a", "c"))
  # The one-variable prefix reuses the top single-variable fit.
  expect_identical(lapply(seen, `[[`, "cols"),
                   list("b", "a", "c", c("b", "a"), c("b", "a", "c")))
  expect_identical(res$selected, NULL)
  expect_identical(res$forward$score_sd[1], res$ranking$score_sd[1])
  expect_identical(nrow(res$analysis$warnings), 6L)
  expect_identical(res$analysis$warnings$stage[4], "forward")
  expect_identical(res$analysis$warnings$step[4], 1L)
  for (call in seen) {
    expect_identical(call$ps, d$ps)
    expect_identical(call$n, nrow(d))
    expect_identical(call$folds, seen[[1]]$folds)
  }
})

test_that("PS-route folds are stratified by arm and sparse outcome class", {
  skip_if_not_installed("personalized")
  d <- hte_select_data()
  # Three positives inside one fold of an unstratified assignment.
  old <- withr::with_seed(123L, sample(rep(1:3, length.out = nrow(d))))
  d$binary <- 0L
  d$binary[which(old == 1L)[1:3]] <- 1L
  res <- get_hte_select(d, "z", c("x", "group"), surv = "binary", ps_var = "ps",
                        fit_args = list(nfolds = 3L))
  folds <- res$analysis$foldid
  expect_identical(as.vector(table(folds[d$binary == 1L])), rep(1L, 3L))
  expect_true(all(abs(table(folds, d$z) - nrow(d) / 6) <= 1))
  d$binary <- 0L
  d$binary[1] <- 1L
  expect_error(get_hte_select(d, "z", "x", surv = "binary", ps_var = "ps",
                              fit_args = list(nfolds = 3L)), "training fold")
  d$DSS <- 0L
  d$DSS[1] <- 1L
  expect_error(get_hte_select(d, "z", "x", ps_var = "ps",
                              fit_args = list(nfolds = 3L)), "training fold")
})

test_that("factor blocks and non-syntactic names retain fixed encoding", {
  skip_if_not_installed("personalized")
  d <- hte_select_data()
  d[["stage / group"]] <- ordered(d$group, levels = c("C", "B", "A"))
  seen <- list()
  local_mocked_bindings(.hte_select_fit = function(x, ...) {
    seen[[length(seen) + 1L]] <<- x
    list(statistics = list(score_sd = 0, score_mean_abs = 0, score_median = 0,
                           score_mean = 0),
         warnings = character())
  })
  res <- get_hte_select(d, "z", c("stage / group", "x"), surv = "y", ps_var = "ps")
  expect_identical(ncol(seen[[1]]), 2L)
  expect_equal(unname(seen[[1]]), unname(1 * cbind(d$group == "B", d$group == "A")))
  expect_equal(seen[[1]], seen[[3]][, 1:2, drop = FALSE])
  expect_identical(res$ranking$variable, c("stage / group", "x"))
  expect_equal(res$ranking$score_sd, c(0, 0))
})

test_that("RNG is restored on success and failure and scores reproduce", {
  skip_if_not_installed("personalized")
  withr::local_seed(99)
  before <- .Random.seed
  first <- hte_select_call()
  expect_identical(.Random.seed, before)
  second <- hte_select_call()
  expect_equal(first$ranking, second$ranking)
  expect_equal(first$forward, second$forward)
  local_mocked_bindings(.hte_select_fit = function(...) {
    runif(1)
    stop("forced failure")
  })
  expect_error(hte_select_call(), "forced failure")
  expect_identical(.Random.seed, before)
  rm(".Random.seed", envir = globalenv())
  expect_error(hte_select_call(), "forced failure")
  expect_identical(exists(".Random.seed", globalenv(), inherits = FALSE), FALSE)
})

test_that("invalid inputs fail before fitting", {
  d <- hte_select_data()
  expect_error(get_hte_select(d, "z", c("x", "x"), ps_var = "ps"), "distinct")
  expect_error(get_hte_select(d, "z", "x"), "exactly one")
  expect_error(get_hte_select(d, "z", "x", match_var = "pair", ps_var = "ps"), "exactly one")
  expect_error(get_hte_select(d, "z", "ps", ps_var = "ps"), "distinct")
  expect_error(get_hte_select(d, "z", "DSS", ps_var = "ps"), "distinct")
  expect_error(get_hte_select(d, "z", "x", surv = FALSE, ps_var = "ps"), "competing risks")
  expect_error(hte_select_call(n_select = 0), "n_select")
  expect_error(hte_select_call(n_select = 1.5), "n_select")
  expect_error(hte_select_call(n_select = 3), "n_select")
  expect_error(get_hte_select(d, "z", "x", ps_var = "ps", fit_args = list(alpha = 0)), "unknown")
  expect_error(get_hte_select(d, "z", "x", ps_var = "ps", fit_args = list(3)), "named")
  expect_error(get_hte_select(d, "z", "x", ps_var = "ps",
                              fit_args = list(nfolds = 3, nfolds = 4)), "duplicated")
  expect_error(get_hte_select(d, "z", "x", ps_var = "ps",
                              fit_args = list(standardize = NULL)), "standardize")
  expect_error(get_hte_select(d, "z", "x", ps_var = "ps", fit_args = list(nfolds = 2)), "nfolds")
  bad <- d; bad$x[1] <- NA
  expect_error(hte_select_call(bad), "missing or non-finite.*x")
  bad <- d; bad$ps[1] <- Inf
  expect_error(hte_select_call(bad), "missing or non-finite.*ps")
  bad <- d; bad$ps[1] <- 0
  expect_error(hte_select_call(bad), "strictly between")
  bad <- d; bad$x <- 1
  expect_error(hte_select_call(bad), "Constant candidate")
  bad <- d; bad$z[1] <- 1
  expect_error(get_hte_select(bad, "z", "x", match_var = "pair"), "Each matched pair")
  bad <- d; bad$time[1] <- 0
  expect_error(get_hte_select(bad, "z", "x", ps_var = "ps"), "positive numeric")
  bad <- d; bad$y <- 1
  expect_error(hte_select_call(bad), "nonconstant")
})

test_that("backend errors identify the affected model and warnings are retained", {
  skip_if_not_installed("personalized")
  withr::local_seed(1)
  local_mocked_bindings(fit.subgroup = function(...) stop("backend problem"),
                        .package = "personalized")
  expect_error(hte_select_call(), "HTE fit failed \\(single 1: x\\): backend problem")
})

test_that("the public signature and documentation expose all defaults", {
  expect_identical(names(formals(get_hte_select)),
                   c("data", "cat_var", "candidate_var", "surv", "match_var", "ps_var",
                     "n_select", "imp_metric", "sel_metric", "fit_args",
                     "eval_args", "seed", "verbose"))
  expect_identical(eval(formals(get_hte_select)$fit_args),
                   list(nfolds = 10L, standardize = TRUE))
  expect_identical(formals(get_hte_select)$imp_metric, "score_sd")
  expect_identical(formals(get_hte_select)$sel_metric, NULL)
  expect_identical(eval(formals(get_hte_select)$eval_args),
                   list(train_frac = 0.5, adjust_var = NULL, target = "RMST",
                        time = NULL, num.trees = 2000L))
  rd <- tools::parse_Rd(test_path("..", "..", "man", "get_hte_select.Rd"))
  text <- paste(capture.output(tools::Rd2txt(rd)), collapse = "\n")
  expect_match(text, "nfolds = 10L, standardize = TRUE", fixed = TRUE)
  expect_match(text, "nfolds Integer", fixed = TRUE)
  expect_match(text, "standardize Logical", fixed = TRUE)
  expect_match(text, "train_frac Numeric", fixed = TRUE)
  expect_match(text, "adjust_var Character", fixed = TRUE)
  expect_match(text, "num.trees Integer", fixed = TRUE)
  usage <- rd[[which(vapply(rd, function(x) identical(attr(x, "Rd_tag"), "\\usage"), logical(1)))]]
  signature <- parse(text = paste(as.character(usage), collapse = ""))[[1L]]
  expect_identical(eval(signature[["eval_args"]]), eval(formals(get_hte_select)$eval_args))
  expect_identical(eval(signature[["imp_metric"]]), "score_sd")
  expect_identical(signature[["sel_metric"]], quote(NULL))
})

test_that("ranking and selection metrics are independent and manual size wins", {
  skip_if_not_installed("personalized")
  local_mocked_bindings(.hte_select_fit = function(x, ...) {
    v <- if (ncol(x) > 2L) 0.5 else if (colnames(x)[1] == "x") 1 else 2
    list(statistics = list(score_sd = 3 - v, score_iqr = v,
                           score_mean_abs = 1, score_mean = 0, score_median = 0),
         warnings = character())
  })
  res <- hte_select_call(imp_metric = "score_iqr", sel_metric = "score_sd")
  expect_identical(res$ranking$variable, c("group", "x"))
  expect_identical(res$analysis$best_step, 2L)
  expect_identical(res$selected, c("group", "x"))
  manual <- hte_select_call(imp_metric = "score_iqr", sel_metric = "score_sd",
                            n_select = 1L)
  expect_identical(manual$selected, "group")
  expect_identical(manual$analysis$best_step, 2L)
  expect_identical(manual$analysis$selected_step, 1L)
})

test_that("validation routes use held-out predictions and match direct calculations", {
  skip_if_not_installed("personalized")
  skip_if_not_installed("grf")
  d <- hte_select_data(360L)
  x <- model.matrix(~x + group, d)[, -1, drop = FALSE]
  for (route in c("ps", "match")) {
    matched <- route == "match"
    res <- get_hte_select(d, "z", c("x", "group"), surv = "y",
      ps_var = if (!matched) "ps", match_var = if (matched) "pair",
      imp_metric = "autoc", sel_metric = "dr_loss", fit_args = list(nfolds = 3L),
      eval_args = list(num.trees = 100L))
    train <- res$analysis$training_rows
    val <- res$analysis$evaluation_rows
    expect_length(intersect(train, val), 0L)
    expect_equal(sort(c(train, val)), seq_len(nrow(d)))
    if (matched) expect_length(intersect(d$pair[train], d$pair[val]), 0L)
    forest <- grf::causal_forest(x[val, , drop = FALSE], d$y[val], d$z[val],
      W.hat = if (matched) rep(0.5, length(val)) else d$ps[val],
      clusters = if (matched) factor(d$pair[val]), num.trees = 100L, seed = 123L)
    dr <- as.numeric(grf::get_scores(forest))
    cols <- unlist(res$analysis$design_columns[res$ranking$variable[1]], use.names = FALSE)
    args <- list(x = x[train, cols, drop = FALSE], y = d$y[train], trt = d$z[train],
                  loss = "sq_loss_lasso", method = "weighting", nfolds = 3L,
                  standardize = TRUE)
    if (matched) args$match.id <- factor(d$pair[train]) else {
      args$propensity.func <- function(x, trt) d$ps[train]
      args$foldid <- res$analysis$foldid
    }
    fit <- withr::with_seed(123L, suppressWarnings(do.call(personalized::fit.subgroup, args)))
    score <- as.numeric(predict(fit, x[val, cols, drop = FALSE]))
    expected <- c(
      autoc = unname(grf::rank_average_treatment_effect.fit(dr, score, R = 0)$estimate),
      qini = unname(grf::rank_average_treatment_effect.fit(dr, score, target = "QINI", R = 0)$estimate),
      r_loss = mean((d$y[val] - forest$Y.hat - (d$z[val] - forest$W.hat) * 2 * score)^2),
      dr_loss = mean((dr - 2 * score)^2))
    expect_equal(unlist(res$ranking[1, names(expected)]), expected, ignore_attr = TRUE)
    expect_equal(unlist(res$forward[1, names(expected)]), expected, ignore_attr = TRUE)
    expect_identical(res$analysis$best_step, which.min(res$forward$dr_loss))
    expect_equal(res$ranking$n_eval, rep(length(val), 2L))
    expect_identical(res$analysis$eval_args$train_frac, 0.5)
    expect_equal(res$ranking$autoc, sort(res$ranking$autoc, decreasing = TRUE))
    expect_identical(res$forward$variables[[2]], res$ranking$variable)
    expect_silent(ggplot2::ggplot_build(res$plots$combined))
  }
})

test_that("binary and both survival targets support RATE on PS and matched data", {
  skip_if_not_installed("personalized")
  skip_if_not_installed("grf")
  d <- hte_select_data(360L)
  for (outcome in c("binary", "RMST", "survival.probability")) {
    for (matched in c(FALSE, TRUE)) {
      res <- get_hte_select(d, "z", "x", surv = if (outcome == "binary") "binary" else TRUE,
        ps_var = if (!matched) "ps", match_var = if (matched) "pair",
        imp_metric = "qini", sel_metric = "autoc", fit_args = list(nfolds = 3L),
        eval_args = list(num.trees = 100L, time = 0.5,
                         target = if (outcome == "binary") "RMST" else outcome))
      expect_identical(res$selected, "x")
      expect_equal(res$ranking$autoc, res$forward$autoc)
      expect_identical(is.finite(res$ranking$qini), TRUE)
      expect_identical(is.na(res$ranking$r_loss), TRUE)
      expect_identical(is.na(res$ranking$dr_loss), TRUE)
      if (outcome != "binary") expect_identical(res$analysis$evaluation$target, outcome)
    }
  }
})

test_that("evaluation nuisances stay fixed, losses sort ascending and ties select first", {
  skip_if_not_installed("personalized")
  skip_if_not_installed("grf")
  d <- hte_select_data()
  d$confounder <- seq_len(nrow(d))
  evaluated <- 0L
  seen <- list()
  local_mocked_bindings(.hte_select_evaluation = function(x, y, trt, ps, ...) {
    evaluated <<- evaluated + 1L
    expect_identical(colnames(x), c("x", "groupB", "groupC", "confounder"))
    list(dr = seq_along(y), marker = ps)
  }, .hte_select_fit = function(x, newx, evaluation, ...) {
    seen[[length(seen) + 1L]] <<- list(x = x, newx = newx, evaluation = evaluation)
    v <- if (colnames(x)[1] == "x" && ncol(x) == 1) 2 else 1
    list(statistics = list(score_sd = 1, score_iqr = 1, score_mean = 0,
      score_median = 0, score_mean_abs = 1, autoc = -v, qini = -v, r_loss = v, dr_loss = v),
      warnings = character())
  })
  res <- hte_select_call(d, imp_metric = "r_loss", sel_metric = "dr_loss",
                         eval_args = list(adjust_var = "confounder"))
  expect_identical(evaluated, 1L)
  expect_identical(res$ranking$variable, c("group", "x"))
  expect_identical(res$analysis$best_step, 1L)
  expect_identical(res$selected, "group")
  expect_identical(res$analysis$rank_direction, "minimize")
  for (item in seen) {
    expect_identical(item$evaluation, seen[[1]]$evaluation)
    expect_equal(nrow(item$x), nrow(seen[[1]]$x))
    expect_equal(nrow(item$newx), nrow(seen[[1]]$newx))
  }
})

test_that("validation restores RNG, reproduces results and reports failures", {
  skip_if_not_installed("personalized")
  skip_if_not_installed("grf")
  withr::local_seed(91)
  before <- .Random.seed
  call <- function() hte_select_call(sel_metric = "r_loss", eval_args = list(num.trees = 100L))
  first <- call()
  expect_identical(.Random.seed, before)
  second <- call()
  expect_equal(first$ranking, second$ranking)
  expect_equal(first$forward, second$forward)
  expect_identical(first$analysis$training_rows, second$analysis$training_rows)
  local_mocked_bindings(.hte_select_evaluation = function(...) {runif(1); stop("forced evaluation")})
  expect_error(call(), "HTE evaluation failed: forced evaluation")
  expect_identical(.Random.seed, before)
})

test_that("new metric and evaluation inputs reject unsupported configurations", {
  d <- hte_select_data()
  expect_error(hte_select_call(imp_metric = "bad"), "imp_metric")
  expect_error(hte_select_call(sel_metric = c("autoc", "qini")), "sel_metric")
  expect_error(hte_select_call(eval_args = list(bad = 1)), "unknown")
  expect_error(hte_select_call(eval_args = list(0.5)), "named")
  expect_error(hte_select_call(eval_args = list(time = 1, time = 2)), "duplicated")
  expect_error(hte_select_call(eval_args = list(train_frac = 1)), "train_frac")
  expect_error(hte_select_call(eval_args = list(target = "HR")), "target")
  expect_error(hte_select_call(eval_args = list(num.trees = 1)), "num.trees")
  expect_error(hte_select_call(eval_args = list(time = Inf)), "time")
  expect_error(hte_select_call(eval_args = list(adjust_var = "y")), "cannot include")
  bad <- d; bad$confounder <- seq_len(nrow(d)); bad$confounder[1] <- NA
  expect_error(hte_select_call(bad, imp_metric = "autoc",
    eval_args = list(adjust_var = "confounder")), "missing or non-finite.*confounder")
  expect_error(get_hte_select(d, "z", "x", ps_var = "ps", imp_metric = "autoc"), "eval_args\\$time")
  for (outcome in list(TRUE, "binary")) {
    for (metric in c("r_loss", "dr_loss"))
      expect_error(get_hte_select(d, "z", "x", surv = outcome, ps_var = "ps",
                                  sel_metric = metric), "continuous outcome")
  }
  expect_error(hte_select_call(imp_metric = "autoc", eval_args = list(train_frac = 0.01)),
                "Each split")
  expect_error(get_hte_select(d, "z", "x", match_var = "pair", imp_metric = "autoc",
    eval_args = list(time = 1, train_frac = 0.99)), "training pairs")
  expect_error(get_hte_select(d, "z", "x", ps_var = "ps", imp_metric = "autoc",
    eval_args = list(time = max(d$time) + 1)), "follow-up")
})

test_that("metric plots keep negative bars, minimize losses and map the top axis", {
  ranking <- data.frame(rank = 1:3, variable = letters[1:3], autoc = c(-0.1, -0.2, -0.4))
  forward <- data.frame(n_vars = 1:3, dr_loss = c(4, 2, 3))
  p <- causalR:::.hte_select_plots(ranking, forward, 2L, "autoc", "dr_loss")
  built <- ggplot2::ggplot_build(p$combined)
  expect_equal(built$data[[1]]$xmin, ranking$autoc)
  expect_equal(built$data[[2]]$y, 2)
  expect_equal(built$layout$panel_scales_x[[1]]$secondary.axis$trans(built$data[[4]]$x),
                forward$dr_loss)
  expect_match(p$combined$labels$caption, "minimum Validation DR-loss")
  expect_equal(ggplot2::ggplot_build(p$ranking)$data[[1]]$y, ranking$autoc)
})

test_that("combined chart aligns ranking rows and cumulative scores on distinct axes", {
  ranking <- data.frame(rank = 1:4, variable = c("a", "b", "c", "d"),
                         score_sd = c(0.8, 0.6, 0.3, 0.1))
  forward <- data.frame(step = 1:4, n_vars = 1:4, score_sd = c(0.8, 1, 0.9, 1),
                         score_mean_abs = c(0.6, 0.9, 0.8, 1.1), score_median = 0,
                         score_mean = c(-0.6, -0.2, -0.4, 0.3))
  plot <- causalR:::.hte_select_plots(ranking, forward, 3L)$combined
  expect_s3_class(plot, "ggplot")
  built <- ggplot2::ggplot_build(plot)
  bars <- built$data[[1]]
  expect_equal(bars$xmax, ranking$score_sd)
  expect_equal(bars$y, 4:1)
  expect_equal(built$data[[2]]$y, 1)
  expect_equal(built$data[[3]]$yintercept, 1)
  path <- built$data[[4]]
  expect_equal(path$y, 4:1)
  secondary <- built$layout$panel_scales_x[[1]]$secondary.axis
  expect_equal(secondary$trans(path$x), forward$score_mean)
  expect_equal(built$data[[6]]$y, 1)
  expect_equal(built$data[[6]]$x, path$x[4])
  expect_equal(built$data[[7]]$yintercept, 2)
  expect_match(secondary$name, "Mean benefit")
})

test_that("combined chart marks the first tied maximum without a selected cutoff", {
  for (n in c(1L, 4L)) {
    ranking <- data.frame(rank = seq_len(n), variable = letters[seq_len(n)], score_sd = 0)
    forward <- data.frame(step = seq_len(n), n_vars = seq_len(n),
                           score_sd = 0, score_mean_abs = 0, score_median = 0, score_mean = 0)
    plot <- causalR:::.hte_select_plots(ranking, forward, NULL)$combined
    expect_s3_class(plot, "ggplot")
    built <- ggplot2::ggplot_build(plot)
    expect_equal(built$data[[1]]$xmax, rep(0, n))
    expect_equal(built$data[[2]]$y, n)
    expect_equal(length(built$data), if (n == 1L) 5L else 6L)
    expect_equal(built$layout$panel_scales_x[[1]]$secondary.axis$trans(
      built$data[[length(built$data) - 1L]]$x), rep(0, n))
  }
})
