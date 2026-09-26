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
        expected <- c(sd(score), mean(abs(score)), median(score), mean(score))
        expect_equal(unname(unlist(res$forward[k, c("score_sd", "score_mean_abs",
                                                     "score_median", "score_mean")])), expected)
        if (k == 1L)
          expect_equal(unname(unlist(res$ranking[1, c("score_sd", "score_mean_abs",
                                                       "score_median", "score_mean")])), expected)
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
  expect_identical(lapply(seen, `[[`, "cols"),
                   list("b", "a", "c", "b", c("b", "a"), c("b", "a", "c")))
  expect_identical(res$selected, NULL)
  expect_identical(nrow(res$analysis$warnings), 6L)
  for (call in seen) {
    expect_identical(call$ps, d$ps)
    expect_identical(call$n, nrow(d))
    expect_identical(call$folds, seen[[1]]$folds)
  }
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
  expect_equal(seen[[1]], seen[[4]][, 1:2, drop = FALSE])
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
                     "n_select", "fit_args", "seed", "verbose"))
  expect_identical(eval(formals(get_hte_select)$fit_args),
                   list(nfolds = 10L, standardize = TRUE))
  rd <- tools::parse_Rd(test_path("..", "..", "man", "get_hte_select.Rd"))
  text <- paste(capture.output(tools::Rd2txt(rd)), collapse = "\n")
  expect_match(text, "nfolds = 10L, standardize = TRUE", fixed = TRUE)
  expect_match(text, "nfolds Integer", fixed = TRUE)
  expect_match(text, "standardize Logical", fixed = TRUE)
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
