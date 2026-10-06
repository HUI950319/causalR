# =============================================================================
# hte-explain.R -- partial dependence, ALE and SHAP of a get_hte() forest
# =============================================================================
#
# Architecture:
#
#   L1  get_hte_pdp()         partial dependence, one covariate at a time
#   L1  get_hte_ale()         accumulated local effects
#   L1  get_hte_shp()         SHAP values: TreeSHAP of an xgboost surrogate of
#                             the forest, or Kernel SHAP of the forest itself
#   L2  .hte_predict_set()    forest CATE with covariates replaced per row
#   L2  .hte_pdp_mean()       forest CATE averaged over patients per block of
#                             settings, streamed in batches; the engine of
#                             get_hte_pdp() and of .hte_pdp() in hte-plt.R
#   L2  .hte_pdp_grid()       grid values or levels of one covariate
#   L2  .hte_explain_size()   patients that fit the time budget
#   L2  .hte_explain_cost()   conservative predict-time model of a forest
#   L2  .hte_explain_rows()   evenly spaced rows, as .hte_pdp() picks them
#   L2  .hte_explain_vars()   `x_var` resolution, as in plt_hte_dep()
#   L2  .hte_explain_check()  checks shared by the three
#   L2  .hte_levels()         levels of a categorical covariate, in code order
#
# All three explain one function: the full forest's CATE on the "diff" scale
# of get_hte(), predicted at new covariate values -- not the out-of-bag CATE
# in x$data$.cate, which is no function of the covariates alone. The causal
# forest is never refitted; only get_hte_shp(method = "surrogate") fits a new
# model, an xgboost imitation of the forest's CATE, and explains that.
# =============================================================================


# ---- L2 helpers ------------------------------------------------------------

# Design cells per predict() call.
.HTE_EXPLAIN_CELLS <- 1e6

# Conservative grf prediction time, the upper envelope of timings taken on
# 2026-10-06 on a 24-core machine, partly under load, over 15 get_hte()
# forests (n 1,000 to 100,000; 10 or 30 covariates; 200 or 2000 trees): a
# predict() call costs `a` seconds whatever its size and `b` more per row.
# Both grow with the forest's size, n x trees.
#' @keywords internal
#' @noRd
.hte_explain_cost <- function(fit) {
  trees <- fit[["_num_trees"]]
  size  <- nrow(fit$X.orig) * trees
  list(a     = 0.05 * size / 1e6,
       b     = 21.5e-9 * trees * (1 + size / 1e7),
       batch = max(1, floor(.HTE_EXPLAIN_CELLS / ncol(fit$X.orig))))
}

# Patients to explain and the estimated seconds. Each patient costs `rows`
# batched predictions, `calls` predict() calls of its own (Kernel SHAP cannot
# batch patients) and `extra` seconds; `fixed` seconds come once. Without
# `max_n` the count is the largest the cost model fits into 90% of
# `time_budget`, and at least 1.
#' @keywords internal
#' @noRd
.hte_explain_size <- function(x, max_n, time_budget, rows, calls = 0,
                              extra = 0, fixed = 0) {
  n_all <- nrow(x$fit$X.orig)
  k     <- .hte_explain_cost(x$fit)
  per   <- rows * (k$b + k$a / k$batch) + calls * k$a + extra
  if (is.null(max_n)) {
    n <- max(1, min(n_all, floor(max(0, 0.9 * time_budget - fixed - k$a) / per)))
    if (n < n_all)
      cli::cli_inform(c("i" = paste(
        "Explaining {n} of {n_all} patients to stay within",
        "{.arg time_budget} = {time_budget} s (cost-model estimate);",
        "set {.arg max_n} to change it.")))
  } else {
    n <- min(n_all, max_n)
  }
  list(n = n, est = fixed + k$a + n * per)
}

#' @keywords internal
#' @noRd
.hte_explain_rows <- function(n_all, n)
  if (n >= n_all) seq_len(n_all) else unique(round(seq(1, n_all, length.out = n)))

# Levels of a categorical covariate in the order get_hte() coded them, or the
# sorted values of a numeric covariate with at most 5 of them.
#' @keywords internal
#' @noRd
.hte_levels <- function(x, v) {
  xv <- x$data[[v]]
  if (is.numeric(xv)) sort(unique(xv)) else levels(droplevels(as.factor(xv)))
}

#' @keywords internal
#' @noRd
.hte_explain_check <- function(x, max_n, time_budget, verbose) {
  if (!inherits(x, "hte_res") || is.null(x$fit))
    stop("`x` must be an `hte_res` object from get_hte().", call. = FALSE)
  if (!requireNamespace("grf", quietly = TRUE))
    stop("Package 'grf' is required to predict from the forest.", call. = FALSE)
  if (!is.null(max_n) && (!is.numeric(max_n) || length(max_n) != 1L ||
                          is.na(max_n) || max_n < 1 ||
                          (is.finite(max_n) && max_n != floor(max_n))))
    stop("`max_n` must be `NULL`, a positive whole number or Inf.", call. = FALSE)
  if (!is.numeric(time_budget) || length(time_budget) != 1L ||
      !is.finite(time_budget) || time_budget <= 0)
    stop("`time_budget` must be a single positive number of seconds.",
         call. = FALSE)
  if (!isTRUE(verbose) && !isFALSE(verbose))
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  invisible(NULL)
}

# Covariates in importance order, selected by `x_var` as in plt_hte_dep().
#' @keywords internal
#' @noRd
.hte_explain_vars <- function(x, x_var) {
  vars_all <- x$importance$variable
  num      <- vapply(x$data[vars_all], .hte_is_num, logical(1L))
  vars <- if (is.null(x_var)) {
    vars_all
  } else if (identical(x_var, "fct")) {
    vars_all[!num]
  } else if (identical(x_var, "num")) {
    vars_all[num]
  } else {
    if (!is.character(x_var) || anyNA(x_var))
      stop("`x_var` must be `NULL`, covariate names, \"fct\" or \"num\".",
           call. = FALSE)
    miss <- setdiff(x_var, vars_all)
    if (length(miss))
      stop(sprintf("`x_var` names no covariate of `x`: %s. Covariates are %s.",
                   paste0("`", miss, "`", collapse = ", "),
                   paste0("`", vars_all, "`", collapse = ", ")), call. = FALSE)
    unique(x_var)
  }
  if (!length(vars))
    stop("No covariate of `x` matches `x_var`.", call. = FALSE)
  vars
}

# The forest's CATE for X.orig[row, ] with covariate var[i] set to value[i]:
# the value itself for a numeric covariate, its position in .hte_levels() for
# a categorical one -- the integer code, or which one-hot column to set.
# var[i] = NA keeps the row as observed. `var` and `value` may also be lists
# of such vectors, each setting one more covariate per row. Rows run in
# batches of at most .HTE_EXPLAIN_CELLS design cells.
#' @keywords internal
#' @noRd
.hte_predict_set <- function(x, row, var, value) {
  if (!is.list(var)) {
    var   <- list(var)
    value <- list(value)
  }
  fit    <- x$fit
  X      <- fit$X.orig
  a      <- attr(x, "analysis")
  src    <- a$covariates[attr(X, "assign")]
  onehot <- identical(a$factor_encoding, "onehot")
  batch  <- .hte_explain_cost(fit)$batch
  out    <- numeric(length(row))
  for (start in seq(1, length(row), by = batch)) {
    i  <- seq(start, min(length(row), start + batch - 1))
    Xk <- X[row[i], , drop = FALSE]
    for (j in seq_along(var)) {
      vi <- var[[j]][i]
      for (v in unique(vi[!is.na(vi)])) {
        s    <- which(vi == v)
        cols <- which(src == v)
        val  <- value[[j]][i][s]
        if (onehot && !is.numeric(x$data[[v]])) {
          Xk[s, cols] <- 0
          Xk[cbind(s, cols[val])] <- 1
        } else {
          Xk[s, cols] <- val
        }
      }
    }
    out[i] <- as.numeric(stats::predict(fit, Xk)$predictions)
  }
  out
}

# Partial dependence: the forest's CATE averaged over the patients `rows` for
# each of B blocks, block b setting covariate var[[j]][b] to value[[j]][b]
# for every j, as .hte_predict_set() takes them. Patients carry the forest's
# observation weights (.hte_weights()), as in get_hte()'s estimates. The
# B x length(rows) jobs are enumerated batch by batch, so memory stays at one
# batch however many patients and blocks there are.
#' @keywords internal
#' @noRd
.hte_pdp_mean <- function(x, rows, var, value) {
  if (!is.list(var)) {
    var   <- list(var)
    value <- list(value)
  }
  m     <- length(rows)
  wr    <- .hte_weights(x$fit)[rows]
  B     <- length(var[[1L]])
  step  <- .hte_explain_cost(x$fit)$batch
  total <- B * as.double(m)
  est   <- numeric(B)
  for (start in seq(1, total, by = step)) {
    idx  <- seq(start, min(total, start + step - 1)) - 1
    k    <- as.integer(idx %/% m) + 1L
    i    <- idx %% m + 1
    pred <- .hte_predict_set(x, rows[i], lapply(var, `[`, k),
                             lapply(value, `[`, k))
    id   <- sort(unique(k))
    est[id] <- est[id] + rowsum(pred * wr[i], k, reorder = TRUE)[, 1L]
  }
  est / sum(wr)
}

# Values a partial dependence sets covariate v to: the quantiles of a
# continuous covariate at grid_n evenly spaced probabilities from 0 to 1,
# tied ones merged, the levels of any other (.hte_levels()). Evenly spaced
# values over the range put 11 of 21 points above the 99th percentile of a
# log-normal covariate, where the forest only extrapolates.
#' @keywords internal
#' @noRd
.hte_pdp_grid <- function(x, v, grid_n) {
  xv <- x$data[[v]]
  if (.hte_is_num(xv))
    unique(stats::quantile(xv, seq(0, 1, length.out = grid_n), na.rm = TRUE,
                           names = FALSE))
  else .hte_levels(x, v)
}


# ---- L1 get_hte_pdp() ------------------------------------------------------

#' Partial dependence of a get_hte() forest
#'
#' Averages the forest's CATE over a set of patients with one covariate at a
#' time set to each grid value or level, the other covariates as observed
#' (Friedman 2001). Nothing is refitted: the forest in `x` predicts at the
#' modified covariates.
#'
#' The explained function is the full forest's CATE on the `"diff"` scale of
#' [get_hte()], the same as the `"pdp"` layer of [plt_hte_dep()]: the S(t) or
#' RMST difference for a survival outcome, the risk difference for a binary
#' one and the mean difference otherwise. A continuous covariate (more than 5
#' distinct values) is set to its quantiles at `grid_n` evenly spaced
#' probabilities from 0 to 1, so the grid follows the data from its minimum to
#' its maximum and a long tail gets few points; a categorical one, or a
#' numeric one with at most 5 values, gets each level. Patients count with the
#' forest's observation weights -- `grf_args$sample.weights`, or equal
#' cluster weights under `equalize.cluster.weights = TRUE` -- as in the
#' estimates of [get_hte()]. When covariates are
#' correlated, the partial dependence also
#' averages over combinations the data never show; [get_hte_ale()] stays
#' within the data.
#'
#' @section Time budget:
#' Every grid point costs one forest prediction per patient, so the run time
#' grows with the number of patients, covariates and grid points, and with the
#' size of the forest (patients times trees). With `max_n = NULL` the patients
#' are cut, evenly spaced through the data, to the number that a conservative
#' cost model -- fitted to timings on a 24-core machine -- puts within
#' `time_budget` seconds, and a message reports the cut. The count depends on
#' the forest and the arguments only, so the same call gives the same result
#' on any machine, while the run time scales with its speed.
#' `attr(res, "analysis")` records the patients used, the estimated and the
#' elapsed seconds.
#'
#' @param x An `hte_res` object from [get_hte()].
#' @param x_var Covariates to explain: `NULL` (default) for every covariate in
#'   the order of `x$importance`, a character vector of covariate names, or
#'   `"fct"` / `"num"` for only the categorical / continuous ones, as in
#'   [plt_hte_dep()].
#' @param grid_n Finite whole number of at least 2, default `21`. Number of
#'   quantile grid points for each continuous covariate; tied quantiles
#'   merge, so a covariate with many ties gets fewer.
#' @param max_n `NULL` (default) to let `time_budget` set the number of
#'   patients averaged over, or a positive whole number or `Inf` to set it
#'   directly. The patients are evenly spaced through the data, so the result
#'   does not depend on the random seed.
#' @param time_budget Positive number of seconds, default `20`, that the cost
#'   model aims for when `max_n = NULL`. Ignored otherwise.
#' @param verbose `TRUE` reports the estimated and the elapsed seconds.
#'   Default `FALSE`.
#'
#' @return A tibble with one row per covariate and grid point: `variable`,
#'   `value` (the grid value of a numeric covariate, `NA` otherwise), `level`
#'   (the level of a categorical covariate, `NA` otherwise) and `estimate`,
#'   the mean CATE. `attr(res, "analysis")` holds `method`, `n` (patients
#'   averaged over), `n_total`, `grid_n`, `time_budget`, `est_sec` and
#'   `elapsed_sec`.
#'
#' @references
#' Friedman JH (2001). Greedy function approximation: a gradient boosting
#' machine. \emph{Annals of Statistics} 29(5):1189-1232.
#'
#' @seealso [get_hte_ale()] for accumulated local effects; [get_hte_shp()]
#'   for SHAP values; [plt_hte_dep()] for the plotted partial dependence.
#'
#' @examplesIf requireNamespace("grf", quietly = TRUE)
#' \donttest{
#' set.seed(20261006)
#' n <- 600
#' d <- data.frame(age = rnorm(n, 60, 10), x2 = rnorm(n),
#'                 stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
#' d$z <- rbinom(n, 1, 0.5)
#' d$y <- d$x2 + d$z * (0.5 + 0.05 * (d$age - 60)) + rnorm(n)
#' res <- get_hte(d, cat_var = "z", adj_var = c("age", "x2", "stage"),
#'                surv = "y", grf_args = list(num.trees = 500, seed = 1))
#' get_hte_pdp(res, x_var = c("age", "stage"))
#' }
#'
#' @export
get_hte_pdp <- function(x,
                        x_var       = NULL,
                        grid_n      = 21,
                        max_n       = NULL,
                        time_budget = 20,
                        verbose     = FALSE) {
  t0 <- proc.time()[[3L]]
  .hte_explain_check(x, max_n, time_budget, verbose)
  if (!is.numeric(grid_n) || length(grid_n) != 1L || !is.finite(grid_n) ||
      grid_n < 2 || grid_n != floor(grid_n))
    stop("`grid_n` must be a finite whole number of at least 2.", call. = FALSE)
  vars <- .hte_explain_vars(x, x_var)

  d    <- x$data
  grid <- lapply(stats::setNames(vars, vars), function(v)
    .hte_pdp_grid(x, v, grid_n))
  size  <- lengths(grid)
  num   <- vapply(d[vars], is.numeric, logical(1L))
  n_all <- nrow(x$fit$X.orig)
  sz    <- .hte_explain_size(x, max_n, time_budget, rows = sum(size))
  rows  <- .hte_explain_rows(n_all, sz$n)
  m     <- length(rows)

  # One block per covariate and grid point, covariate by covariate.
  code <- unlist(lapply(vars, function(v)
    if (num[[v]]) grid[[v]] else seq_along(grid[[v]])), use.names = FALSE)

  out <- tibble::tibble(
    variable = rep(vars, size),
    value    = unlist(lapply(vars, function(v)
      if (num[[v]]) grid[[v]] else rep(NA_real_, size[[v]])), use.names = FALSE),
    level    = unlist(lapply(vars, function(v)
      if (num[[v]]) rep(NA_character_, size[[v]]) else grid[[v]]),
      use.names = FALSE),
    estimate = .hte_pdp_mean(x, rows, rep(vars, size), code))
  elapsed <- proc.time()[[3L]] - t0
  if (verbose)
    cli::cli_inform("Partial dependence over {m} patient{?s}: estimated {round(sz$est, 1)} s, took {round(elapsed, 1)} s.")
  attr(out, "analysis") <- list(method = "pdp", n = m, n_total = n_all,
                                grid_n = grid_n, time_budget = time_budget,
                                est_sec = sz$est, elapsed_sec = elapsed)
  out
}


# ---- L1 get_hte_ale() ------------------------------------------------------

#' Accumulated local effects of a get_hte() forest
#'
#' Accumulates how the forest's CATE changes across small steps of one
#' covariate, among the patients whose value lies in each step, the other
#' covariates as observed (Apley & Zhu 2020). Unlike the partial dependence
#' of [get_hte_pdp()], no patient is moved far from its own covariate values,
#' so correlated covariates do not produce combinations the data never show.
#' Nothing is refitted.
#'
#' The explained function is the full forest's CATE on the `"diff"` scale of
#' [get_hte()]. A continuous covariate (more than 5 distinct values) is cut at
#' `n_bins + 1` quantiles of the patients used; each patient is predicted at
#' both ends of its interval, the differences are averaged within intervals
#' and summed from the lowest end. A categorical covariate, or a numeric one
#' with at most 5 values, steps from level to level, each step averaging the
#' change over the patients of both neighbouring levels. A numeric covariate,
#' a logical one and an ordered factor step in their own order. The order of
#' an unordered factor or a character covariate is arbitrary and the ALE
#' depends on it, so with 3 or more levels they step in the order of Apley &
#' Zhu (2020): one-dimensional scaling places the levels on a line by how far
#' apart the other covariates lie between them -- the Kolmogorov-Smirnov
#' distance of a numeric covariate, half the L1 distance of the level shares
#' of a categorical one, summed over all patients -- so each step joins
#' similar patients. The rows of a categorical covariate follow its steps. A
#' level none of the explained patients has -- likely for a rare level once
#' the patients are cut -- is skipped: the steps join the levels on either
#' side of it, and it gets `ale = NA`, `n = 0`.
#' The curve is centred to a weighted mean of zero, so
#' it shows how the CATE varies with the covariate, not its level. The
#' averages and the centring weight patients as [get_hte_pdp()] does; `n`
#' counts them. A patient missing the covariate is left out of that
#' covariate.
#'
#' @inheritSection get_hte_pdp Time budget
#'
#' @inheritParams get_hte_pdp
#' @param n_bins Finite whole number of at least 1, default `20`. Number of
#'   quantile intervals for each continuous covariate; tied quantiles merge,
#'   so a covariate with many ties gets fewer.
#' @param max_n `NULL` (default) to let `time_budget` set the number of
#'   patients, or a positive whole number or `Inf` to set it directly. The
#'   patients are evenly spaced through the data, so the result does not
#'   depend on the random seed.
#'
#' @return A tibble with one row per covariate and interval end or level:
#'   `variable`, `value` (the interval end of a continuous covariate, or the
#'   value of a numeric one with at most 5 values; `NA` otherwise), `level`
#'   (the level of a categorical covariate, `NA` otherwise), `ale` (the
#'   centred accumulated effect) and `n` (patients in the interval ending at
#'   `value`, `NA` for the lowest end; patients at the level for a categorical
#'   covariate, whose rows follow its steps). `attr(res, "analysis")` holds
#'   `method`, `n`, `n_total`,
#'   `n_bins`, `time_budget`, `est_sec` and `elapsed_sec`.
#'
#' @references
#' Apley DW, Zhu J (2020). Visualizing the effects of predictor variables in
#' black box supervised learning models. \emph{Journal of the Royal
#' Statistical Society Series B} 82(4):1059-1086.
#'
#' @seealso [get_hte_pdp()] for partial dependence; [get_hte_shp()] for SHAP
#'   values.
#'
#' @examplesIf requireNamespace("grf", quietly = TRUE)
#' \donttest{
#' set.seed(20261006)
#' n <- 600
#' d <- data.frame(age = rnorm(n, 60, 10), x2 = rnorm(n),
#'                 stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
#' d$z <- rbinom(n, 1, 0.5)
#' d$y <- d$x2 + d$z * (0.5 + 0.05 * (d$age - 60)) + rnorm(n)
#' res <- get_hte(d, cat_var = "z", adj_var = c("age", "x2", "stage"),
#'                surv = "y", grf_args = list(num.trees = 500, seed = 1))
#' get_hte_ale(res, x_var = "age", n_bins = 10)
#' }
#'
#' @export
get_hte_ale <- function(x,
                        x_var       = NULL,
                        n_bins      = 20,
                        max_n       = NULL,
                        time_budget = 20,
                        verbose     = FALSE) {
  t0 <- proc.time()[[3L]]
  .hte_explain_check(x, max_n, time_budget, verbose)
  if (!is.numeric(n_bins) || length(n_bins) != 1L || !is.finite(n_bins) ||
      n_bins < 1 || n_bins != floor(n_bins))
    stop("`n_bins` must be a finite whole number of at least 1.", call. = FALSE)
  vars <- .hte_explain_vars(x, x_var)

  d     <- x$data
  cont  <- vapply(d[vars], .hte_is_num, logical(1L))
  n_all <- nrow(x$fit$X.orig)
  # Two predictions per patient and covariate, plus the patient as observed
  # when a categorical covariate needs it.
  sz   <- .hte_explain_size(x, max_n, time_budget,
                            rows = 2 * length(vars) + any(!cont))
  rows <- .hte_explain_rows(n_all, sz$n)
  m    <- length(rows)
  d0   <- d[rows, , drop = FALSE]
  wr   <- .hte_weights(x$fit)[rows]

  # The order an unordered factor with 3 or more levels steps through, as in
  # Apley & Zhu (2020) and ALEPlot: one-dimensional scaling places the levels
  # on a line by how far apart the other covariates lie between them -- the
  # Kolmogorov-Smirnov distance of a numeric covariate at 100 quantiles, half
  # the L1 distance of the level shares of a categorical one, summed -- so
  # each step joins similar patients. Computed on all patients; other
  # covariates keep their own order. Returns indices into `lev`.
  covars  <- attr(x, "analysis")$covariates
  path_of <- function(v, lev) {
    K  <- length(lev)
    xv <- d[[v]]
    if (K < 3L || is.numeric(xv) || is.logical(xv) || is.ordered(xv))
      return(seq_len(K))
    g <- factor(as.character(xv), levels = lev)
    D <- matrix(0, K, K)
    for (u in setdiff(covars, v)) {
      xu <- d[[u]]
      ok <- !is.na(g) & !is.na(xu)
      if (!any(ok)) next
      if (is.numeric(xu)) {
        q  <- stats::quantile(xu[ok], seq(0, 1, length.out = 100), names = FALSE)
        Fk <- lapply(lev, function(l) {
          s <- xu[ok][g[ok] == l]
          if (length(s)) stats::ecdf(s)(q)
        })
      } else {
        tab <- table(g[ok], as.character(xu[ok]))
        Fk  <- lapply(seq_len(K), function(k)
          if (sum(tab[k, ])) tab[k, ] / sum(tab[k, ]))
      }
      if (any(vapply(Fk, is.null, logical(1L)))) next
      for (a in seq_len(K - 1L)) for (b in (a + 1L):K) {
        e <- abs(Fk[[a]] - Fk[[b]])
        D[a, b] <- D[b, a] <- D[a, b] + if (is.numeric(xu)) max(e) else sum(e) / 2
      }
    }
    if (all(D == 0)) return(seq_len(K))
    ord <- order(stats::cmdscale(D, k = 1L)[, 1L])
    # the sign of the scaling is arbitrary; start from the earlier end level
    if (ord[K] < ord[1L]) rev(ord) else ord
  }

  spec <- lapply(vars, function(v) {
    xv <- d0[[v]]
    if (cont[[v]]) {
      z <- unique(stats::quantile(xv, seq(0, 1, length.out = n_bins + 1),
                                  type = 1, na.rm = TRUE, names = FALSE))
      i <- if (length(z) > 1L) which(!is.na(xv)) else integer(0)
      j <- pmax(1L, findInterval(xv[i], z, left.open = TRUE))
      list(z = z, i = i, j = j, row = c(rows[i], rows[i]),
           value = c(z[j], z[j + 1L]))
    } else {
      # pos: each patient's level along the path; k: along the levels kept,
      # those some explained patient of positive weight has. The steps join
      # the kept levels, so a level missing from the explained patients
      # neither breaks the chain nor gets an estimate.
      lev  <- .hte_levels(x, v)
      path <- path_of(v, lev)
      pos  <- match(match(if (is.numeric(xv)) xv else as.character(xv), lev),
                    path)
      wl   <- vapply(seq_along(lev), function(q) sum(wr[which(pos == q)]),
                     numeric(1L))
      keep <- which(wl > 0)
      k    <- match(pos, keep)
      up   <- which(!is.na(k) & k < length(keep))
      dn   <- which(!is.na(k) & k > 1L)
      code <- if (is.numeric(xv)) function(q) lev[path[keep[q]]] else
        function(q) path[keep[q]]
      list(lev = lev[path], pos = pos, keep = keep, k = k, up = up, dn = dn,
           row = c(rows[up], rows[dn]),
           value = c(code(k[up] + 1L), code(k[dn] - 1L)))
    }
  })
  n_job <- vapply(spec, function(s) length(s$row), integer(1L))
  base  <- any(!cont)
  pred  <- .hte_predict_set(
    x,
    c(if (base) rows, unlist(lapply(spec, `[[`, "row"))),
    c(if (base) rep(NA_character_, m), rep(vars, n_job)),
    c(if (base) rep(NA_real_, m), unlist(lapply(spec, `[[`, "value"))))
  f0   <- if (base) pred[seq_len(m)]
  pred <- split(pred[seq_along(pred) > base * m],
                factor(rep(vars, n_job), levels = vars))

  out <- lapply(seq_along(vars), function(q) {
    v <- vars[q]; s <- spec[[q]]; p <- pred[[v]]
    if (cont[[v]]) {
      K  <- length(s$z) - 1L
      if (K < 1L)
        return(tibble::tibble(variable = v, value = s$z, level = NA_character_,
                              ale = 0, n = NA_integer_))
      h  <- length(s$i)
      wi <- wr[s$i]
      dk <- wi * (p[h + seq_len(h)] - p[seq_len(h)])
      nk <- tabulate(s$j, K)
      wk <- vapply(seq_len(K), function(b) sum(wi[s$j == b]), numeric(1L))
      st <- vapply(seq_len(K), function(b)
        if (wk[b] > 0) sum(dk[s$j == b]) / wk[b] else 0, numeric(1L))
      A  <- c(0, cumsum(st))
      A  <- A - sum(wk * (A[-1L] + A[-(K + 1L)]) / 2) / sum(wk)
      tibble::tibble(variable = v, value = s$z, level = NA_character_,
                     ale = A, n = c(NA_integer_, nk))
    } else {
      K   <- length(s$keep)
      nu  <- length(s$up)
      fu  <- p[seq_len(nu)]
      fd  <- p[nu + seq_along(s$dn)]
      ale <- rep(NA_real_, length(s$lev))
      if (K) {
        wk <- vapply(seq_len(K), function(q2) sum(wr[which(s$k == q2)]),
                     numeric(1L))
        gu <- wr[s$up] * (fu - f0[s$up])
        gd <- wr[s$dn] * (f0[s$dn] - fd)
        st <- vapply(seq_len(K - 1L), function(q2)
          (sum(gu[s$k[s$up] == q2]) + sum(gd[s$k[s$dn] == q2 + 1L])) /
            (wk[q2] + wk[q2 + 1L]), numeric(1L))
        A  <- c(0, cumsum(st))
        ale[s$keep] <- A - sum(wk * A) / sum(wk)
      }
      num <- is.numeric(d[[v]])
      tibble::tibble(variable = v,
                     value = if (num) s$lev else NA_real_,
                     level = if (num) NA_character_ else s$lev,
                     ale = ale,
                     n = tabulate(s$pos[!is.na(s$pos)], length(s$lev)))
    }
  })
  out <- tibble::as_tibble(do.call(rbind, out))
  elapsed <- proc.time()[[3L]] - t0
  if (verbose)
    cli::cli_inform("ALE over {m} patient{?s}: estimated {round(sz$est, 1)} s, took {round(elapsed, 1)} s.")
  attr(out, "analysis") <- list(method = "ale", n = m, n_total = n_all,
                                n_bins = n_bins, time_budget = time_budget,
                                est_sec = sz$est, elapsed_sec = elapsed)
  out
}


# ---- L1 get_hte_shp() ------------------------------------------------------

#' SHAP values of a get_hte() forest
#'
#' Splits each patient's CATE, minus a baseline, into additive contributions
#' of the covariates (Lundberg & Lee 2017). grf has no TreeSHAP: a forest's
#' prediction solves a weighted estimating equation rather than averaging its
#' trees, so exact tree algorithms do not apply. Two routes are offered.
#' Neither refits the causal forest or uses the outcomes, but the surrogate
#' route fits a new model and explains that model instead of the forest:
#' \describe{
#'   \item{`"surrogate"` (default)}{A new xgboost model is fitted that
#'     imitates the forest: it regresses the forest's CATE on the covariates,
#'     using 80% of the explained patients, and its exact TreeSHAP values
#'     (Lundberg et al. 2020) explain every explained patient. The SHAP values
#'     are exact for this surrogate and only approximate for the forest, so
#'     `attr(res, "analysis")$r2` reports how much
#'     of the forest's CATE it reproduces on the held-out 20%, and a warning
#'     follows below 0.9; with fewer than 10 patients nothing is held out,
#'     `r2` is `NA` and a warning says so. TreeSHAP conditions along the
#'     tree paths, so with
#'     correlated covariates its credit can differ from Kernel SHAP's.}
#'   \item{`"kernel"`}{Kernel SHAP of the forest itself (Covert & Lee 2021),
#'     through [kernelshap::kernelshap()], against `bg_n` background patients:
#'     exact for up to 8 covariates, a hybrid of exact and sampled coalitions
#'     beyond. A patient needs thousands of forest predictions in several
#'     calls, so the time budget admits few patients from a large forest --
#'     2 from a 10,000-patient, 2000-tree forest with 10 covariates -- and a
#'     warning follows below 10.}
#' }
#' Both work on the original covariates: a factor gets one SHAP value, its
#' one-hot columns summed, and the returned feature values keep its levels.
#' With the forest's observation weights (see [get_hte_pdp()]) the surrogate
#' is fitted, and its R^2 computed, with them, and Kernel SHAP weights its
#' background patients; summaries of the result, such as
#' `shapviz::sv_importance()`, still average the explained patients equally.
#'
#' The explained function is the full forest's CATE on the `"diff"` scale of
#' [get_hte()], so a SHAP value is a contribution to the S(t) or RMST
#' difference, the risk difference or the mean difference. A large value marks
#' a covariate the estimated effect varies with -- an effect modifier of the
#' model, or a proxy for one -- not a cause of the outcome.
#'
#' @inheritSection get_hte_pdp Time budget
#'
#' @inheritParams get_hte_pdp
#' @param method `"surrogate"` (default) or `"kernel"`; see Description.
#' @param max_n `NULL` (default) to let `time_budget` set the number of
#'   patients explained, or a positive whole number or `Inf` to set it
#'   directly. The patients are evenly spaced through the data.
#' @param bg_n Positive whole number, default `50`. Background patients of
#'   `method = "kernel"`, evenly spaced through the data; the run time grows
#'   in proportion. Only used by `"kernel"`.
#' @param surrogate_args Named list for `method = "surrogate"`, passed to
#'   [xgboost::xgb.train()]:
#'   \describe{
#'     \item{`nrounds`}{Whole number of at least 1, default `300`. Boosting
#'       rounds.}
#'     \item{`max_depth`}{Whole number of at least 1, default `6`. Tree
#'       depth.}
#'     \item{`learning_rate`}{Number in (0, 1], default `0.1`.}
#'   }
#'   Only used by `"surrogate"`.
#' @param seed Integer seed, default `123`, for the held-out split of
#'   `"surrogate"` and the sampled coalitions of `"kernel"`. The global random
#'   number state is restored on exit.
#'
#' @return A [shapviz::shapviz()] object: the SHAP matrix `S` (patients by
#'   covariates), the original covariate values `X` and the `baseline`, ready
#'   for the `shapviz::sv_*()` plots and `MLR::plt_shp_*()`.
#'   `attr(res, "analysis")` holds `method`, `n` (patients explained),
#'   `n_total`, `bg_n` (kernel), `r2` (surrogate), `surrogate_args`,
#'   `time_budget`, `est_sec` and `elapsed_sec`.
#'
#' @references
#' Lundberg SM, Lee SI (2017). A unified approach to interpreting model
#' predictions. \emph{Advances in Neural Information Processing Systems} 30.
#'
#' Lundberg SM, Erion G, Chen H, et al. (2020). From local explanations to
#' global understanding with explainable AI for trees. \emph{Nature Machine
#' Intelligence} 2:56-67.
#'
#' Covert I, Lee SI (2021). Improving KernelSHAP: practical Shapley value
#' estimation using linear regression. \emph{Proceedings of AISTATS} 130:
#' 3457-3465.
#'
#' Svensson D, Hermansson E, Nikolaou N, Sechidis K, Lipkovich I (2026).
#' Overview and practical recommendations on using Shapley values for
#' identifying predictive biomarkers via CATE modeling. \emph{Statistics in
#' Medicine} 45:e70375.
#'
#' @seealso [get_hte_pdp()] and [get_hte_ale()] for covariate-wise curves;
#'   [get_hte()] for `$importance`.
#'
#' @examplesIf requireNamespace("grf", quietly = TRUE) && requireNamespace("xgboost", quietly = TRUE) && requireNamespace("shapviz", quietly = TRUE)
#' \donttest{
#' set.seed(20261006)
#' n <- 600
#' d <- data.frame(age = rnorm(n, 60, 10), x2 = rnorm(n),
#'                 stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
#' d$z <- rbinom(n, 1, 0.5)
#' d$y <- d$x2 + d$z * (0.5 + 0.05 * (d$age - 60)) + rnorm(n)
#' res <- get_hte(d, cat_var = "z", adj_var = c("age", "x2", "stage"),
#'                surv = "y", grf_args = list(num.trees = 500, seed = 1))
#' sv <- get_hte_shp(res)
#' attr(sv, "analysis")$r2
#' shapviz::sv_importance(sv)
#' }
#'
#' @export
get_hte_shp <- function(x,
                        method         = c("surrogate", "kernel"),
                        max_n          = NULL,
                        bg_n           = 50,
                        surrogate_args = list(nrounds = 300, max_depth = 6,
                                              learning_rate = 0.1),
                        time_budget    = 20,
                        seed           = 123,
                        verbose        = FALSE) {
  t0 <- proc.time()[[3L]]
  .hte_explain_check(x, max_n, time_budget, verbose)
  method <- match.arg(method)
  if (method == "surrogate" && !missing(bg_n))
    stop("`bg_n` only applies to method = \"kernel\".", call. = FALSE)
  if (method == "kernel" && !missing(surrogate_args))
    stop("`surrogate_args` only applies to method = \"surrogate\".",
         call. = FALSE)
  if (!is.numeric(bg_n) || length(bg_n) != 1L || !is.finite(bg_n) ||
      bg_n < 1 || bg_n != floor(bg_n))
    stop("`bg_n` must be a positive whole number.", call. = FALSE)
  sa <- .merge_named_arg(surrogate_args,
                         list(nrounds = 300, max_depth = 6, learning_rate = 0.1),
                         "surrogate_args")
  for (f in c("nrounds", "max_depth"))
    if (!is.numeric(sa[[f]]) || length(sa[[f]]) != 1L || !is.finite(sa[[f]]) ||
        sa[[f]] < 1 || sa[[f]] != floor(sa[[f]]))
      stop(sprintf("`surrogate_args$%s` must be a whole number of at least 1.", f),
           call. = FALSE)
  if (!is.numeric(sa$learning_rate) || length(sa$learning_rate) != 1L ||
      !is.finite(sa$learning_rate) || sa$learning_rate <= 0 ||
      sa$learning_rate > 1)
    stop("`surrogate_args$learning_rate` must be a number in (0, 1].",
         call. = FALSE)
  seed <- .hte_select_count(seed, "seed", 0, .Machine$integer.max - 1e5)
  need <- c("shapviz", if (method == "surrogate") "xgboost" else "kernelshap")
  for (pkg in need)
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required for method = \"%s\".", pkg, method),
           call. = FALSE)

  a      <- attr(x, "analysis")
  covars <- a$covariates
  X      <- x$fit$X.orig
  d      <- x$data
  n_all  <- nrow(X)
  p      <- length(covars)
  w_all  <- .hte_weights(x$fit)

  genv <- globalenv()
  old_seed <- if (exists(".Random.seed", envir = genv, inherits = FALSE))
    get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = genv)
    else if (exists(".Random.seed", envir = genv, inherits = FALSE))
      rm(".Random.seed", envir = genv)
  }, add = TRUE)
  set.seed(seed)

  r2 <- NULL
  if (method == "surrogate") {
    # Per patient: one forest prediction, the TreeSHAP row and its share of
    # the xgboost fit, timed as the forest's predictions were.
    sz   <- .hte_explain_size(x, max_n, time_budget, rows = 1,
                              extra = 2.5e-4 + 0.8 * 2e-5, fixed = 3)
    rows <- .hte_explain_rows(n_all, sz$n)
    m    <- length(rows)
    X0   <- X[rows, , drop = FALSE]
    w0   <- w_all[rows]
    tau  <- .hte_predict_set(x, rows, rep(NA_character_, m), rep(NA_real_, m))
    ho   <- if (m >= 10L) sample.int(m, floor(0.2 * m)) else integer(0)
    tr   <- setdiff(seq_len(m), ho)
    bst  <- xgboost::xgb.train(
      params = list(objective = "reg:squarederror", max_depth = sa$max_depth,
                    learning_rate = sa$learning_rate, seed = seed),
      data = xgboost::xgb.DMatrix(X0[tr, , drop = FALSE], label = tau[tr],
                                  weight = w0[tr]),
      nrounds = sa$nrounds, verbose = 0)
    wh <- w0[ho]
    ss <- if (length(ho))
      sum(wh * (tau[ho] - sum(wh * tau[ho]) / sum(wh))^2) else 0
    r2 <- if (ss > 0) {
      fit_ho <- stats::predict(bst, X0[ho, , drop = FALSE])
      1 - sum(wh * (fit_ho - tau[ho])^2) / ss
    } else NA_real_
    phi <- stats::predict(bst, X0, predcontrib = TRUE)
    src <- covars[attr(X, "assign")]
    S   <- vapply(covars, function(v)
      rowSums(phi[, which(src == v), drop = FALSE]), numeric(m))
    S   <- matrix(S, nrow = m, dimnames = list(NULL, covars))
    sv  <- shapviz::shapviz(S, X = d[rows, covars, drop = FALSE],
                            baseline = unname(phi[1L, ncol(phi)]))
    if (!length(ho))
      cli::cli_warn(paste(
        "Only {m} patient{?s} {?is/are} explained, too few to hold out 20%",
        "and check the surrogate (R{.sup 2} is NA); raise {.arg max_n} or",
        "{.arg time_budget}."))
    if (!is.na(r2) && r2 < 0.9)
      cli::cli_warn(paste(
        "The surrogate reproduces only R{.sup 2} = {round(r2, 3)} of the",
        "forest's CATE on held-out patients; its SHAP values may misattribute.",
        "Raise {.arg surrogate_args} or use {.code method = \"kernel\"}."))
  } else {
    # kernelshap's own defaults, pinned so the cost model counts the same
    # coalitions: all of them up to 8 covariates; beyond, the exact sizes
    # 1 (and 2 up to 16 covariates) and their complements plus 2p sampled
    # per iteration, budgeted at 4 iterations.
    exact <- p <= 8L
    hd    <- 1L + (p %in% 4:16)
    coal  <- if (exact) 2^p - 2 else 2 * sum(choose(p, seq_len(hd))) + 4 * 2 * p
    bg_n  <- min(bg_n, n_all)
    sz    <- .hte_explain_size(x, max_n, time_budget, rows = bg_n * coal,
                               calls = if (exact) 2 else 5)
    rows  <- .hte_explain_rows(n_all, sz$n)
    m     <- length(rows)
    if (is.null(max_n) && m < 10L)
      cli::cli_warn(paste(
        "Only {m} patient{?s} fit{?s/} the time budget for Kernel SHAP of this",
        "forest; {.code method = \"surrogate\"} explains every patient."))
    src    <- covars[attr(X, "assign")]
    onehot <- identical(a$factor_encoding, "onehot")
    cat_v  <- covars[!vapply(d[covars], is.numeric, logical(1L))]
    lev    <- lapply(stats::setNames(cat_v, cat_v), function(v) .hte_levels(x, v))
    # Original covariates -> the forest's design matrix, coded as get_hte()
    # coded them; a missing value stays missing.
    pred_df <- function(object, newdata) {
      Xn <- matrix(0, nrow(newdata), ncol(X))
      for (v in covars) {
        cols <- which(src == v)
        val  <- newdata[[v]]
        if (v %in% cat_v) {
          val <- match(as.character(val), lev[[v]])
          if (onehot) {
            ok <- which(!is.na(val))
            Xn[cbind(ok, cols[val[ok]])] <- 1
            Xn[is.na(val), cols] <- NA
            next
          }
        }
        Xn[, cols] <- val
      }
      as.numeric(stats::predict(object, Xn)$predictions)
    }
    bg <- .hte_explain_rows(n_all, bg_n)
    ks <- kernelshap::kernelshap(
      x$fit, X = d[rows, covars, drop = FALSE],
      bg_X = d[bg, covars, drop = FALSE],
      pred_fun = pred_df, feature_names = covars,
      bg_w = if (any(w_all != w_all[1L])) w_all[bg], exact = exact,
      hybrid_degree = hd, m = 2L * p, verbose = FALSE)
    sv <- shapviz::shapviz(ks)
  }

  elapsed <- proc.time()[[3L]] - t0
  if (verbose)
    cli::cli_inform("SHAP ({method}) of {m} patient{?s}: estimated {round(sz$est, 1)} s, took {round(elapsed, 1)} s.")
  attr(sv, "analysis") <- list(
    method = method, n = m, n_total = n_all,
    bg_n = if (method == "kernel") bg_n, r2 = r2,
    surrogate_args = if (method == "surrogate") sa,
    time_budget = time_budget, est_sec = sz$est, elapsed_sec = elapsed)
  sv
}
