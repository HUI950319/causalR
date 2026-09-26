# Fixed-order benefit-score screening; deliberately separate from grf hte_res.

.hte_select_count <- function(x, name, lower, upper) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < lower || x > upper)
    stop(sprintf("`%s` must be an integer between %s and %s.",
                 name, lower, upper), call. = FALSE)
  as.integer(x)
}

.hte_select_fit <- function(x, y, trt, ps, match_id, foldid, loss,
                            fit_args, seed, context) {
  warnings <- character()
  # Each call starts at the same seed. The public function restores the RNG.
  set.seed(seed)
  args <- c(list(x = x, y = y, trt = trt, loss = loss,
                 method = "weighting", match.id = match_id,
                 cutpoint = 0, larger.outcome.better = TRUE), fit_args)
  if (!is.null(ps)) {
    args$propensity.func <- function(x, trt) ps
    args$foldid <- foldid
  }
  stats <- tryCatch(withCallingHandlers({
    fit <- do.call(personalized::fit.subgroup, args)
    scores <- as.numeric(stats::predict(fit, newx = x, type = "benefit.score"))
    if (length(scores) != nrow(x) || any(!is.finite(scores)))
      stop("Backend returned missing, non-finite or incorrectly sized scores.")
    values <- list(score_sd = stats::sd(scores),
                   score_mean_abs = mean(abs(scores)),
                   score_median = stats::median(scores),
                   score_mean = mean(scores))
    if (any(!is.finite(unlist(values))))
      stop("Score summaries are non-finite.")
    values
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    stop(sprintf("HTE fit failed (%s): %s", context, conditionMessage(e)),
         call. = FALSE)
  })
  list(statistics = stats, warnings = warnings)
}

.hte_select_plots <- function(ranking, forward, n_select) {
  bar_data <- ranking
  bar_data$variable <- factor(bar_data$variable,
                              levels = rev(ranking$variable))
  bar <- ggplot2::ggplot(bar_data, ggplot2::aes(x = variable, y = score_sd)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Benefit-score SD") +
    ggplot2::theme_minimal()
  metrics <- c("score_sd", "score_mean_abs", "score_median")
  curve_data <- do.call(rbind, lapply(metrics, function(metric) {
    data.frame(n_vars = forward$n_vars, statistic = metric,
               value = forward[[metric]])
  }))
  curve_data$statistic <- factor(curve_data$statistic, levels = metrics)
  step_breaks <- if (nrow(forward) <= 10L) seq_len(nrow(forward)) else
    pretty(c(1L, nrow(forward)), n = 6L)
  step_breaks <- sort(unique(c(1L, nrow(forward),
                                step_breaks[step_breaks >= 1 & step_breaks <= nrow(forward)])))
  curve <- ggplot2::ggplot(curve_data,
                           ggplot2::aes(x = n_vars, y = value)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~statistic, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = step_breaks) +
    ggplot2::labs(x = "Number of variables (fixed ranking)",
                   y = "Benefit-score statistic") +
    ggplot2::theme_minimal()
  if (nrow(forward) > 1L) curve <- curve + ggplot2::geom_line()
  if (!is.null(n_select))
    curve <- curve + ggplot2::geom_vline(xintercept = n_select,
                                         linetype = "dashed")
  # One row per ranked candidate; the line follows cumulative-model order.
  combined_data <- ranking
  combined_data$position <- nrow(ranking) + 1L - ranking$rank
  combined_data$score_mean <- forward$score_mean
  primary_max <- max(ranking$score_sd)
  if (primary_max == 0) primary_max <- 1
  score_limits <- range(forward$score_mean)
  score_span <- diff(score_limits)
  # Keep the secondary-axis transformation invertible for constant scores.
  padding <- if (score_span > 0) score_span * 0.05 else
    max(abs(score_limits), 1) * 0.05
  score_lower <- score_limits[1] - padding
  score_width <- score_span + 2 * padding
  combined_data$curve_x <- (combined_data$score_mean - score_lower) /
    score_width * primary_max
  combined <- ggplot2::ggplot(combined_data,
                               ggplot2::aes(x = score_sd, y = position)) +
    ggplot2::geom_col(ggplot2::aes(fill = rank), orientation = "y", width = 0.82) +
    ggplot2::scale_fill_gradient(low = "#548A9A", high = "#9AC3AD", guide = "none")
  peak <- which.max(forward$score_mean)
  peak_row <- combined_data[peak, , drop = FALSE]
  combined <- combined +
    ggplot2::geom_col(data = peak_row, fill = "#E97997",
                       orientation = "y", width = 0.82) +
    ggplot2::geom_hline(yintercept = peak_row$position,
                         linetype = "dashed", colour = "#777777", linewidth = 0.4)
  if (nrow(ranking) > 1L)
    combined <- combined + ggplot2::geom_path(
      ggplot2::aes(x = curve_x, group = 1), colour = "#414141", linewidth = 0.55)
  combined <- combined + ggplot2::geom_point(
    ggplot2::aes(x = curve_x), colour = "#414141", size = 1.6)
  combined <- combined + ggplot2::geom_point(
    data = peak_row, ggplot2::aes(x = curve_x), colour = "#D54B77", size = 2.5)
  caption <- sprintf("Red: maximum mean score at %d variables (first maximum if tied)", peak)
  if (!is.null(n_select) && n_select != peak) {
    combined <- combined + ggplot2::geom_hline(
      yintercept = nrow(ranking) + 1L - n_select,
      colour = "#397DAB", linetype = "dotted", linewidth = 0.55)
    caption <- paste0(caption, "\nBlue dotted line: specified n_select = ", n_select)
  }
  combined <- combined +
    ggplot2::scale_y_continuous(breaks = combined_data$position,
                                 labels = combined_data$variable,
                                 expand = ggplot2::expansion(add = 0.7)) +
    ggplot2::scale_x_continuous(
      name = "Single-variable benefit-score SD",
      limits = c(0, primary_max), expand = ggplot2::expansion(mult = c(0, 0.025)),
      sec.axis = ggplot2::sec_axis(
        transform = ~ . / primary_max * score_width + score_lower,
        name = "Mean benefit score (cumulative model)")) +
    ggplot2::labs(y = NULL, caption = caption) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(colour = "#333333"),
                   axis.title.x.top = ggplot2::element_text(margin = ggplot2::margin(b = 9)),
                   axis.title.x.bottom = ggplot2::element_text(margin = ggplot2::margin(t = 9)),
                   plot.caption = ggplot2::element_text(hjust = 0, colour = "#666666"))
  list(ranking = bar, forward = curve, combined = combined)
}

#' Screen HTE variables by fixed-order benefit-score accumulation
#'
#' Fits a separate personalized subgroup model for each candidate, ranks the
#' candidates by the standard deviation of their fitted benefit scores, then
#' fits every prefix of that fixed ranking. Returns descriptive statistics and
#' plots for choosing a model size manually, without retaining fitted models.
#'
#' @param data Data frame containing all named columns. Used columns must have
#'   no missing or non-finite values. No rows are dropped or imputed.
#' @param cat_var Single binary treatment column, following [get_hte()]: 0/1,
#'   logical, or a two-level factor/character column. The second level is the
#'   treated arm (character levels are sorted).
#' @param candidate_var Nonempty character vector of distinct candidate
#'   effect-modifier columns. Numeric, logical, factor and character columns
#'   are supported; constant columns are rejected. Input order breaks ranking
#'   ties. Treatment, outcome, matching and propensity columns are excluded.
#' @param surv Outcome selector, following [get_hte()]. `TRUE` uses `time`
#'   (positive numeric follow-up) and `DSS` (0/1 event indicator). A single
#'   column name selects a numeric continuous or 0/1 binary outcome; logical
#'   outcomes are also treated as binary. `FALSE` (competing risks) is not
#'   supported. Constant outcomes are rejected.
#' @param match_var Single column identifying existing 1:1 matched pairs;
#'   each pair must have exactly one treated and one control observation.
#'   Passed to `personalized` as `match.id`. Supply exactly one of `match_var`
#'   and `ps_var`; the default `NULL` means this route is not used.
#' @param ps_var Single numeric column containing precomputed probabilities
#'   of the treated arm, strictly between 0 and 1. These probabilities remain
#'   fixed across models; they are not weights. Default `NULL`. No propensity
#'   model is fitted and no matching is performed inside this function.
#' @param n_select Optional integer from 1 to the number of candidates.
#'   Returns the first N ranked variables and marks N on the forward and
#'   combined plots (a blue dotted line when N differs from the mean-score
#'   peak). The combined plot's red marker always denotes the mean-score peak.
#'   Default `NULL` makes no selection. The complete path is always computed.
#' @param fit_args Named list of fitting settings. Partial overrides are
#'   supported; unknown, duplicate and unnamed entries are rejected.
#'   \describe{
#'     \item{nfolds}{Integer, default 10, at least 3 and no greater than the
#'       number of rows (PS route) or matched pairs (matching route). Controls
#'       internal cross-validation for the LASSO penalty.}
#'     \item{standardize}{Logical, default `TRUE`. Standardize design columns
#'       inside the backend's penalized fit. Applies to all outcome types.}
#'   }
#' @param seed Nonnegative integer, default 123. All model fits start from
#'   this seed. The caller's random-number state is restored, also on error.
#' @param verbose Logical, default `FALSE`. Report variables and accumulation
#'   steps as they are fitted.
#'
#' @details
#' The backend is [personalized::fit.subgroup()] with `method = "weighting"`.
#' Losses are `cox_loss_lasso` for survival, `sq_loss_lasso` for continuous
#' outcomes and `logistic_loss_lasso` for binary outcomes. Scores come from
#' `predict(..., type = "benefit.score")`, using the backend's `lambda.min`.
#' Positive scores favour the treated arm for longer survival or larger
#' outcomes (binary outcome 1); no treatment recommendation is returned.
#'
#' Categorical variables use treatment contrasts with the first level as
#' reference. Their indicator columns enter together as one original variable.
#' Encoding, rows and propensity scores are fixed across all models. The PS
#' route shares one balanced random fold assignment. The matched route uses
#' the backend's pair-level folds; internal retries can change these folds, so
#' identical final folds across matched models are not guaranteed.
#'
#' Ranking uses only the unweighted sample standard deviation of training
#' benefit scores. Mean absolute score, median and signed mean are also reported.
#' These are descriptive score summaries, not HRs, absolute CATE estimates,
#' formal variable-importance tests or validated predictive performance.
#' Penalty cross-validation does not validate the entire screening procedure.
#' No automatic maximum-SD, elbow or stopping rule is applied. The combined
#' plot highlights the first maximum signed mean score for visual comparison,
#' without setting `selected` or changing the full path. The forward pass never
#' reorders remaining candidates or removes earlier variables.
#'
#' A penalized model producing constant scores is retained with SD zero.
#' Backend warnings are recorded in `analysis$warnings`; errors stop with
#' variable/step context. This object is not an `hte_res` and is not input to
#' the grf-specific HTE plotting or RATE functions.
#'
#' @return A plain named list:
#' \describe{
#'   \item{ranking}{Data frame with `rank`, `variable`, `n`, `score_sd`,
#'     `score_mean_abs`, `score_median` and `score_mean` (signed mean), sorted
#'     by decreasing `score_sd`.}
#'   \item{forward}{Data frame with `step`, `added_variable`, `n_vars`, a
#'     `variables` list column, `n` and the same four score summaries.}
#'   \item{selected}{Character vector of the first `n_select` variables, or
#'     `NULL`. Included variables can still have zero LASSO coefficients.}
#'   \item{plots}{Named list `ranking`, `forward` and `combined` of ggplot
#'     objects, not printed or saved automatically. Forward panels use separate
#'     y scales. The combined plot aligns ranked bars (bottom axis: single-model
#'     score SD) with cumulative-model signed mean scores (top axis). The
#'     pink/red bar and point with a horizontal dashed line mark the first
#'     step attaining the maximum mean; they do not imply a validated optimum.
#'     The secondary axis uses an invertible linear transformation for display only;
#'     the two quantities do not share a numerical scale.}
#'   \item{analysis}{Outcome and treatment mapping, loss, adjustment settings,
#'     sample size, design-column mapping, seed, fit settings, PS fold IDs
#'     (`NULL` for matching), dependency versions and a warning data frame.}
#' }
#' No fitted model, individual benefit scores or input data are retained.
#'
#' @examplesIf requireNamespace("personalized", quietly = TRUE)
#' \donttest{
#' # Simulated survival data: 40 candidates, 8 true effect modifiers.
#' set.seed(20260926)
#' n <- 600L
#' candidates <- c(sprintf("Modifier_%02d", 1:8),
#'                 sprintf("Prognostic_%02d", 1:8), sprintf("Noise_%02d", 1:24))
#' d <- as.data.frame(matrix(rnorm(n * length(candidates)), nrow = n))
#' names(d) <- candidates
#' d$z <- sample(rep(0:1, length.out = n))
#' d$ps <- 0.5
#' benefit <- as.vector(as.matrix(d[candidates[1:8]]) %*%
#'                      c(0.8, 0.6, 0.45, 0.35, 0.25, 0.18, 0.12, 0.08))
#' baseline <- as.vector(as.matrix(d[candidates[9:16]]) %*% seq(0.3, 0.05, length.out = 8))
#' event <- rexp(n, rate = 0.02 * exp(baseline - (2 * d$z - 1) * benefit / 2))
#' censor <- rexp(n, rate = 0.008)
#' d$time <- pmin(event, censor)
#' d$DSS <- as.integer(event <= censor)
#' # The plot marks the mean-score peak, but selected remains NULL.
#' ans <- get_hte_select(d, "z", candidates, ps_var = "ps",
#'                       fit_args = list(nfolds = 3L))
#' ans$ranking
#' ans$forward
#' ans$plots$combined
#' }
#' @export
get_hte_select <- function(data, cat_var, candidate_var, surv = TRUE,
                           match_var = NULL, ps_var = NULL, n_select = NULL,
                           fit_args = list(nfolds = 10L, standardize = TRUE),
                           seed = 123L, verbose = FALSE) {
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data)))
    stop("`data` must be a non-empty data frame with unique column names.",
         call. = FALSE)
  cat_var <- .sens_check_col(cat_var, data, "cat_var", n = 1L)
  candidate_var <- .sens_check_col(candidate_var, data, "candidate_var")
  if (is.null(cat_var) || !length(candidate_var) || anyDuplicated(candidate_var))
    stop("Supply one `cat_var` and nonempty, distinct `candidate_var` names.",
         call. = FALSE)
  if (is.null(match_var) == is.null(ps_var))
    stop("Supply exactly one of `match_var` and `ps_var`.", call. = FALSE)
  match_var <- .sens_check_col(match_var, data, "match_var", n = 1L)
  ps_var <- .sens_check_col(ps_var, data, "ps_var", n = 1L)
  if (isTRUE(surv)) {
    outcome <- .sens_check_col(c("time", "DSS"), data, "surv = TRUE")
  } else {
    if (!is.character(surv) || length(surv) != 1L)
      stop("`surv` must be TRUE or a single outcome column name; competing risks are not supported.",
           call. = FALSE)
    outcome <- .sens_check_col(surv, data, "surv", n = 1L)
  }
  roles <- c(cat_var, outcome, match_var, ps_var)
  if (anyDuplicated(roles) || any(candidate_var %in% roles))
    stop("Treatment, outcome, adjustment and candidate columns must be distinct.",
         call. = FALSE)
  used <- c(roles, candidate_var)
  valid_type <- vapply(data[used], function(x) {
    is.null(dim(x)) && (is.factor(x) || (!is.object(x) &&
      (is.numeric(x) || is.logical(x) || is.character(x))))
  }, logical(1))
  if (any(!valid_type))
    stop(sprintf("Unsupported column type(s): %s.",
                 paste(used[!valid_type], collapse = ", ")), call. = FALSE)
  bad <- vapply(data[used], function(x) {
    anyNA(x) || (is.numeric(x) && any(!is.finite(x)))
  }, logical(1))
  if (any(bad))
    stop(sprintf("Used columns must have no missing or non-finite values: %s.",
                 paste(used[bad], collapse = ", ")), call. = FALSE)
  constant <- vapply(data[candidate_var], function(x) length(unique(x)) < 2L,
                     logical(1))
  if (any(constant))
    stop(sprintf("Constant candidate variable(s): %s.",
                 paste(candidate_var[constant], collapse = ", ")), call. = FALSE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  seed <- .hte_select_count(seed, "seed", 0, .Machine$integer.max)
  if (!is.null(n_select))
    n_select <- .hte_select_count(n_select, "n_select", 1, length(candidate_var))
  fit_args <- .merge_named_arg(fit_args,
                               list(nfolds = 10L, standardize = TRUE), "fit_args")
  if (!is.logical(fit_args$standardize) || length(fit_args$standardize) != 1L ||
      is.na(fit_args$standardize))
    stop("`fit_args$standardize` must be TRUE or FALSE.", call. = FALSE)
  tz <- .psw_treat(data[[cat_var]], cat_var, arg = "cat_var")
  trt <- tz$z
  if (any(tabulate(trt + 1L, nbins = 2L) < 2L))
    stop("Both treatment arms need at least two observations.", call. = FALSE)
  ps <- match_id <- NULL
  if (!is.null(match_var)) {
    match_id <- droplevels(as.factor(data[[match_var]]))
    pair_counts <- table(match_id, factor(trt, levels = 0:1))
    if (any(pair_counts != 1L))
      stop("Each matched pair must contain exactly one treated and one control observation.",
           call. = FALSE)
  } else {
    ps <- data[[ps_var]]
    if (!is.numeric(ps) || any(ps <= 0 | ps >= 1))
      stop("`ps_var` must contain numeric probabilities strictly between 0 and 1.",
           call. = FALSE)
  }
  fit_args$nfolds <- .hte_select_count(fit_args$nfolds, "fit_args$nfolds", 3,
                                      if (is.null(match_id)) nrow(data) else nlevels(match_id))
  if (isTRUE(surv)) {
    if (!is.numeric(data$time) || any(data$time <= 0) ||
        !(is.numeric(data$DSS) || is.logical(data$DSS)) ||
        any(!data$DSS %in% 0:1) || !any(data$DSS == 1))
      stop("Survival requires positive numeric `time` and 0/1 `DSS` with observed events.",
           call. = FALSE)
    type <- "survival"
    y <- NULL
  } else {
    y <- data[[outcome]]
    if (!(is.numeric(y) || is.logical(y)) || length(unique(y)) < 2L)
      stop("The outcome must be numeric or logical and nonconstant.", call. = FALSE)
    y <- as.numeric(y)
    type <- if (all(y %in% 0:1)) "binary" else "continuous"
  }
  if (!requireNamespace("personalized", quietly = TRUE))
    stop("Install the optional package `personalized` to use `get_hte_select()`.",
         call. = FALSE)
  if (isTRUE(surv)) {
    if (!requireNamespace("survival", quietly = TRUE))
      stop("Install the optional package `survival` for survival outcomes.", call. = FALSE)
    y <- survival::Surv(data$time, as.integer(data$DSS))
  }
  loss <- switch(type, survival = "cox_loss_lasso", continuous = "sq_loss_lasso",
                  binary = "logistic_loss_lasso")
  encoded <- data[candidate_var]
  factors <- candidate_var[!vapply(encoded, is.numeric, logical(1))]
  encoded[factors] <- lapply(encoded[factors], function(x) droplevels(as.factor(x)))
  contrasts <- if (length(factors)) lapply(encoded[factors], function(x) {
    stats::contr.treatment(levels(x), base = 1L)
  }) else NULL
  full_x <- stats::model.matrix(~ ., data = encoded, contrasts.arg = contrasts)
  assignment <- attr(full_x, "assign")
  x <- full_x[, assignment != 0L, drop = FALSE]
  assignment <- assignment[assignment != 0L]
  colnames(x) <- make.unique(colnames(x))
  columns <- stats::setNames(lapply(seq_along(candidate_var), function(i) {
    which(assignment == i)
  }), candidate_var)

  genv <- globalenv()
  old_seed <- if (exists(".Random.seed", envir = genv, inherits = FALSE))
    get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = genv)
    else if (exists(".Random.seed", envir = genv, inherits = FALSE))
      rm(".Random.seed", envir = genv)
  }, add = TRUE)
  set.seed(seed)
  foldid <- if (is.null(match_id))
    sample(rep(seq_len(fit_args$nfolds), length.out = nrow(data))) else NULL
  warnings <- list()
  fit_subset <- function(vars, stage, step) {
    context <- sprintf("%s %d: %s", stage, step, paste(vars, collapse = ", "))
    if (verbose) cli::cli_inform("{context}")
    indices <- unlist(columns[vars], use.names = FALSE)
    fit <- .hte_select_fit(x[, indices, drop = FALSE], y, trt, ps, match_id,
                           foldid, loss, fit_args, seed, context)
    if (length(fit$warnings))
      warnings[[length(warnings) + 1L]] <<- data.frame(
        stage = stage, step = step, variables = paste(vars, collapse = ", "),
        message = fit$warnings)
    fit$statistics
  }
  single <- lapply(seq_along(candidate_var), function(i) {
    data.frame(variable = candidate_var[i], n = nrow(data),
               fit_subset(candidate_var[i], "single", i))
  })
  ranking <- do.call(rbind, single)
  ranking <- ranking[order(-ranking$score_sd, seq_len(nrow(ranking))), , drop = FALSE]
  rownames(ranking) <- NULL
  ranking <- data.frame(rank = seq_len(nrow(ranking)), ranking)
  steps <- lapply(seq_len(nrow(ranking)), function(k) {
    vars <- ranking$variable[seq_len(k)]
    data.frame(step = k, added_variable = vars[k], n_vars = k,
               variables = I(list(vars)), n = nrow(data),
               fit_subset(vars, "forward", k))
  })
  forward <- do.call(rbind, steps)
  warning_table <- if (length(warnings)) do.call(rbind, warnings) else
    data.frame(stage = character(), step = integer(), variables = character(),
               message = character())
  versions <- vapply(c("personalized", "glmnet", "survival"), function(pkg) {
    as.character(utils::packageVersion(pkg))
  }, character(1))
  list(ranking = ranking, forward = forward,
       selected = if (!is.null(n_select)) ranking$variable[seq_len(n_select)] else NULL,
       plots = .hte_select_plots(ranking, forward, n_select),
       analysis = list(backend = "personalized", outcome_type = type,
                       outcome = outcome, cat_var = cat_var,
                       treatment_mapping = unique(data.frame(
                         original = as.character(data[[cat_var]]), encoded = trt)),
                       loss = loss, method = "weighting", n = nrow(data),
                       candidate_var = candidate_var, match_var = match_var,
                       ps_var = ps_var, n_select = n_select,
                       design_columns = lapply(columns, function(idx) colnames(x)[idx]),
                       seed = seed, fit_args = fit_args, foldid = foldid,
                       score_source = "training", versions = versions,
                       warnings = warning_table))
}
