# Fixed-order benefit-score screening; deliberately separate from grf hte_res.

.hte_select_count <- function(x, name, lower, upper) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < lower || x > upper)
    stop(sprintf("`%s` must be an integer between %s and %s.",
                 name, lower, upper), call. = FALSE)
  as.integer(x)
}

.hte_select_metrics <- function(score, evaluation, continuous) {
  rate <- lapply(c("AUTOC", "QINI"), function(target) {
    unname(grf::rank_average_treatment_effect.fit(
      DR.scores = evaluation$dr, priorities = score, target = target,
      R = 0, clusters = evaluation$clusters)$estimate)
  })
  # The squared-loss weighting model uses treatment coded -1/+1, so its
  # fitted benefit score is half the mean treatment contrast. Cox and
  # logistic link scores do not have this outcome-scale interpretation.
  tau <- 2 * score
  list(autoc = rate[[1]], qini = rate[[2]],
       r_loss = if (continuous)
         mean((evaluation$y_residual - evaluation$w_residual * tau)^2) else NA_real_,
       dr_loss = if (continuous) mean((evaluation$dr - tau)^2) else NA_real_)
}

.hte_select_direction <- function(metric) {
  if (metric %in% c("r_loss", "dr_loss")) 1 else -1
}

.hte_select_evaluation <- function(x, y, trt, ps, match_id, type, eval_args, seed) {
  args <- list(X = x, Y = if (type == "survival") y[, 1] else y,
               W = trt, W.hat = ps, clusters = match_id,
               num.trees = eval_args$num.trees, seed = seed)
  if (type == "survival") {
    beyond <- if (eval_args$target == "RMST") y[, 1] >= eval_args$time else
      y[, 1] > eval_args$time
    if (any(vapply(0:1, function(a) !any(beyond & trt == a), logical(1))))
      stop("Each evaluation arm needs follow-up through the target time (beyond it for survival.probability).",
           call. = FALSE)
    if (!any(y[, 2] == 1 & y[, 1] <= eval_args$time))
      stop("The evaluation sample needs observed events by `eval_args$time`.", call. = FALSE)
    args <- c(args, list(D = y[, 2], horizon = eval_args$time, target = eval_args$target))
  }
  forest <- do.call(if (type == "survival") grf::causal_survival_forest else
                     grf::causal_forest, args)
  result <- list(dr = as.numeric(grf::get_scores(forest)), clusters = match_id)
  if (type == "continuous") {
    result$y_residual <- y - forest$Y.hat
    result$w_residual <- trt - ps
  }
  if (any(!is.finite(unlist(result[c("dr", "y_residual", "w_residual")]))))
    stop("Evaluation returned non-finite DR scores or residuals; check sample size, overlap and censoring support.",
         call. = FALSE)
  result
}

.hte_select_fit <- function(x, y, trt, ps, match_id, foldid, loss,
                            fit_args, seed, context, newx = NULL, evaluation = NULL) {
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
                   score_mean = mean(scores),
                   score_iqr = stats::IQR(scores))
    if (any(!is.finite(unlist(values))))
      stop("Score summaries are non-finite.")
    if (!is.null(evaluation)) {
      prediction <- as.numeric(stats::predict(fit, newx = newx, type = "benefit.score"))
      if (length(prediction) != nrow(newx) || any(!is.finite(prediction)))
        stop("Backend returned invalid evaluation scores.")
      metrics <- .hte_select_metrics(prediction, evaluation, loss == "sq_loss_lasso")
      required <- if (loss == "sq_loss_lasso") names(metrics) else c("autoc", "qini")
      if (any(!is.finite(unlist(metrics[required]))))
        stop("Validation metrics are non-finite.")
      values <- c(values, metrics)
    }
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

.hte_select_plots <- function(ranking, forward, n_select,
                              imp_metric = "score_sd", sel_metric = NULL,
                              best_step = NULL) {
  labels <- c(score_sd = "Benefit-score SD", score_iqr = "Benefit-score IQR",
               score_mean_abs = "Mean absolute benefit score",
               score_median = "Median benefit score", score_mean = "Mean benefit score",
               autoc = "Validation AUTOC", qini = "Validation QINI",
               r_loss = "Validation R-loss", dr_loss = "Validation DR-loss")
  line_metric <- if (is.null(sel_metric)) "score_mean" else sel_metric
  bar_data <- ranking
  bar_data$value <- ranking[[imp_metric]]
  bar_data$variable <- factor(bar_data$variable,
                              levels = rev(ranking$variable))
  bar <- ggplot2::ggplot(bar_data, ggplot2::aes(x = variable, y = value)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = labels[[imp_metric]]) +
    ggplot2::theme_minimal()
  metrics <- intersect(unique(c("score_sd", "score_mean_abs", "score_median", "score_iqr",
                                 sel_metric, "autoc", "qini", "r_loss", "dr_loss")),
                       names(forward))
  metrics <- metrics[vapply(forward[metrics], function(v) all(is.finite(v)), logical(1))]
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
                   y = "Metric value") +
    ggplot2::theme_minimal()
  if (nrow(forward) > 1L) curve <- curve + ggplot2::geom_line()
  if (!is.null(n_select))
    curve <- curve + ggplot2::geom_vline(xintercept = n_select,
                                         linetype = "dashed")
  # One row per ranked candidate; the line follows cumulative-model order.
  combined_data <- ranking
  combined_data$position <- nrow(ranking) + 1L - ranking$rank
  combined_data$value <- ranking[[imp_metric]]
  primary_limits <- range(c(0, combined_data$value))
  if (diff(primary_limits) == 0) primary_limits <- c(0, 1)
  primary_lower <- primary_limits[1]
  primary_width <- diff(primary_limits)
  score_limits <- range(forward[[line_metric]])
  score_span <- diff(score_limits)
  # Keep the secondary-axis transformation invertible for constant scores.
  padding <- if (score_span > 0) score_span * 0.05 else
    max(abs(score_limits), 1) * 0.05
  score_lower <- score_limits[1] - padding
  score_width <- score_span + 2 * padding
  combined_data$curve_x <- (forward[[line_metric]] - score_lower) /
    score_width * primary_width + primary_lower
  combined <- ggplot2::ggplot(combined_data,
                               ggplot2::aes(x = value, y = position)) +
    ggplot2::geom_col(ggplot2::aes(fill = rank), orientation = "y", width = 0.82) +
    ggplot2::scale_fill_gradient(low = "#548A9A", high = "#9AC3AD", guide = "none")
  # A supplied best step of 0 means no step beat the constant-effect baseline.
  peak <- if (is.null(best_step))
    which.min(.hte_select_direction(line_metric) * forward[[line_metric]]) else best_step
  peak_row <- combined_data[peak, , drop = FALSE]
  if (peak > 0L)
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
  optimum <- if (.hte_select_direction(line_metric) == 1) "minimum" else "maximum"
  caption <- if (peak > 0L) {
    combined <- combined + ggplot2::geom_point(
      data = peak_row, ggplot2::aes(x = curve_x), colour = "#D54B77", size = 2.5)
    sprintf("Red: %s %s at %d variables (first optimum if tied)",
            optimum, labels[[line_metric]], peak)
  } else sprintf("No cumulative model beats the constant-effect baseline for %s",
                 labels[[line_metric]])
  if (!is.null(n_select) && n_select != peak) {
    combined <- combined + ggplot2::geom_hline(
      yintercept = nrow(ranking) + 1L - n_select,
      colour = "#397DAB", linetype = "dotted", linewidth = 0.55)
    caption <- paste0(caption, "\nBlue dotted line: selected variable count = ", n_select)
  }
  combined <- combined +
    ggplot2::scale_y_continuous(breaks = combined_data$position,
                                 labels = combined_data$variable,
                                 expand = ggplot2::expansion(add = 0.7)) +
    ggplot2::scale_x_continuous(
      name = paste("Single-variable", labels[[imp_metric]]),
      limits = primary_limits, expand = ggplot2::expansion(mult = c(0, 0.025)),
      sec.axis = ggplot2::sec_axis(
        transform = ~ (. - primary_lower) / primary_width * score_width + score_lower,
        name = paste(labels[[line_metric]], "(cumulative model)"))) +
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
#' candidates using `imp_metric`, then fits every prefix of that fixed ranking
#' (the one-variable prefix reuses the top single-variable fit).
#' Returns score summaries and optional validation metrics, with manual or
#' metric-based model-size selection. No fitted models are retained.
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
#'   Overrides the size chosen by `sel_metric`, while its best step is still
#'   reported. Marks the returned size on the plots. Default `NULL` uses
#'   `sel_metric`, or makes no selection when that is also `NULL`.
#'   The complete path is always computed.
#' @param imp_metric Single metric name for single-variable ranking. Default
#'   `"score_sd"` preserves the original decreasing-SD ranking. Choices are
#'   `"score_sd"`, `"score_iqr"`, `"score_mean_abs"`, `"score_median"`,
#'   `"score_mean"`, `"autoc"`, `"qini"`, `"r_loss"` and `"dr_loss"`.
#'   Losses are minimized; all other metrics are maximized. Ties preserve
#'   candidate input order. The four validation metrics activate `eval_args`.
#' @param sel_metric Single metric name with the same choices and directions
#'   as `imp_metric`, independently applied to the cumulative models. Default
#'   `NULL` preserves manual selection. Otherwise choose the best eligible
#'   step, taking the smallest number of variables on exact ties. A step is
#'   eligible only when its training scores are nonconstant (`score_sd > 0`)
#'   and, where a constant-effect baseline exists, it strictly beats that
#'   baseline: 0 for `"autoc"`, `"qini"`, `"score_sd"` and `"score_iqr"`; the
#'   loss of the best constant effect on the evaluation split for `"r_loss"`
#'   and `"dr_loss"`. Location metrics (`"score_mean_abs"`, `"score_median"`,
#'   `"score_mean"`) have no such baseline. Without an eligible step no
#'   variable is selected. This is a minimal null exit, not a significance
#'   test: a small positive AUTOC can still be noise. The combined plot uses
#'   this metric for its line and red optimum marker; with `NULL` it continues
#'   to display the signed mean-score maximum without automatic selection.
#' @param fit_args Named list of fitting settings. Partial overrides are
#'   supported; unknown, duplicate and unnamed entries are rejected.
#'   \describe{
#'     \item{nfolds}{Integer, default 10, at least 3 and no greater than the
#'       number of rows (PS route) or matched pairs (matching route). Controls
#'       internal cross-validation for the LASSO penalty.}
#'     \item{standardize}{Logical, default `TRUE`. Standardize design columns
#'       inside the backend's penalized fit. Applies to all outcome types.}
#'   }
#' @param eval_args Named list for validation, activated when either metric is
#'   `"autoc"`, `"qini"`, `"r_loss"` or `"dr_loss"`. Partial named overrides
#'   are supported; unknown, duplicate and unnamed entries are rejected.
#'   \describe{
#'     \item{train_frac}{Numeric strictly between 0 and 1, default 0.5. Share
#'       used to train every personalized model. Splits are stratified by arm
#'       on the PS route and by whole pairs on the matching route. Training
#'       must retain enough rows or pairs for `fit_args$nfolds`.}
#'     \item{adjust_var}{Character vector of additional adjustment columns,
#'       default `NULL`. The evaluation forest always uses all candidates
#'       plus these columns, fixed across all models. Include necessary
#'       confounders here; effect-modifier selection must not remove them.
#'       Same column-type, missingness and nonconstant requirements as candidates.}
#'     \item{target}{Survival estimand: `"RMST"` (default), restricted mean
#'       survival-time difference, or `"survival.probability"`, survival
#'       probability difference. Used only for survival validation.}
#'     \item{time}{Positive finite numeric horizon in the same units as
#'       `data$time`. Default `NULL`; required for survival validation.
#'       Both evaluation arms must have follow-up through this horizon
#'       (strictly beyond for survival probability) and the evaluation sample
#'       must have observed events by it. Unused for other outcomes.}
#'     \item{num.trees}{Integer at least 2, default 2000. Number of trees for
#'       the fixed GRF evaluation model. Small values are useful for smoke
#'       tests but may give unstable or non-finite evaluation scores.}
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
#' route shares one random fold assignment, dealt within treatment arm by
#' event status (survival) or outcome class (binary) strata, and stops early
#' unless every training fold keeps both arms and some events (both outcome
#' classes for binary outcomes). The matched route uses
#' the backend's pair-level folds; internal retries can change these folds, so
#' identical final folds across matched models are not guaranteed.
#'
#' Descriptive metrics use unweighted training benefit scores: sample SD,
#' mean absolute score, median, signed mean and interquartile range
#' (`score_iqr`, using R's default type-7 quantiles).
#' These are descriptive score summaries, not HRs, absolute CATE estimates,
#' formal variable-importance tests or validated predictive performance.
#' Penalty cross-validation does not validate the entire screening procedure.
#' The default applies no automatic size rule. Explicit `sel_metric` applies
#' a maximum (or minimum loss) rule; maximizing descriptive score spread or
#' location alone is not evidence of predictive accuracy. The forward pass
#' never reorders remaining candidates or removes earlier variables.
#'
#' With validation enabled, every personalized model is trained on the same
#' training split and predicts the same held-out split. One evaluation forest
#' is fitted on the held-out split: [grf::causal_forest()] for continuous/binary
#' outcomes, or [grf::causal_survival_forest()] for censored survival outcomes.
#' It uses fixed adjustment variables, supplied PS (or 0.5 for 1:1 matching),
#' and pair clusters when present. [grf::get_scores()] supplies common doubly
#' robust (DR) evaluation scores using out-of-bag nuisance predictions. The
#' selected HTE learner remains `personalized`; GRF supplies evaluation only.
#' Neither matching nor a DR score guarantees removal of unmeasured confounding.
#'
#' Both `autoc` and `qini` are RATE point estimates from
#' [grf::rank_average_treatment_effect.fit()] with `R = 0`, without bootstrap
#' standard errors. They evaluate benefit-score rankings, need no conversion
#' to an absolute treatment effect and can be negative. QINI here is GRF's
#' rank-weighted effect, not a normalized uplift coefficient. Tied scores use
#' the backend's tie handling; a constant priority has zero RATE.
#'
#' For continuous outcomes, weighting with squared loss codes treatment as
#' -1/+1; its population-optimal benefit score is half the mean contrast.
#' Loss evaluation therefore uses `tau = 2 * benefit.score` and computes
#' `r_loss = mean((Y - m_hat - (A - e_hat) * tau)^2)` and
#' `dr_loss = mean((DR_score - tau)^2)` on the held-out split. The fixed
#' evaluation forest supplies `m_hat` and DR scores; `e_hat` is the fixed PS.
#' These are surrogate losses, not observed individual-effect errors. Raw
#' Cox/logistic benefit scores are not RMST/risk differences, so requesting
#' either loss for survival or binary outcomes raises an error.
#'
#' Validation is activated only by the two metric parameters. When active,
#' both RATE metrics are reported, and both losses for continuous outcomes;
#' unavailable or unrequested validation metrics are `NA`. Descriptive
#' summaries still use the training split. When inactive, all rows are used
#' for training and all validation columns are `NA`.
#'
#' This held-out split is used for variable ranking and/or size tuning, so
#' the winning metric is not an unbiased final performance estimate. Use
#' independent test data or outer resampling to assess the full selection
#' procedure. The function neither performs that assessment nor refits a
#' selected model to all data. Supplied PS should also be constructed without
#' outcome leakage. Selecting a different metric can activate splitting and
#' thus change the training sample as well as the ranking criterion.
#'
#' A penalized model producing constant scores is retained with SD zero.
#' Backend warnings are recorded in `analysis$warnings`; errors stop with
#' variable/step context. This object is not an `hte_res` and is not input to
#' the grf-specific HTE plotting or RATE functions.
#'
#' @return A plain named list:
#' \describe{
#'   \item{ranking}{Data frame with `rank`, `variable`, `n` (training rows),
#'     `n_eval` (validation rows), `score_sd`,
#'     `score_mean_abs`, `score_median`, `score_mean` (signed mean) and
#'     `score_iqr` (interquartile range), `autoc`, `qini`, `r_loss` and `dr_loss`,
#'     sorted by `imp_metric` in its optimization direction.}
#'   \item{forward}{Data frame with `step`, `added_variable`, `n_vars`, a
#'     `variables` list column, `n`, `n_eval` and the same metrics.}
#'   \item{selected}{Character vector of the selected prefix, or `NULL` when
#'     both selection controls are `NULL`. `character(0)` when `sel_metric`
#'     finds no eligible step and `n_select` is `NULL`. Included variables can
#'     still have zero LASSO coefficients.}
#'   \item{plots}{Named list `ranking`, `forward` and `combined` of ggplot
#'     objects, not printed or saved automatically. Forward panels use separate
#'     y scales and omit unavailable metrics. The combined plot aligns bars
#'     of `imp_metric` (bottom axis) with cumulative `sel_metric` values
#'     (top axis; signed mean score when `sel_metric = NULL`). The pink/red
#'     bar and point mark the first optimum (the best eligible step when
#'     `sel_metric` is set; omitted when none is eligible); a different manually selected
#'     size has a blue dotted line. These markers do not prove generalization.
#'     The secondary axis uses an invertible linear transformation for display only;
#'     the two quantities do not share a numerical scale.}
#'   \item{analysis}{Outcome and treatment mapping, loss, adjustment settings,
#'     sample size, design-column mapping, seed, fit/evaluation settings,
#'     original row indices of the fixed split, PS training-fold IDs (`NULL`
#'     for matching), metric directions, `best_step` (0 when no step is
#'     eligible), `sel_baseline` (`NA` without a baseline), `selected_step`, evaluator
#'     metadata, dependency versions and a warning data frame. `n_select`
#'     retains the caller's manual value; `selected_step` is the returned size.}
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
#' # Optional validation-based selection, using the same 40 candidates.
#' if (requireNamespace("grf", quietly = TRUE)) {
#'   validated <- get_hte_select(d, "z", candidates, ps_var = "ps",
#'     imp_metric = "score_iqr", sel_metric = "autoc",
#'     fit_args = list(nfolds = 3L),
#'     eval_args = list(time = 24, num.trees = 500L))
#'   validated$selected
#'   validated$plots$combined
#' }
#' }
#' @export
get_hte_select <- function(data, cat_var, candidate_var, surv = TRUE,
                           match_var = NULL, ps_var = NULL, n_select = NULL,
                           imp_metric = "score_sd", sel_metric = NULL,
                           fit_args = list(nfolds = 10L, standardize = TRUE),
                           eval_args = list(train_frac = 0.5, adjust_var = NULL,
                                            target = "RMST", time = NULL,
                                            num.trees = 2000L),
                           seed = 123L, verbose = FALSE) {
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data)))
    stop("`data` must be a non-empty data frame with unique column names.",
         call. = FALSE)
  cat_var <- .sens_check_col(cat_var, data, "cat_var", n = 1L)
  candidate_var <- .sens_check_col(candidate_var, data, "candidate_var")
  if (is.null(cat_var) || !length(candidate_var) || anyDuplicated(candidate_var))
    stop("Supply one `cat_var` and nonempty, distinct `candidate_var` names.",
         call. = FALSE)
  available <- c("score_sd", "score_iqr", "score_mean_abs", "score_median", "score_mean",
                  "autoc", "qini", "r_loss", "dr_loss")
  for (name in c("imp_metric", "sel_metric")) {
    metric <- get(name)
    if (name == "sel_metric" && is.null(metric)) next
    if (!is.character(metric) || length(metric) != 1L || is.na(metric) ||
        !metric %in% available)
      stop(sprintf("`%s` must be one of: %s.", name, paste(available, collapse = ", ")),
           call. = FALSE)
  }
  use_eval <- any(c(imp_metric, sel_metric) %in% c("autoc", "qini", "r_loss", "dr_loss"))
  eval_args <- .merge_named_arg(eval_args,
    list(train_frac = 0.5, adjust_var = NULL, target = "RMST", time = NULL,
         num.trees = 2000L), "eval_args")
  fraction <- eval_args$train_frac
  if (!is.numeric(fraction) || length(fraction) != 1L || is.na(fraction) ||
      !is.finite(fraction) || fraction <= 0 || fraction >= 1)
    stop("`eval_args$train_frac` must be a number strictly between 0 and 1.", call. = FALSE)
  eval_args$num.trees <- .hte_select_count(eval_args$num.trees, "eval_args$num.trees",
                                          2L, .Machine$integer.max)
  if (!is.character(eval_args$target) || length(eval_args$target) != 1L ||
      is.na(eval_args$target) || !eval_args$target %in% c("RMST", "survival.probability"))
    stop("`eval_args$target` must be 'RMST' or 'survival.probability'.", call. = FALSE)
  if (!is.null(eval_args$time) && (!is.numeric(eval_args$time) ||
      length(eval_args$time) != 1L || is.na(eval_args$time) ||
      !is.finite(eval_args$time) || eval_args$time <= 0))
    stop("`eval_args$time` must be NULL or a positive finite number.", call. = FALSE)
  adjust_var <- .sens_check_col(eval_args$adjust_var, data, "eval_args$adjust_var")
  if (anyDuplicated(adjust_var))
    stop("`eval_args$adjust_var` must contain distinct names.", call. = FALSE)
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
  if (any(adjust_var %in% roles))
    stop("`eval_args$adjust_var` cannot include treatment, outcome, PS or matching columns.",
         call. = FALSE)
  model_vars <- unique(c(candidate_var, if (use_eval) adjust_var))
  used <- c(roles, model_vars)
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
  constant <- vapply(data[model_vars], function(x) length(unique(x)) < 2L,
                     logical(1))
  if (any(constant))
    stop(sprintf("Constant candidate or adjustment variable(s): %s.",
                 paste(model_vars[constant], collapse = ", ")), call. = FALSE)
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
  if (type != "continuous" && any(c(imp_metric, sel_metric) %in% c("r_loss", "dr_loss")))
    stop("`r_loss` and `dr_loss` currently require a continuous outcome; Cox/logistic benefit scores are not outcome-scale effects.",
         call. = FALSE)
  if (use_eval && type == "survival" && is.null(eval_args$time))
    stop("Specify `eval_args$time` for survival validation metrics.", call. = FALSE)
  if (use_eval && !requireNamespace("grf", quietly = TRUE))
    stop("Install the optional package `grf` for validation metrics.", call. = FALSE)
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
  encoded <- data[model_vars]
  factors <- model_vars[!vapply(encoded, is.numeric, logical(1))]
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
  train <- seq_len(nrow(data))
  validation <- integer()
  evaluation <- NULL
  warnings <- list()
  if (use_eval) {
    if (!is.null(match_id)) {
      ids <- levels(match_id)
      k <- floor(length(ids) * fraction)
      if (k < fit_args$nfolds || length(ids) - k < 2L)
        stop("The split needs at least `fit_args$nfolds` training pairs and two evaluation pairs.",
             call. = FALSE)
      train <- which(match_id %in% ids[sample.int(length(ids), k)])
    } else {
      # Stratify the one fixed split by arm; never repeatedly search for a
      # split which improves a metric.
      train <- sort(unlist(lapply(0:1, function(a) {
        rows <- which(trt == a)
        k <- floor(length(rows) * fraction)
        if (k < 2L || length(rows) - k < 2L)
          stop("Each split needs at least two observations per treatment arm.", call. = FALSE)
        rows[sample.int(length(rows), k)]
      }), use.names = FALSE))
      if (length(train) < fit_args$nfolds)
        stop("The training split has fewer rows than `fit_args$nfolds`.", call. = FALSE)
    }
    validation <- setdiff(seq_len(nrow(data)), train)
    if (type == "survival" && !any(y[train, 2] == 1))
      stop("The training split has no observed events.", call. = FALSE)
    if (type != "survival" && length(unique(y[train])) < 2L)
      stop("The outcome is constant in the training split.", call. = FALSE)
    constant_train <- vapply(data[train, candidate_var, drop = FALSE],
                              function(v) length(unique(v)) < 2L, logical(1))
    if (any(constant_train))
      stop(sprintf("Constant candidate in the training split: %s.",
                   paste(candidate_var[constant_train], collapse = ", ")), call. = FALSE)
    evaluation <- tryCatch(withCallingHandlers(
      .hte_select_evaluation(x[validation, , drop = FALSE],
        if (type == "survival") y[validation, , drop = FALSE] else y[validation],
        trt[validation], if (is.null(ps)) rep(0.5, length(validation)) else ps[validation],
        if (!is.null(match_id)) droplevels(match_id[validation]), type, eval_args, seed),
      warning = function(w) {
        warnings[[length(warnings) + 1L]] <<- data.frame(
          stage = "evaluation", step = NA_integer_, variables = paste(model_vars, collapse = ", "),
          message = conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) stop(paste("HTE evaluation failed:", conditionMessage(e)),
                                    call. = FALSE))
  }
  set.seed(seed)
  foldid <- NULL
  if (is.null(match_id)) {
    # Deal folds within arm x event (or binary outcome) strata so that sparse
    # events cannot collapse into one fold and leave a training fold without them.
    outcome_class <- if (type == "survival") y[train, 2] else
      if (type == "binary") y[train] else rep(0, length(train))
    stratum <- interaction(trt[train], outcome_class, drop = TRUE)
    foldid <- integer(length(train))
    foldid[order(stratum, stats::runif(length(train)))] <-
      sample(fit_args$nfolds)[rep(seq_len(fit_args$nfolds), length.out = length(train))]
    usable <- vapply(seq_len(fit_args$nfolds), function(k) {
      keep <- foldid != k
      all(0:1 %in% trt[train][keep]) && switch(type, continuous = TRUE,
        survival = any(outcome_class[keep] == 1), binary = all(0:1 %in% outcome_class[keep]))
    }, logical(1))
    if (!all(usable))
      stop("Every training fold needs both arms and, for survival or binary outcomes, observed events or both outcome classes; reduce `fit_args$nfolds` or add data.",
           call. = FALSE)
  }
  # `fit` reuses an identical earlier fit; its warnings are recorded again.
  fit_subset <- function(vars, stage, step, fit = NULL) {
    if (is.null(fit)) {
      context <- sprintf("%s %d: %s", stage, step, paste(vars, collapse = ", "))
      if (verbose) cli::cli_inform("{context}")
      indices <- unlist(columns[vars], use.names = FALSE)
      args <- list(x = x[train, indices, drop = FALSE],
                    y = if (type == "survival") y[train, , drop = FALSE] else y[train],
                    trt = trt[train], ps = if (!is.null(ps)) ps[train],
                    match_id = if (!is.null(match_id)) droplevels(match_id[train]),
                    foldid = foldid, loss = loss, fit_args = fit_args, seed = seed, context = context)
      if (use_eval) args <- c(args, list(newx = x[validation, indices, drop = FALSE],
                                        evaluation = evaluation))
      fit <- do.call(.hte_select_fit, args)
      if (!use_eval) fit$statistics <- c(fit$statistics,
        list(autoc = NA_real_, qini = NA_real_, r_loss = NA_real_, dr_loss = NA_real_))
    }
    if (length(fit$warnings))
      warnings[[length(warnings) + 1L]] <<- data.frame(
        stage = stage, step = step, variables = paste(vars, collapse = ", "),
        message = fit$warnings)
    fit
  }
  single <- lapply(seq_along(candidate_var), function(i) {
    fit_subset(candidate_var[i], "single", i)
  })
  ranking <- do.call(rbind, lapply(seq_along(candidate_var), function(i) {
    data.frame(variable = candidate_var[i], n = length(train), n_eval = length(validation),
               single[[i]]$statistics)
  }))
  ranking <- ranking[order(.hte_select_direction(imp_metric) * ranking[[imp_metric]],
                            seq_len(nrow(ranking))), , drop = FALSE]
  rownames(ranking) <- NULL
  ranking <- data.frame(rank = seq_len(nrow(ranking)), ranking)
  steps <- lapply(seq_len(nrow(ranking)), function(k) {
    vars <- ranking$variable[seq_len(k)]
    # The one-variable prefix is the top single-variable model.
    fit <- fit_subset(vars, "forward", k,
                      if (k == 1L) single[[match(vars, candidate_var)]])
    data.frame(step = k, added_variable = vars[k], n_vars = k,
               variables = I(list(vars)), n = length(train), n_eval = length(validation),
               fit$statistics)
  })
  forward <- do.call(rbind, steps)
  best_step <- sel_baseline <- NULL
  if (!is.null(sel_metric)) {
    # Constant-effect reference: a constant priority has zero RATE and zero
    # spread; losses use the best constant effect on the evaluation split.
    sel_baseline <- switch(sel_metric,
      autoc = , qini = , score_sd = , score_iqr = 0,
      r_loss = mean((evaluation$y_residual - evaluation$w_residual *
        sum(evaluation$y_residual * evaluation$w_residual) /
        sum(evaluation$w_residual^2))^2),
      dr_loss = mean((evaluation$dr - mean(evaluation$dr))^2),
      NA_real_)
    direction <- .hte_select_direction(sel_metric)
    value <- direction * forward[[sel_metric]]
    eligible <- forward$score_sd > 0 &
      (is.na(sel_baseline) | value < direction * sel_baseline)
    best_step <- if (any(eligible))
      which(eligible)[which.min(value[eligible])] else 0L
  }
  selected_step <- if (!is.null(n_select)) n_select else best_step
  warning_table <- if (length(warnings)) do.call(rbind, warnings) else
    data.frame(stage = character(), step = integer(), variables = character(),
               message = character())
  versions <- vapply(c("personalized", "glmnet", "survival", if (use_eval) "grf"), function(pkg) {
    as.character(utils::packageVersion(pkg))
  }, character(1))
  list(ranking = ranking, forward = forward,
       selected = if (!is.null(selected_step)) ranking$variable[seq_len(selected_step)] else NULL,
       plots = .hte_select_plots(ranking, forward, selected_step, imp_metric, sel_metric,
                                 best_step),
       analysis = list(backend = "personalized", outcome_type = type,
                       outcome = outcome, cat_var = cat_var,
                       treatment_mapping = unique(data.frame(
                         original = as.character(data[[cat_var]]), encoded = trt)),
                       loss = loss, method = "weighting", n = nrow(data),
                       candidate_var = candidate_var, match_var = match_var,
                       ps_var = ps_var, n_select = n_select,
                       imp_metric = imp_metric, sel_metric = sel_metric,
                       best_step = best_step, sel_baseline = sel_baseline,
                       selected_step = selected_step,
                       rank_direction = if (.hte_select_direction(imp_metric) == 1) "minimize" else "maximize",
                       select_direction = if (is.null(sel_metric)) NULL else
                         if (.hte_select_direction(sel_metric) == 1) "minimize" else "maximize",
                       design_columns = lapply(columns, function(idx) colnames(x)[idx]),
                       seed = seed, fit_args = fit_args, foldid = foldid,
                       eval_args = eval_args, training_rows = train,
                       evaluation_rows = validation,
                       evaluation = if (use_eval) list(backend = "grf", n = length(validation),
                         adjust_var = model_vars, target = if (type == "survival") eval_args$target else
                           if (type == "binary") "risk.difference" else "mean.difference",
                         time = if (type == "survival") eval_args$time else NULL,
                         propensity = if (is.null(ps)) "matched 1:1: 0.5" else "fixed ps_var",
                         nuisance_prediction = "out-of-bag", metric_source = "validation",
                         effect_conversion = if (type == "continuous") "2 * benefit.score" else NULL) else NULL,
                       score_source = "training", versions = versions,
                       warnings = warning_table))
}
