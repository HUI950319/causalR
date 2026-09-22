# =============================================================================
# sens-get.R -- unified sensitivity analysis for unmeasured confounding
# =============================================================================
#
# Architecture (3 layers + shared utilities in utils-sens.R):
#
#   L1  get_sens(data, treat, outcome, ..., method)
#         |
#         +-- L2 backend adapters (one external package each)
#         |     .sens_fit_lm    sensemakr        partial R2 / RV / OVB bounds
#         |     .sens_fit_cox   tipr             tipping point + E-value
#         |     .sens_fit_dml   dml.sensemakr    DML confounding bounds
#         |     .sens_fit_iv    iv.sensemakr     Anderson-Rubin IV bounds
#         |
#         +-- L3 helpers
#               .sens_stats_row   the 18-column standardised $stats row
#               .sens_bounds_tbl  the 12-column standardised $bounds table
#               .sens_dml_bounds  binds $bounds and $bench.table from dml
#               .sens_plt_spec    plt_sens() arguments echoed by print()
#
#   print.sens_res reports both scales and the matching plt_sens() call.
# =============================================================================

.SENS_BENCH_DEFAULTS <- list(
  k_treat     = 1,
  k_out       = NULL,
  bound       = NULL,
  bound_label = "Manual bound")

.SENS_EVALUE_DEFAULTS <- list(
  rare           = FALSE,
  confounder     = c("continuous", "binary"),
  smd            = 1,
  exposed_prev   = NULL,
  unexposed_prev = NULL)

.SENS_DML_DEFAULTS <- list(
  model          = c("plm", "npm"),
  target         = "ate",
  reg            = "ranger",
  cf_folds       = 5L,
  cf_reps        = 1L,
  cf_seed        = NULL,
  ps_trim        = 0.01,
  dirty_tuning   = TRUE,
  rho2           = 1,
  combine_method = "median")

.SENS_IV_DEFAULTS <- list(
  parm = c("iv", "fs", "rf"),
  min  = TRUE)

# Which arguments each method actually consumes. Anything supplied outside its
# own method is a hard error: silently ignoring it hides the most common
# mistake, a backend list that never took effect.
.SENS_APPLIES <- list(
  time        = "cox",
  instrument  = "iv",
  bench_var   = c("lm", "dml", "iv"),
  bench_args  = c("lm", "dml", "iv"),
  evalue_args = "cox",
  dml_args    = "dml",
  iv_args     = "iv")


# ---- L3 standardised output ------------------------------------------------

#' @keywords internal
#' @noRd
.sens_stats_row <- function(method, estimand, term, scale, estimate,
                            std.error = NA_real_, statistic = NA_real_,
                            conf.low = NA_real_, conf.high = NA_real_,
                            dof = NA_real_, rv_q = NA_real_, rv_qa = NA_real_,
                            xrv_qa = NA_real_, r2yd_x = NA_real_,
                            evalue_point = NA_real_, evalue_ci = NA_real_,
                            tip_effect = NA_real_, tip_n = NA_real_) {
  tibble::tibble(
    method = method, estimand = estimand, term = term, scale = scale,
    estimate = as.numeric(estimate), std.error = as.numeric(std.error),
    statistic = as.numeric(statistic), conf.low = as.numeric(conf.low),
    conf.high = as.numeric(conf.high), dof = as.numeric(dof),
    rv_q = as.numeric(rv_q), rv_qa = as.numeric(rv_qa),
    xrv_qa = as.numeric(xrv_qa), r2yd_x = as.numeric(r2yd_x),
    evalue_point = as.numeric(evalue_point), evalue_ci = as.numeric(evalue_ci),
    tip_effect = as.numeric(tip_effect), tip_n = as.numeric(tip_n))
}

#' @keywords internal
#' @noRd
.sens_bounds_tbl <- function(method, estimand, bound_label, r2_treat, r2_out,
                             adj_estimate = NA_real_, adj_low = NA_real_,
                             adj_high = NA_real_, adj_conf.low = NA_real_,
                             adj_conf.high = NA_real_, adj_se = NA_real_,
                             adj_statistic = NA_real_) {
  tibble::tibble(
    method = method, estimand = estimand, bound_label = bound_label,
    r2_treat = as.numeric(r2_treat), r2_out = as.numeric(r2_out),
    adj_estimate = as.numeric(adj_estimate), adj_low = as.numeric(adj_low),
    adj_high = as.numeric(adj_high), adj_conf.low = as.numeric(adj_conf.low),
    adj_conf.high = as.numeric(adj_conf.high), adj_se = as.numeric(adj_se),
    adj_statistic = as.numeric(adj_statistic))
}


# ---- L2 backend adapters ---------------------------------------------------

#' @keywords internal
#' @noRd
.sens_fit_lm <- function(data, treat, outcome, adj_var, bench_var,
                         q, alpha, bench_args) {
  if (!requireNamespace("sensemakr", quietly = TRUE))
    stop("Package 'sensemakr' is required for get_sens(method = \"lm\")",
         call. = FALSE)

  form <- stats::reformulate(.sens_quote_names(c(treat, adj_var)),
                              response = as.name(outcome))
  fit <- stats::lm(form, data = data)
  tcol <- .sens_one_term(fit, treat)

  # A factor benchmark expands to several dummies; sensemakr takes those as one
  # named group, which keeps a k-fold benchmark on the variable rather than on
  # an arbitrary level.
  bench <- NULL
  if (!is.null(bench_var)) {
    bench <- .sens_coef_terms(fit, bench_var)
    for (v in bench_var) {
      if (!length(bench[[v]]))
        stop(sprintf("`bench_var` '%s' has no coefficient in the fitted model.",
                     v), call. = FALSE)
    }
    if (all(lengths(bench) == 1L)) bench <- unlist(bench, use.names = FALSE)
  }

  bnd <- bench_args$bound
  s <- sensemakr::sensemakr(
    model                = fit,
    treatment            = tcol,
    benchmark_covariates = bench,
    kd                   = bench_args$k_treat,
    ky                   = if (is.null(bench_args$k_out)) bench_args$k_treat
                           else bench_args$k_out,
    q                    = q,
    alpha                = alpha,
    r2dz.x               = if (is.null(bnd)) NULL else bnd[1L],
    r2yz.dx              = if (is.null(bnd)) NULL else bnd[2L],
    bound_label          = bench_args$bound_label)

  st <- s$sensitivity_stats
  if (!all(c("estimate", "se", "rv_q", "rv_qa", "r2yd.x", "dof") %in% names(st)))
    stop("Unexpected `sensitivity_stats` layout from sensemakr ",
         utils::packageVersion("sensemakr"), ".", call. = FALSE)
  ci <- stats::confint(fit, tcol, level = 1 - alpha)

  list(
    fit = fit, sens = s,
    stats = .sens_stats_row(
      "lm", "beta", tcol, "difference",
      estimate = st$estimate, std.error = st$se, statistic = st$t_statistic,
      conf.low = ci[1L, 1L], conf.high = ci[1L, 2L], dof = st$dof,
      rv_q = st$rv_q, rv_qa = st$rv_qa,
      xrv_qa = as.numeric(sensemakr::xrv(fit, covariates = tcol,
                                         q = q, alpha = alpha)),
      r2yd_x = st$r2yd.x),
    bounds = if (is.null(s$bounds) || !nrow(s$bounds)) NULL else
      .sens_bounds_tbl(
        "lm", "beta", s$bounds$bound_label, s$bounds$r2dz.x, s$bounds$r2yz.dx,
        adj_estimate  = s$bounds$adjusted_estimate,
        adj_low       = s$bounds$adjusted_estimate,
        adj_high      = s$bounds$adjusted_estimate,
        adj_conf.low  = s$bounds$adjusted_lower_CI,
        adj_conf.high = s$bounds$adjusted_upper_CI,
        adj_se        = s$bounds$adjusted_se,
        adj_statistic = s$bounds$adjusted_t))
}

#' @keywords internal
#' @noRd
.sens_fit_cox <- function(data, treat, outcome, time, adj_var,
                          alpha, evalue_args) {
  for (pkg in c("survival", "tipr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required for get_sens(method = \"cox\")",
                   pkg), call. = FALSE)
  }

  response <- substitute(survival::Surv(TIME, OUTCOME),
                         list(TIME = as.name(time), OUTCOME = as.name(outcome)))
  form <- stats::reformulate(.sens_quote_names(c(treat, adj_var)),
                              response = response)
  fit <- survival::coxph(form, data = data)
  tcol <- .sens_one_term(fit, treat)

  sm <- summary(fit, conf.int = 1 - alpha)
  hr <- unname(sm$conf.int[tcol, 1L])
  lo <- unname(sm$conf.int[tcol, 3L])
  hi <- unname(sm$conf.int[tcol, 4L])

  rare <- isTRUE(evalue_args$rare)
  rr  <- .sens_hr_to_rr(hr, rare)
  rlo <- .sens_hr_to_rr(lo, rare)
  rhi <- .sens_hr_to_rr(hi, rare)

  # E-value of the interval is driven by the confidence limit closest to the
  # null; an interval covering the null needs no confounding at all, so its
  # E-value is 1 by definition.
  ci_lim <- if (rlo <= 1 && rhi >= 1) NA_real_ else if (rr < 1) rhi else rlo
  evalue_ci <- if (is.na(ci_lim)) 1 else tipr::e_value(ci_lim)

  # tipr's `hr_correction` and the `rare` switch are the same VanderWeele-Ding
  # transform seen from opposite sides, so one argument drives both.
  hrc  <- !rare
  near <- if (hr < 1) hi else lo
  tp <- if (identical(evalue_args$confounder, "binary")) {
    if (is.null(evalue_args$exposed_prev) || is.null(evalue_args$unexposed_prev))
      stop("`evalue_args$exposed_prev` and `$unexposed_prev` are required when `confounder = \"binary\"`.",
           call. = FALSE)
    tipr::tip_hr_with_binary(
      near,
      exposed_confounder_prev   = evalue_args$exposed_prev,
      unexposed_confounder_prev = evalue_args$unexposed_prev,
      verbose = FALSE, hr_correction = hrc)
  } else {
    tipr::tip_hr(near, exposure_confounder_effect = evalue_args$smd,
                 verbose = FALSE, hr_correction = hrc)
  }

  list(
    fit = fit, sens = NULL,
    stats = .sens_stats_row(
      "cox", "HR", tcol, "HR",
      estimate = hr,
      std.error = sm$coefficients[tcol, "se(coef)"],   # log-HR scale
      statistic = sm$coefficients[tcol, "z"],
      conf.low = lo, conf.high = hi,
      evalue_point = tipr::e_value(rr),
      evalue_ci    = evalue_ci,
      tip_effect   = tp$confounder_outcome_effect[1L],
      tip_n        = tp$n_unmeasured_confounders[1L]),
    bounds = NULL)
}

#' @keywords internal
#' @noRd
.sens_dml_bounds <- function(s, target) {
  parts <- list()
  if (!is.null(s$bounds) && nrow(s$bounds))
    parts[[length(parts) + 1L]] <- .sens_bounds_tbl(
      "dml", s$bounds$target, s$bounds$bound.label,
      s$bounds$cf.d, s$bounds$cf.y,
      adj_low = s$bounds$theta.minus, adj_high = s$bounds$theta.plus,
      adj_conf.low = s$bounds$lwr, adj_conf.high = s$bounds$upr)
  if (!is.null(s$bench.table) && nrow(s$bench.table))
    parts[[length(parts) + 1L]] <- .sens_bounds_tbl(
      "dml", s$bench.table$target, s$bench.table$bound.label,
      s$bench.table$cf.d, s$bench.table$cf.y,
      adj_low = s$bench.table$theta.minus, adj_high = s$bench.table$theta.plus,
      adj_conf.low = s$bench.table$lwr, adj_conf.high = s$bench.table$upr)
  if (!length(parts)) return(NULL)
  out <- do.call(rbind, parts)
  out[out$estimand %in% target, , drop = FALSE]
}

#' @keywords internal
#' @noRd
.sens_fit_dml <- function(data, treat, outcome, adj_var, bench_var,
                          q, alpha, bench_args, dml_args, verbose) {
  if (!requireNamespace("dml.sensemakr", quietly = TRUE))
    stop("Package 'dml.sensemakr' is required for get_sens(method = \"dml\")",
         call. = FALSE)

  xm <- .sens_model_matrix(data, adj_var)
  d  <- data[[treat]]
  if (!is.numeric(d)) {
    lv <- levels(factor(d))
    if (length(lv) != 2L)
      stop("`treat` must be numeric or two-level for method = \"dml\".",
           call. = FALSE)
    d <- as.numeric(factor(d, levels = lv)) - 1
  }

  fit <- dml.sensemakr::dml(
    y = data[[outcome]], d = d, x = xm,
    model        = dml_args$model,
    target       = dml_args$target,
    cf.folds     = dml_args$cf_folds,
    cf.reps      = dml_args$cf_reps,
    cf.seed      = dml_args$cf_seed,
    ps.trim      = dml_args$ps_trim,
    reg          = dml_args$reg,
    dirty.tuning = dml_args$dirty_tuning,
    verbose      = isTRUE(verbose),
    warnings     = FALSE)

  cm  <- dml_args$combine_method
  est <- stats::coef(fit, combine.method = cm)
  ses <- dml.sensemakr::se(fit, combine.method = cm)
  ci  <- stats::confint(fit, level = 1 - alpha, combine.method = cm)

  # dml.sensemakr parameterises the null as `theta`, not as a fraction `q` of
  # the estimate; iv.sensemakr encodes exactly this mapping internally.
  theta <- (1 - q) * unname(est[[1L]])

  bnd <- bench_args$bound
  call_sens <- function() dml.sensemakr::sensemakr(
    fit,
    benchmark_covariates = bench_var,
    cf.y        = if (is.null(bnd)) NULL else bnd[2L],
    cf.d        = if (is.null(bnd)) NULL else bnd[1L],
    rho2        = dml_args$rho2,
    kD          = bench_args$k_treat,
    kY          = if (is.null(bench_args$k_out)) bench_args$k_treat
                  else bench_args$k_out,
    bound_label = bench_args$bound_label,
    theta       = theta, alpha = alpha)
  s <- if (isTRUE(verbose)) call_sens() else {
    utils::capture.output(out <- call_sens()); out
  }

  rv <- s$sensitivity_stats
  if (!is.matrix(rv) || !all(c("rv", "rva") %in% colnames(rv)))
    stop("Unexpected `sensitivity_stats` layout from dml.sensemakr ",
         utils::packageVersion("dml.sensemakr"), ".", call. = FALSE)
  tg <- rownames(rv)
  xrv <- dml.sensemakr::extreme_robustness_value(fit, theta = theta, alpha = alpha)

  list(
    fit = fit, sens = s,
    stats = .sens_stats_row(
      "dml", tg, treat, "difference",
      estimate = est[tg], std.error = ses[tg],
      statistic = est[tg] / ses[tg],
      conf.low = ci[tg, 1L], conf.high = ci[tg, 2L],
      rv_q = rv[tg, "rv"], rv_qa = rv[tg, "rva"],
      xrv_qa = xrv[tg]),
    bounds = .sens_dml_bounds(s, tg))
}

#' @keywords internal
#' @noRd
.sens_fit_iv <- function(data, treat, outcome, instrument, adj_var, bench_var,
                         q, alpha, bench_args, iv_args) {
  if (!requireNamespace("iv.sensemakr", quietly = TRUE))
    stop("Package 'iv.sensemakr' is required for get_sens(method = \"iv\")",
         call. = FALSE)

  xm  <- .sens_model_matrix(data, adj_var)
  # iv.sensemakr uses the same names for data columns and lm coefficients.
  # Syntactic matrix names keep those identical, including after benchmarking.
  original_names <- colnames(xm)
  colnames(xm) <- make.names(original_names, unique = TRUE)
  bench <- bench_var
  mapped <- match(bench_var, original_names)
  bench[!is.na(mapped)] <- colnames(xm)[mapped[!is.na(mapped)]]
  fit <- iv.sensemakr::iv_fit(y = data[[outcome]], d = data[[treat]],
                              z = data[[instrument]],
                              x = if (ncol(xm)) xm else NULL,
                              h0 = 0, alpha = alpha)

  bnd <- bench_args$bound
  s <- iv.sensemakr::sensemakr(
    fit,
    benchmark_covariates = bench,
    kz          = bench_args$k_treat,
    ky          = if (is.null(bench_args$k_out)) bench_args$k_treat
                  else bench_args$k_out,
    r2zw.x      = if (is.null(bnd)) NULL else bnd[1L],
    r2y0w.zx    = if (is.null(bnd)) NULL else bnd[2L],
    bound_label = bench_args$bound_label,
    q = q, alpha = alpha, min = isTRUE(iv_args$min))

  parm <- intersect(iv_args$parm, names(s$sensitivity_stats))
  if (!length(parm))
    stop(sprintf("None of `iv_args$parm` is reported by iv.sensemakr; available: %s.",
                 paste0("\"", names(s$sensitivity_stats), "\"", collapse = ", ")),
         call. = FALSE)
  st <- do.call(rbind, s$sensitivity_stats[parm])

  bnd_tbl <- NULL
  b <- s$bounds[[parm[1L]]]
  if (!is.null(b) && nrow(b)) {
    # The outcome-side column is named per estimand (r2y0w.zx / r2dw.zx /
    # r2yw.zx); only its position is stable, so take it positionally.
    bnd_tbl <- .sens_bounds_tbl(
      "iv", parm[1L], b$bound_label, b$r2zw.x, b[[3L]],
      adj_conf.low = b$lwr, adj_conf.high = b$upr)
  }

  list(
    fit = fit, sens = s,
    stats = .sens_stats_row(
      "iv", parm, instrument, "difference",
      estimate = st$estimate, statistic = st$t.value,
      conf.low = st$lwr, conf.high = st$upr, dof = st$dof,
      rv_qa = st$rv_qa, xrv_qa = st$xrv_qa),
    bounds = bnd_tbl)
}


# ---- L1 public entry point -------------------------------------------------

#' Sensitivity analysis for unmeasured confounding
#'
#' Single entry point that routes an observational estimate to one of four
#' omitted-variable-bias backends and returns their results in one shared
#' shape, so that a linear model, a Cox model, a double machine learning fit
#' and an instrumental-variable fit can be reported side by side.
#'
#' @param data A data frame holding every column named below.
#' @param treat,outcome Length-1 character. The treatment and outcome columns.
#'   `treat` must resolve to exactly one model coefficient: numeric, or a
#'   two-level factor.
#' @param time Length-1 character. Follow-up time column. Required for, and
#'   only accepted by, `method = "cox"`.
#' @param instrument Length-1 character. Instrument column. Required for, and
#'   only accepted by, `method = "iv"`.
#' @param adj_var Character vector of measured covariates to adjust for, or
#'   `NULL`.
#' @param bench_var Character vector of covariates used to calibrate how strong
#'   an unmeasured confounder would have to be, expressed as multiples of an
#'   observed covariate. Must be a subset of `adj_var`. Not used by
#'   `method = "cox"`.
#' @param method Backend selector: `"lm"` (sensemakr partial \eqn{R^2}),
#'   `"cox"` (tipr tipping point and E-value), `"dml"` (dml.sensemakr) or
#'   `"iv"` (iv.sensemakr).
#' @param q Fraction of the estimate a confounder would have to explain away
#'   for the result to be considered overturned. `q = 1` asks what it takes to
#'   reach the null. Must be 1 for `method = "cox"`.
#' @param conf_level Confidence level. The backends receive `1 - conf_level`
#'   as their `alpha`; there is deliberately no separate `alpha` argument.
#' @param bench_args Named list controlling benchmark and manual confounding
#'   scenarios, used by `method` `"lm"`, `"dml"` and `"iv"`.
#'   \describe{
#'     \item{`k_treat`}{Positive numeric vector. How many times as strongly as
#'       `bench_var` the confounder is related to the treatment (or, for
#'       `"iv"`, to the instrument). Default `1`.}
#'     \item{`k_out`}{Positive numeric vector or `NULL`. Same multiple on the
#'       outcome side. `NULL` (default) reuses `k_treat`.}
#'     \item{`bound`}{Length-2 numeric in (0, 1) giving one manual scenario as
#'       `c(treatment-side R2, outcome-side R2)`, or `NULL` (default). Mapped
#'       to `r2dz.x`/`r2yz.dx` for `"lm"`, `cf.d`/`cf.y` for `"dml"` and
#'       `r2zw.x`/`r2y0w.zx` for `"iv"`.}
#'     \item{`bound_label`}{Length-1 character labelling that manual scenario.
#'       Default `"Manual bound"`.}
#'   }
#' @param evalue_args Named list for `method = "cox"` only.
#'   \describe{
#'     \item{`rare`}{Logical, default `FALSE`. `TRUE` treats the outcome as
#'       rare, so the hazard ratio is used directly as a risk ratio. `FALSE`
#'       applies the VanderWeele-Ding square-root transform. The same switch
#'       drives `tipr`'s `hr_correction`, which it sets to `!rare`.}
#'     \item{`confounder`}{`"continuous"` (default) or `"binary"`, selecting
#'       [tipr::tip_hr()] or [tipr::tip_hr_with_binary()].}
#'     \item{`smd`}{Positive numeric, default `1`. Standardised mean difference
#'       of a continuous confounder between treatment groups.}
#'     \item{`exposed_prev`,`unexposed_prev`}{Numeric in (0, 1) or `NULL`
#'       (default). Confounder prevalence in each arm; both are required when
#'       `confounder = "binary"`.}
#'   }
#' @param dml_args Named list for `method = "dml"` only. Fields map to
#'   [dml.sensemakr::dml()] and its `sensemakr()` method: `model`
#'   (`"plm"`/`"npm"`), `target` (`"ate"`, `"att"`, `"atu"`), `reg` (caret
#'   learner, default `"ranger"`), `cf_folds` (default `5L`), `cf_reps`
#'   (`1L`), `cf_seed` (`NULL`), `ps_trim` (`0.01`), `dirty_tuning` (`TRUE`),
#'   `rho2` (`1`, how adversarial the confounder is) and `combine_method`
#'   (`"median"` or `"mean"`).
#' @param iv_args Named list for `method = "iv"` only, with `parm` (any of
#'   `"iv"`, `"fs"`, `"rf"`; all three by default) and `min` (logical,
#'   default `TRUE`, passed straight to [iv.sensemakr::sensemakr()]).
#' @param verbose Logical. `TRUE` prints backend progress and reports dropped
#'   incomplete rows. Default `FALSE`.
#'
#' @section Robustness values versus E-values:
#' The two headline quantities are on different scales and must not be
#' compared. A robustness value (`rv_q`, `rv_qa`, `xrv_qa`) is a partial
#' \eqn{R^2} in \eqn{[0, 1)}: the share of residual variation in both
#' treatment and outcome that a confounder would have to explain. An E-value
#' (`evalue_point`, `evalue_ci`) is a risk-ratio multiplier in
#' \eqn{[1, \infty)}: how strongly a confounder would have to be associated
#' with both exposure and outcome. There is no monotone mapping between them,
#' so they are kept in separate columns and each is `NA` outside its own
#' backend.
#'
#' @section Choosing rare:
#' `evalue_args$rare` decides whether the hazard ratio is treated as a risk
#' ratio. The default `FALSE` applies the VanderWeele-Ding transform
#' \eqn{RR = (1 - 0.5^{\sqrt{HR}}) / (1 - 0.5^{\sqrt{1/HR}})}, which shrinks
#' the ratio towards the null and therefore reports a smaller, more
#' conservative E-value. This differs from `tipr`'s own default, which leaves
#' the hazard ratio untransformed. Set `rare = TRUE` only when the outcome is
#' genuinely rare over follow-up. The E-values match
#' `EValue::evalues.HR()` exactly on both settings; `EValue` itself is not
#' used because its `estimate` methods clash with `lava`, which
#' `dml.sensemakr` loads.
#'
#' @return An object of class `sens_res`: a list of
#'   \describe{
#'     \item{`stats`}{One tibble row per estimand, always with the same 18
#'       columns: `method`, `estimand`, `term`, `scale`, `estimate`,
#'       `std.error`, `statistic`, `conf.low`, `conf.high`, `dof`, `rv_q`,
#'       `rv_qa`, `xrv_qa`, `r2yd_x`, `evalue_point`, `evalue_ci`,
#'       `tip_effect`, `tip_n`. Columns a backend does not produce are `NA`:
#'       the robustness and partial-\eqn{R^2} columns for `"cox"`, the E-value
#'       and tipping-point columns for the other three, `std.error` for
#'       `"iv"`, and `rv_q` for `"iv"` (upstream reports only the
#'       alpha-adjusted value). For `"cox"`, `estimate` and the confidence
#'       limits are on the hazard-ratio scale while `std.error` is on the
#'       log-hazard scale.}
#'     \item{`bounds`}{Tibble of confounding scenarios, or `NULL` when none
#'       was requested and for `method = "cox"`. Columns `method`, `estimand`,
#'       `bound_label`, `r2_treat`, `r2_out`, `adj_estimate`, `adj_low`,
#'       `adj_high`, `adj_conf.low`, `adj_conf.high`, `adj_se`,
#'       `adj_statistic`. DML and IV bounds are intervals rather than point
#'       adjustments, so `adj_estimate`, `adj_se` and `adj_statistic` are `NA`
#'       there.}
#'     \item{`fit`}{The fitted model: `lm`, `coxph`, `dml` or `iv_fit`.}
#'     \item{`sens`}{The backend sensitivity object, or `NULL` for `"cox"`.}
#'   }
#'   Analysis metadata is attached as `attr(x, "analysis")`; `n` counts rows
#'   used in the fitted model after missing-value exclusion (subjects,
#'   including censored subjects, for Cox models). Note that dplyr verbs drop
#'   the `method` attribute carried by `stats` and `bounds`.
#'   Column names are treated literally, including spaces and punctuation.
#'   The IV backend requires syntactic covariate names; it uses
#'   [make.names()] internally, so its model and benchmark labels use those
#'   names while the analysis metadata retains the supplied names.
#'
#' @seealso [plt_sens()] for the matching plots.
#'
#' @examplesIf requireNamespace("sensemakr", quietly = TRUE)
#' \donttest{
#' # Linear outcome: partial R2 and robustness values
#' get_sens(sensemakr::darfur,
#'          treat = "directlyharmed", outcome = "peacefactor",
#'          adj_var = c("age", "farmer_dar", "herder_dar", "pastvoted",
#'                      "hhsize_darfur", "female", "village"),
#'          bench_var = "female", method = "lm",
#'          bench_args = list(k_treat = 1:3))
#' }
#'
#' @examplesIf requireNamespace("survival", quietly = TRUE) && requireNamespace("tipr", quietly = TRUE)
#' \donttest{
#' # Cox / HR: tipping point and E-value
#' lung <- stats::na.omit(survival::lung[, c("time", "status", "sex", "age")])
#' lung$status <- lung$status - 1L
#' get_sens(lung, treat = "sex", outcome = "status", time = "time",
#'          adj_var = "age", method = "cox")
#' }
#'
#' @export
get_sens <- function(data,
                     treat,
                     outcome,
                     time        = NULL,
                     instrument  = NULL,
                     adj_var     = NULL,
                     bench_var   = NULL,
                     method      = c("lm", "cox", "dml", "iv"),
                     q           = 1,
                     conf_level  = 0.95,
                     bench_args  = list(k_treat     = 1,
                                        k_out       = NULL,
                                        bound       = NULL,
                                        bound_label = "Manual bound"),
                     evalue_args = list(rare           = FALSE,
                                        confounder     = c("continuous", "binary"),
                                        smd            = 1,
                                        exposed_prev   = NULL,
                                        unexposed_prev = NULL),
                     dml_args    = list(model          = c("plm", "npm"),
                                        target         = "ate",
                                        reg            = "ranger",
                                        cf_folds       = 5L,
                                        cf_reps        = 1L,
                                        cf_seed        = NULL,
                                        ps_trim        = 0.01,
                                        dirty_tuning   = TRUE,
                                        rho2           = 1,
                                        combine_method = "median"),
                     iv_args     = list(parm = c("iv", "fs", "rf"),
                                        min  = TRUE),
                     verbose     = FALSE) {

  method <- match.arg(method)

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  if (!is.numeric(q) || length(q) != 1L || is.na(q) || q < 0 || q > 1)
    stop("`q` must be a single number between 0 and 1.", call. = FALSE)
  if (identical(method, "cox") && !isTRUE(all.equal(q, 1)))
    stop("`q` must be 1 for method = \"cox\"; tipping points and E-values are defined against the null only.",
         call. = FALSE)
  alpha <- 1 - conf_level

  supplied <- c(time        = !missing(time),
                instrument  = !missing(instrument),
                bench_var   = !missing(bench_var),
                bench_args  = !missing(bench_args),
                evalue_args = !missing(evalue_args),
                dml_args    = !missing(dml_args),
                iv_args     = !missing(iv_args))
  for (nm in names(supplied)) {
    if (supplied[[nm]] && !method %in% .SENS_APPLIES[[nm]])
      stop(sprintf("`%s` does not apply to method = \"%s\"; it is only used by method %s.",
                   nm, method,
                   paste0("\"", .SENS_APPLIES[[nm]], "\"", collapse = " / ")),
           call. = FALSE)
  }

  if (!is.data.frame(data) || !nrow(data))
    stop("`data` must be a non-empty data frame.", call. = FALSE)
  treat      <- .sens_check_col(treat, data, "treat", n = 1L)
  outcome    <- .sens_check_col(outcome, data, "outcome", n = 1L)
  time       <- .sens_check_col(time, data, "time", n = 1L)
  instrument <- .sens_check_col(instrument, data, "instrument", n = 1L)
  adj_var    <- .sens_check_col(adj_var, data, "adj_var")
  bench_var  <- .sens_check_col(bench_var, data, "bench_var")
  if (!is.null(bench_var) && !all(bench_var %in% adj_var))
    stop(sprintf("`bench_var` must be a subset of `adj_var`; %s not adjusted for.",
                 paste0("`", setdiff(bench_var, adj_var), "`", collapse = ", ")),
         call. = FALSE)
  if (identical(method, "cox") && is.null(time))
    stop("`time` is required for method = \"cox\".", call. = FALSE)
  if (identical(method, "iv") && is.null(instrument))
    stop("`instrument` is required for method = \"iv\".", call. = FALSE)

  bench_args  <- .merge_named_arg(bench_args,  .SENS_BENCH_DEFAULTS,  "bench_args")
  evalue_args <- .merge_named_arg(evalue_args, .SENS_EVALUE_DEFAULTS, "evalue_args")
  dml_args    <- .merge_named_arg(dml_args,    .SENS_DML_DEFAULTS,    "dml_args")
  iv_args     <- .merge_named_arg(iv_args,     .SENS_IV_DEFAULTS,     "iv_args")
  evalue_args$confounder <- match.arg(evalue_args$confounder,
                                      c("continuous", "binary"))
  dml_args$model  <- match.arg(dml_args$model, c("plm", "npm"))
  dml_args$target <- match.arg(dml_args$target, c("ate", "att", "atu"),
                               several.ok = TRUE)
  iv_args$parm    <- match.arg(iv_args$parm, c("iv", "fs", "rf"),
                               several.ok = TRUE)
  if (!is.null(bench_args$bound) &&
      (!is.numeric(bench_args$bound) || length(bench_args$bound) != 2L ||
       any(bench_args$bound <= 0) || any(bench_args$bound >= 1)))
    stop("`bench_args$bound` must be `NULL` or two numbers in (0, 1).",
         call. = FALSE)

  used <- c(treat, outcome, time, instrument, adj_var)
  if (method %in% c("dml", "iv")) data <- .sens_complete(data, used, verbose)

  res <- switch(
    method,
    lm  = .sens_fit_lm(data, treat, outcome, adj_var, bench_var,
                       q, alpha, bench_args),
    cox = .sens_fit_cox(data, treat, outcome, time, adj_var,
                        alpha, evalue_args),
    dml = .sens_fit_dml(data, treat, outcome, adj_var, bench_var,
                        q, alpha, bench_args, dml_args, verbose),
    iv  = .sens_fit_iv(data, treat, outcome, instrument, adj_var, bench_var,
                       q, alpha, bench_args, iv_args),
    stop(sprintf("Unsupported method: '%s'", method), call. = FALSE))

  backend <- c(lm = "sensemakr", cox = "tipr",
               dml = "dml.sensemakr", iv = "iv.sensemakr")[[method]]
  attr(res$stats, "method") <- method
  if (!is.null(res$bounds)) attr(res$bounds, "method") <- method

  structure(
    list(stats = res$stats, bounds = res$bounds, fit = res$fit, sens = res$sens),
    class = c("sens_res", "list"),
    analysis = list(
      method = method, backend = backend,
      backend_version = as.character(utils::packageVersion(backend)),
      treat = treat, outcome = outcome, time = time, instrument = instrument,
      adj_var = adj_var, bench_var = bench_var, bench_args = bench_args,
      q = q, conf_level = conf_level, alpha = alpha,
      n = switch(method, lm = stats::nobs(res$fit), cox = res$fit$n, nrow(data)),
      evalue_args = if (identical(method, "cox")) evalue_args else NULL,
      rho2 = if (identical(method, "dml")) dml_args$rho2 else NULL,
      seed = if (identical(method, "dml")) dml_args$cf_seed else NULL,
      call = match.call()))
}


# ---- L3 print --------------------------------------------------------------

# Window that comfortably contains both the robustness value and every
# benchmark point, so the echoed `lim` really shows the contour of interest.
#' @keywords internal
#' @noRd
.sens_plt_spec <- function(x) {
  a <- attr(x, "analysis")
  if (identical(a$method, "cox")) {
    # The tipping point can sit either side of 1 depending on the direction of
    # the effect, so the window has to span both rather than start at 1.
    rng <- range(c(1, x$stats$tip_effect), na.rm = TRUE)
    pad <- diff(rng) * 0.25
    if (!is.finite(pad) || pad <= 0) pad <- 0.25
    return(list(type = "tip",
                lim = c(max(0.05, floor((rng[1L] - pad) * 20) / 20),
                        ceiling((rng[2L] + pad) * 20) / 20),
                estimand = NULL))
  }
  # Per axis, not one square window: the treatment-side and outcome-side
  # partial R2 routinely differ by an order of magnitude (an IV benchmark sits
  # at 0.007 against 0.22), and a single window then hides the whole contour
  # against one axis. The robustness value enters both, since it is the point
  # the reader is looking for.
  rv <- c(x$stats$rv_q, x$stats$rv_qa)
  axis_top <- function(v) {
    v <- c(v, rv)
    v <- v[is.finite(v)]
    if (!length(v) || max(v) <= 0) return(0.15)
    min(0.9, ceiling(max(v) * 1.25 * 1000) / 1000)
  }
  list(type = "contour",
       lim = c(axis_top(x$bounds$r2_treat), axis_top(x$bounds$r2_out)),
       estimand = if (identical(a$method, "lm")) NULL else x$stats$estimand[1L])
}

#' @export
#' @noRd
print.sens_res <- function(x, ...) {
  a <- attr(x, "analysis")
  cat(sprintf("<sens_res> method = \"%s\" (%s %s), n = %d, conf_level = %s, q = %s\n",
              a$method, a$backend, a$backend_version, a$n,
              format(a$conf_level), format(a$q)))
  cat(sprintf("  treat = %s, outcome = %s%s%s\n",
              a$treat, a$outcome,
              if (is.null(a$time)) "" else paste0(", time = ", a$time),
              if (is.null(a$instrument)) "" else
                paste0(", instrument = ", a$instrument)))
  cat("\n")
  print(x$stats)
  if (!is.null(x$bounds)) {
    cat("\nConfounding scenarios:\n")
    print(x$bounds)
  }
  cat("\n# RV / xRV are partial R2 in [0, 1); E-values are risk-ratio",
      "multipliers in [1, Inf). Not comparable.\n")
  spec <- .sens_plt_spec(x)
  cat("# plt_sens: type = \"", spec$type, "\", lim = ",
      .sens_fmt_vec(spec$lim),
      if (is.null(spec$estimand)) "" else
        paste0(", estimand = \"", spec$estimand, "\""),
      "\n", sep = "")
  invisible(x)
}
