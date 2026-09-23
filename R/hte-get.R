# =============================================================================
# hte-get.R -- heterogeneous treatment effects from causal forests
# =============================================================================
#
# Architecture:
#
#   L1  get_hte()          validate, fit one grf forest, assemble the result
#   L2  .hte_arm_scores()  arm-specific AIPW scores backed out of the forest
#   L2  .hte_estimate()    one row per estimand x measure for a set of units
#   L2  .hte_dr_var()      doubly robust CATE along one covariate + p_het
#                          (shared with plt_hte_dep())
#   L3  print.hte_res()
#
# =============================================================================

# Estimand -> grf `target.sample`, in the order rows are reported.
.HTE_ESTIMANDS <- c(ATE = "all", ATT = "treated", ATC = "control",
                    ATO = "overlap")
.HTE_MEASURES  <- c("diff", "ratio", "OR")
# Natural-spline df of the doubly robust curve for a continuous covariate.
.HTE_SPLINE_DF <- 3L


# ---- L2 scores and estimates -----------------------------------------------

# grf keeps m(x) = E[Y | X] (`Y.hat`), e(x) (`W.hat`) and the out-of-bag CATE,
# and mu_w(x) = m(x) + (w - e(x)) tau(x) recovers both arms -- the identity grf
# itself uses for its ATT / ATC estimators. The residual Y - mu_W(X) is the
# plain outcome residual for causal_forest; causal_survival_forest keeps only
# its censoring-adjusted form, `_psi$numerator` = (W - e) (Y* - m), whose
# `denominator` is (W - e)^2. Either way g1 - g0 reproduces grf::get_scores()
# to machine precision, so the ratio measures and grf's difference share one
# set of scores.
#' @keywords internal
#' @noRd
.hte_arm_scores <- function(fit) {
  tau <- as.numeric(stats::predict(fit)$predictions)
  e   <- fit$W.hat
  wc  <- fit$W.orig - e
  r   <- if (inherits(fit, "causal_survival_forest"))
    fit[["_psi"]]$numerator / wc - wc * tau
  else fit$Y.orig - fit$Y.hat - wc * tau
  list(tau = tau,
       g1  = fit$Y.hat + (1 - e) * tau + fit$W.orig / e * r,
       g0  = fit$Y.hat - e * tau + (1 - fit$W.orig) / (1 - e) * r)
}

# `diff` is grf's own average_treatment_effect(), which also honours clusters
# and sample weights. `ratio` and `OR` average the arm scores and take a Wald
# interval on the log scale, the standard error coming from the influence
# function (delta method). `event_risk` turns S(t) into the event risk
# 1 - S(t) first, so a ratio below 1 favours treatment as a hazard ratio does.
# `beyond` flags the survival patients still followed at `time`; where an arm
# has none, grf extrapolates to a near-null effect, so the rows are NA.
#' @keywords internal
#' @noRd
.hte_estimate <- function(fit, s, idx, grid, event_risk, z, label,
                          beyond = NULL) {
  w <- fit$W.orig[idx]
  one <- function(estimand, measure) {
    if (measure == "diff") {
      a  <- grf::average_treatment_effect(
        fit, target.sample = .HTE_ESTIMANDS[[estimand]], subset = which(idx))
      est <- a[["estimate"]]
      se  <- a[["std.err"]]
      return(c(est, se, est - z * se, est + z * se,
               2 * stats::pnorm(-abs(est / se))))
    }
    g1 <- s$g1[idx]
    g0 <- s$g0[idx]
    if (event_risk) {
      g1 <- 1 - g1
      g0 <- 1 - g0
    }
    a <- mean(g1)
    b <- mean(g0)
    if (!(a > 0 && b > 0 && (measure == "ratio" || (a < 1 && b < 1)))) {
      warning(sprintf("%s: `%s` is undefined because an arm mean is out of range; set to NA.",
                      label, measure), call. = FALSE)
      return(rep(NA_real_, 5L))
    }
    if (measure == "ratio") {
      th  <- log(a / b)
      inf <- (g1 - a) / a - (g0 - b) / b
    } else {
      th  <- stats::qlogis(a) - stats::qlogis(b)
      inf <- (g1 - a) / (a * (1 - a)) - (g0 - b) / (b * (1 - b))
    }
    se <- stats::sd(inf) / sqrt(length(inf))
    c(exp(th), se, exp(th - z * se), exp(th + z * se),
      2 * stats::pnorm(-abs(th / se)))
  }
  vals <- if (sum(w == 1) < 2L || sum(w == 0) < 2L) {
    warning(sprintf("%s: fewer than two units in one arm; estimates set to NA.",
                    label), call. = FALSE)
    matrix(NA_real_, nrow(grid), 5L)
  } else if (!is.null(beyond) &&
             (!any(beyond[idx] & w == 1) || !any(beyond[idx] & w == 0))) {
    warning(sprintf("%s: no patient in one arm is followed past `time`; estimates set to NA.",
                    label), call. = FALSE)
    matrix(NA_real_, nrow(grid), 5L)
  } else {
    t(vapply(seq_len(nrow(grid)),
             function(i) one(grid$estimand[i], grid$measure[i]), numeric(5L)))
  }
  data.frame(estimand = grid$estimand, measure = grid$measure,
             estimate = vals[, 1L], std.error = vals[, 2L],
             conf.low = vals[, 3L], conf.high = vals[, 4L],
             p.value = vals[, 5L], stringsAsFactors = FALSE)
}


# Numeric covariates with more than 5 distinct values are continuous, the
# threshold get_hte() already applies to `sub_var`.
#' @keywords internal
#' @noRd
.hte_is_num <- function(x) is.numeric(x) && length(unique(x[!is.na(x)])) > 5L

# Doubly robust CATE along one covariate, from the AIPW scores in
# `data$.dr_score`. A categorical covariate gets the score mean per level --
# grf::average_treatment_effect(subset = level) -- and a K - 1 df Wald test
# that the levels are equal, the subgroup p_inter of get_hte(). A continuous
# one gets a natural spline with HC3 errors and a joint Wald test of the
# spline terms. Either way `p_het` tests whether the CATE varies with it.
#' @keywords internal
#' @noRd
.hte_dr_var <- function(data, v, w, z, spline_df) {
  x <- data[[v]]
  s <- data$.dr_score
  if (.hte_is_num(x)) {
    fit <- stats::lm(s ~ splines::ns(x, df = spline_df))
    V   <- sandwich::vcovHC(fit, type = "HC3")
    b   <- stats::coef(fit)[-1L]
    g   <- data.frame(x = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE),
                              length.out = 100L))
    M   <- stats::model.matrix(stats::delete.response(stats::terms(fit)), g)
    est <- drop(M %*% stats::coef(fit))
    se  <- sqrt(rowSums((M %*% V) * M))
    return(list(
      type  = "continuous", df = length(b),
      p_het = stats::pchisq(drop(b %*% solve(V[-1L, -1L], b)), length(b),
                            lower.tail = FALSE),
      curve = data.frame(x = g$x, estimate = est, conf.low = est - z * se,
                         conf.high = est + z * se)))
  }
  g  <- droplevels(as.factor(x))
  lv <- do.call(rbind, lapply(levels(g), function(l) {
    i  <- which(g == l)
    nt <- sum(w[i])
    # the same two-per-arm floor as the get_hte() subgroup rows
    ok  <- nt >= 2 && length(i) - nt >= 2
    est <- if (ok) mean(s[i]) else NA_real_
    se  <- if (ok) stats::sd(s[i]) / sqrt(length(i)) else NA_real_
    data.frame(level = l, n = length(i), n_treat = nt, estimate = est,
               std.error = se, conf.low = est - z * se,
               conf.high = est + z * se, stringsAsFactors = FALSE)
  }))
  ok <- is.finite(lv$std.error) & lv$std.error > 0
  if (sum(ok) < 2L)
    return(list(type = "categorical", df = NA_integer_, p_het = NA_real_,
                levels = lv))
  wt <- 1 / lv$std.error[ok]^2
  th <- lv$estimate[ok]
  list(type = "categorical", df = sum(ok) - 1L,
       p_het = stats::pchisq(sum(wt * (th - sum(wt * th) / sum(wt))^2),
                             sum(ok) - 1L, lower.tail = FALSE),
       levels = lv)
}


# ---- L1 public entry point -------------------------------------------------

#' Heterogeneous treatment effects with causal forests
#'
#' Fits one generalized random forest for a binary exposure and returns, in
#' one object, the average treatment effect on an absolute or relative scale,
#' subgroup averages of the conditional average treatment effect (CATE), the
#' per-row CATE and doubly robust scores for drawing CATE curves, and the
#' fitted forest. Arguments follow [RegR::get_eff()] where the two overlap.
#'
#' @param data A data frame holding every column named below.
#' @param cat_var Length-1 character. The binary exposure column: `0`/`1`,
#'   logical, or a two-level factor or character column whose **second** level
#'   is the treated arm, as in [get_PSW()]. Unlike [RegR::get_eff()], exactly
#'   one exposure is analysed.
#' @param sub_var Character vector of categorical subgroup columns, or `NULL`
#'   (default). Each level gets one `$subgroup` row per estimand and measure.
#'   Numeric columns with more than 5 distinct values are rejected; cut them
#'   into groups first. Columns missing from `adj_var` are added to the forest
#'   covariates with a message, because a subgroup estimate is only guaranteed
#'   for variables the forest conditions on. A patient whose value is missing
#'   is left out of that column's rows.
#' @param adj_var Character vector of covariates the forest conditions on, or
#'   `NULL`. Factor, character and logical columns enter as one indicator
#'   column per level, with no reference level dropped, so a tree can split
#'   any single level off from the rest. Missing values are kept and left to
#'   grf, which splits on missingness; a missing factor value leaves all of
#'   that factor's indicator columns missing. Other column types, such as
#'   dates, are rejected: convert them to numbers first.
#' @param surv Outcome selector, following [RegR::get_eff()]:
#'   \itemize{
#'     \item `TRUE` (default): survival outcome in the fixed columns `time`
#'       (follow-up) and `DSS` (event, 0/1), fitted with
#'       [grf::causal_survival_forest()].
#'     \item A single column name: binary (coded 0/1) or continuous outcome,
#'       fitted with [grf::causal_forest()].
#'     \item `FALSE` (competing risks in `get_eff()`) is rejected: grf has no
#'       competing-risk forest.
#'   }
#' @param method Backend. Only `"grf"` is available.
#' @param estimand Character vector, any of `"ATE"` (default), `"ATT"`,
#'   `"ATC"`, `"ATO"`, passed to grf as `target.sample` `"all"`, `"treated"`,
#'   `"control"` and `"overlap"`. Survival outcomes support `"ATE"` only.
#' @param measure Character vector, any of `"diff"` (default), `"ratio"`,
#'   `"OR"`; see the Effect measures section. `"ratio"` and `"OR"` are
#'   computed for the ATE only.
#' @param conf_level Confidence level of the Wald intervals. Default `0.95`.
#' @param time Single positive number: the time point \eqn{t} of a survival
#'   outcome, passed to grf as `horizon`. It defines the estimand and is fixed
#'   when the forest is grown, so another time point needs another call.
#'   Default `120`, as in [RegR::get_eff()]. Only accepted with `surv = TRUE`.
#'   Both arms need patients still followed beyond `time` (up to it for
#'   RMST): with none in an arm the call stops, because grf would return a
#'   near-null effect without warning, and with fewer than 10 it warns.
#' @param grf_args Named list forwarded to [grf::causal_forest()] or
#'   [grf::causal_survival_forest()], for example `num.trees`, `seed`,
#'   `tune.parameters` or a known propensity `W.hat` (as in a trial). `X`,
#'   `Y`, `W`, `D` and `horizon` are set by `get_hte()`. For survival outcomes
#'   `target` defaults to `"survival.probability"` (grf's own default is
#'   `"RMST"`); `target = "RMST"` switches to restricted mean survival time up
#'   to `time`. Per-row fields (`W.hat`, `Y.hat`, `sample.weights`,
#'   `clusters`) may be given for every row of `data` -- rows dropped for a
#'   missing `cat_var` or outcome are dropped from them too -- or for the
#'   rows analysed only. `clusters` and `sample.weights` are only supported
#'   with `measure = "diff"`.
#' @param verbose Logical. `TRUE` reports how many rows were dropped for a
#'   missing `cat_var` or outcome. Default `FALSE`.
#'
#' @section Effect measures:
#' Every measure compares the mean outcome of the two arms:
#'
#' | Outcome | `"diff"` | `"ratio"` | `"OR"` |
#' |:--|:--|:--|:--|
#' | continuous | mean difference | ratio of means (needs positive means) | -- |
#' | binary (0/1) | risk difference | risk ratio | odds ratio |
#' | survival, `target = "survival.probability"` | \eqn{S_1(t) - S_0(t)} | event-risk ratio \eqn{(1 - S_1(t)) / (1 - S_0(t))} | odds ratio of the event by \eqn{t} |
#' | survival, `target = "RMST"` | RMST difference | RMST ratio | -- |
#'
#' `"diff"` is [grf::average_treatment_effect()] itself. `"ratio"` and `"OR"`
#' come from arm-specific AIPW scores backed out of the forest, whose
#' difference reproduces [grf::get_scores()] exactly. Their intervals and
#' p-values are Wald tests on the log scale with an influence-function
#' (delta-method) standard error, and `std.error` is reported on that log
#' scale. For a survival probability the relative measures compare the event
#' risk \eqn{1 - S(t)}, so a value below 1 favours treatment, as a hazard
#' ratio would, while `"diff"` stays on the survival scale of
#' [RegR::get_eff()]. A hazard ratio itself is not available: a causal forest
#' estimates contrasts of mean outcomes, which a hazard ratio is not.
#'
#' Every measure needs overlap: if the propensity (`W.hat`, estimated or
#' supplied) is exactly 0 or 1 for any patient, the doubly robust scores are
#' undefined and the call stops. Propensities at or beyond 0.05 and 0.95 give
#' one warning with their range, which grf itself would repeat for every
#' estimand and subgroup.
#'
#' Requested combinations that are not available -- `"OR"` for a continuous
#' outcome or RMST, a relative measure for anything but the ATE, a survival
#' estimand other than the ATE -- are skipped with a message; the call stops
#' only when none remains.
#'
#' @section Subgroups:
#' `estimate` is the doubly robust effect within the subgroup:
#' [grf::average_treatment_effect()] with `subset` for `"diff"`, and the arm
#' scores averaged over the subgroup for `"ratio"` and `"OR"`. `cate_mean` is
#' the estimand-weighted mean of the out-of-bag CATE in the subgroup, given
#' for `"diff"` rows only and only as a description: forest estimates are
#' shrunk towards the overall mean, so it carries no interval. `p_inter` tests
#' whether the subgroup estimates of one `sub_var` are equal (Wald
#' chi-square with K - 1 degrees of freedom, on the log scale for `"ratio"`
#' and `"OR"`). A subgroup in which an arm has no patient followed beyond
#' `time` gets `NA`, like one with fewer than two patients in an arm.
#'
#' @section CATE curves:
#' `$data` supports two univariate curves over a covariate `x`. Regressing
#' `.dr_score` on a smooth function of `x`, for example
#' `lm(.dr_score ~ splines::ns(x, 4), data = res$data)` with heteroskedasticity
#' robust standard errors, estimates \eqn{E[\tau(X) \mid x]} with valid
#' pointwise intervals. Plotting `.cate` against `x` shows the forest's own
#' out-of-bag estimates, whose smoother bands are not confidence intervals.
#'
#' @return An object of class `hte_res`: a list of
#'   \describe{
#'     \item{`stats`}{Tibble with one row per estimand and measure: `method`,
#'       `estimand`, `measure`, `estimate`, `std.error`, `conf.low`,
#'       `conf.high`, `p.value`, `n`, `n_treat`.}
#'     \item{`subgroup`}{Tibble with one row per `sub_var` level, estimand and
#'       measure: `sub_var`, `level`, `estimand`, `measure`, `n`, `n_treat`,
#'       `estimate`, `std.error`, `conf.low`, `conf.high`, `p.value`,
#'       `cate_mean`, `p_inter`. `NULL` without `sub_var`.}
#'     \item{`importance`}{Tibble with one row per covariate, sorted by
#'       `importance`: `variable`, `importance` ([grf::variable_importance()]
#'       summed over the covariate's design columns, so the column sums to 1)
#'       `n_col` (how many design columns the covariate occupies), `df` and
#'       `p_het`. The first two columns are what `MLR::plt_bar_per()` reads,
#'       so `plt_bar_per(res$importance)` plots it directly. `importance` is a
#'       depth-weighted split frequency, not a test: covariates with more
#'       columns or more distinct values score higher even without any effect
#'       modification. `p_het` is the test: a Wald test, on the AIPW scores,
#'       that the CATE does not vary with the covariate -- equal level means
#'       for a categorical covariate (`df` = levels - 1, the `p_inter` of a
#'       `sub_var`), and zero natural-spline terms (`df` = 3, HC3 errors) for
#'       a numeric one with more than 5 distinct values. Levels with fewer
#'       than two patients in either arm are left out. [plt_hte_dep()] draws
#'       the same estimates.}
#'     \item{`data`}{The rows analysed plus `.cate`, the out-of-bag
#'       CATE, and `.dr_score`, the AIPW score (equal to
#'       [grf::get_scores()]). Both are on the `"diff"` scale whatever
#'       `measure` is.}
#'     \item{`fit`}{The grf forest.}
#'   }
#'   Analysis metadata is attached as `attr(x, "analysis")`, including the
#'   covariates actually used, the treated level, the propensity range
#'   (`ps_range`) and the seed of the forest.
#'   Rows missing `cat_var` or the outcome are dropped. Missing covariates are
#'   left to grf, and each subgroup row, `p_het` and [plt_hte_dep()] panel
#'   uses the patients whose value of that covariate is observed.
#'
#' @references
#' Wager S, Athey S (2018). Estimation and inference of heterogeneous
#' treatment effects using random forests. \emph{Journal of the American
#' Statistical Association} 113(523):1228-1242.
#'
#' Athey S, Tibshirani J, Wager S (2019). Generalized random forests.
#' \emph{The Annals of Statistics} 47(2):1148-1178.
#'
#' Cui Y, Kosorok MR, Sverdrup E, Wager S, Zhu R (2023). Estimating
#' heterogeneous treatment effects with right-censored data via causal
#' survival forests. \emph{Journal of the Royal Statistical Society Series B}
#' 85(2):179-211.
#'
#' @seealso [RegR::get_eff()] for regression-based effects and hazard ratios;
#'   [get_PSW()] for the propensity score and weighting diagnostics.
#'
#' @examplesIf requireNamespace("grf", quietly = TRUE)
#' \donttest{
#' set.seed(20260923)
#' n <- 800
#' d <- data.frame(age = rnorm(n, 60, 10), x2 = rnorm(n),
#'                 sex = factor(sample(c("F", "M"), n, replace = TRUE)))
#' d$z <- rbinom(n, 1, plogis(0.3 * d$x2))
#' d$y <- rbinom(n, 1, plogis(-1 + 0.5 * d$x2 + d$z * (0.3 + 0.8 * (d$sex == "M"))))
#'
#' # Binary outcome: risk difference, risk ratio and odds ratio, by sex
#' res <- get_hte(d, cat_var = "z", sub_var = "sex", adj_var = c("age", "x2"),
#'                surv = "y", measure = c("diff", "ratio", "OR"),
#'                grf_args = list(num.trees = 500, seed = 1))
#' res
#'
#' # CATE curve over age: AIPW scores smoothed (valid intervals) and the
#' # out-of-bag CATE (descriptive)
#' curve_fit <- lm(.dr_score ~ poly(age, 3), data = res$data)
#' plot(res$data$age, res$data$.cate, xlab = "age", ylab = "CATE")
#'
#' # Survival: S(t) difference and event-risk ratio at t = 60
#' ev   <- rexp(n, 0.02 * exp(0.3 * d$x2 - 0.5 * d$z))
#' cens <- pmin(rexp(n, 0.01), 120)
#' d$time <- pmin(ev, cens)
#' d$DSS  <- as.integer(ev <= cens)
#' get_hte(d, cat_var = "z", adj_var = c("age", "x2", "sex"), surv = TRUE,
#'         time = 60, measure = c("diff", "ratio"),
#'         grf_args = list(num.trees = 500, seed = 1))
#' }
#'
#' @export
get_hte <- function(data,
                    cat_var,
                    sub_var    = NULL,
                    adj_var    = NULL,
                    surv       = TRUE,
                    method     = "grf",
                    estimand   = "ATE",
                    measure    = "diff",
                    conf_level = 0.95,
                    time       = 120,
                    grf_args   = list(),
                    verbose    = FALSE) {

  method <- match.arg(method, "grf")
  if (!requireNamespace("grf", quietly = TRUE))
    stop("Package 'grf' is required for method = \"grf\".", call. = FALSE)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  if (!is.character(estimand) || !length(estimand))
    stop("`estimand` must be a character vector.", call. = FALSE)
  est_all  <- names(.HTE_ESTIMANDS)
  estimand <- est_all[est_all %in%
                        match.arg(toupper(estimand), est_all, several.ok = TRUE)]
  measure  <- .HTE_MEASURES[.HTE_MEASURES %in%
                              match.arg(measure, .HTE_MEASURES, several.ok = TRUE)]

  if (!is.data.frame(data) || !nrow(data))
    stop("`data` must be a non-empty data frame.", call. = FALSE)
  cat_var <- .sens_check_col(cat_var, data, "cat_var", n = 1L)
  sub_var <- setdiff(.sens_check_col(sub_var, data, "sub_var"), cat_var)
  adj_var <- setdiff(.sens_check_col(adj_var, data, "adj_var"), cat_var)
  if (!length(sub_var)) sub_var <- NULL
  if (!length(adj_var)) adj_var <- NULL

  # ---- Outcome: surv = TRUE reads the fixed get_eff() columns -------------
  if (isFALSE(surv))
    stop("`surv = FALSE` (competing risks) is not supported: grf has no competing-risk forest.",
         call. = FALSE)
  is_surv <- isTRUE(surv)
  if (is_surv) {
    outcome <- c("time", "DSS")
    if (!all(outcome %in% names(data)))
      stop("`surv = TRUE` requires the columns `time` and `DSS`.", call. = FALSE)
    if (!is.numeric(time) || length(time) != 1L || is.na(time) || time <= 0)
      stop("`time` must be a single positive number.", call. = FALSE)
  } else {
    if (!is.character(surv) || length(surv) != 1L || is.na(surv))
      stop("`surv` must be TRUE (columns `time` / `DSS`) or a single outcome column name.",
           call. = FALSE)
    outcome <- .sens_check_col(surv, data, "surv", n = 1L)
    if (!missing(time))
      stop("`time` only applies to survival outcomes (`surv = TRUE`).",
           call. = FALSE)
  }

  # ---- Covariates: subgroup variables must enter the forest ---------------
  covars <- unique(c(adj_var, sub_var))
  if (!length(covars))
    stop("Supply at least one covariate through `adj_var` or `sub_var`.",
         call. = FALSE)
  if (any(covars %in% outcome))
    stop(sprintf("`adj_var` / `sub_var` must not include the outcome column(s) %s.",
                 paste0("`", intersect(covars, outcome), "`", collapse = ", ")),
         call. = FALSE)
  # Anything else -- a Date is not is.numeric() -- would be one-hot encoded
  # into one column per distinct value.
  odd <- covars[!vapply(data[covars], function(x) is.numeric(x) ||
                          is.factor(x) || is.character(x) || is.logical(x),
                        logical(1L))]
  if (length(odd))
    stop(sprintf("Covariate column(s) %s must be numeric, factor, character or logical; convert a date, for example, to years since a reference date.",
                 paste0("`", odd, "`", collapse = ", ")), call. = FALSE)
  for (v in sub_var) {
    x <- data[[v]]
    if (is.numeric(x) && length(unique(x[!is.na(x)])) > 5L)
      stop(sprintf("`sub_var` column `%s` is numeric with more than 5 distinct values; cut it into groups first.",
                   v), call. = FALSE)
  }
  added <- setdiff(sub_var, adj_var)
  if (length(added))
    cli::cli_inform(c("i" = paste(
      "Added {.field {added}} to the forest covariates: subgroup estimates",
      "are only guaranteed for variables the forest conditions on.")))

  # Only the exposure and the outcome must be complete. grf splits on missing
  # covariates itself, so dropping those rows -- or letting a sub_var shrink
  # the overall sample -- would only throw patients away.
  keep <- stats::complete.cases(data[c(cat_var, outcome)])
  data <- .sens_complete(data, c(cat_var, outcome), verbose)
  empty <- covars[vapply(data[covars], function(x) all(is.na(x)), logical(1L))]
  if (length(empty))
    stop(sprintf("Covariate column(s) %s have no observed value.",
                 paste0("`", empty, "`", collapse = ", ")), call. = FALSE)
  tz <- .psw_treat(data[[cat_var]], cat_var, arg = "cat_var")
  W  <- tz$z
  if (sum(W == 1L) < 2L || sum(W == 0L) < 2L)
    stop("Both arms of `cat_var` need at least two complete rows.",
         call. = FALSE)

  if (is_surv) {
    Y <- data[["time"]]
    D <- data[["DSS"]]
    if (is.logical(D)) D <- as.integer(D)
    if (!is.numeric(Y) || !is.numeric(D) || !all(D %in% c(0, 1)))
      stop("`surv = TRUE` needs a numeric `time` column and a `DSS` column coded 0/1.",
           call. = FALSE)
    type <- "survival"
  } else {
    Y <- data[[outcome]]
    if (is.logical(Y)) Y <- as.integer(Y)
    if (!is.numeric(Y))
      stop(sprintf("Outcome column `%s` must be numeric or logical; code a binary outcome as 0/1.",
                   outcome), call. = FALSE)
    u <- unique(Y)
    if (length(u) == 2L && !all(u %in% c(0, 1)))
      stop(sprintf("Outcome column `%s` has two values; code a binary outcome as 0/1.",
                   outcome), call. = FALSE)
    type <- if (all(u %in% c(0, 1))) "binary" else "continuous"
  }

  # ---- Backend arguments ---------------------------------------------------
  fun <- if (is_surv) grf::causal_survival_forest else grf::causal_forest
  # Name the argument that sets a reserved field, rather than calling it
  # "unknown" and listing every other grf argument.
  fixed <- c(X = "`adj_var` / `sub_var`", Y = "`surv`", W = "`cat_var`",
             D = "`surv = TRUE` (column `DSS`)", horizon = "`time`")
  hit <- intersect(names(grf_args), names(fixed))
  if (length(hit))
    stop(sprintf("`grf_args` cannot set %s: get_hte() sets %s.",
                 paste0("`", hit, "`", collapse = ", "),
                 paste(sprintf("`%s` from %s", hit, fixed[hit]),
                       collapse = " and ")), call. = FALSE)
  grf_args <- .merge_named_arg(
    grf_args, list(), "grf_args",
    allowed = setdiff(names(formals(fun)), c("X", "Y", "W", "D", "horizon")))
  # A per-row field given for every row of `data` follows the rows kept above;
  # one already matching the rows analysed passes unchanged.
  for (f in intersect(c("W.hat", "Y.hat", "sample.weights", "clusters"),
                      names(grf_args)))
    if (!all(keep) && length(grf_args[[f]]) == length(keep))
      grf_args[[f]] <- grf_args[[f]][keep]
  target <- NULL
  if (is_surv) {
    target <- match.arg(
      if (is.null(grf_args$target)) "survival.probability" else grf_args$target,
      c("survival.probability", "RMST"))
    grf_args$target <- target
  }

  # ---- Follow-up past the time point ----------------------------------------
  # grf identifies S(t) or RMST(t) only from patients still followed at `time`
  # (its source uses Y > horizon for a survival probability, Y >= horizon for
  # RMST). With none left in an arm it returns a near-null effect without a
  # warning -- measured 0.018 against a true 0.154 -- so stop here instead.
  beyond <- NULL
  if (is_surv) {
    beyond  <- if (target == "RMST") Y >= time else Y > time
    past    <- if (target == "RMST") "up to" else "beyond"
    arms    <- c(1L, 0L)
    arm_lab <- c(tz$treated,
                 setdiff(levels(factor(data[[cat_var]])), tz$treated)[1L])
    n_after <- vapply(arms, function(a) sum(beyond & W == a), integer(1L))
    if (any(n_after == 0L)) {
      k <- which(n_after == 0L)[1L]
      stop(sprintf("No patient with `%s` = %s is followed %s `time` = %s (longest follow-up %s), so the effect at that time is not identified; choose a smaller `time`.",
                   cat_var, arm_lab[k], past, format(time),
                   format(max(Y[W == arms[k]]))), call. = FALSE)
    }
    if (any(n_after < 10L))
      warning(sprintf("Few patients are followed %s `time` = %s: %s. The estimate at that time rests on them.",
                      past, format(time),
                      paste(sprintf("%d with `%s` = %s", n_after, cat_var,
                                    arm_lab), collapse = ", ")),
              call. = FALSE)
  }

  # ---- Requested estimand x measure grid -----------------------------------
  grid <- expand.grid(estimand = estimand, measure = measure,
                      stringsAsFactors = FALSE)
  why <- vapply(seq_len(nrow(grid)), function(i) {
    e <- grid$estimand[i]
    m <- grid$measure[i]
    if (is_surv && e != "ATE") "grf's survival forest estimates the ATE only"
    else if (m != "diff" && e != "ATE") "ratio and OR are available for the ATE only"
    else if (m == "OR" && (type == "continuous" || identical(target, "RMST")))
      "OR needs an outcome probability"
    else NA_character_
  }, character(1L))
  skipped <- sprintf("%s x %s: %s", grid$estimand, grid$measure, why)[!is.na(why)]
  if (all(!is.na(why)))
    stop(paste0("No requested estimand / measure combination is available:\n",
                paste0("* ", skipped, collapse = "\n")), call. = FALSE)
  if (length(skipped))
    cli::cli_inform(c("i" = "Skipped {length(skipped)} combination{?s}:",
                      stats::setNames(skipped, rep("*", length(skipped)))))
  grid <- grid[is.na(why), , drop = FALSE]
  if (any(grid$measure != "diff") &&
      any(c("clusters", "sample.weights") %in% names(grf_args)))
    stop("`ratio` and `OR` do not support `grf_args$clusters` or `grf_args$sample.weights`; use measure = \"diff\".",
         call. = FALSE)

  # ---- Fit and estimate ----------------------------------------------------
  X   <- .sens_model_matrix(data, covars, one_hot = TRUE)
  fit <- do.call(fun, c(list(X = X, Y = Y, W = W),
                        if (is_surv) list(D = D, horizon = time),
                        grf_args))
  # A propensity of exactly 0 or 1 -- a regression forest reaches it when a
  # covariate region holds one arm only -- divides the AIPW score by zero:
  # grf's own estimate turns NaN and the ratio scores fail. This also checks
  # a W.hat passed through grf_args.
  bad <- !is.finite(fit$W.hat) | fit$W.hat <= 0 | fit$W.hat >= 1
  if (any(bad))
    stop(sprintf("The estimated propensity of `%s` is exactly 0 or 1 for %d patient%s, so the doubly robust estimates are undefined. Restrict the data to the region of overlap (for example the patients get_PSW(trim_args = list(method = \"cr\")) keeps), drop covariates that fully determine `%s`, or pass a bounded `W.hat` in `grf_args`.",
                 cat_var, sum(bad), if (sum(bad) == 1L) "" else "s",
                 cat_var), call. = FALSE)
  s  <- .hte_arm_scores(fit)
  z  <- stats::qnorm(1 - (1 - conf_level) / 2)
  event_risk <- identical(target, "survival.probability")

  # grf repeats its overlap warning in every average_treatment_effect() call,
  # once per estimand and subgroup (18 times in one call, measured), so the
  # estimates below keep it quiet and one warning covers the whole sample.
  ps_rng    <- range(fit$W.hat)
  ps_warned <- FALSE
  quiet_ps  <- function(expr) withCallingHandlers(expr, warning = function(w) {
    if (startsWith(conditionMessage(w), "Estimated treatment propensities")) {
      ps_warned <<- TRUE
      invokeRestart("muffleWarning")
    }
  })

  overall <- quiet_ps(.hte_estimate(fit, s, rep(TRUE, length(W)), grid,
                                    event_risk, z, "Overall", beyond))
  stats_tbl <- tibble::as_tibble(data.frame(
    method = method, overall, n = length(W), n_treat = sum(W),
    stringsAsFactors = FALSE))

  sub_tbl <- NULL
  if (length(sub_var)) {
    # Plug-in weights matching each estimand, for the descriptive cate_mean.
    h <- list(ATE = rep(1, length(W)), ATT = W, ATC = 1 - W,
              ATO = fit$W.hat * (1 - fit$W.hat))
    sub_tbl <- quiet_ps(do.call(rbind, lapply(sub_var, function(v) {
      g <- droplevels(as.factor(data[[v]]))
      rows <- do.call(rbind, lapply(levels(g), function(lv) {
        idx <- g %in% lv                 # FALSE where `v` is missing
        est <- .hte_estimate(fit, s, idx, grid, event_risk, z,
                             sprintf("%s = %s", v, lv), beyond)
        cate_mean <- vapply(seq_len(nrow(est)), function(i) {
          if (est$measure[i] != "diff") return(NA_real_)
          stats::weighted.mean(s$tau[idx], h[[est$estimand[i]]][idx])
        }, numeric(1L))
        data.frame(sub_var = v, level = lv, est[c("estimand", "measure")],
                   n = sum(idx), n_treat = sum(W[idx]),
                   est[c("estimate", "std.error", "conf.low", "conf.high",
                         "p.value")],
                   cate_mean = cate_mean, stringsAsFactors = FALSE)
      }))
      # Equal subgroup effects: Wald chi-square on the analysis scale.
      rows$p_inter <- NA_real_
      key <- paste(rows$estimand, rows$measure)
      for (k in unique(key)) {
        j  <- which(key == k)
        th <- if (rows$measure[j[1L]] == "diff") rows$estimate[j]
              else log(rows$estimate[j])
        se <- rows$std.error[j]
        ok <- is.finite(th) & is.finite(se) & se > 0
        if (sum(ok) >= 2L) {
          wt <- 1 / se[ok]^2
          q  <- sum(wt * (th[ok] - sum(wt * th[ok]) / sum(wt))^2)
          rows$p_inter[j] <- stats::pchisq(q, df = sum(ok) - 1L,
                                           lower.tail = FALSE)
        }
      }
      rows
    })))
    rownames(sub_tbl) <- NULL
    sub_tbl <- tibble::as_tibble(sub_tbl)
  }
  if (ps_warned)
    warning(sprintf("Estimated propensities of `%s` range from %.3f to %.3f; grf flags values at or beyond 0.05 and 0.95, where effects are poorly identified. %s",
                    cat_var, ps_rng[1L], ps_rng[2L],
                    if (is_surv) "Trimming to the region of overlap, as get_PSW(trim_args = list(method = \"cr\")) does, is more stable."
                    else "estimand = \"ATO\", or trimming to the region of overlap as get_PSW(trim_args = list(method = \"cr\")) does, is more stable."),
            call. = FALSE)

  data$.cate     <- s$tau
  data$.dr_score <- s$g1 - s$g0

  # grf scores one importance per design column; summing a covariate's columns
  # gives one row per covariate, laid out for MLR::plt_bar_per() (a categorical
  # first column, then a numeric one). p_het is the doubly robust test that
  # plt_hte_dep() prints in its strips.
  src <- covars[attr(X, "assign")]
  vi  <- as.numeric(grf::variable_importance(fit))
  dr  <- lapply(covars, function(v) .hte_dr_var(data, v, W, z, .HTE_SPLINE_DF))
  imp_tbl <- tibble::tibble(
    variable   = covars,
    importance = vapply(covars, function(v) sum(vi[src == v]), numeric(1L),
                        USE.NAMES = FALSE),
    n_col      = tabulate(attr(X, "assign"), nbins = length(covars)),
    df         = vapply(dr, function(r) r$df, integer(1L)),
    p_het      = vapply(dr, function(r) r$p_het, numeric(1L)))
  imp_tbl <- imp_tbl[order(imp_tbl$importance, decreasing = TRUE), ]

  structure(
    list(stats = stats_tbl, subgroup = sub_tbl, importance = imp_tbl,
         data = data, fit = fit),
    class = c("hte_res", "list"),
    analysis = list(
      method = method, backend = "grf",
      backend_version = as.character(utils::packageVersion("grf")),
      forest = class(fit)[1L], outcome_type = type,
      cat_var = cat_var, treated = tz$treated, surv = surv, outcome = outcome,
      sub_var = sub_var, adj_var = adj_var, covariates = covars,
      estimand = estimand, measure = measure, target = target,
      time = if (is_surv) time else NULL, conf_level = conf_level,
      n = length(W), n_treat = sum(W), ps_range = ps_rng,
      seed = fit[["seed"]], call = match.call()))
}


# ---- L3 print --------------------------------------------------------------

#' @export
#' @noRd
print.hte_res <- function(x, ...) {
  a <- attr(x, "analysis")
  cat(sprintf("<hte_res> %s (grf %s), n = %d (treated = %d), conf_level = %s\n",
              a$forest, a$backend_version, a$n, a$n_treat,
              format(a$conf_level)))
  cat(sprintf("  cat_var = %s (treated: %s), outcome = %s%s\n",
              a$cat_var, a$treated, paste(a$outcome, collapse = " / "),
              if (is.null(a$time)) "" else
                sprintf(", target = %s at time = %s", a$target, format(a$time))))
  cat(sprintf("  propensity range %.3f to %.3f\n", a$ps_range[1L],
              a$ps_range[2L]))
  cat("\n")
  print(x$stats)
  if (!is.null(x$subgroup)) {
    cat("\nSubgroups:\n")
    print(x$subgroup)
  }
  cat("\n")
  if (identical(a$target, "survival.probability") &&
      any(x$stats$measure != "diff"))
    cat("# ratio / OR compare the event risk 1 - S(t); diff is S1(t) - S0(t).\n")
  cat("# $data: .cate = out-of-bag CATE, .dr_score = AIPW score (diff scale).\n")
  cat("# $importance: grf split frequency by covariate (not a test) + p_het; plt_bar_per()-ready.\n")
  cat("# plt_hte_dep(x) draws the CATE against each covariate.\n")
  invisible(x)
}
