# =============================================================================
# hte-unihtee.R -- effect-modifier screening by univariate TEM-VIP projections
# =============================================================================
#
# Architecture:
#
#   L1  get_hte_unihtee()  validate, obtain per-row AIPW arm scores from a
#                          causal forest (get_hte()) or from cross-fitted
#                          GLMs, then regress one pseudo-outcome per measure
#                          on each candidate alone
#   L1  plt_hte_unihtee()  bar / volcano of the table, or one panel per
#                          candidate with its projection line
#   L2  .uni_glm_fit()     cross-fitted propensity and per-arm outcome GLMs
#   L2  .uni_project()     the projection on one candidate, shared by the
#                          table and its plotted line
#
# =============================================================================

# Cross-fitted nuisances for method = "glm": a logistic propensity model and
# one outcome GLM per arm, each predicted on the fold it was not fitted on.
# Folds are drawn within each arm so every training set holds both arms; the
# caller sets the seed and restores the RNG.
#' @keywords internal
#' @noRd
.uni_glm_fit <- function(data, covars, W, Y, family, folds, ps_trim) {
  n    <- length(W)
  fold <- integer(n)
  for (a in 0:1) fold[W == a] <- sample(rep_len(seq_len(folds), sum(W == a)))
  rhs <- .sens_quote_names(covars)
  dd  <- data[covars]
  dd$.w <- W
  dd$.y <- Y
  e <- mu1 <- mu0 <- numeric(n)
  for (k in seq_len(folds)) {
    tr <- fold != k
    te <- dd[!tr, , drop = FALSE]
    e[!tr] <- stats::predict(
      stats::glm(stats::reformulate(rhs, ".w"), family = stats::binomial(),
                 data = dd[tr, , drop = FALSE]), te, type = "response")
    mu1[!tr] <- stats::predict(
      stats::glm(stats::reformulate(rhs, ".y"), family = family,
                 data = dd[tr & W == 1L, , drop = FALSE]), te, type = "response")
    mu0[!tr] <- stats::predict(
      stats::glm(stats::reformulate(rhs, ".y"), family = family,
                 data = dd[tr & W == 0L, , drop = FALSE]), te, type = "response")
  }
  list(e = pmin(pmax(e, ps_trim), 1 - ps_trim), mu1 = mu1, mu0 = mu0,
       n_clip = sum(e < ps_trim | e > 1 - ps_trim))
}

# Least-squares projection of a pseudo-outcome `phi` on one candidate `x`
# (numeric, or 0/1 for a binary candidate) on its own scale. The table
# rescales the slope to per SD -- HC3 errors rescale with it -- and the plot
# draws the line, so both rest on the same rows, weights and covariance.
# With `beyond` (survival) a continuous candidate keeps the values both arms
# reach among patients followed past `time`, and a binary one needs both
# levels in both arms there, as .hte_dr_var() does. Returns the fit, its
# covariance, the rows used, the candidate's SD and range, or `why`.
#' @keywords internal
#' @noRd
.uni_project <- function(x, phi, W, wt, cl, beyond, cont) {
  if (!is.null(beyond)) {
    if (cont) {
      x1  <- x[beyond & W == 1 & !is.na(x)]
      x0  <- x[beyond & W == 0 & !is.na(x)]
      lim <- if (length(x1) && length(x0))
        c(max(min(x1), min(x0)), min(max(x1), max(x0))) else c(Inf, -Inf)
      x[!is.na(x) & (x < lim[1L] | x > lim[2L])] <- NA
    } else if (!all(vapply(0:1, function(l)
      any(beyond & W == 1 & x %in% l) && any(beyond & W == 0 & x %in% l),
      logical(1L))))
      return(list(why = "a level has no patient in one arm followed past `time`"))
  }
  i <- which(!is.na(x) & is.finite(phi) & wt > 0)
  enough <- if (cont) length(unique(x[i])) >= 3L else
    all(vapply(0:1, function(l) sum(x[i] == l & W[i] == 1) >= 2 &&
                 sum(x[i] == l & W[i] == 0) >= 2, logical(1L)))
  if (!enough)
    return(list(why = if (cont) "fewer than three distinct values remain"
                else "a level has fewer than two patients in one arm"))
  fit <- stats::lm(p ~ x, weights = w,
                   data = data.frame(p = phi[i], x = x[i], w = wt[i]))
  if (fit$rank < 2L || fit$df.residual <= 0L ||
      (length(cl) && length(unique(cl[i])) < 2L))
    return(list(why = "the regression is rank deficient or has too few clusters"))
  V <- if (length(cl))
    sandwich::vcovCL(fit, cluster = cl[i], type = "HC3") else
    sandwich::vcovHC(fit, type = "HC3")
  list(fit = fit, V = V, i = i, range = range(x[i]),
       sd = if (cont) sqrt(stats::cov.wt(cbind(x[i]), wt[i])$cov[1L, 1L])
            else NA_real_)
}


# ---- L1 public entry point -------------------------------------------------

#' Screen effect modifiers by univariate TEM-VIP projections
#'
#' For each candidate covariate \eqn{W_j}, estimates the treatment effect
#' modifier variable importance parameter (TEM-VIP) of Boileau et al. (2025)
#' -- the parameter of the `unihtee` package -- with a confidence interval
#' and a Benjamini-Hochberg adjusted p-value. It is the least-squares slope of
#' the conditional average treatment effect (CATE) projected on \eqn{W_j}
#' alone,
#' \deqn{\Psi_j = \mathrm{Cov}[\tau(X), W_j] / \mathrm{Var}[W_j],}
#' estimated by regressing cross-fitted AIPW pseudo-outcomes on \eqn{W_j}.
#' The doubly robust scores come from the causal forest of [get_hte()]
#' (`method = "grf"`) or from cross-fitted generalized linear models
#' (`method = "glm"`). The `unihtee` package itself is not used.
#'
#' @param data A data frame holding every column named below.
#' @param cat_var Length-1 character. The binary exposure column, coded as in
#'   [get_hte()]: `0`/`1`, logical, or a two-level factor or character column
#'   whose second level is the treated arm.
#' @param candidate_var Character vector of the candidate modifiers, one table
#'   row each, or `NULL` (default) for every `adj_var` column. A numeric
#'   column with more than two distinct values is a continuous candidate; any
#'   column with exactly two values (numeric, logical, factor or character) is
#'   binary. Categorical columns with three or more levels, constant columns
#'   and `cat_var` itself are rejected.
#' @param adj_var Character vector of covariates the nuisance models condition
#'   on. The models adjust for the union of `adj_var` and `candidate_var`, so
#'   every candidate is adjusted for.
#' @param surv Outcome selector, as in [get_hte()]: `TRUE` (default) for the
#'   survival columns `time` and `DSS` (`method = "grf"` only), or a single
#'   continuous or 0/1 binary outcome column.
#' @param method Where the AIPW scores come from. `"grf"` (default) fits one
#'   [get_hte()] forest and uses its out-of-bag scores; `"glm"` fits a
#'   logistic propensity model and one outcome GLM per arm (Gaussian or
#'   binomial) by `glm_args$folds`-fold cross-fitting, with no forest.
#' @param measure Character vector, any of `"diff"` (default), `"ratio"` and
#'   `"OR"`, on the outcome scales of [get_hte()]: for a survival probability
#'   `"ratio"` and `"OR"` refer to the event risk \eqn{1 - S(t)}. `"OR"` is
#'   skipped, with a message, for a continuous outcome or RMST.
#' @param conf_level Confidence level of the Wald intervals. Default `0.95`.
#' @param time Time point of a survival outcome, as in [get_hte()]. Default
#'   `120`. Only accepted with `surv = TRUE`.
#' @param grf_args Named list forwarded to [get_hte()] and so to grf, for
#'   example `num.trees`. Without `seed`, `seed` below is used.
#' @param glm_args Named list for `method = "glm"`; partial overrides keep the
#'   other defaults.
#'   \describe{
#'     \item{`folds`}{Whole number from 2 to 100, default `5`: cross-fitting
#'       folds, drawn within each arm.}
#'     \item{`ps_trim`}{Number in \eqn{[0, 0.5)}, default `0.01`: propensities
#'       are clipped to \eqn{[ps\_trim, 1 - ps\_trim]}, with a warning that
#'       counts the clipped rows.}
#'   }
#' @param seed Nonnegative whole number, default `123`. It draws the
#'   cross-fitting folds (the global random-number state is restored
#'   afterwards) or seeds the forest when `grf_args` has no `seed`.
#' @param verbose Logical. `TRUE` reports dropped incomplete rows. Default
#'   `FALSE`.
#'
#' @section Estimate and scale:
#' A continuous candidate gets the slope per standard deviation of the
#' candidate (`contrast = "+1 SD"`; weighted by the forest's observation
#' weights when `grf_args` sets `sample.weights`), so candidates are
#' comparable within one table. A binary candidate gets the difference
#' between its second and first level (`contrast`, for example `"M vs F"`;
#' factor levels, otherwise sorted values), which is the same slope per
#' unit of a 0/1 indicator. `"diff"` projects \eqn{\tau(X)} itself. `"ratio"`
#' and `"OR"` project the log risk (or mean) ratio and log odds ratio
#' \eqn{\log \mu_1(X) - \log \mu_0(X)} and \eqn{\mathrm{logit}\,\mu_1(X) -
#' \mathrm{logit}\,\mu_0(X)}, through the pseudo-outcome
#' \eqn{\log \hat\mu_1 - \log \hat\mu_0 + (\Gamma_1 - \hat\mu_1) / \hat\mu_1 -
#' (\Gamma_0 - \hat\mu_0) / \hat\mu_0} (with \eqn{\hat\mu_a(1 - \hat\mu_a)} for
#' the odds), where \eqn{\Gamma_a} are the arm-specific AIPW scores.
#' As in [get_hte()], `estimate` and the interval are then exponentiated --
#' the factor by which the ratio changes per SD or between the two levels --
#' while `std.error` stays on the log scale. Rows whose arm prediction falls
#' outside the range a log ratio or log odds needs are left out of that
#' measure, with a warning. Survival candidates keep only values (continuous)
#' or require levels (binary) that both arms reach among patients followed
#' past `time`, as `p_het` in [get_hte()] does.
#'
#' Standard errors are HC3 heteroskedasticity-robust (cluster-robust with
#' `grf_args$clusters`) errors of the least-squares slope, so for
#' `method = "grf"` and `"diff"` the row of a continuous candidate equals
#' `grf::best_linear_projection(res$fit, A = x / sd(x))` on the [get_hte()]
#' forest. `p.adj` is Benjamini-Hochberg across the candidates of one
#' measure.
#'
#' @section Reading the table:
#' Each \eqn{\Psi_j} is marginal, not partial: a covariate correlated with a
#' true modifier gets a nonzero slope of its own. It is also linear: a
#' U-shaped modification can have a true slope of zero. [get_hte()] tests the
#' same univariate projection without the linearity -- `$importance$p_het`
#' uses a natural spline with 2 degrees of freedom and a joint test of all
#' levels of a categorical covariate -- and on a [get_hte()] result
#' `plt_hte_dep(res, display = "dr", dr_args = list(spline_df = 1))` draws
#' the straight line whose slope, per unit of the covariate, is the
#' `"diff"` estimate here divided by `sd`. [plt_hte_unihtee()] draws the
#' table and these lines, with a spline beside each line to show what the
#' slope leaves out. The table screens candidates; it
#' does not confirm modifiers, and on observational data its ranking also
#' reflects how well `adj_var` captures confounding.
#'
#' Compared with `unihtee::unihtee()`, the estimand and the one-step
#' estimator are the same, but the nuisance models differ (a forest or
#' per-arm GLMs instead of one main-effects GLM, which makes unihtee's
#' plug-in estimate exactly zero), the standard error is the least-squares
#' influence function rather than one built on that plug-in estimate, and
#' continuous outcomes are not rescaled to \eqn{[0, 1]}, so `"ratio"` is a
#' ratio of means.
#'
#' @return A list of
#'   \describe{
#'     \item{`vip`}{Tibble with one row per candidate and measure, sorted by
#'       `p.value` within each measure: `variable`, `type` (`"continuous"` or
#'       `"binary"`), `contrast`, `n` (rows in the regression), `sd` (the
#'       candidate's SD, `NA` for binary), `measure`, `estimate`,
#'       `std.error`, `conf.low`, `conf.high`, `p.value` and `p.adj`.}
#'     \item{`data`}{The rows analysed plus `.ps` (propensity), `.mu1` and
#'       `.mu0` (predicted outcome under each arm, on the survival scale for
#'       a survival outcome) and one `.score_<measure>` pseudo-outcome per
#'       measure. With `method = "grf"` it also holds the [get_hte()] columns
#'       `.cate` and `.dr_score` (equal to `.score_diff`).}
#'   }
#'   Analysis metadata is attached as `attr(x, "analysis")`, including the
#'   covariates adjusted for and the arguments of the backend.
#'
#' @references
#' Boileau P, Leng N, Hejazi NS, van der Laan M, Dudoit S (2025). A
#' nonparametric framework for treatment effect modifier discovery in high
#' dimensions. \emph{Journal of the Royal Statistical Society Series B}
#' 87(1):157-185.
#'
#' Hines O, Diaz-Ordaz K, Vansteelandt S (2022). Variable importance measures
#' for heterogeneous causal effects. arXiv:2204.06030.
#'
#' @seealso [plt_hte_unihtee()] to plot the result; [get_hte()] for the
#'   forest, its spline test `p_het` and subgroup effects; [plt_hte_dep()]
#'   for the doubly robust curves; [get_hte_select()] for benefit-score
#'   variable selection.
#'
#' @examplesIf requireNamespace("sandwich", quietly = TRUE)
#' set.seed(1)
#' n <- 600
#' d <- data.frame(w1 = rnorm(n), w2 = rnorm(n), w3 = rnorm(n),
#'                 sex = factor(sample(c("F", "M"), n, replace = TRUE)))
#' d$a <- rbinom(n, 1, plogis(0.5 * d$w1))
#' d$y <- d$w1 + d$w2 + d$a * (d$w3 + (d$sex == "M")) + rnorm(n)
#'
#' # Forest-free screening: per-arm GLMs, 5-fold cross-fitting
#' get_hte_unihtee(d, cat_var = "a", adj_var = c("w1", "w2", "w3", "sex"),
#'                 surv = "y", method = "glm")$vip
#'
#' \donttest{
#' if (requireNamespace("grf", quietly = TRUE)) {
#'   # Scores from the get_hte() causal forest
#'   res <- get_hte_unihtee(d, cat_var = "a", adj_var = c("w1", "w2", "w3", "sex"),
#'                          surv = "y", grf_args = list(num.trees = 500))
#'   res$vip
#' }
#' }
#'
#' @export
get_hte_unihtee <- function(data,
                            cat_var,
                            candidate_var = NULL,
                            adj_var    = NULL,
                            surv       = TRUE,
                            method     = c("grf", "glm"),
                            measure    = "diff",
                            conf_level = 0.95,
                            time       = 120,
                            grf_args   = list(),
                            glm_args   = list(),
                            seed       = 123,
                            verbose    = FALSE) {

  method  <- match.arg(method)
  measure <- .HTE_MEASURES[.HTE_MEASURES %in%
                             match.arg(measure, .HTE_MEASURES, several.ok = TRUE)]
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  seed <- .hte_select_count(seed, "seed", 0, .Machine$integer.max)
  if (!requireNamespace("sandwich", quietly = TRUE))
    stop("Package 'sandwich' is required for the robust standard errors.",
         call. = FALSE)
  if (!is.data.frame(data) || !nrow(data))
    stop("`data` must be a non-empty data frame.", call. = FALSE)
  cat_var <- .sens_check_col(cat_var, data, "cat_var", n = 1L)
  adj_var <- setdiff(.sens_check_col(adj_var, data, "adj_var"), cat_var)
  candidate_var <- .sens_check_col(candidate_var, data, "candidate_var")
  if (cat_var %in% candidate_var)
    stop("`candidate_var` must not include `cat_var`.", call. = FALSE)
  if (anyDuplicated(candidate_var))
    stop("`candidate_var` must not contain duplicates.", call. = FALSE)
  if (is.null(candidate_var)) candidate_var <- adj_var
  if (!length(candidate_var))
    stop("Supply the candidate modifiers through `candidate_var` or `adj_var`.",
         call. = FALSE)
  covars <- unique(c(adj_var, candidate_var))
  is_surv <- isTRUE(surv)
  if (!is_surv && !missing(time))
    stop("`time` only applies to survival outcomes (`surv = TRUE`).",
         call. = FALSE)

  # Candidate type and levels are fixed on the full data, before any fit.
  kind <- vapply(candidate_var, function(v) {
    x <- data[[v]]
    if (!(is.numeric(x) || is.factor(x) || is.character(x) || is.logical(x)))
      stop(sprintf("`candidate_var` column `%s` must be numeric, logical, factor or character.",
                   v), call. = FALSE)
    k <- length(unique(x[!is.na(x)]))
    if (k < 2L)
      stop(sprintf("`candidate_var` column `%s` is constant.", v), call. = FALSE)
    if (k == 2L) return("binary")
    if (is.numeric(x)) return("continuous")
    stop(sprintf("`candidate_var` column `%s` has %d levels, but a categorical candidate needs exactly two levels: code the others as 0/1 indicators, or read its joint test in get_hte()$importance$p_het.",
                 v, k), call. = FALSE)
  }, character(1L), USE.NAMES = FALSE)
  levs <- lapply(candidate_var, function(v) {
    x <- data[[v]]
    if (is.factor(x)) levels(droplevels(x)) else
      as.character(sort(unique(x[!is.na(x)])))
  })

  # ---- AIPW arm scores -----------------------------------------------------
  if (method == "grf") {
    if (is.null(grf_args)) grf_args <- list()
    if (is.list(grf_args) && is.null(grf_args[["seed"]])) grf_args$seed <- seed
    res <- if (is_surv)
      get_hte(data, cat_var, adj_var = covars, surv = TRUE, measure = measure,
              conf_level = conf_level, time = time, grf_args = grf_args,
              verbose = verbose) else
      get_hte(data, cat_var, adj_var = covars, surv = surv, measure = measure,
              conf_level = conf_level, grf_args = grf_args, verbose = verbose)
    fit  <- res$fit
    an   <- attr(res, "analysis")
    d    <- res$data
    s    <- .hte_arm_scores(fit)
    W    <- fit$W.orig
    e    <- fit$W.hat
    mu1  <- fit$Y.hat + (1 - e) * s$tau
    mu0  <- fit$Y.hat - e * s$tau
    g1   <- s$g1
    g0   <- s$g0
    wt   <- .hte_weights(fit)
    cl   <- fit$clusters
    type <- an$outcome_type
    target  <- an$target
    treated <- an$treated
    outcome <- an$outcome
    backend_args <- an$grf_args
    beyond <- if (is_surv) {
      if (identical(target, "RMST")) d$time >= time else d$time > time
    }
    # get_hte() has already reported the skipped combination.
    if (type == "continuous" || identical(target, "RMST"))
      measure <- setdiff(measure, "OR")
  } else {
    if (is_surv)
      stop("method = \"glm\" does not handle survival outcomes; use method = \"grf\", which fits a causal survival forest.",
           call. = FALSE)
    if (!is.character(surv) || length(surv) != 1L || is.na(surv))
      stop("`surv` must be TRUE (columns `time` / `DSS`) or a single outcome column name.",
           call. = FALSE)
    outcome <- .sens_check_col(surv, data, "surv", n = 1L)
    if (outcome %in% c(cat_var, covars))
      stop("`adj_var` / `candidate_var` must not include the outcome column.",
           call. = FALSE)
    backend_args <- .merge_named_arg(glm_args, list(folds = 5L, ps_trim = 0.01),
                                     "glm_args")
    backend_args$folds <- .hte_select_count(backend_args$folds, "glm_args$folds",
                                            2, 100)
    pt <- backend_args$ps_trim
    if (!is.numeric(pt) || length(pt) != 1L || is.na(pt) || pt < 0 || pt >= 0.5)
      stop("`glm_args$ps_trim` must be a single number in [0, 0.5).",
           call. = FALSE)
    d  <- .sens_complete(data, c(cat_var, outcome, covars), verbose)
    tz <- .psw_treat(d[[cat_var]], cat_var, arg = "cat_var")
    W  <- tz$z
    treated <- tz$treated
    Y  <- d[[outcome]]
    if (is.logical(Y)) Y <- as.integer(Y)
    if (!is.numeric(Y))
      stop(sprintf("Outcome column `%s` must be numeric or logical; code a binary outcome as 0/1.",
                   outcome), call. = FALSE)
    u <- unique(Y)
    if (length(u) == 2L && !all(u %in% c(0, 1)))
      stop(sprintf("Outcome column `%s` has two values; code a binary outcome as 0/1.",
                   outcome), call. = FALSE)
    type <- if (all(u %in% c(0, 1))) "binary" else "continuous"
    if (min(sum(W == 1L), sum(W == 0L)) < 2L * backend_args$folds)
      stop("Each arm of `cat_var` needs at least 2 x `glm_args$folds` complete rows.",
           call. = FALSE)
    if ("OR" %in% measure && type == "continuous") {
      measure <- setdiff(measure, "OR")
      if (!length(measure))
        stop("`measure = \"OR\"` needs a binary outcome.", call. = FALSE)
      cli::cli_inform(c("i" = "Skipped {.val OR}: it needs an outcome probability."))
    }

    genv <- globalenv()
    old_seed <- if (exists(".Random.seed", envir = genv, inherits = FALSE))
      get(".Random.seed", envir = genv, inherits = FALSE)
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = genv)
      else if (exists(".Random.seed", envir = genv, inherits = FALSE))
        rm(".Random.seed", envir = genv)
    }, add = TRUE)
    set.seed(seed)
    nu <- .uni_glm_fit(d, covars, W, Y,
                       if (type == "binary") stats::binomial() else stats::gaussian(),
                       backend_args$folds, pt)
    if (nu$n_clip)
      warning(sprintf("%d estimated propensit%s of `%s` fell outside [%s, %s] and %s clipped (`glm_args$ps_trim`).",
                      nu$n_clip, if (nu$n_clip == 1L) "y" else "ies", cat_var,
                      format(pt), format(1 - pt),
                      if (nu$n_clip == 1L) "was" else "were"), call. = FALSE)
    e   <- nu$e
    mu1 <- nu$mu1
    mu0 <- nu$mu0
    g1  <- mu1 + W / e * (Y - mu1)
    g0  <- mu0 + (1 - W) / (1 - e) * (Y - mu0)
    wt  <- rep(1, length(W))
    cl  <- NULL
    target <- NULL
    beyond <- NULL
  }

  # ---- Pseudo-outcome per measure --------------------------------------------
  event_risk <- identical(target, "survival.probability")
  scores <- lapply(stats::setNames(measure, measure), function(m) {
    if (m == "diff") return(g1 - g0)
    a1 <- mu1
    a0 <- mu0
    b1 <- g1
    b0 <- g0
    if (event_risk) {
      a1 <- 1 - mu1
      a0 <- 1 - mu0
      b1 <- 1 - g1
      b0 <- 1 - g0
    }
    ok <- is.finite(a1) & is.finite(a0) & a1 > 0 & a0 > 0 &
      (m == "ratio" | (a1 < 1 & a0 < 1))
    if (any(!ok))
      warning(sprintf("`%s`: %d of %d rows have an arm prediction outside the range a %s needs and are left out.",
                      m, sum(!ok), length(ok),
                      if (m == "ratio") "log ratio" else "log odds"),
              call. = FALSE)
    a1[!ok] <- NA_real_
    a0[!ok] <- NA_real_
    if (m == "ratio")
      log(a1 / a0) + (b1 - a1) / a1 - (b0 - a0) / a0
    else stats::qlogis(a1) - stats::qlogis(a0) +
      (b1 - a1) / (a1 * (1 - a1)) - (b0 - a0) / (a0 * (1 - a0))
  })

  # ---- One least-squares projection per candidate x measure ----------------
  # n, sd, and the slope with its HC3 standard error -- per SD for a
  # continuous candidate -- or NA with a warning.
  project <- function(j, phi, m) {
    cont <- kind[j] == "continuous"
    x    <- d[[candidate_var[j]]]
    if (!cont) x <- match(as.character(x), levs[[j]]) - 1L
    pr <- .uni_project(x, phi, W, wt, cl, beyond, cont)
    if (!is.null(pr$why)) {
      warning(sprintf("`%s` (%s): %s; estimate set to NA.", candidate_var[j],
                      m, pr$why), call. = FALSE)
      return(rep(NA_real_, 4L))
    }
    k <- if (cont) pr$sd else 1
    c(length(pr$i), pr$sd, stats::coef(pr$fit)[[2L]] * k,
      sqrt(pr$V[2L, 2L]) * k)
  }

  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  contrast <- ifelse(kind == "continuous", "+1 SD",
                     vapply(levs, function(l) paste(l[2L], "vs", l[1L]),
                            character(1L)))
  vip <- do.call(rbind, lapply(measure, function(m) {
    vals <- vapply(seq_along(candidate_var),
                   function(j) project(j, scores[[m]], m), numeric(4L))
    b  <- vals[3L, ]
    se <- vals[4L, ]
    tr <- if (m == "diff") identity else exp
    p  <- 2 * stats::pnorm(-abs(b / se))
    out <- data.frame(variable = candidate_var, type = kind,
                      contrast = contrast, n = as.integer(vals[1L, ]),
                      sd = vals[2L, ], measure = m, estimate = tr(b),
                      std.error = se, conf.low = tr(b - z * se),
                      conf.high = tr(b + z * se), p.value = p,
                      p.adj = stats::p.adjust(p, method = "BH"),
                      stringsAsFactors = FALSE)
    out[order(out$p.value), , drop = FALSE]
  }))
  rownames(vip) <- NULL

  d$.ps  <- e
  d$.mu1 <- mu1
  d$.mu0 <- mu0
  for (m in measure) d[[paste0(".score_", m)]] <- scores[[m]]

  structure(
    list(vip = tibble::as_tibble(vip), data = d),
    analysis = list(
      method = method, outcome_type = type, cat_var = cat_var,
      treated = treated, surv = surv, candidate_var = candidate_var,
      adj_var = adj_var, covariates = covars, measure = measure,
      target = target, time = if (is_surv) time else NULL, outcome = outcome,
      conf_level = conf_level, n = length(W), n_treat = sum(W),
      ps_range = range(e), seed = seed, backend_args = backend_args,
      candidate_type = stats::setNames(kind, candidate_var),
      candidate_levels = stats::setNames(levs, candidate_var),
      weights = if (any(wt != 1)) wt, clusters = if (length(cl)) cl))
}

# ---- L1 plot -----------------------------------------------------------------

#' Plot an effect-modifier screen from get_hte_unihtee()
#'
#' Draws the table of [get_hte_unihtee()] for one measure as a bar chart of
#' the signed estimates with their confidence intervals (`"bar"`) or as a
#' volcano plot (`"volcano"`), or draws one panel per candidate with the
#' projection behind its row (`"dep"`). Candidates whose Benjamini-Hochberg
#' `p.adj` is below `sig_level` are coloured by the direction of the
#' estimate.
#'
#' @param x The list returned by [get_hte_unihtee()].
#' @param type `"bar"` (default), `"volcano"` or `"dep"`.
#' @param x_var Candidates to draw: `NULL` (default) for every candidate, in
#'   the `p.value` order of the table, or a character vector of candidate
#'   names in the order to draw. Rows whose estimate is `NA` are left out of
#'   `"bar"` and `"volcano"`.
#' @param measure The measure to draw, one of those in `x`; `NULL` (default)
#'   takes the first.
#' @param sig_level Threshold on `p.adj`, strictly between 0 and 1. Default
#'   `0.05`. Candidates below it are red (positive) or blue (negative) and, in
#'   the volcano plot, labelled; the volcano's dashed line marks the largest
#'   raw p-value that passes it.
#' @param dr_args Named list for `type = "dep"` only:
#'   \describe{
#'     \item{`spline_df`}{Finite whole number of at least 1, default `2`, as
#'       in [plt_hte_dep()]: degrees of freedom of the natural spline drawn
#'       beside the line of a continuous candidate with more than 5 distinct
#'       values.}
#'   }
#' @param title Plot title, or `NULL` (default).
#' @param save `NULL` or a list with `filename`, `width` and `height`, passed
#'   to `RegR::save_plt()` for PDF output. `list()` and `NULL` skip saving; a
#'   list naming only the file is completed with this figure's pinned size.
#'   Do not include the `plot` argument; it is supplied internally.
#'
#' @section What the panels show:
#' In `"dep"` the red line and band are the least-squares projection of the
#' pseudo-outcome `.score_<measure>` on the candidate with the confidence
#' band of its HC3 covariance: the slope times `sd` is the table's estimate
#' (on the log scale for `"ratio"` and `"OR"`). A binary candidate gets the
#' fitted mean of each level instead, whose difference is the estimate. The
#' grey dashed curve is a natural spline of the same pseudo-outcome, the
#' `"dr"` layer of [plt_hte_dep()]; where it bends away from the line the
#' slope summarises the modification poorly, and a U shape can hide behind a
#' flat line. The dashed horizontal line is the mean pseudo-outcome, the
#' doubly robust ATE for `"diff"`. All panels share one y range, which covers
#' the lines but not the bands.
#'
#' @return A `ggplot`, or for `"dep"` with several candidates a patchwork
#'   with one panel per candidate, three per row. The pinned size is in
#'   `attr(p, "plot_size")`. If `save` is non-empty, the same plot is also
#'   written to PDF through `RegR::save_plt()`.
#'
#' @seealso [get_hte_unihtee()]; [plt_hte_dep()] for the forest's own CATE.
#'
#' @examplesIf requireNamespace("sandwich", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)
#' set.seed(1)
#' n <- 600
#' d <- data.frame(w1 = rnorm(n), w2 = rnorm(n), w3 = rnorm(n),
#'                 sex = factor(sample(c("F", "M"), n, replace = TRUE)))
#' d$a <- rbinom(n, 1, plogis(0.5 * d$w1))
#' d$y <- d$w1 + d$w2 + d$a * (d$w3 + (d$sex == "M")) + rnorm(n)
#' res <- get_hte_unihtee(d, cat_var = "a", adj_var = c("w1", "w2", "w3", "sex"),
#'                        surv = "y", method = "glm")
#'
#' plt_hte_unihtee(res)                      # signed estimates with CIs
#' plt_hte_unihtee(res, type = "volcano")
#' plt_hte_unihtee(res, type = "dep", x_var = c("w3", "sex", "w2"))
#'
#' @export
plt_hte_unihtee <- function(x,
                            type      = c("bar", "volcano", "dep"),
                            x_var     = NULL,
                            measure   = NULL,
                            sig_level = 0.05,
                            dr_args   = list(spline_df = 2),
                            title     = NULL,
                            save      = list()) {

  a <- attr(x, "analysis")
  if (!is.list(x) || !all(c("vip", "data") %in% names(x)) ||
      is.null(a$candidate_levels))
    stop("`x` must be the list returned by get_hte_unihtee().", call. = FALSE)
  type <- match.arg(type)
  if (type != "dep" && !missing(dr_args))
    stop("`dr_args` only applies to type = \"dep\".", call. = FALSE)
  dr_args <- .merge_named_arg(dr_args, list(spline_df = 2), "dr_args")
  if (!is.numeric(dr_args$spline_df) || length(dr_args$spline_df) != 1L ||
      !is.finite(dr_args$spline_df) || dr_args$spline_df < 1 ||
      dr_args$spline_df != floor(dr_args$spline_df))
    stop("`dr_args$spline_df` must be a finite whole number of at least 1.",
         call. = FALSE)
  if (is.null(measure)) measure <- a$measure[1L]
  if (!is.character(measure) || length(measure) != 1L ||
      !measure %in% a$measure)
    stop(sprintf("`measure` must be one measure of `x`: %s.",
                 paste0("\"", a$measure, "\"", collapse = ", ")), call. = FALSE)
  if (!is.numeric(sig_level) || length(sig_level) != 1L || is.na(sig_level) ||
      sig_level <= 0 || sig_level >= 1)
    stop("`sig_level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  if (!is.null(save) && !is.list(save))
    stop("`save` must be `NULL` or a list.", call. = FALSE)

  v    <- x$vip[x$vip$measure == measure, , drop = FALSE]
  vars <- if (is.null(x_var)) v$variable else {
    if (!is.character(x_var) || anyNA(x_var))
      stop("`x_var` must be `NULL` or candidate names.", call. = FALSE)
    miss <- setdiff(x_var, v$variable)
    if (length(miss))
      stop(sprintf("`x_var` names no candidate of `x`: %s. Candidates are %s.",
                   paste0("`", miss, "`", collapse = ", "),
                   paste0("`", v$variable, "`", collapse = ", ")), call. = FALSE)
    unique(x_var)
  }
  v   <- v[match(vars, v$variable), , drop = FALSE]
  rel <- measure != "diff"

  # the measures as get_hte() names them
  ref   <- setdiff(levels(factor(x$data[[a$cat_var]])), a$treated)[1L]
  scale <- switch(
    a$outcome_type,
    survival   = if (identical(a$target, "RMST"))
      c(diff = "RMST difference", ratio = "RMST ratio")[[measure]] else
      c(diff = sprintf("S(%s) difference", format(a$time)),
        ratio = "event-risk ratio", OR = "event odds ratio")[[measure]],
    binary     = c(diff = "risk difference", ratio = "risk ratio",
                   OR = "odds ratio")[[measure]],
    continuous = c(diff = "mean difference", ratio = "ratio of means")[[measure]])
  conf <- 100 * a$conf_level
  sig  <- factor(ifelse(!is.na(v$p.adj) & v$p.adj < sig_level,
                        ifelse((if (rel) log(v$estimate) else v$estimate) > 0,
                               "Positive", "Negative"), "n.s."),
                 levels = c("Positive", "Negative", "n.s."))
  cols <- c(Positive = "firebrick", Negative = "steelblue", n.s. = "grey70")
  null <- if (rel) 1 else 0

  # ---- Bar and volcano ------------------------------------------------------
  if (type != "dep") {
    ok <- is.finite(v$estimate)
    if (!any(ok))
      stop("No candidate of `x` has an estimate for this measure.", call. = FALSE)
    pd <- data.frame(variable = v$variable,
                     label = paste0(v$variable, " (", v$contrast, ")"),
                     estimate = v$estimate, conf.low = v$conf.low,
                     conf.high = v$conf.high, neglog10p = -log10(v$p.value),
                     sig = sig, stringsAsFactors = FALSE)[ok, , drop = FALSE]
    xlab <- sprintf("%s the %s per SD or between levels",
                    if (rel) "Factor on" else "Change in", scale)
    # About 10 caption characters fit per inch of the 6.5-inch figure, as in
    # plt_hte_dep().
    wrap <- function(s) paste(strwrap(s, width = 65L), collapse = "\n")
    if (type == "bar") {
      pd$label <- factor(pd$label, levels = pd$label[order(
        if (rel) log(pd$estimate) else pd$estimate)])
      p <- ggplot2::ggplot() +
        ggplot2::geom_col(data = pd, ggplot2::aes(x = estimate, y = label,
                                                  fill = sig), width = 0.7) +
        ggplot2::geom_errorbar(data = pd, ggplot2::aes(xmin = conf.low,
                                                       xmax = conf.high,
                                                       y = label),
                               width = 0.25, orientation = "y") +
        ggplot2::geom_vline(xintercept = null, colour = "grey40") +
        ggplot2::scale_fill_manual(values = cols, name = NULL) +
        ggplot2::labs(x = xlab, y = NULL, title = title,
                      caption = wrap(sprintf("Bars: estimate with %g%% CI; coloured: BH p.adj < %g",
                                             conf, sig_level)))
      size <- c(6.5, max(3, 1.2 + 0.35 * nrow(pd)))
    } else {
      pass <- !is.na(v$p.adj[ok]) & v$p.adj[ok] < sig_level
      p <- ggplot2::ggplot() +
        ggplot2::geom_vline(xintercept = null, colour = "grey40")
      if (any(pass))
        p <- p + ggplot2::geom_hline(yintercept = -log10(max(v$p.value[ok][pass])),
                                     linetype = 2, colour = "grey40")
      p <- p +
        ggplot2::geom_point(data = pd, ggplot2::aes(x = estimate, y = neglog10p,
                                                    colour = sig), size = 2.5) +
        ggplot2::geom_text(data = pd[pd$sig != "n.s.", , drop = FALSE],
                           ggplot2::aes(x = estimate, y = neglog10p,
                                        label = variable),
                           vjust = -0.8, size = 3.2) +
        ggplot2::scale_colour_manual(values = cols, name = NULL) +
        # room above the top point for its label
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.12))) +
        ggplot2::labs(x = xlab, y = "-log10(p)", title = title,
                      caption = wrap(sprintf("Coloured and labelled: BH p.adj < %g%s",
                                             sig_level, if (any(pass))
                                               "; dashed: the largest raw p that passes" else "")))
      size <- c(6.5, 5)
    }
    if (rel) p <- p + ggplot2::scale_x_log10()
    p <- p + UtilsR::theme_my(base_rect_size = 1.5)

  # ---- One panel per candidate ---------------------------------------------
  } else {
    if (!requireNamespace("sandwich", quietly = TRUE))
      stop("Package 'sandwich' is required for the projection bands.",
           call. = FALSE)
    if (length(vars) > 1L && !requireNamespace("patchwork", quietly = TRUE))
      stop("Package 'patchwork' is required to combine several panels.",
           call. = FALSE)
    d      <- x$data
    W      <- as.integer(as.character(d[[a$cat_var]]) == a$treated)
    phi    <- d[[paste0(".score_", measure)]]
    wt     <- if (is.null(a$weights)) rep(1, nrow(d)) else a$weights
    beyond <- .hte_beyond(x)
    z      <- stats::qnorm(1 - (1 - a$conf_level) / 2)
    mid    <- stats::weighted.mean(phi, wt, na.rm = TRUE)
    ylab   <- if (rel) sprintf("log %s: %s vs %s", scale, a$treated, ref)
              else .hte_cate_label(x)
    fmt_p  <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "< 0.001"
                          else sprintf("= %.3f", p)
    dd     <- d
    dd$.dr_score <- phi
    has_spline <- FALSE

    panel <- function(j) {
      vn    <- vars[j]
      cont  <- a$candidate_type[[vn]] == "continuous"
      lv    <- a$candidate_levels[[vn]]
      xv    <- d[[vn]]
      if (!cont) xv <- match(as.character(xv), lv) - 1L
      pr    <- .uni_project(xv, phi, W, wt, a$clusters, beyond, cont)
      label <- sprintf("%s: %s %s, p.adj %s", vn, v$contrast[j],
                       format(signif(v$estimate[j], 3)), fmt_p(v$p.adj[j]))
      q  <- ggplot2::ggplot(data.frame(panel = label)) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey75") +
        ggplot2::geom_hline(yintercept = mid, linetype = 2)
      yv <- c(0, mid)
      if (is.null(pr$why)) {
        g   <- if (cont) seq(pr$range[1L], pr$range[2L], length.out = 100L) else 0:1
        M   <- cbind(1, g)
        est <- drop(M %*% stats::coef(pr$fit))
        se  <- sqrt(rowSums((M %*% pr$V) * M))
        ln  <- data.frame(x = if (cont) g else factor(lv, levels = lv),
                          estimate = est, conf.low = est - z * se,
                          conf.high = est + z * se, panel = label)
        yv  <- c(yv, est)
        q <- q + if (cont) list(
          ggplot2::geom_ribbon(data = ln, ggplot2::aes(x = x, ymin = conf.low,
                                                       ymax = conf.high),
                               fill = "firebrick", alpha = 0.15),
          ggplot2::geom_line(data = ln, ggplot2::aes(x = x, y = estimate),
                             colour = "firebrick", linewidth = 0.8)) else
          ggplot2::geom_pointrange(data = ln, ggplot2::aes(x = x, y = estimate,
                                                           ymin = conf.low,
                                                           ymax = conf.high),
                                   colour = "firebrick")
        if (cont && .hte_is_num(d[[vn]])) {
          sp <- .hte_dr_var(dd, vn, W, z, dr_args$spline_df, beyond, wt,
                            a$clusters)$curve
          if (NROW(sp)) {
            sc <- data.frame(x = sp$x, spline = sp$estimate, panel = label)
            q  <- q + ggplot2::geom_line(data = sc, ggplot2::aes(x = x, y = spline),
                                         colour = "grey30", linetype = "longdash",
                                         linewidth = 0.7)
            yv <- c(yv, sc$spline)
            has_spline <<- TRUE
          }
        }
      }
      q <- q + ggplot2::facet_wrap(~panel) +
        ggplot2::labs(x = vn, y = ylab) +
        UtilsR::theme_my(base_rect_size = 1.5)
      list(plot = q, yv = yv)
    }

    built <- lapply(seq_along(vars), panel)
    rng   <- range(unlist(lapply(built, `[[`, "yv")), na.rm = TRUE)
    plots <- lapply(built, function(b) b$plot + ggplot2::coord_cartesian(ylim = rng))
    caption <- paste(c(
      sprintf("red: projection on the candidate with %g%% CI (level means for a binary one)",
              conf),
      if (has_spline) sprintf("grey dashed: natural spline, df = %d",
                              as.integer(dr_args$spline_df)),
      sprintf("dashed: mean pseudo-outcome%s", if (rel) "" else " (ATE)")),
      collapse = "; ")
    caption <- paste0(toupper(substr(caption, 1L, 1L)), substring(caption, 2L))
    ncol <- min(3L, length(plots))
    size <- if (length(plots) == 1L) c(5, 4.2) else
      c(3.4 * ncol + 0.6, 3 * ceiling(length(plots) / ncol) + 0.8)
    caption <- paste(strwrap(caption, width = floor(10 * size[1L])),
                     collapse = "\n")
    p <- if (length(plots) == 1L) {
      plots[[1L]] + ggplot2::labs(title = title, caption = caption)
    } else {
      patchwork::wrap_plots(plots, ncol = ncol) +
        patchwork::plot_layout(axis_titles = "collect") +
        patchwork::plot_annotation(title = title, caption = caption,
                                   theme = UtilsR::theme_my(base_rect_size = 1.5))
    }
  }

  attr(p, "plot_size") <- stats::setNames(size, c("width", "height"))
  if (!is.null(save) && length(save) > 0L) {
    if (!requireNamespace("RegR", quietly = TRUE))
      stop("Package 'RegR' is required for a non-empty `save`.", call. = FALSE)
    if (is.null(save$width))  save$width  <- size[1L]
    if (is.null(save$height)) save$height <- size[2L]
    do.call(RegR::save_plt, c(list(plot = p), save))
  }
  p
}
