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
#   L2  .uni_glm_fit()     cross-fitted propensity and per-arm outcome GLMs
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
#' levels of a categorical covariate -- and
#' `plt_hte_dep(res, display = "dr", dr_args = list(spline_df = 1))` draws
#' the straight line whose slope, per unit of the covariate, is the
#' `"diff"` estimate here divided by `sd`. The table screens candidates; it
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
#' @seealso [get_hte()] for the forest, its spline test `p_het` and subgroup
#'   effects; [plt_hte_dep()] for the doubly robust curves;
#'   [get_hte_select()] for benefit-score variable selection.
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
  # Returns n, sd, slope and its HC3 standard error, or NA with a warning.
  project <- function(j, phi, m) {
    v <- candidate_var[j]
    x <- d[[v]]
    na_row <- function(why) {
      warning(sprintf("`%s` (%s): %s; estimate set to NA.", v, m, why),
              call. = FALSE)
      c(NA_real_, NA_real_, NA_real_, NA_real_)
    }
    cont <- kind[j] == "continuous"
    if (!cont) x <- match(as.character(x), levs[[j]]) - 1L
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
        return(na_row("a level has no patient in one arm followed past `time`"))
    }
    i <- which(!is.na(x) & is.finite(phi) & wt > 0)
    enough <- if (cont) length(unique(x[i])) >= 3L else
      all(vapply(0:1, function(l) sum(x[i] == l & W[i] == 1) >= 2 &&
                   sum(x[i] == l & W[i] == 0) >= 2, logical(1L)))
    if (!enough)
      return(na_row(if (cont) "fewer than three distinct values remain"
                    else "a level has fewer than two patients in one arm"))
    sdv <- if (cont) sqrt(stats::cov.wt(cbind(x[i]), wt[i])$cov[1L, 1L]) else NA_real_
    lm_fit <- stats::lm(p ~ xs, weights = w,
                        data = data.frame(p = phi[i],
                                          xs = if (cont) x[i] / sdv else x[i],
                                          w = wt[i]))
    if (lm_fit$rank < 2L || lm_fit$df.residual <= 0L ||
        (length(cl) && length(unique(cl[i])) < 2L))
      return(na_row("the regression is rank deficient or has too few clusters"))
    V <- if (length(cl))
      sandwich::vcovCL(lm_fit, cluster = cl[i], type = "HC3") else
      sandwich::vcovHC(lm_fit, type = "HC3")
    c(length(i), sdv, stats::coef(lm_fit)[[2L]], sqrt(V[2L, 2L]))
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
      target = target, time = if (is_surv) time else NULL,
      conf_level = conf_level, n = length(W), n_treat = sum(W),
      ps_range = range(e), seed = seed, backend_args = backend_args))
}
