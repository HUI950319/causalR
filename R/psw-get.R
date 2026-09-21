# =============================================================================
# psw-get.R -- propensity score weights for a binary exposure
# =============================================================================
#
# Architecture (3 layers + shared helpers from utils-sens.R):
#
#   L1  get_PSW(data, treat, adj_var, ps, estimand, ...)
#         |
#         +-- L2 pipeline stages, run in this fixed order
#         |     .psw_score     the score: .psw_fit() or the `ps` column
#         |     .psw_fit       WeightIt::weightit(), score only
#         |     .psw_trim      set trimmed units to NA, optionally refit
#         |     .psw_trunc     clamp the score into an interval
#         |     .psw_balance   halfmoon::check_balance() over every weight
#         |
#         +-- L3 helpers
#               .psw_tilt        the tilting function h(e) of one estimand
#               .psw_crump       Crump et al. optimal symmetric trim point
#               .psw_ess         effective sample size of a weight vector
#               .psw_smd         smd_max / smd_over of one weight column
#               .psw_stats_row   the 15-column standardised $stats row
#               .psw_treat       resolve the exposure to a 0/1 integer
#               .psw_plt_spec    plt_PSW() arguments echoed by print()
#
#   print.psw_res reports the diagnostics table and the matching plt_PSW()
#   call. .psw_score, .psw_fit, .psw_balance, .psw_ess, .psw_smd and
#   .psw_treat are shared with get_PSM(), the matching dual in psm-get.R.
#
# Every estimand shares one processed propensity score. Trimming or truncating
# per estimand would leave the weight columns describing different analysis
# populations, which is exactly the side-by-side comparison the function
# exists for.
#
# WeightIt estimates the score and nothing else. The weights themselves are
# the closed-form tilting functions of Li, Morgan & Zaslavsky (2018), one line
# each, because WeightIt cannot produce the entropy weight at all and its
# other weights would then come from a different code path than that one.
# Every weight column is checked against propensity::wt_*() in the tests.
# =============================================================================

.PSW_ESTIMANDS <- c("ATE", "ATT", "ATC", "ATO", "ATM", "EW")

# Stabilisation multiplies by the marginal probability of the observed
# exposure. That only recentres a weight whose tilting function is free of the
# score, i.e. ATE -- which is also the only weight propensity gives a
# `stabilize` argument to.
.PSW_STABILIZE <- "ATE"

# WeightIt methods that actually estimate a propensity score. The balancing
# weight methods ("ebal", "energy", "optweight", ...) return weights without a
# score, so they cannot feed a tilting function.
.PSW_PS_METHODS <- c("glm", "gbm", "cbps", "bart", "super")

.PSW_TRIM_DEFAULTS  <- list(method = "none", lower = NULL, upper = NULL,
                            refit = TRUE)
.PSW_TRUNC_DEFAULTS <- list(method = "none", lower = 0.01, upper = 0.99)

# Arguments the function manages itself and will not forward to the backend.
.PSW_MANAGED <- c("formula", "data", "method", "estimand", "ps", "x", "y")


# ---- L3 helpers ------------------------------------------------------------

# The six tilting functions h(e) of Li, Morgan & Zaslavsky (2018). The weight
# is h(e) / P(Z = z | X), so a treated unit gets h(e)/e and a control unit
# h(e)/(1-e). log1p() keeps the entropy tilt accurate as e approaches 1.
#' @keywords internal
#' @noRd
.psw_tilt <- function(estimand, ps) {
  switch(
    estimand,
    ATE = rep(1, length(ps)),
    ATT = ps,
    ATC = 1 - ps,
    ATO = ps * (1 - ps),
    ATM = pmin(ps, 1 - ps),
    EW  = -(ps * log(ps) + (1 - ps) * log1p(-ps)),
    stop(sprintf("Unsupported estimand: '%s'", estimand), call. = FALSE))
}

# Crump, Hotz, Imbens & Mitnik (2009): the optimal symmetric cut-off is the
# smallest alpha whose retained set satisfies
#   1 / (alpha (1 - alpha))  <=  2 * mean(1 / (e (1 - e)))
# The criterion is evaluated at the observed scores rather than on a grid.
# Sorting by min(e, 1 - e) makes the retained set of every candidate alpha a
# suffix of that order, so each mean is a running mean and the search is one
# pass. The alpha-by-alpha loop this replaces recomputed the mean from
# scratch each time and was quadratic: 420 s at n = 2e5 with poor overlap,
# against 0.05 s here, for the same alpha.
#' @keywords internal
#' @noRd
.psw_crump <- function(ps) {
  m   <- pmin(ps, 1 - ps)
  o   <- order(m)
  ms  <- m[o]
  v   <- 1 / (ps[o] * (1 - ps[o]))
  sfx <- rev(cumsum(rev(v))) / rev(seq_along(v))   # mean(v) over units k..n
  k   <- which(!duplicated(ms) & ms < 0.5)         # first unit of each alpha
  ok  <- 1 / (ms[k] * (1 - ms[k])) <= 2 * sfx[k]
  if (any(ok)) ms[k[which(ok)[1L]]] else 0
}

# Effective sample size (sum w)^2 / sum w^2 of the units that carry weight:
# NA (trimmed) and 0 (unmatched) contribute nothing either way.
#' @keywords internal
#' @noRd
.psw_ess <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (!length(w)) return(NA_real_)
  sum(w)^2 / sum(w^2)
}

# The two balance summaries of one weight column: the largest absolute
# standardised mean difference across adj_var, and how many sit above 0.1,
# the line the love plot draws by default.
#' @keywords internal
#' @noRd
.psw_smd <- function(bal, col) {
  if (is.null(bal)) return(c(NA_real_, NA_real_))
  s <- abs(bal$estimate[bal$metric == "smd" & bal$method == col])
  s <- s[is.finite(s)]
  if (!length(s)) return(c(NA_real_, NA_real_))
  c(max(s), sum(s > 0.1))
}

#' @keywords internal
#' @noRd
.psw_stats_row <- function(estimand, w, z, keep, smd_max = NA_real_,
                           smd_over = NA_real_) {
  wk <- w[keep]
  tibble::tibble(
    estimand  = estimand,
    n         = sum(keep),
    n_treat   = sum(keep & z == 1L),
    n_ctrl    = sum(keep & z == 0L),
    ess       = .psw_ess(wk),
    ess_treat = .psw_ess(w[keep & z == 1L]),
    ess_ctrl  = .psw_ess(w[keep & z == 0L]),
    ess_pct   = .psw_ess(wk) / sum(keep),
    w_min     = min(wk),
    w_max     = max(wk),
    w_mean    = mean(wk),
    w_sd      = stats::sd(wk),
    w_cv      = stats::sd(wk) / mean(wk),
    smd_max   = as.numeric(smd_max),
    smd_over  = as.numeric(smd_over))
}

# Resolve the exposure to a 0/1 integer and name the arm that became 1. A
# two-level factor or character column is read as "the second level is the
# treated arm", which is also what WeightIt scores, so the score and the
# weights cannot disagree about direction. A character column is ordered
# alphabetically, which is the wrong way round for "case" / "control", so the
# label is returned alongside the coding for print() to show.
#' @keywords internal
#' @noRd
.psw_treat <- function(x, nm) {
  if (is.logical(x)) return(list(z = as.integer(x), treated = "TRUE"))
  if (is.factor(x) || is.character(x)) {
    lv <- if (is.factor(x)) levels(droplevels(x)) else sort(unique(x))
    if (length(lv) != 2L)
      stop(sprintf("`treat` column `%s` must have exactly 2 levels; found %d.",
                   nm, length(lv)), call. = FALSE)
    return(list(z = as.integer(match(as.character(x), lv) - 1L),
                treated = lv[[2L]]))
  }
  u <- sort(unique(x))
  if (!all(u %in% c(0, 1)))
    stop(sprintf("`treat` column `%s` must be 0/1, logical, or a two-level factor or character column.",
                 nm), call. = FALSE)
  list(z = as.integer(x), treated = "1")
}


# ---- L2 pipeline stages ----------------------------------------------------

# Every score comes from WeightIt, including the default logistic one: routing
# the backends through one estimator keeps `method` a single switch with no
# special case, and `weightit()` accepts factor and character exposures that
# stats::glm() rejects outright. Only `$ps` is used -- the weights weightit
# computes for its own `estimand` are discarded, since the six tilting
# functions are built here.
#
# `estimand = "ATE"` is fixed because weightit() requires one; for every
# method below it leaves the score itself untouched.
#
# Direction: weightit() scores the *second* level of a factor or character
# exposure, matching .psw_treat() and stats::glm(). Verified equal to
# stats::glm(family = binomial()) to 4.4e-16 on numeric, factor and character
# exposures, so nothing about the weights changed when the glm branch went.
#
# The call is assembled with `data` as a symbol rather than do.call()'d:
# do.call() inlines the function body and the whole data frame into
# `obj$call`, which then makes up most of the object (measured at n = 400:
# 665 KB against 113 KB) and is what any printed call would show. The symbol
# resolves in this frame, which the formula environment keeps reachable.
#' @keywords internal
#' @noRd
.psw_fit <- function(data, treat, adj_var, method, ps_args) {
  form <- stats::reformulate(adj_var, response = treat)
  cl   <- as.call(c(list(quote(WeightIt::weightit)),
                    list(formula = form, data = quote(data), method = method,
                         estimand = "ATE"),
                    ps_args))
  obj  <- eval(cl)
  if (is.null(obj$ps))
    stop(sprintf("WeightIt method \"%s\" returns balancing weights without a propensity score, so it cannot feed a tilting function. Use one of %s.",
                 method,
                 paste0("\"", .PSW_PS_METHODS, "\"", collapse = " / ")),
         call. = FALSE)
  list(ps = unname(as.numeric(obj$ps)), fit = obj)
}

# The score, fitted through WeightIt or taken from the `ps` column, validated
# either way: a supplied score can carry NA or sit on 0 / 1.
#' @keywords internal
#' @noRd
.psw_score <- function(data, treat, adj_var, ps, method, ps_args) {
  f <- if (is.null(ps)) .psw_fit(data, treat, adj_var, method, ps_args)
       else list(ps = as.numeric(data[[ps]]), fit = NULL)
  if (anyNA(f$ps) || any(f$ps <= 0 | f$ps >= 1))
    stop("The propensity score must be non-missing and strictly between 0 and 1.",
         call. = FALSE)
  f
}

# Trimming removes units from the analysis population but not rows from the
# data: a trimmed unit's score becomes NA, so every weight column stays
# aligned with the input and downstream code sees one NA pattern. Same
# convention as propensity::ps_trim().
#' @keywords internal
#' @noRd
.psw_trim <- function(ps, ta) {
  if (identical(ta$method, "none")) return(list(ps = ps, bounds = NULL))

  b <- switch(
    ta$method,
    ps   = c(if (is.null(ta$lower)) 0.1 else ta$lower,
             if (is.null(ta$upper)) 0.9 else ta$upper),
    pctl = unname(stats::quantile(
      ps,
      c(if (is.null(ta$lower)) 0.01 else ta$lower,
        if (is.null(ta$upper)) 0.99 else ta$upper))),
    cr   = {
      a <- .psw_crump(ps)
      c(a, 1 - a)
    })

  if (!is.numeric(b) || length(b) != 2L || anyNA(b) || b[1L] >= b[2L])
    stop("`trim_args` must give a lower bound strictly below the upper bound.",
         call. = FALSE)
  ps[ps < b[1L] | ps > b[2L]] <- NA_real_
  list(ps = ps, bounds = b)
}

#' @keywords internal
#' @noRd
.psw_trunc <- function(ps, ua) {
  if (identical(ua$method, "none")) return(list(ps = ps, bounds = NULL))

  b <- switch(
    ua$method,
    ps   = c(ua$lower, ua$upper),
    pctl = unname(stats::quantile(ps, c(ua$lower, ua$upper), na.rm = TRUE)))

  if (!is.numeric(b) || length(b) != 2L || anyNA(b) || b[1L] >= b[2L])
    stop("`trunc_args` must give a lower bound strictly below the upper bound.",
         call. = FALSE)
  list(ps = pmin(pmax(ps, b[1L]), b[2L]), bounds = b)
}

# halfmoon returns a long table whose `method` column holds the weight column
# name, or "observed" for the unweighted comparison.
#
# Two things about the call are load-bearing. `.exposure` is captured with
# rlang::enquo() and resolved by name, so it has to reach check_balance() as a
# literal string rather than as the symbol `treat`; do.call() inlines the
# values and is also the only form tidyselect takes for `.vars` / `.weights`
# without a deprecation warning. And get_PSW() drops the trimmed rows through
# `keep`: a weight column with any NA makes check_balance() return NA for
# that whole weight -- silently, and even with na.rm = TRUE. Balance belongs
# to the analysis population in any case, which is exactly the retained set.
# get_PSM() passes the whole cohort; its call site says why.
#' @keywords internal
#' @noRd
.psw_balance <- function(data, adj_var, treat, wcols, keep = NULL,
                         caller = "get_PSW") {
  if (!requireNamespace("halfmoon", quietly = TRUE))
    stop(sprintf("Package 'halfmoon' is required for %s(balance = TRUE); use balance = FALSE to skip the balance table.",
                 caller), call. = FALSE)
  if (!is.null(keep)) data <- data[keep, , drop = FALSE]
  do.call(halfmoon::check_balance,
          list(.data     = data,
               .vars     = adj_var,
               .exposure = treat,
               .weights  = wcols,
               .metrics  = "smd",
               na.rm     = TRUE))
}

# The love plot is the figure to reach for first, but it needs the balance
# table; without one the effective sample sizes are all there is to show.
#' @keywords internal
#' @noRd
.psw_plt_spec <- function(x) {
  list(type = if (is.null(x$balance)) "ess" else "love")
}


# ---- L1 public entry point -------------------------------------------------

#' Propensity score weights for a binary exposure
#'
#' Single entry point that turns one propensity score into any combination of
#' the six tilting-function weights -- inverse probability, SMR, overlap,
#' matching and entropy weights -- and reports their effective sample size,
#' weight distribution and covariate balance side by side, so the weighting
#' scheme can be chosen on evidence rather than habit.
#'
#' `get_PSW()` constructs and diagnoses weights; it does not estimate a
#' treatment effect. Feed `result$data` and the weight column of your choice
#' to an outcome model, for example
#' `survival::coxph(Surv(t, d) ~ z, data = result$data, weights = w_ato, robust = TRUE)`.
#' Note that `RegR::get_eff_cat()` is **not** a downstream consumer of these
#' columns: it refits its own propensity model internally, so trimming,
#' truncation and refitting done here would not carry over to it.
#'
#' @section Weighting schemes:
#' With \eqn{e = P(Z = 1 \mid X)} and tilting function \eqn{h(e)}, the weight
#' is \eqn{h(e) / P(Z = z \mid X)}: \eqn{h(e)/e} for a treated unit and
#' \eqn{h(e)/(1-e)} for a control unit (Li, Morgan and Zaslavsky, 2018).
#'
#' \describe{
#'   \item{`"ATE"`}{\eqn{h = 1}; weights \eqn{1/e} and \eqn{1/(1-e)}. Inverse
#'     probability of treatment weighting, IPTW; `PSweight`'s `"IPW"`.}
#'   \item{`"ATT"`}{\eqn{h = e}; weights \eqn{1} and \eqn{e/(1-e)}. SMR
#'     weighting; `PSweight`'s `"treated"`.}
#'   \item{`"ATC"`}{\eqn{h = 1-e}; weights \eqn{(1-e)/e} and \eqn{1}. Also
#'     written ATU.}
#'   \item{`"ATO"`}{\eqn{h = e(1-e)}; weights \eqn{1-e} and \eqn{e}. Overlap
#'     weights, OW; `PSweight`'s `"overlap"`. Balances the covariate means of
#'     the two arms exactly when the score is a logistic fit on those
#'     covariates, and minimises the variance of the weighted contrast.}
#'   \item{`"ATM"`}{\eqn{h = \min(e, 1-e)}. Matching weights; `PSweight`'s
#'     `"matching"`.}
#'   \item{`"EW"`}{\eqn{h = -(e \ln e + (1-e) \ln(1-e))}. Entropy weights,
#'     also written ATEN; `PSweight`'s `"entropy"`.}
#' }
#'
#' Only those six names are accepted; the aliases are listed for orientation,
#' not as arguments.
#'
#' @section Agreement with WeightIt:
#' For the default `method = "glm"` the score does not depend on the estimand,
#' so every weight here is bit-identical to `WeightIt::weightit()`'s own: ATE,
#' ATT, ATC, ATO and ATM all agree to 1.8e-15, in both arms and without any
#' rescaling. `"EW"` has no WeightIt counterpart.
#'
#' The score-adaptive backends are a different matter. `"cbps"` fits the score
#' so as to balance for a *particular* estimand, and `"gbm"` tunes its trees
#' against an estimand-specific criterion, so their `$ps` changes with
#' `estimand` (measured here: 3.8e-02 for cbps, 4.6e-01 for gbm). `get_PSW()`
#' fits the score once, with `estimand = "ATE"`, and shares it across every
#' tilting function -- weight columns built on six different scores could not
#' be compared, and `"EW"` has no score of its own to fit against. So with
#' these backends `get_PSW(method = "cbps")` is deliberately *not* the same
#' estimator as `WeightIt::weightit(method = "cbps", estimand = "ATO")`.
#'
#' `"gbm"` is also stochastic: two identical calls differ unless you
#' [set.seed()] first.
#'
#' @section Order of operations:
#' The score is fitted (or taken from `ps`), then trimmed, then optionally
#' refitted on the retained units, then truncated, and only then turned into
#' weights; stabilisation is applied last. Every estimand shares that one
#' processed score, because weight columns built on different analysis
#' populations could not be compared.
#'
#' @param data A data frame holding every column named below.
#' @param treat Length-1 character. The binary exposure column: `0`/`1`,
#'   logical, or a two-level factor or character column whose **second** level
#'   is the treated one. A character column is ordered alphabetically, so
#'   `"case"` / `"control"` would make `"control"` the treated arm; convert
#'   such a column to a factor with the control level first. The arm actually
#'   taken as treated is recorded in `attr(x, "analysis")$treated` and shown
#'   by `print()`. The column is returned in `$data` exactly as supplied: the
#'   0/1 coding is used internally, so neither the direction of the score nor
#'   the balance table depends on the labels.
#' @param adj_var Character vector of covariates the propensity model adjusts
#'   for. Required unless `ps` is supplied, and required either way when
#'   `balance = TRUE`.
#' @param ps Length-1 character or `NULL` (default). Column holding an
#'   already-estimated propensity score, strictly inside (0, 1). Supplying it
#'   skips the modelling step, so `method` and `ps_args` no longer apply.
#' @param estimand Character vector, any of `"ATE"`, `"ATT"`, `"ATC"`,
#'   `"ATO"`, `"ATM"`, `"EW"`. All six by default; each produces one weight
#'   column and one `$stats` row.
#' @param method Propensity model backend, passed to [WeightIt::weightit()],
#'   from which only the score is taken. `"glm"` (default) is logistic
#'   regression and matches `stats::glm(family = binomial())` exactly; the
#'   others are `"gbm"`, `"cbps"`, `"bart"` and `"super"`; `"super"` needs its
#'   learners named through `ps_args`, as in
#'   `ps_args = list(SL.library = c("SL.glm", "SL.mean"))`. The
#'   balancing-weight methods (`"ebal"`, `"energy"`, `"optweight"`) return no
#'   propensity score and are rejected.
#'
#'   Every `method` works with every `estimand`, including combinations
#'   `WeightIt::weightit()` refuses -- `"cbps"` with `"ATM"`, and any method
#'   with `"EW"` -- because the score is taken from WeightIt but the weight is
#'   built here.
#' @param stabilize Logical, default `FALSE`. Multiply the weight by the
#'   marginal probability of the observed exposure, which recentres it on 1.
#'   Defined for `"ATE"` only; other weights are left untouched. `TRUE` when
#'   no selected estimand can use it is an error, not a silent no-op.
#' @param trim_args Named list controlling which units leave the analysis
#'   population. Trimmed units keep their row and take `NA` in the score and
#'   in every weight column.
#'   \describe{
#'     \item{`method`}{`"none"` (default), `"ps"` (absolute score bounds,
#'       defaulting to 0.1 and 0.9), `"pctl"` (score quantiles, defaulting to
#'       0.01 and 0.99), or `"cr"` (the Crump et al. 2009 optimal symmetric
#'       cut-off computed from the data, which takes no `lower` or `upper`
#'       and rejects them rather than ignoring them).}
#'     \item{`lower`,`upper`}{Numeric bounds or quantile probabilities, or
#'       `NULL` (default) for the per-method defaults above.}
#'     \item{`refit`}{Logical, default `TRUE`. Re-estimate the propensity
#'       model on the retained units, which is what trimming is meant to be
#'       followed by. Checked only when trimming is actually requested, and an
#'       error when the score came from `ps` and there is no model to refit.
#'       The refitted score is not trimmed again, so a retained unit can end
#'       up with a score outside the window: `.trimmed` records the first
#'       pass, and `print()` shows the window next to the final range.}
#'   }
#' @param trunc_args Named list controlling score truncation, which keeps every
#'   unit but pulls extreme scores in. `method` is `"none"` (default), `"ps"`
#'   (clamp to `lower`, `upper`) or `"pctl"` (clamp to those quantiles);
#'   `lower` and `upper` default to 0.01 and 0.99.
#' @param balance Logical, default `TRUE`. Compute the standardised mean
#'   differences of `adj_var` under every weight with
#'   `halfmoon::check_balance()`.
#' @param ps_args Named list forwarded to [WeightIt::weightit()], for example
#'   `list(link = "probit")` for `method = "glm"`. `formula`, `data`,
#'   `method` and `estimand` are managed here and rejected if supplied.
#' @param verbose Logical, default `FALSE`. Report how many units trimming
#'   removed and how many rows were dropped as incomplete.
#'
#' @return An object of class `psw_res`: a list of
#'   \describe{
#'     \item{`data`}{The input data plus `ps` (the processed score, `NA` where
#'       trimmed), `.trimmed` (logical) and one weight column per estimand,
#'       named `w_ate`, `w_att`, `w_atc`, `w_ato`, `w_atm`, `w_ew`. The `w_`
#'       prefix is deliberate: it is what
#'       `halfmoon::plot_ess(.weights = starts_with("w_"))` selects on.}
#'     \item{`stats`}{One tibble row per estimand, always with the same 15
#'       columns: `estimand`, `n`, `n_treat`, `n_ctrl`, `ess`, `ess_treat`,
#'       `ess_ctrl`, `ess_pct`, `w_min`, `w_max`, `w_mean`, `w_sd`, `w_cv`,
#'       `smd_max`, `smd_over`. The effective sample size is
#'       \eqn{(\sum w)^2 / \sum w^2}; `smd_max` is the largest absolute
#'       standardised mean difference across `adj_var`, and `smd_over` counts
#'       those above 0.1. Both are `NA` when `balance = FALSE`.}
#'     \item{`balance`}{The `halfmoon::check_balance()` long table
#'       (`variable`, `group_level`, `method`, `metric`, `estimate`), whose
#'       `method` column holds the weight column name or `"observed"`. It
#'       covers the retained units only, because a weight column containing
#'       any `NA` makes `check_balance()` report `NA` for that whole weight.
#'       `NULL` when `balance = FALSE`. The standardised difference is
#'       halfmoon's, through `smd::smd()`: the weighted mean difference over
#'       the pooled standard deviation of the two arms. `cobalt::bal.tab()`
#'       picks its denominator from the estimand (the treated arm's standard
#'       deviation for ATT, for instance), so its numbers can differ
#'       slightly.}
#'     \item{`fit`}{The `weightit` object the score came from, or `NULL` when
#'       the score was supplied through `ps`. Only its `$ps` is used; the
#'       weights it carries are its own, not the ones in `$data`.}
#'   }
#'   Analysis metadata is attached as `attr(x, "analysis")`.
#'
#' @references
#' Li F, Morgan KL, Zaslavsky AM (2018). Balancing covariates via propensity
#' score weighting. \emph{Journal of the American Statistical Association}
#' 113(521):390-400.
#'
#' Crump RK, Hotz VJ, Imbens GW, Mitnik OA (2009). Dealing with limited
#' overlap in estimation of average treatment effects. \emph{Biometrika}
#' 96(1):187-199.
#'
#' @seealso [plt_PSW()] for the matching plots; [get_PSM()] for the matching
#'   dual; [get_sens()] for sensitivity to unmeasured confounding.
#'
#'   One name differs between the two propensity functions. Here `method`
#'   selects the **propensity model**, because the weighting scheme is what
#'   `estimand` selects. In [get_PSM()] `method` selects the **matching
#'   algorithm**, which is MatchIt's own meaning, and the propensity model
#'   moves to `ps_method`.
#'
#' @examples
#' set.seed(20260921)
#' n <- 400
#' d <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1, 0.4), x3 = runif(n))
#' d$z <- rbinom(n, 1, plogis(-0.3 + 0.8 * d$x1 - 0.6 * d$x2 + 1.1 * d$x3))
#'
#' res <- get_PSW(d, treat = "z", adj_var = c("x1", "x2", "x3"),
#'                balance = FALSE)
#' res
#'
#' # The overlap weight keeps the most information, the ATE weight the least
#' res$stats[, c("estimand", "ess", "ess_pct", "w_max")]
#'
#' @examplesIf requireNamespace("halfmoon", quietly = TRUE)
#' \donttest{
#' set.seed(20260921)
#' n <- 400
#' d <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1, 0.4), x3 = runif(n))
#' d$z <- rbinom(n, 1, plogis(-0.3 + 0.8 * d$x1 - 0.6 * d$x2 + 1.1 * d$x3))
#'
#' # Overlap and matching weights only, with balance diagnostics and
#' # percentile truncation of the score
#' get_PSW(d, treat = "z", adj_var = c("x1", "x2", "x3"),
#'         estimand   = c("ATO", "ATM"),
#'         trunc_args = list(method = "pctl", lower = 0.01, upper = 0.99))
#' }
#'
#' @export
get_PSW <- function(data,
                    treat,
                    adj_var    = NULL,
                    ps         = NULL,
                    estimand   = c("ATE", "ATT", "ATC", "ATO", "ATM", "EW"),
                    method     = "glm",
                    stabilize  = FALSE,
                    trim_args  = list(method = "none", lower = NULL,
                                      upper = NULL, refit = TRUE),
                    trunc_args = list(method = "none", lower = 0.01,
                                      upper = 0.99),
                    balance    = TRUE,
                    ps_args    = list(),
                    verbose    = FALSE) {

  if (!is.character(estimand) || !length(estimand))
    stop("`estimand` must be a character vector.", call. = FALSE)
  estimand <- match.arg(toupper(estimand), .PSW_ESTIMANDS, several.ok = TRUE)
  estimand <- .PSW_ESTIMANDS[.PSW_ESTIMANDS %in% estimand]

  if (!is.data.frame(data) || !nrow(data))
    stop("`data` must be a non-empty data frame.", call. = FALSE)
  for (nm in c("balance", "stabilize", "verbose")) {
    v <- get(nm)
    if (!is.logical(v) || length(v) != 1L || is.na(v))
      stop(sprintf("`%s` must be TRUE or FALSE.", nm), call. = FALSE)
  }

  treat   <- .sens_check_col(treat, data, "treat", n = 1L)
  adj_var <- .sens_check_col(adj_var, data, "adj_var")
  ps      <- .sens_check_col(ps, data, "ps", n = 1L)
  if (is.null(ps) && is.null(adj_var))
    stop("Supply `adj_var` to fit a propensity model, or `ps` to use an existing score.",
         call. = FALSE)
  if (balance && is.null(adj_var))
    stop("`adj_var` is required for the balance table; supply it, or use balance = FALSE.",
         call. = FALSE)

  # Modelling arguments cannot take effect once the score is given. Only an
  # explicitly supplied value is an error: the `method` default has to stay
  # harmless, or every call carrying a `ps` would fail.
  if (!is.null(ps)) {
    given <- c(method = !missing(method), ps_args = !missing(ps_args))
    if (any(given))
      stop(sprintf("`%s` does not apply when `ps` is supplied; the score is taken as given and no model is fitted.",
                   paste(names(given)[given], collapse = "` / `")),
           call. = FALSE)
  } else {
    if (!is.character(method) || length(method) != 1L || is.na(method))
      stop("`method` must be a single string.", call. = FALSE)
    if (!method %in% .PSW_PS_METHODS)
      stop(sprintf("`method` must be one of %s; got \"%s\".",
                   paste0("\"", .PSW_PS_METHODS, "\"", collapse = " / "),
                   method), call. = FALSE)
  }

  trim_args  <- .merge_named_arg(trim_args,  .PSW_TRIM_DEFAULTS,  "trim_args")
  trunc_args <- .merge_named_arg(trunc_args, .PSW_TRUNC_DEFAULTS, "trunc_args")
  trim_args$method  <- match.arg(trim_args$method, c("none", "ps", "pctl", "cr"))
  trunc_args$method <- match.arg(trunc_args$method, c("none", "ps", "pctl"))
  if (!is.list(ps_args) || (length(ps_args) && is.null(names(ps_args))))
    stop("`ps_args` must be a named list.", call. = FALSE)
  bad <- intersect(names(ps_args), .PSW_MANAGED)
  if (length(bad))
    stop(sprintf("`ps_args` may not set %s; %s managed by get_PSW().",
                 paste0("`", bad, "`", collapse = ", "),
                 if (length(bad) > 1L) "they are" else "it is"), call. = FALSE)

  if (isTRUE(stabilize) && !any(estimand %in% .PSW_STABILIZE))
    stop(sprintf("`stabilize = TRUE` has no effect on %s; it is defined for %s only.",
                 paste0("\"", estimand, "\"", collapse = " / "),
                 paste0("\"", .PSW_STABILIZE, "\"", collapse = " / ")),
         call. = FALSE)
  if (!identical(trim_args$method, "none") && isTRUE(trim_args$refit) &&
      !is.null(ps))
    stop("`trim_args$refit = TRUE` needs a propensity model to refit, but the score came from `ps`. Use refit = FALSE.",
         call. = FALSE)
  if (identical(trim_args$method, "cr") &&
      (!is.null(trim_args$lower) || !is.null(trim_args$upper)))
    stop("`trim_args$lower` / `upper` do not apply to method = \"cr\", which computes its own symmetric cut-off. Drop them, or use method = \"ps\" / \"pctl\".",
         call. = FALSE)

  # Trimming keeps its rows, so incomplete cases are the only rows dropped.
  used <- c(treat, adj_var, ps)
  data <- .sens_complete(data, used, verbose)

  wcols <- paste0("w_", tolower(estimand))
  clash <- intersect(c("ps", ".trimmed", wcols), names(data))
  if (length(clash))
    stop(sprintf("`data` already has column(s) %s, which get_PSW() writes. Rename them first.",
                 paste0("`", clash, "`", collapse = ", ")), call. = FALSE)

  tz <- .psw_treat(data[[treat]], treat)
  z  <- tz$z
  if (length(unique(z)) != 2L)
    stop(sprintf("`treat` column `%s` has only one arm after dropping incomplete rows.",
                 treat), call. = FALSE)
  # The model and the balance table see the 0/1 coding, so the direction of
  # the score and the signs in the balance table cannot depend on how the arms
  # were labelled. The column goes back into the returned data as it came in.
  treat_col     <- data[[treat]]
  data[[treat]] <- z

  f  <- .psw_score(data, treat, adj_var, ps, method, ps_args)
  e  <- f$ps
  ft <- f$fit

  tr      <- .psw_trim(e, trim_args)
  e       <- tr$ps
  trimmed <- is.na(e)
  if (all(trimmed))
    stop("Trimming removed every unit; widen `trim_args`.", call. = FALSE)
  if (isTRUE(verbose) && any(trimmed))
    cli::cli_inform(c("i" = paste0(
      "Trimmed {sum(trimmed)} unit{?s} outside [",
      format(round(tr$bounds[1L], 4)), ", ",
      format(round(tr$bounds[2L], 4)), "].")))

  # Trimming redefines the analysis population, so the score that weights it
  # should be estimated on that population rather than carry information from
  # the units just removed.
  if (any(trimmed) && isTRUE(trim_args$refit)) {
    f <- .psw_fit(data[!trimmed, , drop = FALSE], treat, adj_var, method,
                  ps_args)
    e[!trimmed] <- f$ps
    ft <- f$fit
  }

  un <- .psw_trunc(e, trunc_args)
  e  <- un$ps

  keep <- !trimmed
  p1   <- mean(z[keep])
  den  <- z * e + (1 - z) * (1 - e)

  data[["ps"]]       <- e
  data[[".trimmed"]] <- trimmed
  for (i in seq_along(estimand)) {
    w <- .psw_tilt(estimand[[i]], e) / den
    if (isTRUE(stabilize) && estimand[[i]] %in% .PSW_STABILIZE)
      w <- w * (z * p1 + (1 - z) * (1 - p1))
    data[[wcols[[i]]]] <- w
  }

  bal <- if (balance) .psw_balance(data, adj_var, treat, wcols, keep) else NULL
  data[[treat]] <- treat_col

  st <- do.call(rbind, lapply(seq_along(estimand), function(i) {
    m <- .psw_smd(bal, wcols[[i]])
    .psw_stats_row(estimand[[i]], data[[wcols[[i]]]], z, keep,
                   smd_max = m[[1L]], smd_over = m[[2L]])
  }))

  structure(
    list(data = data, stats = st, balance = bal, fit = ft),
    class = c("psw_res", "list"),
    analysis = list(
      treat = treat, treated = tz$treated, adj_var = adj_var, ps = ps,
      estimand = estimand, wcols = wcols,
      method = if (is.null(ps)) method else NA_character_,
      stabilize = stabilize,
      trim = trim_args, trim_bounds = tr$bounds,
      trunc = trunc_args, trunc_bounds = un$bounds,
      n = nrow(data), n_trimmed = sum(trimmed),
      ps_range = range(e, na.rm = TRUE),
      call = match.call()))
}


# ---- L3 print --------------------------------------------------------------

#' @export
#' @noRd
print.psw_res <- function(x, ...) {
  a <- attr(x, "analysis")
  cat(sprintf("<psw_res> n = %d (%d trimmed), treat = %s (treated = %s), score %s\n",
              a$n, a$n_trimmed, a$treat, a$treated,
              if (is.na(a$method)) paste0("supplied (", a$ps, ")")
              else paste0("from method = \"", a$method, "\"")))
  # The processing steps come first and the range last, because a refit can
  # move retained scores outside the trimming window.
  steps <- c(
    if (!is.null(a$trim_bounds))
      paste0("trim ", a$trim$method, " ", .sens_fmt_vec(a$trim_bounds),
             if (isTRUE(a$trim$refit)) " + refit" else ""),
    if (!is.null(a$trunc_bounds))
      paste0("truncate ", a$trunc$method, " ", .sens_fmt_vec(a$trunc_bounds)),
    if (isTRUE(a$stabilize)) "stabilized")
  cat(sprintf("  %sscore range %s\n",
              if (length(steps)) paste0(paste(steps, collapse = ", "),
                                        "; final ") else "",
              .sens_fmt_vec(a$ps_range)))
  cat("\n")
  print(x$stats)
  cat("\n# ESS is (sum w)^2 / sum w^2.",
      "smd_max / smd_over are NA unless balance = TRUE.\n")
  cat("# plt_PSW: type = \"", .psw_plt_spec(x)$type, "\"\n", sep = "")
  invisible(x)
}
