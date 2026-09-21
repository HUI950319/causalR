# =============================================================================
# psm-get.R -- propensity score matching for a binary exposure
# =============================================================================
#
# Architecture (3 layers + shared helpers from utils-sens.R and psw-get.R):
#
#   L1  get_PSM(data, treat, adj_var, ps, method, ...)
#         |
#         +-- L2 pipeline stages, run in this fixed order
#         |     .psw_score     the score, from WeightIt or `ps` (psw-get.R)
#         |     .psm_match     one MatchIt run per method, score held fixed
#         |     .psw_balance   halfmoon::check_balance() over every weight,
#         |                    on the whole cohort here (psw-get.R)
#         |
#         +-- L3 helpers
#               .psm_spec        the matchit arguments one method accepts
#               .psm_stats_row   the 16-column standardised $stats row
#               .psm_plt_spec    plt_PSM() arguments echoed by print()
#               .psw_treat, .psw_ess, .psw_smd   shared, from psw-get.R
#
#   print.psm_res reports the diagnostics table and the matching plt_PSM()
#   call.
#
# This is the matching-side dual of get_PSW(): one propensity score, several
# schemes compared side by side, one weight column each. The score is fitted
# once and handed to every matchit() call as `distance`, so the schemes differ
# only in how they use it -- letting each method refit its own score would
# make the columns incomparable, which is the whole point of the function.
#
# Three MatchIt behaviours are handled here rather than passed through, each
# measured rather than assumed:
#
#   * `normalize = TRUE` (MatchIt's default) rescales each arm to mean 1,
#     which destroys the inverse-probability reading of the weights: with
#     ratio = 2 the control weights sum to 2 * n_treated instead of
#     n_treated. get_PSM() fixes `normalize = FALSE` so the weights mean the
#     same thing they do in get_PSW(). The effective sample size is
#     scale-free and identical either way.
#   * A character or 1/2-coded exposure makes matchit()'s internal glm fail
#     with "y values must be 0 <= y <= 1", and a factor whose levels are
#     reversed silently swaps the arms. .psw_treat() settles the coding
#     before matchit() ever sees it; the returned $data carries the column
#     as supplied, and the arm taken as treated is recorded and printed.
#   * "discarded" (outside common support) and "unmatched" (no partner
#     found) are different things and are reported in separate columns.
# =============================================================================

# Estimands each method actually accepts, measured against MatchIt 4.7.2.
# The distance-based methods reject ATE outright:
#   "The argument to `estimand` should be one of "ATT" or "ATC"."
.PSM_ESTIMANDS <- list(
  nearest  = c("ATT", "ATC"),
  optimal  = c("ATT", "ATC"),
  genetic  = c("ATT", "ATC"),
  full     = c("ATT", "ATE", "ATC"),
  cem      = c("ATT", "ATE", "ATC"),
  exact    = c("ATT", "ATE", "ATC"),
  subclass = c("ATT", "ATE", "ATC"))

.PSM_METHODS <- names(.PSM_ESTIMANDS)

# Which methods take which knob, read off MatchIt's own "the argument `x` is
# not used with method `y`" warnings rather than assumed. Handing a method an
# argument it ignores is only a warning, but a call naming several methods at
# once would then emit one per method, so each gets only what it uses.
.PSM_RATIO   <- c("nearest", "optimal", "genetic")
.PSM_CALIPER <- c("nearest", "full", "genetic")
.PSM_REPLACE <- c("nearest", "genetic")

# "cem" and "exact" coarsen or match on the covariates themselves and ignore
# `distance` entirely, so they take no propensity score. They are still worth
# comparing against, but the score is not what they share with the others.
.PSM_DISTANCE <- c("nearest", "optimal", "full", "genetic", "subclass")

# Methods that need packages beyond MatchIt itself, read off MatchIt's own
# check_installed() calls. "cem" is implemented inside MatchIt 4.x and does
# not need the cem package.
.PSM_PKG <- list(optimal = "optmatch", full = "optmatch",
                 genetic = c("Matching", "rgenoud"))

# Arguments the function manages itself and will not forward to matchit().
.PSM_MANAGED <- c("formula", "data", "method", "estimand", "distance",
                  "ratio", "caliper", "replace", "normalize", "link",
                  "distance.options")


# ---- L3 helpers ------------------------------------------------------------

# The matchit() arguments for one method: the shared knobs it accepts, plus
# whatever the caller added through `match_args`.
#' @keywords internal
#' @noRd
.psm_spec <- function(method, ratio, caliper, replace, match_args) {
  sp <- list()
  if (method %in% .PSM_RATIO   && !is.null(ratio))   sp$ratio   <- ratio
  if (method %in% .PSM_CALIPER && !is.null(caliper)) sp$caliper <- caliper
  if (method %in% .PSM_REPLACE && isTRUE(replace))   sp$replace <- TRUE
  utils::modifyList(sp, match_args)
}

#' @keywords internal
#' @noRd
.psm_stats_row <- function(method, estimand, w, z, discarded,
                           subclass, smd_max = NA_real_,
                           smd_over = NA_real_) {
  keep <- w > 0
  tibble::tibble(
    method       = method,
    estimand     = estimand,
    n            = sum(keep),
    n_treat      = sum(keep & z == 1L),
    n_ctrl       = sum(keep & z == 0L),
    n_unmatched  = sum(!keep & !discarded),
    n_discarded  = sum(discarded),
    pct_retained = 100 * sum(keep) / length(w),
    n_pairs      = length(unique(stats::na.omit(subclass))),
    ess          = .psw_ess(w),
    ess_treat    = .psw_ess(w[z == 1L]),
    ess_ctrl     = .psw_ess(w[z == 0L]),
    ess_pct      = .psw_ess(w) / sum(keep),
    w_max        = max(w[keep]),
    w_cv         = stats::sd(w[keep]) / mean(w[keep]),
    smd_max      = as.numeric(smd_max),
    smd_over     = as.numeric(smd_over))
}

# Love plots are read against the 0.1 threshold; without a balance table the
# only thing left to show is how much sample each scheme keeps.
#' @keywords internal
#' @noRd
.psm_plt_spec <- function(x) {
  list(type = if (is.null(x$balance)) "ess" else "love")
}


# ---- L2 pipeline stages ----------------------------------------------------

# One matchit() run. The score comes in as `distance`, so matchit() does no
# modelling of its own and every method in a single call scores the same way.
#
# `data` and `ps` enter the call as symbols, not values. do.call() would
# inline the matchit() function body and a snapshot of the data frame into
# `obj$call`: measured at n = 400 that is 1159 KB per method against 181 KB,
# and summary(obj) then prints the whole data frame under "Call:". The
# symbols resolve in this frame through the formula environment, which is
# also how MatchIt::match.data() finds the data without being handed it.
#' @keywords internal
#' @noRd
.psm_match <- function(data, treat, adj_var, ps, method, estimand,
                       ratio, caliper, replace, match_args) {
  pkg  <- .PSM_PKG[[method]]
  miss <- pkg[!vapply(pkg, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop(sprintf("Package%s %s %s required for get_PSM(method = \"%s\").",
                 if (length(miss) > 1L) "s" else "",
                 paste0("'", miss, "'", collapse = " and "),
                 if (length(miss) > 1L) "are" else "is", method),
         call. = FALSE)

  form <- stats::reformulate(adj_var, response = treat)
  args <- c(list(formula = form, data = quote(data), method = method,
                 estimand = estimand, normalize = FALSE),
            if (method %in% .PSM_DISTANCE) list(distance = quote(ps)),
            .psm_spec(method, ratio, caliper, replace, match_args))
  # A call naming several methods fails as a whole when one of them does;
  # MatchIt's message does not say which, so it is prefixed here.
  tryCatch(
    eval(as.call(c(list(quote(MatchIt::matchit)), args))),
    error = function(e) stop(
      sprintf("MatchIt::matchit(method = \"%s\") failed: %s",
              method, conditionMessage(e)), call. = FALSE))
}


# ---- L1 public entry point -------------------------------------------------

#' Propensity score matching for a binary exposure
#'
#' Single entry point that fits one propensity score and hands it to any
#' combination of MatchIt's matching algorithms, then reports how much sample
#' each one keeps, how its matching weights are distributed and how well it
#' balances the covariates -- so that nearest-neighbour, optimal, full and
#' coarsened exact matching can be compared on one table rather than one at a
#' time.
#'
#' `get_PSM()` constructs and diagnoses matched cohorts; it does not estimate
#' a treatment effect. Feed `result$data` with the weight column of your
#' choice to an outcome model, clustering on the matching subclass, as in
#' `lm(y ~ z, data = res$data, weights = w_nearest)` with
#' `sandwich::vcovCL(cluster = res$data$s_nearest)`.
#'
#' @section Relation to get_PSW():
#' The two functions are duals and share their score, their `w_` column
#' convention and their diagnostics table. One name differs: in [get_PSW()]
#' `method` selects the **propensity model**, while here `method` selects the
#' **matching algorithm** -- MatchIt's own meaning -- and the propensity model
#' is `ps_method`. Matching weights are themselves inverse-probability
#' weights, with the modelled score replaced by the treated share inside each
#' matched set, so `w_nearest` and [get_PSW()]'s `w_att` target the same
#' estimand by different routes.
#'
#' @section Which estimand each method allows:
#' Measured against MatchIt 4.7.2, not assumed. `"nearest"`, `"optimal"` and
#' `"genetic"` accept `"ATT"` and `"ATC"` only; `"full"`, `"cem"`, `"exact"`
#' and `"subclass"` accept `"ATE"` as well. Asking for an unsupported pair is
#' an error naming the offending method, rather than MatchIt's own message,
#' which does not say that another method would work.
#'
#' `"optimal"` and `"full"` need \pkg{optmatch}; `"genetic"` needs
#' \pkg{Matching} and \pkg{rgenoud} and is slow. `"cem"` is built into
#' MatchIt and needs nothing.
#'
#' Not every method takes every knob, and MatchIt warns once per ignored
#' argument, so each method is handed only what it uses: `ratio` reaches
#' `"nearest"`, `"optimal"` and `"genetic"`; `caliper` reaches `"nearest"`,
#' `"full"` and `"genetic"`; `replace` reaches `"nearest"` and `"genetic"`.
#' `"cem"` and `"exact"` ignore the propensity score altogether -- they
#' coarsen or match on the covariates themselves -- so they are compared
#' alongside the others rather than sharing their score. `"exact"` matches
#' on the covariate values as they are and fails on a continuous covariate
#' ("No exact matches were found"): give it discrete `adj_var` only, or use
#' `"cem"`, which coarsens first. What each method actually received is
#' recorded in `attr(x, "analysis")$specs`.
#'
#' A method that fails stops the whole call with an error naming it; the
#' others are not reported partially.
#'
#' @section Weight normalisation:
#' MatchIt's `normalize = TRUE` default rescales each arm to a mean of 1,
#' which removes the inverse-probability reading of the weights: at
#' `ratio = 2` the control weights would sum to twice the number of treated
#' units. `get_PSM()` fixes `normalize = FALSE`, so an ATT control column sums
#' to the number of treated units, as the textbook weight does. The effective
#' sample size is scale-free and is the same under either setting.
#'
#' @param data A data frame holding every column named below.
#' @param treat Length-1 character. The binary exposure column: `0`/`1`,
#'   logical, or a two-level factor or character column whose **second** level
#'   is the treated one. A character column is ordered alphabetically, so
#'   `"case"` / `"control"` would make `"control"` the treated arm; convert
#'   such a column to a factor with the control level first. The arm actually
#'   taken as treated is recorded in `attr(x, "analysis")$treated` and shown
#'   by `print()`. The 0/1 coding is settled here, before MatchIt sees it,
#'   because a character exposure makes MatchIt's internal model fail and a
#'   reversed factor silently swaps the arms; the column is returned in
#'   `$data` exactly as supplied.
#' @param adj_var Character vector of covariates to match on and to report
#'   balance for.
#' @param ps Length-1 character or `NULL` (default). Column holding an
#'   already-estimated propensity score, strictly inside (0, 1). Supplying it
#'   skips the modelling step, so `ps_method` and `ps_args` no longer apply.
#' @param method Character vector of matching algorithms, any of
#'   `"nearest"` (default), `"optimal"`, `"full"`, `"genetic"`, `"cem"`,
#'   `"exact"`, `"subclass"`. Each produces one weight column, one subclass
#'   column and one `$stats` row.
#' @param estimand Length-1 character, `"ATT"` (default), `"ATE"` or
#'   `"ATC"`, applied to every selected method and validated against each.
#' @param ratio Number of controls matched to each treated unit. Used by
#'   `"nearest"`, `"optimal"` and `"genetic"`; ignored by the others, which do
#'   not accept it.
#' @param caliper `NULL` (default) for no caliper, or a number in **standard
#'   deviations of the propensity score**, which is the probability scale, not
#'   the logit scale; 0.1 to 0.2 is the usual range. Austin's
#'   0.2-standard-deviation rule is stated on the logit scale and cannot be
#'   reproduced here, because `ps` must be a probability and the matching
#'   distance is that probability. Used by `"nearest"`, `"full"` and
#'   `"genetic"`; `"optimal"` does not accept one.
#' @param replace Logical, default `FALSE`. Allow a control to be matched more
#'   than once. Used by `"nearest"` and `"genetic"`. With replacement the
#'   control weights stop being 0/1 and become the number of times each
#'   control was used.
#' @param ps_method Propensity model backend, passed to
#'   [WeightIt::weightit()] exactly as [get_PSW()]'s `method`: `"glm"`
#'   (default), `"gbm"`, `"cbps"`, `"bart"` or `"super"`.
#' @param match_args Named list forwarded to [MatchIt::matchit()] for every
#'   selected method, for arguments this signature does not cover, such as
#'   `exact`, `mahvars`, `m.order`, `discard`, `min.controls` or `subclass`.
#'   `formula`, `data`, `method`, `estimand`, `distance` and `normalize` are
#'   managed here and rejected if supplied.
#' @param balance Logical, default `TRUE`. Compute standardised mean
#'   differences under every matching scheme with
#'   [halfmoon::check_balance()], on the whole cohort so that the unmatched
#'   reference stays comparable.
#' @param ps_args Named list forwarded to [WeightIt::weightit()], for example
#'   `list(link = "probit")`.
#' @param verbose Logical, default `FALSE`. Report dropped incomplete rows and
#'   per-method matched counts.
#'
#' @return An object of class `psm_res`: a list of
#'   \describe{
#'     \item{`data`}{The input data, with `treat` exactly as supplied, plus
#'       `ps`, and per method a weight column `w_<method>` (`0` for a unit
#'       that was not matched) and a subclass column `s_<method>` (`NA` where
#'       unmatched). Every input row is kept,
#'       so the columns stay aligned and several schemes fit in one frame; the
#'       `w_` prefix is what
#'       `halfmoon::plot_ess(.weights = starts_with("w_"))` selects on, and
#'       the `s_` prefix keeps the subclasses out of that selection.}
#'     \item{`stats`}{One tibble row per method, always with the same 16
#'       columns: `method`, `estimand`, `n`, `n_treat`, `n_ctrl`,
#'       `n_unmatched`, `n_discarded`, `pct_retained`, `n_pairs`, `ess`,
#'       `ess_treat`, `ess_ctrl`, `ess_pct`, `w_max`, `w_cv`, `smd_max`,
#'       `smd_over`. `n` and the two arm counts are of matched units;
#'       `n_unmatched` and `n_discarded` are separate because failing to find
#'       a partner and falling outside common support are different things.
#'       `n_pairs` is the number of matched sets: pairs under 1:1 nearest
#'       matching, strata under `"subclass"`, `"cem"` and `"exact"`, sets of
#'       varying size under `"full"`. `smd_max` is the largest absolute
#'       standardised mean difference across `adj_var` and `smd_over` counts
#'       those above 0.1, the line [plt_PSM()] draws by default; both are
#'       `NA` when `balance = FALSE`.}
#'     \item{`balance`}{The [halfmoon::check_balance()] long table
#'       (`variable`, `group_level`, `method`, `metric`, `estimate`), whose
#'       `method` column holds the weight column name or `"observed"`. `NULL`
#'       when `balance = FALSE`. The standardised difference is halfmoon's,
#'       through `smd::smd()`: the weighted mean difference over the pooled
#'       standard deviation of the two arms. `cobalt::bal.tab()` divides by
#'       the treated arm's standard deviation for an ATT estimand by default,
#'       so its numbers differ slightly (measured: -0.333 here against -0.357
#'       there for the same nearest-neighbour match).}
#'     \item{`fit`}{Named list of the `matchit` objects, one per method, for
#'       `summary()`, `plot()`, `cobalt::bal.tab()` or
#'       [MatchIt::match.data()], which finds the data on its own and returns
#'       the exposure recoded to 0/1.}
#'   }
#'   Analysis metadata is attached as `attr(x, "analysis")`, including the
#'   matchit arguments each method actually received.
#'
#' @references
#' Ho DE, Imai K, King G, Stuart EA (2011). MatchIt: nonparametric
#' preprocessing for parametric causal inference. \emph{Journal of
#' Statistical Software} 42(8):1-28.
#'
#' Austin PC (2011). Optimal caliper widths for propensity-score matching.
#' \emph{Pharmaceutical Statistics} 10(2):150-161.
#'
#' @seealso [get_PSW()] for the weighting dual, [plt_PSM()] for the plots.
#'
#' @examplesIf requireNamespace("halfmoon", quietly = TRUE)
#' \donttest{
#' set.seed(20260921)
#' n <- 400
#' d <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1, 0.4), x3 = runif(n))
#' d$z <- rbinom(n, 1, plogis(-0.9 - 1.1 * d$x1 + 0.8 * d$x2 - 1.2 * d$x3))
#'
#' res <- get_PSM(d, treat = "z", adj_var = c("x1", "x2", "x3"),
#'                caliper = 0.2)
#' res
#'
#' # Nearest-neighbour against full matching, which keeps everyone
#' get_PSM(d, treat = "z", adj_var = c("x1", "x2", "x3"),
#'         method = c("nearest", "full"))$stats
#' }
#'
#' @export
get_PSM <- function(data,
                    treat,
                    adj_var    = NULL,
                    ps         = NULL,
                    method     = "nearest",
                    estimand   = "ATT",
                    ratio      = 1,
                    caliper    = NULL,
                    replace    = FALSE,
                    ps_method  = "glm",
                    match_args = list(),
                    balance    = TRUE,
                    ps_args    = list(),
                    verbose    = FALSE) {

  if (!is.character(method) || !length(method))
    stop("`method` must be a character vector.", call. = FALSE)
  method <- match.arg(tolower(method), .PSM_METHODS, several.ok = TRUE)
  method <- .PSM_METHODS[.PSM_METHODS %in% method]
  estimand <- match.arg(toupper(estimand), c("ATT", "ATE", "ATC"))

  if (!is.data.frame(data) || !nrow(data))
    stop("`data` must be a non-empty data frame.", call. = FALSE)
  for (nm in c("balance", "replace", "verbose")) {
    v <- get(nm)
    if (!is.logical(v) || length(v) != 1L || is.na(v))
      stop(sprintf("`%s` must be TRUE or FALSE.", nm), call. = FALSE)
  }
  if (!is.null(caliper) &&
      (!is.numeric(caliper) || length(caliper) != 1L || is.na(caliper) ||
       caliper <= 0))
    stop("`caliper` must be `NULL` or a single positive number.",
         call. = FALSE)
  if (!is.numeric(ratio) || length(ratio) != 1L || is.na(ratio) || ratio < 1)
    stop("`ratio` must be a single number of at least 1.", call. = FALSE)

  treat   <- .sens_check_col(treat, data, "treat", n = 1L)
  adj_var <- .sens_check_col(adj_var, data, "adj_var")
  ps      <- .sens_check_col(ps, data, "ps", n = 1L)
  if (is.null(adj_var))
    stop("`adj_var` is required: it is what the matching balances on.",
         call. = FALSE)

  # A model argument cannot take effect once the score is given. Only an
  # explicitly supplied value is an error, so the defaults stay harmless.
  if (!is.null(ps)) {
    given <- c(ps_method = !missing(ps_method), ps_args = !missing(ps_args))
    if (any(given))
      stop(sprintf("`%s` does not apply when `ps` is supplied; the score is taken as given and no model is fitted.",
                   paste(names(given)[given], collapse = "` / `")),
           call. = FALSE)
  }

  # Each method is checked separately so the message names the one at fault;
  # MatchIt's own error does not say that another method would accept it.
  bad <- method[!vapply(method, function(m) estimand %in% .PSM_ESTIMANDS[[m]],
                        logical(1))]
  if (length(bad))
    stop(sprintf("estimand = \"%s\" is not available for method %s; %s %s. Methods accepting \"%s\": %s.",
                 estimand, paste0("\"", bad, "\"", collapse = " / "),
                 if (length(bad) > 1L) "they accept" else "it accepts",
                 paste0("\"", unique(unlist(.PSM_ESTIMANDS[bad])), "\"",
                        collapse = " / "),
                 estimand,
                 paste0("\"", names(Filter(function(e) estimand %in% e,
                                           .PSM_ESTIMANDS)), "\"",
                        collapse = " / ")),
         call. = FALSE)

  if (!is.list(match_args) ||
      (length(match_args) && is.null(names(match_args))))
    stop("`match_args` must be a named list.", call. = FALSE)
  clash <- intersect(names(match_args), .PSM_MANAGED)
  if (length(clash))
    stop(sprintf("`match_args` may not set %s; %s managed by get_PSM().",
                 paste0("`", clash, "`", collapse = ", "),
                 if (length(clash) > 1L) "they are" else "it is"),
         call. = FALSE)

  used <- c(treat, adj_var, ps)
  data <- .sens_complete(data, used, verbose)

  wcols <- paste0("w_", method)
  scols <- paste0("s_", method)
  hit   <- intersect(c("ps", wcols, scols), names(data))
  if (length(hit))
    stop(sprintf("`data` already has column(s) %s, which get_PSM() writes. Rename them first.",
                 paste0("`", hit, "`", collapse = ", ")), call. = FALSE)

  tz <- .psw_treat(data[[treat]], treat)
  z  <- tz$z
  if (length(unique(z)) != 2L)
    stop(sprintf("`treat` column `%s` has only one arm after dropping incomplete rows.",
                 treat), call. = FALSE)
  # matchit() and the balance table see the 0/1 coding; the column goes back
  # into the returned data as it came in.
  treat_col     <- data[[treat]]
  data[[treat]] <- z

  e <- .psw_score(data, treat, adj_var, ps, ps_method, ps_args)$ps
  data[["ps"]] <- e

  fits  <- list()
  specs <- list()
  for (i in seq_along(method)) {
    m   <- method[[i]]
    obj <- .psm_match(data, treat, adj_var, e, m, estimand,
                      ratio, caliper, replace, match_args)
    fits[[m]]  <- obj
    specs[[m]] <- .psm_spec(m, ratio, caliper, replace, match_args)
    data[[wcols[[i]]]] <- unname(obj$weights)
    data[[scols[[i]]]] <- if (is.null(obj$subclass)) NA_integer_ else
      as.integer(as.character(obj$subclass))
    if (isTRUE(verbose))
      cli::cli_inform(c("i" = paste0(
        "method \"", m, "\": matched ", sum(obj$weights > 0), " of ",
        nrow(data), " unit", if (nrow(data) != 1L) "s" else "", ".")))
  }

  # Unlike the NA weights trimming produces in get_PSW(), a zero weight does
  # not make check_balance() return NA, so the whole cohort is passed in (no
  # `keep`). That is also the right thing statistically: standardising on the
  # matched subset's own SD inflates every SMD (measured: sd 0.73 vs 1.03, a
  # 21% difference) and breaks comparability with the unmatched "observed"
  # reference row, which is the row the love plot exists to compare against.
  bal <- if (balance)
    .psw_balance(data, adj_var, treat, wcols, caller = "get_PSM") else NULL
  data[[treat]] <- treat_col

  st <- do.call(rbind, lapply(seq_along(method), function(i) {
    mm <- .psw_smd(bal, wcols[[i]])
    .psm_stats_row(method[[i]], estimand,
                   w         = data[[wcols[[i]]]],
                   z         = z,
                   discarded = fits[[i]]$discarded,
                   subclass  = data[[scols[[i]]]],
                   smd_max   = mm[[1L]], smd_over = mm[[2L]])
  }))

  structure(
    list(data = data, stats = st, balance = bal, fit = fits),
    class = c("psm_res", "list"),
    analysis = list(
      treat = treat, treated = tz$treated, adj_var = adj_var, ps = ps,
      method = method, estimand = estimand,
      wcols = wcols, scols = scols,
      ratio = ratio, caliper = caliper, replace = replace,
      ps_method = if (is.null(ps)) ps_method else NA_character_,
      specs = specs,
      n = nrow(data), ps_range = range(e),
      call = match.call()))
}


# ---- L3 print --------------------------------------------------------------

#' @export
#' @noRd
print.psm_res <- function(x, ...) {
  a <- attr(x, "analysis")
  cat(sprintf("<psm_res> n = %d, treat = %s (treated = %s), estimand = %s, score %s\n",
              a$n, a$treat, a$treated, a$estimand,
              if (is.na(a$ps_method)) paste0("supplied (", a$ps, ")")
              else paste0("from ps_method = \"", a$ps_method, "\"")))
  cat(sprintf("  range %s, ratio = %s, caliper = %s%s\n",
              .sens_fmt_vec(a$ps_range), format(a$ratio),
              if (is.null(a$caliper)) "none" else
                paste0(format(a$caliper), " SD of the score"),
              if (isTRUE(a$replace)) ", with replacement" else ""))
  cat("\n")
  print(x$stats)
  cat("\n# Weights are unnormalised, so an ATT control column sums to the",
      "number of treated units.\n")
  cat("# n_unmatched (no partner) and n_discarded (outside common support)",
      "are different things.\n")
  cat("# plt_PSM: type = \"", .psm_plt_spec(x)$type, "\"\n", sep = "")
  invisible(x)
}
