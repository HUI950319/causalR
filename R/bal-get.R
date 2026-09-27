# =============================================================================
# bal-get.R -- covariate balance of matching and weighting schemes, one plot
# =============================================================================
#
# Architecture (2 layers):
#
#   L1  get_bal(data, treat, adj_var, methods, threshold, save)
#         |
#         +-- .bal_specs          shorthand or named list -> validated specs
#         +-- get_PSM / get_PSW   one call per scheme, weight column only
#         +-- cobalt::bal.tab / cobalt::love.plot   one shared denominator
#         +-- .psw_save           pin the size, save through RegR::save_plt()
#
# get_PSM() and get_PSW() report balance through halfmoon, which standardises
# every weight by its own weighted SD. Across schemes that target different
# populations the denominator then moves with the scheme, and part of each
# difference in SMD is a change of scale rather than of means. Here every
# scheme is divided by the unadjusted pooled SD (cobalt's s.d.denom =
# "pooled"), as RegR::get_psm_iptw() does for its three-series plot, so the
# numbers deliberately differ from plt_PSM() / plt_PSW().
# =============================================================================

# Shorthand -> legend label. The weighting shorthands are get_PSW() estimands;
# "PSM" is the common clinical default, nearest 1:1 ATT with a 0.2 caliper.
.BAL_LABELS <- c(PSM = "PS matching (nearest, ATT)",
                 ATE = "IPTW (ATE)",
                 ATT = "SMR weighting (ATT)",
                 ATC = "SMR weighting (ATC)",
                 ATO = "Overlap weighting (ATO)",
                 ATM = "Matching weighting (ATM)",
                 EW  = "Entropy weighting (EW)")

# The first three are RegR's Original / PSM / IPTW shapes; the rest are the
# remaining unfilled shapes, one per colour of UtilsR::pal_lancet (15).
.BAL_SHAPES <- c(17, 16, 15, 18, 8, 3, 4, 7, 9, 10, 11, 12, 13, 14, 6)

# Fields get_bal() supplies to every scheme itself.
.BAL_MANAGED <- c("data", "treat", "adj_var", "balance")


# Turn `methods` into a named list of specs, each a named list with `design`
# plus get_PSM() / get_PSW() arguments. Everything that can be checked
# without fitting is checked here, so a typo in the seventh scheme does not
# wait for the first six to be fitted.
#' @keywords internal
#' @noRd
.bal_specs <- function(methods) {
  if (is.character(methods)) {
    if (!length(methods) || anyNA(methods))
      stop("`methods` must be a non-empty character vector or a named list.",
           call. = FALSE)
    key <- toupper(methods)
    bad <- methods[!key %in% names(.BAL_LABELS)]
    if (length(bad))
      stop(sprintf("Unknown shorthand(s) %s in `methods`; use any of %s, or a named list of scheme specifications.",
                   paste0("\"", bad, "\"", collapse = ", "),
                   paste0("\"", names(.BAL_LABELS), "\"", collapse = ", ")),
           call. = FALSE)
    specs <- lapply(key, function(k)
      if (k == "PSM")
        list(design = "matching", method = "nearest", estimand = "ATT",
             ratio = 1, caliper = 0.2)
      else list(design = "weighting", method = "glm", estimand = k))
    names(specs) <- unname(.BAL_LABELS[key])
  } else if (is.list(methods) && !is.data.frame(methods)) {
    specs <- methods
    nms   <- names(specs)
    if (!length(specs) || is.null(nms) || anyNA(nms) || any(!nzchar(nms)))
      stop("`methods` must be a character vector of shorthands or a fully named list; the names are the legend labels.",
           call. = FALSE)
  } else {
    stop("`methods` must be a character vector of shorthands or a fully named list; the names are the legend labels.",
         call. = FALSE)
  }

  labs <- names(specs)
  dups <- unique(labs[duplicated(c("Unadjusted", labs))[-1L]])
  if (length(dups))
    stop(sprintf("`methods` has duplicated label(s) %s; each label names one legend entry and one weight column, and \"Unadjusted\" is taken by the reference.",
                 paste0("`", dups, "`", collapse = ", ")), call. = FALSE)
  if (length(specs) > length(.BAL_SHAPES) - 1L)
    stop(sprintf("`methods` can hold at most %d schemes, one colour each next to the unadjusted sample; got %d.",
                 length(.BAL_SHAPES) - 1L, length(specs)), call. = FALSE)

  for (nm in labs) {
    s <- specs[[nm]]
    if (!is.list(s) || is.data.frame(s))
      stop(sprintf("Scheme `%s` must be a named list.", nm), call. = FALSE)
    f <- names(s)
    if (length(s) && (is.null(f) || anyNA(f) || any(!nzchar(f)) ||
                      anyDuplicated(f)))
      stop(sprintf("Scheme `%s` must be a list with unique, non-empty field names.",
                   nm), call. = FALSE)
    if (is.null(s$design))
      stop(sprintf("Scheme `%s` needs a `design` field: \"matching\" or \"weighting\".",
                   nm), call. = FALSE)
    if (!is.character(s$design) || length(s$design) != 1L ||
        !s$design %in% c("matching", "weighting"))
      stop(sprintf("Scheme `%s`: `design` must be \"matching\" or \"weighting\".",
                   nm), call. = FALSE)

    fn   <- if (s$design == "matching") "get_PSM" else "get_PSW"
    args <- setdiff(f, "design")
    hit  <- intersect(args, .BAL_MANAGED)
    if (length(hit))
      stop(sprintf("Scheme `%s`: %s %s set by get_bal() for every scheme.",
                   nm, paste0("`", hit, "`", collapse = ", "),
                   if (length(hit) > 1L) "are" else "is"), call. = FALSE)
    if ("matchit_args" %in% args) {
      ma  <- s$matchit_args
      eg  <- if (is.list(ma) && length(ma) && !is.null(names(ma)))
        paste(sprintf("%s = %s", names(ma), vapply(ma, deparse1, "")),
              collapse = ", ") else "ratio = 1, caliper = 0.2"
      stop(sprintf("Scheme `%s`: `matchit_args` is not a field. Write get_PSM()'s own arguments straight into the list (%s), and any other MatchIt argument in `match_args`.",
                   nm, eg), call. = FALSE)
    }
    allowed <- c("design", setdiff(names(formals(fn)), .BAL_MANAGED))
    unknown <- setdiff(args, allowed)
    if (length(unknown))
      stop(sprintf("Scheme `%s` contains unknown field(s) %s. Allowed fields for design = \"%s\" are: %s.",
                   nm, paste0("`", unknown, "`", collapse = ", "), s$design,
                   paste0("`", allowed, "`", collapse = ", ")), call. = FALSE)
    # one scheme is one weight column, so the vector-valued arguments of
    # get_PSM() / get_PSW() are restricted to a single value here
    if (s$design == "matching" && !is.null(s$method) && length(s$method) != 1L)
      stop(sprintf("Scheme `%s`: name one matching `method`; write each algorithm as its own scheme.",
                   nm), call. = FALSE)
    if (s$design == "weighting" &&
        (!is.character(s$estimand) || length(s$estimand) != 1L))
      stop(sprintf("Scheme `%s`: name exactly one `estimand`; write each estimand as its own scheme.",
                   nm), call. = FALSE)
  }
  specs
}


#' Compare covariate balance across matching and weighting schemes
#'
#' Fits any mix of propensity score matching ([get_PSM()]) and weighting
#' ([get_PSW()]) schemes on the same data and draws their covariate balance on
#' a single love plot, so the schemes can be compared on one scale.
#'
#' @section Specifying schemes:
#' `methods` takes either shorthands or a named list.
#'
#' Shorthands expand to these schemes and legend labels:
#' \describe{
#'   \item{`"PSM"`}{`get_PSM(method = "nearest", estimand = "ATT", ratio = 1,
#'     caliper = 0.2)`, "PS matching (nearest, ATT)". The caliper is in
#'     standard deviations of the score on the probability scale, as in
#'     [get_PSM()].}
#'   \item{`"ATE"`}{`get_PSW(estimand = "ATE")` with a logistic score,
#'     "IPTW (ATE)".}
#'   \item{`"ATT"`, `"ATC"`}{"SMR weighting (ATT)" / "(ATC)".}
#'   \item{`"ATO"`, `"ATM"`, `"EW"`}{"Overlap weighting (ATO)",
#'     "Matching weighting (ATM)", "Entropy weighting (EW)".}
#' }
#'
#' A named list gives one scheme per element, any number of them and in any
#' mix. The name is the legend label; the element is a list with
#' `design = "matching"` or `"weighting"` and any arguments of [get_PSM()] or
#' [get_PSW()] respectively, under their own names:
#'
#' ```
#' list(
#'   `PSM 1:2, caliper 0.1` = list(design = "matching", method = "nearest",
#'                                 estimand = "ATT", ratio = 2, caliper = 0.1),
#'   `Stabilized IPTW (ATE)` = list(design = "weighting", method = "glm",
#'                                  estimand = "ATE", stabilize = TRUE))
#' ```
#'
#' `ratio`, `caliper` and `replace` are fields of their own; other MatchIt
#' arguments go in `match_args`. Each scheme names one matching `method` or
#' one weighting `estimand`. `data`, `treat`, `adj_var` and `balance` are set
#' here and rejected in a scheme. Omitted fields take the defaults of
#' [get_PSM()] / [get_PSW()]; note that in a scheme `method` means the
#' matching algorithm for `"matching"` and the propensity model for
#' `"weighting"`, exactly as in those two functions.
#'
#' @section Balance measure:
#' Every scheme, and the unadjusted sample, is summarised by the weighted
#' mean difference between arms divided by the **unadjusted pooled** standard
#' deviation, \eqn{\sqrt{(s_1^2 + s_0^2)/2}}, computed with
#' `cobalt::bal.tab(s.d.denom = "pooled")`; binary covariates use
#' \eqn{p(1-p)} in place of \eqn{s^2}. A common denominator keeps the schemes
#' on one scale even though they target different populations, so these
#' numbers differ from [plt_PSM()] / [plt_PSW()], which standardise each
#' weight by its own weighted standard deviation.
#'
#' Two consequences are exact rather than approximate. Stabilising the ATE
#' weight multiplies it by a constant within each arm, so its balance is
#' identical to the unstabilised weight's. And overlap weights from a
#' logistic score on the same covariates balance every covariate mean
#' exactly.
#'
#' @param data A data frame holding every column named below.
#' @param treat Length-1 character. The binary exposure column, coded as in
#'   [get_PSW()]: `0`/`1`, logical, or a two-level factor or character column
#'   whose second level is the treated one.
#' @param adj_var Character vector of covariates every scheme adjusts for and
#'   the plot reports. A two-level factor or character covariate is shown
#'   under its own name; a factor with more levels gets one row per level.
#' @param methods Character vector of shorthands, any of `"PSM"`, `"ATE"`,
#'   `"ATT"`, `"ATC"`, `"ATO"`, `"ATM"`, `"EW"` (all seven by default), or a
#'   named list of scheme specifications; see *Specifying schemes*. At most
#'   14 schemes.
#' @param threshold Numeric, default `0.1`. Where the red reference line is
#'   drawn.
#' @param save `NULL` or a list with `filename`, and optionally `width` and
#'   `height`, passed to `RegR::save_plt()`. `NULL` and `list()` both mean no
#'   file is written; the defaults come from the plot's own pinned size.
#'
#' @return A list of
#'   \describe{
#'     \item{`plt`}{The love plot, a `ggplot` carrying its pinned size in
#'       `attr(., "plot_size")`.}
#'     \item{`balance`}{A data frame with `variable`, `method` (the scheme
#'       label, or `"Unadjusted"`) and `smd`, the signed standardised mean
#'       difference, treated minus control.}
#'     \item{`data`}{The complete cases of `treat`, `adj_var` and any `ps`
#'       column a scheme names, with one weight column per scheme named by
#'       its label. Unmatched and trimmed units have weight `0`, so every
#'       column describes the same rows.}
#'   }
#'
#' @seealso [get_PSM()] and [get_PSW()] for one design at a time, with
#'   effective sample sizes and weight diagnostics.
#'
#' @examples
#' set.seed(20260927)
#' n <- 400
#' d <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1, 0.4), x3 = runif(n))
#' d$z <- rbinom(n, 1, plogis(-0.9 - 1.1 * d$x1 + 0.8 * d$x2 - 1.2 * d$x3))
#'
#' res <- get_bal(d, treat = "z", adj_var = c("x1", "x2", "x3"),
#'                methods = c("PSM", "ATE", "ATO"))
#' res$plt
#'
#' \donttest{
#' # Several matching schemes next to weighting, each with its own settings
#' methods <- list(
#'   `PSM 1:1, caliper 0.2` = list(design = "matching", method = "nearest",
#'                                 estimand = "ATT", ratio = 1, caliper = 0.2),
#'   `PSM 1:2, caliper 0.1` = list(design = "matching", method = "nearest",
#'                                 estimand = "ATT", ratio = 2, caliper = 0.1),
#'   `Stabilized IPTW (ATE)` = list(design = "weighting", method = "glm",
#'                                  estimand = "ATE", stabilize = TRUE),
#'   `Overlap weighting (ATO)` = list(design = "weighting", method = "glm",
#'                                    estimand = "ATO"))
#' get_bal(d, treat = "z", adj_var = c("x1", "x2", "x3"),
#'         methods = methods)$balance
#' }
#'
#' @export
get_bal <- function(data,
                    treat,
                    adj_var,
                    methods   = c("PSM", "ATE", "ATT", "ATC", "ATO", "ATM",
                                  "EW"),
                    threshold = 0.1,
                    save      = list()) {

  if (!is.data.frame(data) || !nrow(data))
    stop("`data` must be a non-empty data frame.", call. = FALSE)
  treat   <- .sens_check_col(treat, data, "treat", n = 1L)
  adj_var <- .sens_check_col(adj_var, data, "adj_var")
  if (is.null(adj_var))
    stop("`adj_var` is required: it is what every scheme adjusts for and the plot reports.",
         call. = FALSE)
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold))
    stop("`threshold` must be a single number.", call. = FALSE)
  specs <- .bal_specs(methods)
  labs  <- names(specs)

  # Incomplete rows go once, here: left to get_PSM() / get_PSW(), each scheme
  # could drop a different set and the columns would describe different
  # people.
  ps_cols <- unique(unlist(lapply(specs, `[[`, "ps")))
  ps_cols <- .sens_check_col(ps_cols, data, "ps")
  data    <- .sens_complete(data, c(treat, adj_var, ps_cols))
  hit <- intersect(labs, names(data))
  if (length(hit))
    stop(sprintf("Column(s) %s already exist in `data`; get_bal() writes one weight column per scheme label. Rename them, or relabel the scheme.",
                 paste0("`", hit, "`", collapse = ", ")), call. = FALSE)

  w <- lapply(labs, function(nm) {
    s  <- specs[[nm]]
    fn <- if (s$design == "matching") get_PSM else get_PSW
    r  <- tryCatch(
      do.call(fn, c(list(data = data, treat = treat, adj_var = adj_var,
                         balance = FALSE),
                    s[names(s) != "design"])),
      error = function(e)
        stop(sprintf("Scheme `%s`: %s", nm, conditionMessage(e)),
             call. = FALSE))
    # NA marks a trimmed unit; like an unmatched one it is out of the sample
    x <- r$data[[attr(r, "analysis")$wcols]]
    x[is.na(x)] <- 0
    x
  })
  names(w) <- labs
  w <- as.data.frame(w, check.names = FALSE)

  z  <- .psw_treat(data[[treat]], treat)$z
  bt <- cobalt::bal.tab(
    data[adj_var], treat = z, weights = w,
    method = unname(vapply(specs, `[[`, "", "design")),
    s.d.denom = "pooled", binary = "std", continuous = "std", un = TRUE)

  # cobalt writes a two-level factor as "<var>_<second level>"; the plot and
  # the table name it after the variable, as RegR::get_ps() does.
  two <- adj_var[vapply(data[adj_var], function(v)
    (is.factor(v) || is.character(v)) && length(unique(v)) == 2L,
    logical(1))]
  map <- stats::setNames(two, vapply(two, function(v)
    paste0(v, "_", levels(factor(data[[v]]))[2L]), ""))
  map <- map[names(map) %in% rownames(bt$Balance)]

  B    <- bt$Balance
  vars <- rownames(B)
  vars[vars %in% names(map)] <- unname(map[vars[vars %in% names(map)]])
  bal <- data.frame(
    variable = rep(vars, length(labs) + 1L),
    method   = rep(c("Unadjusted", labs), each = nrow(B)),
    # "Diff.Un" then one column per scheme, in order; cobalt calls a lone
    # scheme "Diff.Adj" rather than by its name
    smd      = unlist(B[grep("^Diff\\.", names(B))], use.names = FALSE),
    stringsAsFactors = FALSE)

  k <- length(labs) + 1L
  p <- cobalt::love.plot(
    bt, stats = "mean.diffs", abs = TRUE, var.order = "unadjusted",
    line = TRUE, thresholds = c(m = threshold),
    colors = UtilsR::pal_lancet[seq_len(k)], shapes = .BAL_SHAPES[seq_len(k)],
    size = 3.5, drop.distance = TRUE, labels = FALSE,
    sample.names = c("Unadjusted", labs),
    var.names = if (length(map)) map, title = NULL)
  p <- p + UtilsR::theme_my(base_size = 16)
  p <- UtilsR::fmt_ref(p, x = threshold, color = "red")
  # RegR's framed inside legend, but anchored bottom-right: rows are sorted
  # by unadjusted SMD, so that corner is the one the reference line leaves
  # empty, and with many schemes RegR's c(0.8, 0.4) sits on top of it
  p <- UtilsR::fmt_legend(
    p, legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
    legend.background = ggplot2::element_rect(
      fill = "white", colour = "#D6D6D6", linewidth = 1))
  p <- p + ggplot2::labs(x = "Absolute standardized mean difference",
                         y = "Covariates")

  data[labs] <- w
  plt <- .psw_save(p, c(8, max(5, 1.5 + 0.45 * nrow(B))), save)
  list(plt = plt, balance = bal, data = data)
}
