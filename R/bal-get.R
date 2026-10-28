# =============================================================================
# bal-get.R -- covariate balance of matching and weighting schemes, one plot
# =============================================================================
#
# Architecture (2 layers):
#
#   L1  get_bal(data, treat, adj_var, methods, cat_smd, tbl, cores, var_names,
#               love_args, save_plt, save_tbl)
#         |
#         +-- .bal_specs          shorthand or named list -> validated specs
#         +-- get_PSM / get_PSW   one call per matching scheme, one per group
#         |                       of weighting schemes differing in estimand
#         +-- cobalt::bal.tab / cobalt::love.plot   one shared denominator
#         +-- gtsummary::tbl_merge   tbl = TRUE, one svysummary per scheme,
#         |                          built on `cores` workers
#         +-- .psw_save           pin the size, save through RegR::save_plt()
#         +-- RegR::save_tb       save_tbl, whose fields are .BAL_TB_SAVE
#
# Every scheme is divided by the unadjusted pooled SD of the complete cases
# (cobalt's s.d.denom = "pooled"), as RegR::get_psm_iptw() does for its
# three-series plot, so that schemes targeting different populations stay on
# one scale. get_PSM() and get_PSW() report balance through halfmoon, whose
# smd::smd() also divides by the unweighted pooled SD, but with the variance
# over n rather than n - 1 and, for get_PSW(), on the rows left after
# trimming. A factor of three or more levels is one Yang & Dalton (2012) row
# by default (cat_smd = "overall"), spliced into cobalt's table in place of
# its level rows; halfmoon and cobalt give one row per level.
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

# love_args defaults. NULL colours / shapes fall back to UtilsR::pal_lancet /
# .BAL_SHAPES once the number of schemes is known; a NULL justification is
# left to UtilsR::fmt_legend(), which derives it from the position.
.BAL_LOVE_DEFAULTS <- list(threshold = 0.1, colors = NULL, shapes = NULL,
                           size = 3.5, line = TRUE, var_order = "unadjusted",
                           base_size = 16, ref_color = "red",
                           legend_position = c(0.99, 0.02),
                           legend_justification = NULL)

# `save_tbl` fields, the arguments of RegR::save_tb() other than `data`.
.BAL_TB_SAVE <- c("path", "title", "note", "header_value", "header_colwidths",
                  "line_spacing")


# Turn `methods` into a named list of specs, each a named list with `design`
# plus get_PSM() / get_PSW() arguments. Everything that can be checked
# without fitting is checked here, so a typo in the seventh scheme does not
# wait for the first six to be fitted.
#' @keywords internal
#' @noRd
.bal_specs <- function(methods) {
  shorthand <- function(x) {
    key <- toupper(x)
    bad <- x[!key %in% names(.BAL_LABELS)]
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
    specs
  }
  is_short <- function(x) is.character(x) && length(x) && !anyNA(x)

  if (is.character(methods)) {
    if (!is_short(methods))
      stop("`methods` must be a non-empty character vector or a named list.",
           call. = FALSE)
    specs <- shorthand(methods)
  } else if (is.list(methods) && !is.data.frame(methods)) {
    # a named element is a specification labelled by its name, an unnamed
    # one a shorthand that brings its own label
    nms <- names(methods)
    if (is.null(nms)) nms <- rep("", length(methods))
    nms[is.na(nms)] <- ""
    short <- !nzchar(nms) & vapply(methods, is_short, NA)
    if (!length(methods) || any(!nzchar(nms) & !short))
      stop("`methods` must be a character vector of shorthands or a list whose elements are named lists of scheme specifications, the names being the legend labels, or unnamed shorthands.",
           call. = FALSE)
    specs <- unlist(lapply(seq_along(methods), function(i)
      if (short[i]) shorthand(methods[[i]]) else methods[i]),
      recursive = FALSE)
  } else {
    stop("`methods` must be a character vector of shorthands or a list whose elements are named lists of scheme specifications, the names being the legend labels, or unnamed shorthands.",
         call. = FALSE)
  }

  labs <- names(specs)
  dups <- unique(labs[duplicated(c("Unadjusted", labs))[-1L]])
  if (length(dups))
    stop(sprintf("`methods` has duplicated label(s) %s; each label names one legend entry and one weight column, and \"Unadjusted\" is taken by the reference.",
                 paste0("`", dups, "`", collapse = ", ")), call. = FALSE)
  if ("Un" %in% labs)
    stop("`methods` label `Un` is reserved: cobalt calls the unadjusted column \"Diff.Un\", which a scheme of that name would overwrite. Relabel the scheme.",
         call. = FALSE)
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
#' An unnamed element of that list is a shorthand, expanded and labelled as
#' above, so tuned schemes and defaults can share one call:
#'
#' ```
#' list(`PSM 1:2` = list(design = "matching", ratio = 2), "ATE", "ATO")
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
#' on one scale even though they target different populations. [plt_PSM()] /
#' [plt_PSW()] standardise by the same unweighted pooled standard deviation
#' (through halfmoon and `smd::smd()`), so for continuous and binary
#' covariates the two agree, except that halfmoon divides the variance by
#' \eqn{n} rather than \eqn{n - 1} and that [get_PSW()] computes it on the
#' rows left after trimming. They report each level of a factor, as
#' `cat_smd = "level"` does here.
#'
#' A factor or character covariate with three or more levels is, under
#' `cat_smd = "overall"`, summarised by the multivariate SMD of Yang & Dalton
#' (2012), \eqn{\sqrt{D^\top S^{+} D}}: \eqn{D} is the vector of weighted
#' level proportions, treated minus control, and \eqn{S} the mean of the two
#' arms' unweighted multinomial covariance matrices
#' \eqn{\mathrm{diag}(p) - p p^\top}, with \eqn{S^{+}} its pseudo-inverse. It
#' has no sign. If the covariance rank is below the number of levels minus
#' one, the overall SMD is unavailable (`NA`) with a warning; use
#' `cat_smd = "level"` to inspect the individual levels. Otherwise it is
#' never smaller than the largest absolute SMD of a
#' single level, so a threshold applied to it is the stricter one. It is the
#' number `gtsummary::add_difference(test = ~ "smd")` and
#' `halfmoon::tidy_smd()` report for the variable.
#'
#' Two consequences are exact rather than approximate. Stabilising the ATE
#' weight multiplies it by a constant within each arm, so its balance is
#' identical to the unstabilised weight's. And overlap weights from a
#' logistic score on the same covariates balance every covariate mean
#' exactly.
#'
#' @param data A data frame holding every column named below. Analysis
#'   columns stored as `integer64` are converted through bit64 within the
#'   exact double range [-2^53, 2^53]. Labelled storage is unwrapped; SPSS
#'   user-defined missing values and explicit missing factor levels become
#'   `NA` before complete-case selection. Matrix, list, complex and S4
#'   analysis columns require explicit conversion. Grouped/rowwise frames
#'   and data.tables are analysed and returned as ordinary data frames.
#' @param treat Length-1 character. The binary exposure column, coded as in
#'   [get_PSW()]: `0`/`1`, logical, or a two-level factor or character column
#'   whose second level is the treated one.
#' @param adj_var Character vector of covariates every scheme adjusts for and
#'   the plot reports. A two-level factor or character covariate is shown
#'   under its own name; a factor with more levels as set by `cat_smd`.
#' @param methods Character vector of shorthands, any of `"PSM"`, `"ATE"`,
#'   `"ATT"`, `"ATC"`, `"ATO"`, `"ATM"`, `"EW"` (all seven by default), or a
#'   list of scheme specifications named by their legend labels, which may
#'   also hold unnamed shorthands; see *Specifying schemes*. At most 14
#'   schemes.
#' @param cat_smd How a factor or character covariate with three or more
#'   levels is summarised. `"overall"` (default): one row under the variable's
#'   name, the unsigned multivariate SMD of Yang & Dalton (2012); see
#'   *Balance measure*. `"level"`: one signed row per level, named
#'   `<variable>_<level>`, as cobalt and [plt_PSM()] / [plt_PSW()] report it.
#'   Two-level and numeric covariates are the same either way.
#' @param tbl `FALSE` (default) or `TRUE` to also build `$tbl`, which takes
#'   far longer than the balance itself. Needs gtsummary and survey. A
#'   non-empty `save_tbl` builds it too.
#' @param cores `NULL` (default) or a positive whole number: how many worker
#'   processes build `$tbl`, whose tables, one for the unadjusted sample and
#'   one per scheme, are independent. `NULL` takes the number of tables, at
#'   most 8 and at most the logical cores less two; `1` builds them one after
#'   another. Workers are forked through [parallel::mclapply()] on Linux and
#'   macOS (WSL included) and started as a PSOCK cluster on Windows, which
#'   costs a few seconds. The result does not depend on `cores`, and it is
#'   ignored when no table is built.
#' @param var_names Display labels keyed by column name, a named list or
#'   character vector such as `c(bmi = "Body mass index")`, used by the plot
#'   and `$tbl`; `character()` (default) gives none of your own. When RegR is
#'   installed they are merged into `RegR::name_map_seer`, your labels
#'   winning; a name given more than once keeps its first label, as in the
#'   `colname_map` of `RegR::get_tb_gtsummary()`. Names not in `adj_var` are
#'   ignored, and a covariate in neither
#'   keeps its column name. `NULL` keeps every column name. In the plot a
#'   factor's label is applied to each of its level rows; `$balance` always
#'   keeps the column names.
#' @param love_args Named list of love plot settings; `list()` (default)
#'   keeps every default. Unknown or duplicated fields are an error.
#'   \describe{
#'     \item{`threshold`}{Single number, default `0.1`. Where the reference
#'       line is drawn.}
#'     \item{`colors`}{Character vector, or `NULL` (default) for
#'       `UtilsR::pal_lancet`. The first colour is the unadjusted sample's,
#'       then one per scheme; at least that many are required.}
#'     \item{`shapes`}{Point shapes, numeric or character, in the same order;
#'       `NULL` (default) for `17, 16, 15, 18, 8, ...`.}
#'     \item{`size`}{Single positive number, default `3.5`. Point size.}
#'     \item{`line`}{`TRUE` (default) or `FALSE`. Connect each scheme's
#'       points.}
#'     \item{`var_order`}{Row order: `"unadjusted"` (default, largest
#'       unadjusted imbalance on top), `"alphabetical"`, a scheme label to
#'       sort by that scheme, or `NULL` for the order of `adj_var`.}
#'     \item{`base_size`}{Single positive number, default `16`, passed to
#'       [UtilsR::theme_my()].}
#'     \item{`ref_color`}{Single string, default `"red"`. Colour of the
#'       reference line.}
#'     \item{`legend_position`}{Two numbers in panel coordinates for a legend
#'       inside the plot, default `c(0.99, 0.02)` (the corner the unadjusted
#'       ordering leaves empty), or one string: `"right"`, `"bottom"`,
#'       `"none"`, or [UtilsR::fmt_legend()]'s corners `"br"`, `"bl"`,
#'       `"tr"`, `"tl"`. A `"top"` or `"bottom"` legend is laid out in two
#'       columns.}
#'     \item{`legend_justification`}{Two numbers or one string, or `NULL`
#'       (default) to let [UtilsR::fmt_legend()] derive it from the
#'       position: the nearest corner for an inside legend, `c(1, 0)` for the
#'       default, and the centre of the edge for an outside one.}
#'   }
#' @param save_plt `NULL` or a named list saving the love plot as a PDF
#'   through `RegR::save_plt()`: `filename`, and optionally `width` and
#'   `height`, which default to the plot's own pinned size. `NULL` and
#'   `list()` both mean no file is written. Unknown or duplicated fields are
#'   an error.
#' @param save_tbl `NULL` or a named list saving `$tbl` as a Word file through
#'   `RegR::save_tb()`: `path` (the output directory), `title` (also the file
#'   name), `note`, `header_value`, `header_colwidths` and `line_spacing`.
#'   `NULL` and `list()` both mean no file is written. A non-empty list builds
#'   the table even when `tbl = FALSE`. Unknown or duplicated fields are an
#'   error. Neither save changes the return value.
#'
#' @return A list of class `c("bal_res", "list")`, printed as the largest
#'   SMD of each method, of
#'   \describe{
#'     \item{`plt`}{The love plot, a `ggplot` carrying its pinned size in
#'       `attr(., "plot_size")`.}
#'     \item{`balance`}{A data frame with `variable`, `method` (the scheme
#'       label, or `"Unadjusted"`) and `smd`, the signed standardised mean
#'       difference, treated minus control; the `cat_smd = "overall"` row of a
#'       factor with three or more levels is unsigned.}
#'     \item{`data`}{The complete cases of `treat`, `adj_var` and any `ps`
#'       column a scheme names, with one weight column per scheme named by
#'       its label. Unmatched and trimmed units have weight `0`, so every
#'       column describes the same rows.}
#'     \item{`tbl`}{Only with `tbl = TRUE` or a non-empty `save_tbl`: a
#'       gtsummary `tbl_merge` with
#'       one spanner for the unadjusted sample and one per scheme, each
#'       holding both arms' mean (SD) or n (%), an SMD and a p-value. A scheme
#'       is summarised through a survey design on its weight column
#'       (`tbl_svysummary()`), zero weights included. The SMD is gtsummary's
#'       own, `add_difference(test = ~ "smd")` through `smd::smd()`, shown as
#'       an absolute value: a factor's value equals `cat_smd = "overall"`
#'       whatever `cat_smd` is, and a continuous covariate's divides the
#'       variance by \eqn{n}, so it is slightly larger than the plot's. The
#'       p-values are `add_p()`'s defaults. For a continuous covariate that is
#'       a rank test, which compares whole distributions, so it can be small
#'       beside a small SMD; under a matching or trimming scheme it also
#'       ranks the zero-weight units, unlike the chi-squared test.}
#'   }
#'
#' @seealso [get_PSM()] and [get_PSW()] for one design at a time, with
#'   effective sample sizes and weight diagnostics.
#'
#' @references
#' Yang D, Dalton JE. A unified approach to measuring the effect size between
#' two groups using SAS. SAS Global Forum 2012; Paper 335-2012.
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
#'
#' # Restyle: stricter threshold, legend below, readable covariate names
#' get_bal(d, treat = "z", adj_var = c("x1", "x2", "x3"),
#'         methods = c("PSM", "ATO"),
#'         var_names = c(x1 = "Age", x2 = "Male", x3 = "Frailty"),
#'         love_args = list(threshold = 0.05, legend_position = "bottom"))$plt
#'
#' # Table 1 per scheme: both arms, SMD and p-value under one spanner each
#' if (requireNamespace("gtsummary", quietly = TRUE) &&
#'     requireNamespace("survey", quietly = TRUE))
#'   get_bal(d, treat = "z", adj_var = c("x1", "x2", "x3"),
#'           methods = c("PSM", "ATO"), tbl = TRUE)$tbl
#' }
#'
#' @export
get_bal <- function(data,
                    treat,
                    adj_var,
                    methods   = c("PSM", "ATE", "ATT", "ATC", "ATO", "ATM",
                                  "EW"),
                    cat_smd   = c("overall", "level"),
                    tbl       = FALSE,
                    cores     = NULL,
                    var_names = character(),
                    love_args = list(),
                    save_plt  = list(),
                    save_tbl  = list()) {

  if (!is.data.frame(data) || !nrow(data))
    stop("`data` must be a non-empty data frame.", call. = FALSE)
  cat_smd <- match.arg(cat_smd)
  if (!is.logical(tbl) || length(tbl) != 1L || is.na(tbl))
    stop("`tbl` must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(cores) &&
      !(is.numeric(cores) && length(cores) == 1L && !is.na(cores) &&
        cores >= 1 && cores == round(cores)))
    stop("`cores` must be NULL or a positive whole number.", call. = FALSE)
  save_plt <- .merge_named_arg(save_plt, list(), "save_plt",
                               c("filename", "width", "height"))
  save_tbl <- .merge_named_arg(save_tbl, list(), "save_tbl", .BAL_TB_SAVE)
  tbl <- tbl || length(save_tbl) > 0L
  if (tbl) for (pkg in c("gtsummary", "survey"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required for get_bal(tbl = TRUE) and `save_tbl`.",
                   pkg), call. = FALSE)
  treat   <- .sens_check_col(treat, data, "treat", n = 1L)
  adj_var <- .sens_check_col(adj_var, data, "adj_var")
  if (is.null(adj_var))
    stop("`adj_var` is required: it is what every scheme adjusts for and the plot reports.",
         call. = FALSE)
  specs <- .bal_specs(methods)
  labs  <- names(specs)
  k     <- length(labs) + 1L          # the unadjusted sample, then each scheme

  la  <- .merge_named_arg(love_args, .BAL_LOVE_DEFAULTS, "love_args")
  bad <- function(field, msg)
    stop(sprintf("`love_args$%s` %s", field, msg), call. = FALSE)
  one_num <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x)
  place   <- function(x) (is.character(x) && length(x) == 1L && !is.na(x)) ||
    (is.numeric(x) && length(x) == 2L && !anyNA(x))
  if (!one_num(la$threshold)) bad("threshold", "must be a single number.")
  for (f in c("size", "base_size"))
    if (!one_num(la[[f]]) || la[[f]] <= 0)
      bad(f, "must be a single positive number.")
  if (!is.logical(la$line) || length(la$line) != 1L || is.na(la$line))
    bad("line", "must be TRUE or FALSE.")
  if (!is.character(la$ref_color) || length(la$ref_color) != 1L)
    bad("ref_color", "must be a single colour string.")
  if (is.null(la$colors)) la$colors <- UtilsR::pal_lancet
  if (is.null(la$shapes)) la$shapes <- .BAL_SHAPES
  for (f in c("colors", "shapes"))
    if (length(la[[f]]) < k)
      bad(f, sprintf("needs at least %d values, one for the unadjusted sample and one per scheme; got %d.",
                     k, length(la[[f]])))
  if (!is.null(la$var_order) &&
      !(is.character(la$var_order) && length(la$var_order) == 1L &&
        la$var_order %in% c("unadjusted", "alphabetical", labs)))
    bad("var_order", "must be NULL, \"unadjusted\", \"alphabetical\" or one scheme label.")
  if (!place(la$legend_position))
    bad("legend_position", "must be one string such as \"right\" or \"none\", or two numbers.")
  if (!is.null(la$legend_justification) && !place(la$legend_justification))
    bad("legend_justification", "must be NULL, one string, or two numbers.")
  # var_names: the caller's labels first, then RegR::name_map_seer, as
  # RegR::get_tb_gtsummary() merges them, the first of a repeated name
  # winning; only adj_var entries are kept.
  vn <- NULL
  if (!is.null(var_names)) {
    one_lab <- function(x) is.character(x) && length(x) == 1L && !is.na(x)
    if (!(is.list(var_names) || is.character(var_names)) ||
        (length(var_names) &&
         (is.null(names(var_names)) || any(!nzchar(names(var_names))) ||
          !all(vapply(var_names, one_lab, logical(1))))))
      stop("`var_names` must be NULL or a named list or character vector of single labels.",
           call. = FALSE)
    vn <- c(as.list(var_names),
            if (requireNamespace("RegR", quietly = TRUE)) RegR::name_map_seer)
    vn <- unlist(vn[!duplicated(names(vn)) & names(vn) %in% adj_var])
    if (!length(vn)) vn <- NULL
  }
  # cobalt keys plot rows by label, so two covariates under one label would
  # be drawn as a single row
  disp <- adj_var
  lbd  <- adj_var %in% names(vn)
  disp[lbd] <- vn[adj_var[lbd]]
  clash <- disp %in% disp[duplicated(disp)]
  if (any(clash)) {
    grp <- split(adj_var[clash], disp[clash])
    stop(sprintf("`var_names` gives covariates the same display label, which would share one plot row: %s. Give them distinct labels in `var_names`.",
                 paste(sprintf("%s -> \"%s\"", vapply(grp, function(v)
                   paste0("`", v, "`", collapse = ", "), ""), names(grp)),
                   collapse = "; ")), call. = FALSE)
  }

  # Incomplete rows go once, here: left to get_PSM() / get_PSW(), each scheme
  # could drop a different set and the columns would describe different
  # people.
  ps_cols <- unique(unlist(lapply(specs, `[[`, "ps")))
  ps_cols <- .sens_check_col(ps_cols, data, "ps")
  if (inherits(data, c("data.table", "grouped_df", "rowwise_df")))
    data <- as.data.frame(data)
  # Normalise only analysis columns; complete.cases() ignores class-specific
  # missingness and reads integer64's bit storage as ordinary doubles.
  for (v in unique(c(treat, adj_var, ps_cols))) {
    x <- data[[v]]
    if (!is.null(dim(x)) || is.list(x) || isS4(x))
      stop(sprintf("Column `%s` must be an atomic vector, not a matrix, list or S4 column.",
                   v), call. = FALSE)
    if (is.complex(x))
      stop(sprintf("Column `%s` must not contain complex values.", v), call. = FALSE)
    if (inherits(x, "integer64")) {
      if (!requireNamespace("bit64", quietly = TRUE))
        stop("Package 'bit64' is required for integer64 columns.", call. = FALSE)
      if (any(x < -(2^53) | x > 2^53, na.rm = TRUE))
        stop(sprintf("Column `%s` has integer64 values outside the exact double range [-2^53, 2^53].",
                     v), call. = FALSE)
      x <- as.double(x)
    } else if (inherits(x, "haven_labelled")) {
      spss <- inherits(x, "haven_labelled_spss")
      values <- attr(x, "na_values", exact = TRUE)
      range <- attr(x, "na_range", exact = TRUE)
      label <- attr(x, "label", exact = TRUE)
      attributes(x) <- NULL
      if (spss) {
        missing <- is.na(x) | x %in% values
        if (length(range) == 2L)
          missing <- missing | (x >= range[1L] & x <= range[2L])
        x[missing] <- NA
      }
      attr(x, "label") <- label
    } else if (is.factor(x) && anyNA(levels(x))) {
      levels(x) <- levels(x)
    }
    if (is.numeric(x) && any(is.infinite(x)))
      stop(sprintf("Column `%s` must not contain infinite values; use NA for missingness.",
                   v), call. = FALSE)
    data[[v]] <- x
  }
  data    <- .sens_complete(data, c(treat, adj_var, ps_cols))
  hit <- intersect(labs, names(data))
  if (length(hit))
    stop(sprintf("Column(s) %s already exist in `data`; get_bal() writes one weight column per scheme label. Rename them, or relabel the scheme.",
                 paste0("`", hit, "`", collapse = ", ")), call. = FALSE)

  # Weighting schemes that differ only in `estimand` share one score, trimming
  # and truncation, so each such group is one get_PSW() call, whose weights
  # are those of separate calls; a matching scheme is fitted on its own.
  fit <- function(nms, fn, args) tryCatch(
    do.call(fn, c(list(data = data, treat = treat, adj_var = adj_var,
                       balance = FALSE), args)),
    error = function(e)
      stop(sprintf("Scheme %s: %s", paste0("`", nms, "`", collapse = ", "),
                   conditionMessage(e)), call. = FALSE))
  key <- vapply(labs, function(nm) {
    s <- specs[[nm]]
    if (s$design == "matching") return(paste0("matching:", nm))
    s <- s[setdiff(names(s), c("design", "estimand"))]
    paste(deparse(s[order(names(s))]), collapse = "")
  }, "")
  w <- list()
  for (g in unique(key)) {
    nms <- labs[key == g]
    s   <- specs[[nms[1L]]]
    if (s$design == "matching") {
      r  <- fit(nms, get_PSM, s[names(s) != "design"])
      wc <- stats::setNames(attr(r, "analysis")$wcols, nms)
    } else {
      est <- toupper(vapply(specs[nms], `[[`, "", "estimand"))
      r   <- fit(nms, get_PSW, c(s[setdiff(names(s), c("design", "estimand"))],
                                 list(estimand = unique(est))))
      wc  <- stats::setNames(paste0("w_", tolower(est)), nms)
    }
    # NA marks a trimmed unit; like an unmatched one it is out of the sample
    for (nm in nms) {
      x <- r$data[[wc[[nm]]]]
      x[is.na(x)] <- 0
      w[[nm]] <- x
    }
  }
  w <- as.data.frame(w[labs], check.names = FALSE)

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

  # cat_smd = "overall": one Yang & Dalton row replaces a factor's level rows,
  # D weighted by each column's scheme, S from the unweighted arm proportions
  # as every other row's denominator; the pseudo-inverse of S over all k
  # levels equals the inverse over k - 1, whichever level is dropped.
  B     <- bt$Balance
  dcols <- grep("^Diff\\.", names(B))
  multi <- if (cat_smd == "overall") adj_var[vapply(data[adj_var], function(v)
    (is.factor(v) || is.character(v)) && length(unique(v)) > 2L, logical(1))]
  for (v in multi) {
    x    <- factor(data[[v]])
    prop <- function(s, wt)
      vapply(levels(x), function(l) sum(wt[s] * (x[s] == l)) / sum(wt[s]), 0)
    covm <- function(s) {
      q <- prop(s, rep(1, length(x)))
      diag(q) - outer(q, q)
    }
    e  <- svd((covm(z == 1) + covm(z == 0)) / 2)
    ok <- e$d > sqrt(.Machine$double.eps) * e$d[1L]
    singular <- sum(ok) < nlevels(x) - 1L
    if (singular)
      warning(sprintf("Covariate `%s` has a singular pooled covariance beyond the redundant level; overall SMD is NA. Use cat_smd = \"level\" to inspect levels.",
                      v), call. = FALSE)
    Si <- e$v[, ok, drop = FALSE] %*% (t(e$u[, ok, drop = FALSE]) / e$d[ok])
    at  <- which(rownames(B) %in% paste0(v, "_", levels(x)))
    row <- B[at[1L], , drop = FALSE]
    row[dcols] <- vapply(c(list(rep(1, length(x))), w), function(wt) {
      if (singular) return(NA_real_)
      D <- prop(z == 1, wt) - prop(z == 0, wt)
      sqrt(max(0, drop(D %*% Si %*% D)))
    }, 0)
    rownames(row) <- v
    keep <- setdiff(seq_len(nrow(B)), at)
    B <- rbind(B[keep[keep < at[1L]], , drop = FALSE], row,
               B[keep[keep > at[1L]], , drop = FALSE])
  }
  bt$Balance <- B

  vars <- rownames(B)
  vars[vars %in% names(map)] <- unname(map[vars[vars %in% names(map)]])
  bal <- data.frame(
    variable = rep(vars, length(labs) + 1L),
    method   = rep(c("Unadjusted", labs), each = nrow(B)),
    # "Diff.Un" then one column per scheme, in order; cobalt calls a lone
    # scheme "Diff.Adj" rather than by its name
    smd      = unlist(B[grep("^Diff\\.", names(B))], use.names = FALSE),
    stringsAsFactors = FALSE)

  # Plot labels: a two-level factor under its variable's label; any other
  # var_names entry goes to cobalt, which labels each level of a factor.
  lab_of <- function(v) if (v %in% names(vn)) vn[[v]] else v
  # cobalt reads a factor's name in var.names as the label of its level rows,
  # so a cat_smd = "overall" row is relabelled by renaming it instead.
  own <- intersect(multi, names(vn))
  if (length(own))
    rownames(bt$Balance)[match(own, rownames(bt$Balance))] <- unname(vn[own])
  pn <- c(vapply(map, lab_of, ""), vn[setdiff(names(vn), c(map, own))])

  # cobalt refuses an entirely unavailable statistic. Keep the variable
  # labels on an empty panel, without inventing zero-valued points.
  p <- if (all(is.na(bal$smd))) {
    empty <- data.frame(variable = factor(rownames(bt$Balance),
                        levels = rev(rownames(bt$Balance))), smd = NA_real_)
    ggplot2::ggplot(empty, ggplot2::aes(x = smd, y = variable)) +
      ggplot2::geom_blank(na.rm = TRUE) +
      ggplot2::scale_x_continuous(limits = c(0, max(0.1, la$threshold))) +
      ggplot2::labs(x = "Absolute Standardized Mean Differences", y = NULL,
                    subtitle = "Overall SMD unavailable: singular covariance")
  } else cobalt::love.plot(
    bt, stats = "mean.diffs", abs = TRUE, var.order = la$var_order,
    line = la$line, thresholds = c(m = la$threshold),
    colors = la$colors[seq_len(k)], shapes = la$shapes[seq_len(k)],
    size = la$size, drop.distance = TRUE, drop.missing = FALSE, labels = FALSE,
    sample.names = c("Unadjusted", labs),
    var.names = if (length(pn)) pn, title = NULL)
  p <- p + UtilsR::theme_my(base_size = la$base_size)
  p <- UtilsR::fmt_ref(p, x = la$threshold, color = la$ref_color)
  # RegR's framed inside legend, but in the bottom-right corner by default:
  # rows are sorted by unadjusted SMD, so that corner is the one the
  # reference line leaves empty, and with many schemes RegR's c(0.8, 0.4)
  # covers it. fmt_legend() derives the justification from the position and
  # overwrites one passed to it, so an explicit one is set afterwards, also
  # as legend.justification.inside: fmt_legend() sets that one for an inside
  # legend, which reads it in preference. Scheme labels are long, so a legend
  # above or below the panel takes two columns rather than one clipped row.
  p <- UtilsR::fmt_legend(
    p, legend.position = la$legend_position,
    ncol = if (identical(la$legend_position, "top") ||
               identical(la$legend_position, "bottom")) 2,
    legend.background = ggplot2::element_rect(
      fill = "white", colour = "#D6D6D6", linewidth = 1))
  if (!is.null(la$legend_justification))
    p <- p + ggplot2::theme(
      legend.justification = la$legend_justification,
      legend.justification.inside = la$legend_justification)
  p <- p + ggplot2::labs(x = "Absolute standardized mean difference",
                         y = "Covariates")

  # tbl: gtsummary's own summaries, SMD and p-value per sample, merged. The
  # zero weights of unmatched and trimmed units stay in each design, so the
  # SMD's unweighted denominator is the same complete cases as the plot's.
  tb <- NULL
  if (tbl) {
    one <- function(wt) {
      st <- list(gtsummary::all_continuous() ~ "{mean} ({sd})")
      t <- if (is.null(wt))
        gtsummary::tbl_summary(dd, by = gtsummary::all_of(treat), statistic = st,
                               label = lb)
      else
        gtsummary::tbl_svysummary(
          survey::svydesign(ids = ~1, weights = wt, data = dd),
          by = gtsummary::all_of(treat), statistic = st, label = lb)
      t <- gtsummary::add_difference(t, test = gtsummary::everything() ~ "smd")
      t <- gtsummary::add_p(t, pvalue_fun = gtsummary::label_style_pvalue(digits = 3))
      t <- gtsummary::modify_column_hide(t, "conf.low")
      t <- gtsummary::modify_header(t, estimate = "**SMD**")
      gtsummary::modify_fmt_fun(t, estimate ~ function(x)
        gtsummary::style_number(abs(x), digits = 3))
    }
    # Only what a table needs travels to a PSOCK worker, not this frame; the
    # stats namespace is where gtsummary finds "{mean} ({sd})".
    dd <- data[c(treat, adj_var)]
    if (is.factor(dd[[treat]])) dd[[treat]] <- droplevels(dd[[treat]])
    environment(one) <- list2env(
      list(dd = dd, lb = if (!is.null(vn)) as.list(vn),
           treat = treat), parent = asNamespace("stats"))
    W  <- c(list(NULL), as.list(w))
    dc <- parallel::detectCores()
    nc <- if (is.null(cores)) min(length(W), 8L, max(1L, dc - 2L, na.rm = TRUE))
          else min(as.integer(cores), length(W))
    # survey reads weights summing below n as mis-scaled sampling weights;
    # balancing weights (ATO, matching zeros) do that by design.
    tabs <- withCallingHandlers({
      if (nc < 2L) {
        lapply(W, one)
      } else if (.Platform$OS.type == "windows") {
        cl <- parallel::makeCluster(nc)
        tryCatch(parallel::parLapply(cl, W, one),
                 finally = parallel::stopCluster(cl))
      } else {
        r   <- parallel::mclapply(W, one, mc.cores = nc)
        bad <- vapply(r, function(x) is.null(x) || inherits(x, "try-error"),
                      logical(1))
        if (any(bad))
          stop(sprintf("Building the table of %s failed: %s",
                       c("Unadjusted", labs)[which(bad)[1L]],
                       if (is.null(r[[which(bad)[1L]]])) "the worker exited"
                       else conditionMessage(attr(r[[which(bad)[1L]]], "condition"))),
               call. = FALSE)
        r
      }
    }, warning = function(cnd)
      if (grepl("Sample size greater than population size",
                conditionMessage(cnd), fixed = TRUE))
        invokeRestart("muffleWarning"))
    tb <- gtsummary::tbl_merge(tabs, tab_spanner = paste0("**", c("Unadjusted", labs), "**"))
  }

  data[labs] <- w
  plt <- .psw_save(p, c(8, max(5, 1.5 + 0.45 * nrow(B))), save_plt)
  if (length(save_tbl)) {
    if (!requireNamespace("RegR", quietly = TRUE))
      stop("Package 'RegR' is required for a non-empty `save_tbl`.", call. = FALSE)
    do.call(RegR::save_tb, c(list(data = gtsummary::as_flex_table(tb)), save_tbl))
  }
  out <- list(plt = plt, balance = bal, data = data)
  if (tbl) out$tbl <- tb
  class(out) <- c("bal_res", "list")
  out
}


# ---- print -------------------------------------------------------------------

#' @export
#' @noRd
print.bal_res <- function(x, ...) {
  b   <- x$balance
  lv  <- unique(b$method)
  m   <- factor(b$method, lv)
  tab <- data.frame(method = lv,
                    max = as.numeric(tapply(abs(b$smd), m, function(x)
                      if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))),
                    over = as.integer(tapply(abs(b$smd) > 0.1, m, function(x)
                      if (all(is.na(x))) NA_integer_ else sum(x, na.rm = TRUE))))
  names(tab) <- c("method", "max |SMD|", "over 0.1")
  k <- sum(lv != "Unadjusted")
  UtilsR::pr_head("bal_res", "Covariate balance",
                  sprintf("%d scheme%s", k, if (k == 1L) "" else "s"),
                  sprintf("n = %d", NROW(x$data)))
  UtilsR::pr_kv(list(covariates = unique(b$variable)))
  UtilsR::pr_section("Largest standardised mean difference")
  over <- !is.na(tab[["max |SMD|"]]) & tab[["max |SMD|"]] > 0.1
  UtilsR::pr_table(tab, highlight = over, n = NULL, digits = 3)
  if (any(over)) UtilsR::pr_note("Highlighted: an |SMD| above 0.1.")
  UtilsR::pr_note("<result>$plt is the love plot, $balance every SMD, $data the ",
                  "weights", if (!is.null(x$tbl)) ", $tbl the table")
  UtilsR::pr_note("Components: {names(x)}; see {.code UtilsR::pr_tree()}.")
  invisible(x)
}
