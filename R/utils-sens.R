# =============================================================================
# utils-sens.R -- shared argument-contract and data helpers for get_sens()
# =============================================================================
#
# Architecture (1 layer, no public API):
#
#   .merge_named_arg    validate + merge a fully named list argument
#   .sens_check_col     validate that a character argument names data columns
#   .sens_complete      drop incomplete rows on the columns a backend uses
#   .sens_model_matrix  numeric design matrix (no intercept) for dml / iv / hte
#   .sens_hr_to_rr      VanderWeele-Ding square-root HR -> RR approximation
#   .sens_fmt_vec       compact vector rendering for print methods
#
# =============================================================================

# Validate a fully named list argument and merge it over `defaults`. Unnamed,
# duplicated and unknown fields are rejected before the merge. Behaviour is
# deliberately identical to `RegR:::.merge_named_arg()`; RegR internals cannot
# be called from here, so the contract is duplicated rather than imported.
#' @keywords internal
#' @noRd
.merge_named_arg <- function(x, defaults, arg_name,
                             allowed = names(defaults)) {
  if (is.null(x) || (is.list(x) && length(x) == 0L)) return(defaults)
  if (!is.list(x) || is.data.frame(x))
    stop(sprintf("`%s` must be `NULL` or a named list.", arg_name),
         call. = FALSE)
  nms <- names(x)
  if (is.null(nms) || anyNA(nms) || any(!nzchar(nms)))
    stop(sprintf("`%s` must be a fully named list.", arg_name), call. = FALSE)
  dups <- unique(nms[duplicated(nms)])
  if (length(dups))
    stop(sprintf("`%s` contains duplicated field(s): %s.", arg_name,
                 paste0("`", dups, "`", collapse = ", ")), call. = FALSE)
  unknown <- setdiff(nms, allowed)
  if (length(unknown))
    stop(sprintf(
      "`%s` contains unknown field(s): %s. Allowed fields are: %s.",
      arg_name, paste0("`", unknown, "`", collapse = ", "),
      paste0("`", allowed, "`", collapse = ", ")), call. = FALSE)
  utils::modifyList(defaults, x, keep.null = TRUE)
}

#' @keywords internal
#' @noRd
.sens_check_col <- function(x, data, arg_name, n = NULL) {
  if (is.null(x)) return(NULL)
  if (!is.character(x) || anyNA(x) || !all(nzchar(x)))
    stop(sprintf("`%s` must be a character vector of column names.", arg_name),
         call. = FALSE)
  if (!is.null(n) && length(x) != n)
    stop(sprintf("`%s` must name exactly %d column(s).", arg_name, n),
         call. = FALSE)
  miss <- setdiff(x, names(data))
  if (length(miss))
    stop(sprintf("`%s` column(s) not found in `data`: %s.", arg_name,
                 paste0("`", miss, "`", collapse = ", ")), call. = FALSE)
  x
}

# The dml and iv backends take plain vectors and a design matrix, so they get
# no model-frame NA handling of their own. Drop incomplete rows up front and
# report how many, rather than letting the backend fail deep inside a fold.
#' @keywords internal
#' @noRd
.sens_complete <- function(data, cols, verbose = FALSE) {
  keep <- stats::complete.cases(data[, cols, drop = FALSE])
  if (!any(keep))
    stop(sprintf("No complete cases remain across %s.",
                 paste0("`", cols, "`", collapse = ", ")), call. = FALSE)
  if (isTRUE(verbose) && any(!keep))
    cli::cli_inform(c("i" = "Dropped {sum(!keep)} incomplete row{?s}."))
  data[keep, , drop = FALSE]
}

#' @keywords internal
#' @noRd
.sens_quote_names <- function(vars) {
  vapply(vars, function(v) deparse1(as.name(v), backtick = TRUE),
         character(1), USE.NAMES = FALSE)
}

#' @keywords internal
#' @noRd
.sens_model_matrix <- function(data, vars, one_hot = FALSE) {
  if (is.null(vars) || !length(vars))
    return(matrix(numeric(0), nrow = nrow(data), ncol = 0L))
  quoted <- .sens_quote_names(vars)
  # Without an intercept only the first factor keeps every level. `one_hot`
  # keeps every level of every factor, which forests (get_hte) need so that
  # any level can be split off in one step and no factor gets an extra column.
  ca <- NULL
  if (one_hot) {
    fac <- vars[!vapply(data[vars], is.numeric, logical(1L))]
    if (length(fac)) {
      data[fac] <- lapply(data[fac], function(x) droplevels(as.factor(x)))
      ca <- lapply(data[fac], stats::contrasts, contrasts = FALSE)
    }
  }
  mm <- stats::model.matrix(stats::reformulate(quoted, intercept = FALSE),
                            data = data, contrasts.arg = ca)
  # Matrix backends select numeric benchmarks by the original column name.
  hit <- match(colnames(mm), quoted)
  colnames(mm)[!is.na(hit)] <- vars[hit[!is.na(hit)]]
  mm
}

# VanderWeele & Ding (2017) approximation of a risk ratio from a hazard ratio
# for a non-rare outcome. With a rare outcome the hazard ratio already
# approximates the risk ratio and no transform is applied.
#' @keywords internal
#' @noRd
.sens_hr_to_rr <- function(hr, rare = FALSE) {
  if (isTRUE(rare)) return(hr)
  (1 - 0.5^sqrt(hr)) / (1 - 0.5^sqrt(1 / hr))
}

#' @keywords internal
#' @noRd
.sens_fmt_vec <- function(x, digits = 3L) {
  if (is.null(x) || !length(x)) return("NULL")
  v <- paste(formatC(x, format = "g", digits = digits), collapse = ", ")
  if (length(x) == 1L) v else paste0("c(", v, ")")
}

# Map variables to their model-matrix columns without matching name prefixes.
#' @keywords internal
#' @noRd
.sens_coef_terms <- function(fit, vars) {
  mm <- stats::model.matrix(fit)
  terms <- attr(stats::terms(fit), "term.labels")
  stats::setNames(lapply(vars, function(v) {
    term <- match(deparse1(as.name(v), backtick = TRUE), terms)
    colnames(mm)[which(attr(mm, "assign") == term)]
  }), vars)
}

# Resolve `treat` to exactly one model coefficient. A multi-level factor cannot
# be handled by any of the four backends, so fail with the matched terms listed
# instead of silently sensitising the first dummy.
#' @keywords internal
#' @noRd
.sens_one_term <- function(fit, treat) {
  hit <- .sens_coef_terms(fit, treat)[[1L]]
  if (length(hit) != 1L)
    stop(sprintf(
      "`treat` must resolve to exactly one model coefficient; %s matched %d (%s). Use a numeric or two-level treatment.",
      paste0("`", treat, "`"), length(hit),
      if (length(hit)) paste(hit, collapse = ", ") else "none"),
      call. = FALSE)
  hit
}
