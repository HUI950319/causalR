# =============================================================================
# psm-plt.R -- plots for a get_PSM() result
# =============================================================================
#
# Architecture (2 layers):
#
#   L1  plt_PSM(x, type, method, ...)
#         |
#         +-- L2 renderers
#               .psm_plt_love     balance across schemes, halfmoon::geom_love()
#               .psm_plt_ess      effective sample size, native ggplot2
#               .psm_plt_weight   matching weights, native ggplot2
#               .psm_plt_ps       score overlap, halfmoon::geom_mirror_histogram()
#
# The four types mirror plt_PSW() so a weighting result and a matching result
# can be read the same way. Two things differ, both because matching drops
# units rather than down-weighting them: the weight figure plots only the
# matched units and names the discarded count in the subtitle, and the score
# histogram compares the matched sample against the whole cohort rather than
# weighted against unweighted.
# =============================================================================

# The weight columns are what the tables key on; the method names are what a
# reader recognises.
#' @keywords internal
#' @noRd
.psm_labels <- function(wcols) {
  stats::setNames(c("Unmatched", sub("^w_", "", wcols)),
                  c("observed", wcols))
}


# ---- L2 renderers ----------------------------------------------------------

#' @keywords internal
#' @noRd
.psm_plt_love <- function(x, wcols, threshold, title) {
  if (!requireNamespace("halfmoon", quietly = TRUE))
    stop("Package 'halfmoon' is required for plt_PSM(type = \"love\")",
         call. = FALSE)
  lab <- .psm_labels(wcols)
  d   <- x$balance[x$balance$metric == "smd" &
                     x$balance$method %in% names(lab), , drop = FALSE]
  d$Matching <- factor(unname(lab[d$method]), levels = unname(lab))
  d$smd      <- abs(d$estimate)

  ggplot2::ggplot(d, ggplot2::aes(x = smd, y = variable,
                                  group = Matching, colour = Matching)) +
    halfmoon::geom_love(vline_xintercept = threshold) +
    ggplot2::labs(
      x = "Absolute standardised mean difference",
      y = NULL, colour = NULL,
      title = if (is.null(title)) "Covariate balance by matching scheme"
              else title) +
    ggplot2::theme(legend.position = "bottom")
}

#' @keywords internal
#' @noRd
.psm_plt_ess <- function(x, wcols, title) {
  a <- attr(x, "analysis")
  d <- x$stats[x$stats$method %in% sub("^w_", "", wcols), , drop = FALSE]
  d <- d[order(d$ess), , drop = FALSE]
  d$method <- factor(d$method, levels = d$method)
  # not `n`: the data frame has an `n` column and aes() would resolve to it
  n_total <- a$n

  ggplot2::ggplot(d, ggplot2::aes(x = ess, y = method)) +
    ggplot2::geom_col(fill = "grey35", width = 0.65) +
    ggplot2::geom_vline(xintercept = n_total, linetype = 2,
                        colour = "grey50") +
    ggplot2::annotate("text", x = n_total, y = Inf,
                      label = sprintf("n = %d", n_total),
                      hjust = 1.05, vjust = 1.6, size = 3.2,
                      colour = "grey35") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.0f matched, ESS %.0f", n, ess)),
      hjust = -0.06, size = 3.2) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(c(0, 0.28))) +
    ggplot2::labs(
      x = "Effective sample size", y = NULL,
      title = if (is.null(title))
        "How much information each matching scheme keeps" else title)
}

#' @keywords internal
#' @noRd
.psm_plt_weight <- function(x, wcols, title) {
  a   <- attr(x, "analysis")
  lab <- .psm_labels(wcols)
  arm <- factor(.psw_treat(x$data[[a$treat]], a$treat)$z, levels = c(0L, 1L),
                labels = c("Control", "Treated"))

  long <- do.call(rbind, lapply(wcols, function(w) {
    keep <- x$data[[w]] > 0
    data.frame(method = unname(lab[[w]]), arm = arm[keep],
               weight = x$data[[w]][keep], stringsAsFactors = FALSE)
  }))
  long$method <- factor(long$method, levels = unname(lab[wcols]))

  dropped <- vapply(wcols, function(w) sum(x$data[[w]] == 0), integer(1))
  sub <- paste0("matched units only; dropped ",
                paste(sprintf("%s %d", unname(lab[wcols]), dropped),
                      collapse = ", "))

  ggplot2::ggplot(long, ggplot2::aes(x = weight, y = arm, fill = arm)) +
    ggplot2::geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.35,
                          width = 0.55, show.legend = FALSE) +
    ggplot2::geom_vline(xintercept = 1, linetype = 3, colour = "grey50") +
    ggplot2::facet_wrap(~ method, ncol = 2, scales = "free_x") +
    ggplot2::labs(
      x = "Matching weight", y = NULL, subtitle = sub,
      title = if (is.null(title)) "Matching weights by arm" else title)
}

#' @keywords internal
#' @noRd
.psm_plt_ps <- function(x, wcols, bins, title) {
  if (!requireNamespace("halfmoon", quietly = TRUE))
    stop("Package 'halfmoon' is required for plt_PSM(type = \"ps\")",
         call. = FALSE)
  a   <- attr(x, "analysis")
  lab <- .psm_labels(wcols)
  arm <- factor(.psw_treat(x$data[[a$treat]], a$treat)$z, levels = c(0L, 1L),
                labels = c("Control", "Treated"))

  long <- do.call(rbind, lapply(wcols, function(w) data.frame(
    method = unname(lab[[w]]), arm = arm, ps = x$data$ps,
    weight = x$data[[w]], stringsAsFactors = FALSE)))
  long$method <- factor(long$method, levels = unname(lab[wcols]))

  ggplot2::ggplot(long, ggplot2::aes(x = ps, group = arm)) +
    halfmoon::geom_mirror_histogram(
      ggplot2::aes(fill = arm), bins = bins, alpha = 0.35) +
    halfmoon::geom_mirror_histogram(
      ggplot2::aes(fill = arm, weight = weight), bins = bins) +
    ggplot2::facet_wrap(~ method, ncol = 2) +
    ggplot2::labs(
      x = "Propensity score", y = "Count", fill = NULL,
      title = if (is.null(title))
        "Score overlap: whole cohort (pale) and matched (solid)" else title) +
    ggplot2::theme(legend.position = "bottom")
}


# ---- L1 public entry point -------------------------------------------------

#' Plots for a propensity score matching result
#'
#' Four views of a [get_PSM()] result, mirroring [plt_PSW()]: how well each
#' matching scheme balances the covariates, how much sample each one keeps,
#' what the matching weights look like, and how much of the score range
#' survives matching.
#'
#' @param x A `psm_res` object from [get_PSM()].
#' @param type Which figure to draw. `"love"` (default when `x` carries a
#'   balance table) overlays the absolute standardised mean differences of
#'   every scheme on the unmatched ones. `"ess"` (the default otherwise)
#'   compares effective sample sizes. `"weight"` shows the matching weights of
#'   the matched units in each arm. `"ps"` mirrors the score histograms of the
#'   whole cohort against the matched sample.
#' @param method Character vector selecting which matching schemes to draw, or
#'   `NULL` (default) for all of those `x` holds. `"love"` always keeps the
#'   unmatched series as the reference.
#' @param threshold Numeric, default `0.1`. Where `"love"` draws its
#'   reference line.
#' @param bins Number of histogram bins for `"ps"`, default `40`.
#' @param title Plot title, or `NULL` (default) for a generated one.
#' @param save `NULL` or a list with `filename`, and optionally `width` and
#'   `height`, passed to `RegR::save_plt()`. `NULL` and `list()` both mean no
#'   file is written; the defaults come from the plot's own pinned size.
#'
#' @return A `ggplot` object, carrying its pinned size in
#'   `attr(., "plot_size")`.
#'
#' @seealso [get_PSM()], and [plt_PSW()] for the weighting dual.
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
#' plt_PSM(res, "love")
#' plt_PSM(res, "ess")
#' }
#'
#' @export
plt_PSM <- function(x,
                    type      = c("love", "ess", "weight", "ps"),
                    method    = NULL,
                    threshold = 0.1,
                    bins      = 40L,
                    title     = NULL,
                    save      = list()) {

  if (!inherits(x, "psm_res"))
    stop("`x` must be a `psm_res` object from get_PSM().", call. = FALSE)
  a    <- attr(x, "analysis")
  spec <- .psm_plt_spec(x)

  # The default type depends on whether balance was computed, so it is
  # resolved before, not by, match.arg().
  if (missing(type)) type <- spec$type
  type <- match.arg(type)
  if (identical(type, "love") && is.null(x$balance))
    stop("type = \"love\" needs the balance table; re-run get_PSM() with balance = TRUE.",
         call. = FALSE)

  if (is.null(method)) {
    wcols <- a$wcols
  } else {
    method <- match.arg(tolower(method), a$method, several.ok = TRUE)
    wcols  <- paste0("w_", a$method[a$method %in% method])
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold))
    stop("`threshold` must be a single number.", call. = FALSE)
  if (!is.numeric(bins) || length(bins) != 1L || is.na(bins) || bins < 2)
    stop("`bins` must be a single number of at least 2.", call. = FALSE)

  p <- switch(
    type,
    love   = .psm_plt_love(x, wcols, threshold, title),
    ess    = .psm_plt_ess(x, wcols, title),
    weight = .psm_plt_weight(x, wcols, title),
    ps     = .psm_plt_ps(x, wcols, bins, title),
    stop(sprintf("Unsupported type: '%s'", type), call. = FALSE))

  plot_size <- switch(type,
                      love   = c(7, 5),
                      ess    = c(7.5, 4.5),
                      weight = c(8, 2 + 1.6 * ceiling(length(wcols) / 2)),
                      ps     = c(8, 2 + 2.2 * ceiling(length(wcols) / 2)))
  attr(p, "plot_size") <- stats::setNames(plot_size, c("width", "height"))

  if (!is.null(save) && !is.list(save))
    stop("`save` must be `NULL` or a list.", call. = FALSE)
  if (!is.null(save) && length(save) > 0L) {
    if (!requireNamespace("RegR", quietly = TRUE))
      stop("Package 'RegR' is required for a non-empty `save`", call. = FALSE)
    if (is.null(save$width))  save$width  <- plot_size[1L]
    if (is.null(save$height)) save$height <- plot_size[2L]
    do.call(RegR::save_plt, c(list(plot = p), save))
  }
  p
}
