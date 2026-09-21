# =============================================================================
# psw-plt.R -- plots for a get_PSW() result
# =============================================================================
#
# Architecture (2 layers):
#
#   L1  plt_PSW(x, type, estimand, ...)
#         |
#         +-- L2 renderers
#         |     .psw_plt_love     balance across weights, halfmoon::geom_love()
#         |     .psw_plt_ess      effective sample size, native ggplot2
#         |     .psw_plt_weight   weight distribution by arm, native ggplot2
#         |     .psw_plt_ps       score overlap, halfmoon::geom_mirror_histogram()
#         |
#         +-- .psw_save   pin the size, save through RegR::save_plt();
#                         shared with plt_PSM()
#
# The love plot and the mirror histogram have established upstream geoms and
# use them. The other two read straight off $stats and $data, which already
# hold exactly what they draw -- recomputing them through halfmoon would
# re-derive the effective sample size from the untrimmed data and disagree
# with the printed table.
# =============================================================================

# The weight columns are what the tables key on; the estimand names are what
# a reader recognises.
#' @keywords internal
#' @noRd
.psw_labels <- function(wcols) {
  stats::setNames(c("Unweighted", toupper(sub("^w_", "", wcols))),
                  c("observed", wcols))
}

# Pin the size on the plot and, when `save` is a non-empty list, hand it to
# RegR::save_plt() with that size as the default. `NULL` and `list()` both
# mean no file is written.
#' @keywords internal
#' @noRd
.psw_save <- function(p, plot_size, save) {
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


# ---- L2 renderers ----------------------------------------------------------

#' @keywords internal
#' @noRd
.psw_plt_love <- function(x, wcols, threshold, title) {
  if (!requireNamespace("halfmoon", quietly = TRUE))
    stop("Package 'halfmoon' is required for plt_PSW(type = \"love\")",
         call. = FALSE)
  lab <- .psw_labels(wcols)
  d   <- x$balance[x$balance$metric == "smd" &
                     x$balance$method %in% names(lab), , drop = FALSE]
  d$Weighting <- factor(unname(lab[d$method]), levels = unname(lab))
  d$smd       <- abs(d$estimate)

  ggplot2::ggplot(d, ggplot2::aes(x = smd, y = variable,
                                  group = Weighting,
                                  colour = Weighting)) +
    halfmoon::geom_love(vline_xintercept = threshold) +
    ggplot2::labs(
      x = "Absolute standardised mean difference",
      y = NULL, colour = NULL,
      title = if (is.null(title)) "Covariate balance by weighting scheme"
              else title) +
    ggplot2::theme(legend.position = "bottom")
}

#' @keywords internal
#' @noRd
.psw_plt_ess <- function(x, wcols, title) {
  lab <- .psw_labels(wcols)
  d   <- x$stats[x$stats$estimand %in% toupper(sub("^w_", "", wcols)), ,
                 drop = FALSE]
  d   <- d[order(d$ess), , drop = FALSE]
  d$estimand <- factor(d$estimand, levels = d$estimand)
  n <- d$n[[1L]]

  ggplot2::ggplot(d, ggplot2::aes(x = ess, y = estimand)) +
    ggplot2::geom_col(fill = "grey35", width = 0.65) +
    ggplot2::geom_vline(xintercept = n, linetype = 2, colour = "grey50") +
    ggplot2::annotate("text", x = n, y = Inf, label = sprintf("n = %d", n),
                      hjust = 1.05, vjust = 1.6, size = 3.2,
                      colour = "grey35") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.0f (%.0f%%)", ess,
                                   100 * ess_pct)),
      hjust = -0.08, size = 3.2) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(c(0, 0.18))) +
    ggplot2::labs(
      x = "Effective sample size", y = NULL,
      title = if (is.null(title))
        "How much information each weighting scheme keeps" else title)
}

#' @keywords internal
#' @noRd
.psw_plt_weight <- function(x, wcols, title) {
  a   <- attr(x, "analysis")
  lab <- .psw_labels(wcols)
  d   <- x$data[!x$data$.trimmed, , drop = FALSE]
  arm <- factor(.psw_treat(d[[a$treat]], a$treat)$z, levels = c(0L, 1L),
                labels = c("Control", "Treated"))

  long <- do.call(rbind, lapply(wcols, function(w) data.frame(
    estimand = unname(lab[[w]]), arm = arm, weight = d[[w]],
    stringsAsFactors = FALSE)))
  long$estimand <- factor(long$estimand, levels = unname(lab[wcols]))

  ggplot2::ggplot(long, ggplot2::aes(x = weight, y = arm,
                                     fill = arm)) +
    ggplot2::geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.35,
                          width = 0.55, show.legend = FALSE) +
    ggplot2::geom_vline(xintercept = 1, linetype = 3, colour = "grey50") +
    ggplot2::facet_wrap(~ estimand, ncol = 2, scales = "free_x") +
    ggplot2::labs(
      x = "Weight", y = NULL,
      title = if (is.null(title))
        "Weight distribution by arm" else title)
}

#' @keywords internal
#' @noRd
.psw_plt_ps <- function(x, wcols, bins, title) {
  if (!requireNamespace("halfmoon", quietly = TRUE))
    stop("Package 'halfmoon' is required for plt_PSW(type = \"ps\")",
         call. = FALSE)
  a   <- attr(x, "analysis")
  lab <- .psw_labels(wcols)
  d   <- x$data[!x$data$.trimmed, , drop = FALSE]
  arm <- factor(.psw_treat(d[[a$treat]], a$treat)$z, levels = c(0L, 1L),
                labels = c("Control", "Treated"))

  long <- do.call(rbind, lapply(wcols, function(w) data.frame(
    estimand = unname(lab[[w]]), arm = arm, ps = d$ps, weight = d[[w]],
    stringsAsFactors = FALSE)))
  long$estimand <- factor(long$estimand, levels = unname(lab[wcols]))

  ggplot2::ggplot(long, ggplot2::aes(x = ps, group = arm)) +
    halfmoon::geom_mirror_histogram(
      ggplot2::aes(fill = arm), bins = bins, alpha = 0.35) +
    halfmoon::geom_mirror_histogram(
      ggplot2::aes(fill = arm, weight = weight), bins = bins) +
    ggplot2::facet_wrap(~ estimand, ncol = 2) +
    ggplot2::labs(
      x = "Propensity score", y = "Count", fill = NULL,
      title = if (is.null(title))
        "Score overlap before (pale) and after (solid) weighting" else title) +
    ggplot2::theme(legend.position = "bottom")
}


# ---- L1 public entry point -------------------------------------------------

#' Plots for a propensity score weighting result
#'
#' Four views of a [get_PSW()] result: how well each weighting scheme balances
#' the covariates, how much information each one costs, what the weights
#' themselves look like, and how much the two arms overlap on the score.
#'
#' @param x A `psw_res` object from [get_PSW()].
#' @param type Which figure to draw. `"love"` (default when `x` carries a
#'   balance table) overlays the absolute standardised mean differences of
#'   every weighting scheme on the unweighted ones. `"ess"` (the default
#'   otherwise) compares effective sample sizes. `"weight"` shows the weight
#'   distribution in each arm. `"ps"` mirrors the score histograms of the two
#'   arms before and after weighting.
#' @param estimand Character vector selecting which estimands to draw, or
#'   `NULL` (default) for all of those `x` holds. `"love"` always keeps the
#'   unweighted series as the reference.
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
#' @seealso [get_PSW()].
#'
#' @examplesIf requireNamespace("halfmoon", quietly = TRUE)
#' \donttest{
#' set.seed(20260921)
#' n <- 400
#' d <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1, 0.4), x3 = runif(n))
#' d$z <- rbinom(n, 1, plogis(-0.3 + 0.8 * d$x1 - 0.6 * d$x2 + 1.1 * d$x3))
#'
#' res <- get_PSW(d, treat = "z", adj_var = c("x1", "x2", "x3"))
#' plt_PSW(res, "love")
#' plt_PSW(res, "ess")
#' }
#'
#' @export
plt_PSW <- function(x,
                    type      = c("love", "ess", "weight", "ps"),
                    estimand  = NULL,
                    threshold = 0.1,
                    bins      = 40L,
                    title     = NULL,
                    save      = list()) {

  if (!inherits(x, "psw_res"))
    stop("`x` must be a `psw_res` object from get_PSW().", call. = FALSE)
  a    <- attr(x, "analysis")
  spec <- .psw_plt_spec(x)

  # The default type depends on whether balance was computed, so it is
  # resolved before, not by, match.arg().
  if (missing(type)) type <- spec$type
  type <- match.arg(type)
  if (identical(type, "love") && is.null(x$balance))
    stop("type = \"love\" needs the balance table; re-run get_PSW() with balance = TRUE.",
         call. = FALSE)

  if (is.null(estimand)) {
    wcols <- a$wcols
  } else {
    estimand <- match.arg(toupper(estimand), a$estimand, several.ok = TRUE)
    wcols <- paste0("w_", tolower(.PSW_ESTIMANDS[.PSW_ESTIMANDS %in% estimand]))
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold))
    stop("`threshold` must be a single number.", call. = FALSE)
  if (!is.numeric(bins) || length(bins) != 1L || is.na(bins) || bins < 2)
    stop("`bins` must be a single number of at least 2.", call. = FALSE)

  p <- switch(
    type,
    love   = .psw_plt_love(x, wcols, threshold, title),
    ess    = .psw_plt_ess(x, wcols, title),
    weight = .psw_plt_weight(x, wcols, title),
    ps     = .psw_plt_ps(x, wcols, bins, title),
    stop(sprintf("Unsupported type: '%s'", type), call. = FALSE))

  plot_size <- switch(type,
                      love   = c(7, 5),
                      ess    = c(7, 4.5),
                      weight = c(8, 2 + 1.6 * ceiling(length(wcols) / 2)),
                      ps     = c(8, 2 + 2.2 * ceiling(length(wcols) / 2)))
  .psw_save(p, plot_size, save)
}
