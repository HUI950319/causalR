# =============================================================================
# hte-save.R -- the heterogeneity figure of a causal forest, ready to save
# =============================================================================
#
# save_hte() stacks a top panel -- the variable importance bars of
# MLR::plt_bar_per() or the CATE waterfall of plt_hte_cate() -- over the
# subgroup forest plot of plt_hte_sub(), and sizes the page from the forest
# plot, which is pinned to the size its text needs.
# =============================================================================


#' Save the heterogeneity figure of a causal forest
#'
#' Stacks a top panel over the subgroup forest plot of [plt_hte_sub()]: the
#' variable importance bars of `x$importance` drawn by `MLR::plt_bar_per()`,
#' or the CATE waterfall of [plt_hte_cate()]. The page takes the width of the
#' forest plot, so it is no wider than the table needs, and is written to PDF
#' through `RegR::save_plt()` when `save` is given.
#'
#' @param x An `hte_res` object from [get_hte()].
#' @param sub_var A single categorical covariate: the subgroups of the forest
#'   plot, and with `panel = "auto"` the variable whose importance picks the
#'   top panel.
#' @param panel The top panel. `"auto"` (default) draws the importance bars
#'   when `sub_var` is among the three most important variables of
#'   `x$importance`, and the CATE waterfall otherwise; `"importance"` and
#'   `"cate"` draw one of them.
#' @param name_map `NULL` (default) or the display names of variables and
#'   levels, passed as `name_map` to [plt_hte_sub()] and as `add_name` to
#'   `MLR::plt_bar_per()`.
#' @param sub_args Named list of arguments of [plt_hte_sub()], other than
#'   `x`, `sub_var`, `name_map` and `save`. A partial list is merged over the
#'   defaults written out in the usage, so only the fields given change;
#'   unknown or duplicated fields are rejected.
#'   \describe{
#'     \item{`show_n`}{`TRUE` (default): the patients per arm.}
#'     \item{`fixed_size`}{`TRUE` (default): the plot is pinned to the size its
#'       text needs, which sets the page width.}
#'     \item{`row_height`}{Height of a table row in inches, default `0.5`.}
#'     \item{`ci_shape`}{Default `c("normal", "summary")`: squares by inverse
#'       variance for the subgroups, a summary diamond for all patients.}
#'     \item{`header_list`, `xlab`}{For a difference of survival
#'       probabilities (the default `measure = "diff"` of a survival forest
#'       of `S(time)`), the count and effect headers and the axis title of
#'       [plt_hte_sub()] break onto two lines, e.g. `"N\n(PORT / no PORT)"`,
#'       `"5-year DSS\ndifference, %"` and `"Difference in 5-year DSS
#'       (%),\nPORT - no PORT"`, so the table is narrower. A `header_list`
#'       given here is merged over these headers; `xlab = NULL` restores the
#'       one-line title. Other measures and outcomes keep the one-line
#'       defaults.}
#'     \item{others}{Any other argument of [plt_hte_sub()], e.g. `box_col`
#'       for the colours or `axis_args = list(xlim = c(-0.15, 0.15))` for an
#'       axis centred on 0.}
#'   }
#' @param imp_args Named list of arguments of `MLR::plt_bar_per()` for the
#'   importance bars, other than `data`, `add_name` and `save`, merged over
#'   the defaults in the usage as `sub_args` is. Fields are checked against
#'   `MLR::plt_bar_per()` only when the importance bars are drawn.
#'   \describe{
#'     \item{`mode`}{`"bar_per2"` (default): bars with two rings of shares.}
#'     \item{`per_position`}{Default `c(0.26, 0.08, 0.93, 0.93)`: left,
#'       bottom, width and height of the rings, as fractions of the panel.}
#'     \item{`labs_args`}{Default `list(x = "Variable importance", y = NULL)`.}
#'     \item{others}{Any other argument of `MLR::plt_bar_per()`, e.g.
#'       `bar_per_palette` or `fmt_args`.}
#'   }
#' @param cate_args Named list of arguments of [plt_hte_cate()] for the CATE
#'   waterfall, other than `x` and `save`; `list()` (default) keeps all of
#'   its defaults. The caption of [plt_hte_cate()] is dropped, to be written
#'   in the figure legend.
#' @param top_height `NULL` (default) or a single positive number: the height
#'   of the top panel in inches. `NULL` gives 4.1 in to the importance bars
#'   and 3.6 in to the CATE waterfall.
#' @param save `NULL` or a list with `filename`, `width` and `height`, passed
#'   to `RegR::save_plt()` for PDF output. `list()` and `NULL` skip saving; a
#'   `width` or `height` left out is taken from `attr(p, "plot_size")`. Do not
#'   include the `plot` argument; it is supplied internally.
#'
#' @section Page size:
#' The width is that of the forest plot, `attr(, "plot_size")` of
#' [plt_hte_sub()], plus 0.5 in for the panel tags and margins under the
#' importance bars, or 1.2 in under the CATE waterfall, whose axis title and
#' tick labels push the forest plot to the right. The height is the top
#' panel, the forest plot and 0.6 in of margins. Each panel has a margin of 20
#' points at the top and left, and the panels are tagged A and B outside their
#' top left corner.
#'
#' @return A \pkg{patchwork} of the two panels, with `attr(p, "plot_size")`
#'   the `c(width, height)` of the page in inches. If `save` is non-empty,
#'   the figure is also written to PDF through `RegR::save_plt()`.
#'
#' @seealso [plt_hte_sub()], [plt_hte_cate()]; `MLR::plt_bar_per()`.
#'
#' @examplesIf requireNamespace("grf", quietly = TRUE) && requireNamespace("forestplot", quietly = TRUE) && requireNamespace("ggplotify", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)
#' \donttest{
#' set.seed(20261115)
#' n <- 600
#' d <- data.frame(age = round(runif(n, 20, 85)),
#'                 sex = factor(sample(c("F", "M"), n, replace = TRUE)),
#'                 x1  = rnorm(n))
#' d$PORT <- factor(sample(c("No", "Yes"), n, replace = TRUE),
#'                  levels = c("No", "Yes"))
#' ev   <- rexp(n, 0.02 * exp(-0.8 * (d$PORT == "Yes") * (d$sex == "M")))
#' cens <- pmin(rexp(n, 0.01), 120)
#' d$time <- pmin(ev, cens)
#' d$DSS  <- as.integer(ev <= cens)
#' res <- get_hte(d, cat_var = "PORT", sub_var = "sex", adj_var = c("age", "x1"),
#'                surv = TRUE, time = 60,
#'                grf_args = list(num.trees = 500, seed = 1))
#'
#' # The CATE waterfall over the forest plot, the axis centred on 0
#' p <- save_hte(res, sub_var = "sex", panel = "cate",
#'               sub_args = list(axis_args = list(xlim = c(-0.3, 0.3))))
#' attr(p, "plot_size")
#'
#' # The importance bars on top, written to PDF at the page size
#' if (requireNamespace("MLR", quietly = TRUE) &&
#'     requireNamespace("RegR", quietly = TRUE))
#'   save_hte(res, sub_var = "sex", panel = "importance",
#'            save = list(filename = file.path(tempdir(), "Supplementary Figure 5")))
#' }
#'
#' @export
save_hte <- function(x,
                     sub_var,
                     panel      = c("auto", "importance", "cate"),
                     name_map   = NULL,
                     sub_args   = list(show_n = TRUE, fixed_size = TRUE,
                                       row_height = 0.5,
                                       ci_shape = c("normal", "summary")),
                     imp_args   = list(mode = "bar_per2",
                                       per_position = c(0.26, 0.08, 0.93, 0.93),
                                       labs_args = list(x = "Variable importance",
                                                        y = NULL)),
                     cate_args  = list(),
                     top_height = NULL,
                     save       = list()) {

  if (!inherits(x, "hte_res"))
    stop("`x` must be an `hte_res` object from get_hte().", call. = FALSE)
  if (!is.character(sub_var) || length(sub_var) != 1L || is.na(sub_var) ||
      !nzchar(sub_var))
    stop("`sub_var` must be a single column name.", call. = FALSE)
  panel <- match.arg(panel)
  if (!is.null(top_height) &&
      (!is.numeric(top_height) || length(top_height) != 1L ||
       !is.finite(top_height) || top_height <= 0))
    stop("`top_height` must be `NULL` or a single positive height in inches.",
         call. = FALSE)
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required for save_hte().", call. = FALSE)

  # A partial list changes only the fields it names; `x`, `sub_var`,
  # `name_map` and `save` belong to save_hte() itself
  usage    <- formals(save_hte)
  sub_args <- .merge_named_arg(
    sub_args, eval(usage$sub_args), "sub_args",
    setdiff(names(formals(plt_hte_sub)), c("x", "sub_var", "name_map", "save", "...")))
  cate_args <- .merge_named_arg(
    cate_args, list(), "cate_args",
    setdiff(names(formals(plt_hte_cate)), c("x", "save")))

  imp <- x$importance
  if (panel == "auto") {
    top3  <- utils::head(imp$variable[order(-imp$importance)], 3L)
    panel <- if (sub_var %in% top3) "importance" else "cate"
  }
  if (panel == "importance") {
    if (!requireNamespace("MLR", quietly = TRUE))
      stop("Package 'MLR' is required for the importance bars of save_hte(); ",
           "install it with remotes::install_github(\"HUI950319/MLR\"), ",
           "or use `panel = \"cate\"`.", call. = FALSE)
    imp_args <- .merge_named_arg(
      imp_args, eval(usage$imp_args), "imp_args",
      setdiff(names(formals(MLR::plt_bar_per)), c("data", "add_name", "save")))
  }

  # Two-line headers and axis title for a difference of survival
  # probabilities, written as plt_hte_sub() writes them on one line
  a       <- attr(x, "analysis")
  measure <- if (is.null(sub_args[["measure"]])) "diff" else sub_args[["measure"]][1L]
  pct     <- sub_args[["axis_args"]][["percent"]]
  pct     <- is.null(pct) || isTRUE(pct)
  if (identical(a$outcome_type, "survival") && !identical(a$target, "RMST") &&
      identical(measure, "diff")) {
    ref  <- setdiff(levels(factor(x$data[[a$cat_var]])), a$treated)[1L]
    arms <- .hte_arm_names(a$cat_var, a$treated, ref)
    t_y  <- format(a$time / 12)
    two  <- list(
      header_list = list(count  = sprintf("N\n(%s / %s)", arms[1L], arms[2L]),
                         effect = sprintf("%s-year DSS\ndifference%s", t_y,
                                          if (pct) ", %" else "")),
      xlab = sprintf("Difference in %s-year DSS%s,\n%s - %s", t_y,
                     if (pct) " (%)" else "", arms[1L], arms[2L]))
    sub_args <- utils::modifyList(two, sub_args, keep.null = TRUE)
  }

  # `x` goes in by name, so the calls do not carry the whole result
  p_sub <- do.call(plt_hte_sub,
                   c(list(quote(x), sub_var = sub_var, name_map = name_map),
                     sub_args))
  if (panel == "importance") {
    p_top <- do.call(MLR::plt_bar_per,
                     c(list(quote(imp), add_name = name_map), imp_args))
    h_top <- 4.1
    w_pad <- 0.5
  } else {
    p_top <- do.call(plt_hte_cate, c(list(quote(x)), cate_args)) +
      ggplot2::labs(caption = NULL)
    h_top <- 3.6
    # the waterfall's axis title and tick labels push the forest plot right
    w_pad <- 1.2
  }
  if (!is.null(top_height)) h_top <- top_height

  sz <- attr(p_sub, "plot_size")
  p  <- (p_top + p_sub +
           patchwork::plot_layout(ncol = 1, heights = c(h_top, sz[[2L]])) &
           ggplot2::theme(plot.margin = ggplot2::margin(20, 5.5, 5.5, 20)))
  p  <- UtilsR::fmt_strip(p, strip = FALSE)
  p  <- UtilsR::fmt_tag(p, label_position = "tl-out")

  size <- c(width = sz[[1L]] + w_pad, height = h_top + sz[[2L]] + 0.6)
  attr(p, "plot_size") <- size
  if (!is.null(save) && length(save) > 0L) {
    if (!requireNamespace("RegR", quietly = TRUE))
      stop("Package 'RegR' is required for a non-empty `save`.", call. = FALSE)
    if (is.null(save[["width"]]))  save$width  <- size[[1L]]
    if (is.null(save[["height"]])) save$height <- size[[2L]]
    do.call(RegR::save_plt, c(list(plot = p), save))
  }
  p
}
