# =============================================================================
# sens-plt.R -- plots for a get_sens() result
# =============================================================================
#
# Architecture (2 layers):
#
#   L1  plt_sens(x, type, ...)
#         |
#         +-- L2 renderers
#               .sens_plt_contour   lm / dml / iv, upstream base graphics
#                                   wrapped with ggplotify::as.ggplot()
#               .sens_plt_extreme   lm only, same wrapping
#               .sens_plt_tip       cox, native ggplot2
#               .sens_plt_evalue    cox, native ggplot2
#
# The three omitted-variable-bias packages each draw their contour with base
# graphics and each use a different parameterisation, bound placement and
# critical-value line. Redrawing them would mean reimplementing three published
# methods, so they are wrapped rather than reproduced. The two Cox plots have
# no upstream figure at all and are drawn natively.
# =============================================================================

.SENS_CONTOUR_DEFAULTS <- list(
  sensitivity_of = c("auto", "estimate", "t-value", "lwr", "upr"),
  n_levels       = 10L,
  grid_n         = 70L,
  round          = 3L)

.SENS_TYPES <- list(
  lm  = c("contour", "extreme"),
  cox = c("tip", "evalue"),
  dml = "contour",
  iv  = "contour")

# Which bound each backend can put on the contour axis.
.SENS_SENS_OF <- list(
  lm  = c("estimate", "t-value", "lwr", "upr"),
  dml = c("lwr", "upr"),
  iv  = c("lwr", "upr", "t-value"))


# ---- L2 renderers ----------------------------------------------------------

# sensemakr, dml.sensemakr and iv.sensemakr all draw with base graphics, and
# gridGraphics cannot emulate labels on contour lines -- wrapping them with
# ggplotify silently drops every level label, which is the readable content of
# the figure. The lm surface can be recomputed exactly from sensemakr's own
# adjusted_* helpers, so it is drawn natively and keeps its labels. dml and iv
# expose no equivalent grid function, so those two are still wrapped.
# ggplotify captures a base plot by replaying it onto a device. With no device
# open -- a script, a test, a render -- R starts the default one and leaves an
# Rplots.pdf behind in the working directory, so a throwaway null device is
# opened first.
#' @keywords internal
#' @noRd
.sens_grab <- function(draw) {
  if (is.null(grDevices::dev.list())) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  ggplotify::as.ggplot(draw)
}

#' @keywords internal
#' @noRd
.sens_contour_lm <- function(x, threshold, lim, ca, alpha) {
  st <- x$sens$sensitivity_stats
  n  <- max(20L, as.integer(ca$grid_n))
  xs <- seq(0, lim[1L], length.out = n)
  ys <- seq(0, lim[2L], length.out = n)
  g  <- expand.grid(r2dz.x = xs, r2yz.dx = ys)

  surface <- function(r2d, r2y) switch(
    ca$sensitivity_of,
    estimate  = sensemakr::adjusted_estimate(st$estimate, st$se, st$dof,
                                             r2d, r2y),
    `t-value` = sensemakr::adjusted_t(st$estimate, st$se, st$dof, r2d, r2y,
                                      h0 = threshold),
    lwr       = sensemakr::adjusted_ci(st$estimate, st$se, st$dof, r2d, r2y,
                                       which = "lwr", alpha = alpha),
    upr       = sensemakr::adjusted_ci(st$estimate, st$se, st$dof, r2d, r2y,
                                       which = "upr", alpha = alpha))
  zm <- matrix(as.numeric(surface(g$r2dz.x, g$r2yz.dx)), nrow = n)

  crit <- if (identical(ca$sensitivity_of, "t-value")) {
    abs(stats::qt(alpha / 2, st$dof - 1))
  } else threshold

  lines_at <- function(levels) {
    cl <- grDevices::contourLines(xs, ys, zm, levels = levels)
    if (!length(cl)) return(NULL)
    do.call(rbind, lapply(seq_along(cl), function(i) data.frame(
      id = i, level = cl[[i]]$level, x = cl[[i]]$x, y = cl[[i]]$y)))
  }
  lv <- pretty(range(zm, finite = TRUE), ca$n_levels)
  lv <- lv[lv > min(zm, na.rm = TRUE) & lv < max(zm, na.rm = TRUE) &
             abs(lv - crit) > .Machine$double.eps^0.5]
  paths <- lines_at(lv)
  red   <- lines_at(crit)

  labs <- NULL
  if (!is.null(paths)) {
    labs <- do.call(rbind, lapply(split(paths, paths$id), function(d) {
      k <- which.max(d$x)
      data.frame(x = d$x[k], y = d$y[k],
                 label = formatC(d$level[k], format = "g",
                                 digits = ca$round))
    }))
  }

  pts <- data.frame(
    x = 0, y = 0, label = sprintf("Unadjusted
(%s)",
                                  formatC(surface(0, 0), format = "g",
                                          digits = ca$round)))
  bnd <- NULL
  if (!is.null(x$bounds) && nrow(x$bounds)) {
    bnd <- data.frame(
      x = x$bounds$r2_treat, y = x$bounds$r2_out,
      label = sprintf("%s
(%s)", x$bounds$bound_label,
                      formatC(surface(x$bounds$r2_treat, x$bounds$r2_out),
                              format = "g", digits = ca$round)))
    bnd <- bnd[bnd$x <= lim[1L] & bnd$y <= lim[2L], , drop = FALSE]
    if (!nrow(bnd)) bnd <- NULL
  }

  ylab <- switch(ca$sensitivity_of,
                 estimate = "estimate", `t-value` = "t-value",
                 lwr = "lower confidence limit",
                 upr = "upper confidence limit")
  p <- ggplot2::ggplot()
  if (!is.null(paths))
    p <- p + ggplot2::geom_path(data = paths,
                                ggplot2::aes(x = x, y = y, group = id),
                                colour = "grey40", linewidth = 0.4)
  if (!is.null(red))
    p <- p + ggplot2::geom_path(data = red,
                                ggplot2::aes(x = x, y = y, group = id),
                                colour = "red", linetype = 2, linewidth = 0.9)
  if (!is.null(labs))
    p <- p + ggplot2::geom_label(data = labs,
                                 ggplot2::aes(x = x, y = y, label = label),
                                 size = 2.9, linewidth = 0, hjust = 1,
                                 fill = "white", alpha = 0.85)
  p +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = x, y = y),
                        shape = 17, size = 2.6) +
    ggplot2::geom_text(data = pts, ggplot2::aes(x = x, y = y, label = label),
                       size = 3, hjust = -0.1, vjust = -0.3) +
    (if (is.null(bnd)) NULL else list(
      ggplot2::geom_point(data = bnd, ggplot2::aes(x = x, y = y),
                          shape = 18, size = 3.4, colour = "red"),
      ggplot2::geom_text(data = bnd,
                         ggplot2::aes(x = x, y = y, label = label),
                         size = 3, hjust = -0.1, vjust = -0.3))) +
    ggplot2::coord_cartesian(xlim = c(0, lim[1L]), ylim = c(0, lim[2L])) +
    ggplot2::labs(
      x = expression(paste("Partial ", R^2, " of confounder(s) with the treatment")),
      y = expression(paste("Partial ", R^2, " of confounder(s) with the outcome")),
      title = sprintf("Sensitivity of the %s to unmeasured confounding", ylab))
}

#' @keywords internal
#' @noRd
.sens_plt_contour <- function(x, estimand, threshold, lim, contour_args) {
  a <- attr(x, "analysis")
  if (identical(a$method, "lm"))
    return(.sens_contour_lm(x, threshold, lim, contour_args, a$alpha))
  so <- contour_args$sensitivity_of
  draw <- switch(
    a$method,
    dml = function() dml.sensemakr::ovb_contour_plot(
      x$fit, parameter = estimand, which.bound = so, level = a$conf_level,
      rho2 = a$rho2, threshold = threshold,
      lim.x = lim[1L], lim.y = lim[2L],
      nlevels = contour_args$n_levels, grid.number = contour_args$grid_n,
      round = contour_args$round),
    iv = function() iv.sensemakr::ovb_contour_plot(
      x$fit, benchmark_covariates = a$bench_var,
      kz = a$bench_args$k_treat,
      ky = if (is.null(a$bench_args$k_out)) a$bench_args$k_treat
           else a$bench_args$k_out,
      sensitivity.of = so, parm = estimand),
    stop(sprintf("Unsupported method: '%s'", a$method), call. = FALSE))
  .sens_grab(draw)
}

#' @keywords internal
#' @noRd
.sens_plt_extreme <- function(x, threshold, extreme_r2) {
  # As with the contour, the sensemakr method supplies its own threshold, so
  # the numeric method is called directly to keep ours.
  .sens_grab(function() {
    st <- x$sens$sensitivity_stats
    b  <- x$sens$bounds
    sensemakr::ovb_extreme_plot(
      estimate = st$estimate, se = st$se, dof = st$dof,
      r2dz.x = b$r2dz.x, r2yz.dx = extreme_r2, threshold = threshold)
  })
}

# tipr rejects a vectorised confounder-outcome effect (`check_gamma` coerces it
# to a single logical), so the curve is built one grid point at a time.
#' @keywords internal
#' @noRd
.sens_adjust_curve <- function(rr, gamma, ea) {
  vapply(gamma, function(g) {
    out <- if (identical(ea$confounder, "binary")) {
      tipr::adjust_rr_with_binary(rr, ea$exposed_prev, ea$unexposed_prev, g,
                                  verbose = FALSE)
    } else {
      tipr::adjust_rr(rr, ea$smd, g, verbose = FALSE)
    }
    as.numeric(out$rr_adjusted[1L])
  }, numeric(1))
}

#' @keywords internal
#' @noRd
.sens_plt_tip <- function(x, lim, title) {
  a  <- attr(x, "analysis")
  ea <- a$evalue_args
  st <- x$stats
  rr_point <- .sens_hr_to_rr(st$estimate, ea$rare)
  near     <- if (st$estimate < 1) st$conf.high else st$conf.low
  rr_near  <- .sens_hr_to_rr(near, ea$rare)

  g <- seq(lim[1L], lim[2L], length.out = 80L)
  g <- g[g > 0]
  d <- rbind(
    data.frame(gamma = g, adjusted = .sens_adjust_curve(rr_point, g, ea),
               which = "Point estimate", stringsAsFactors = FALSE),
    data.frame(gamma = g, adjusted = .sens_adjust_curve(rr_near, g, ea),
               which = "Confidence limit nearest the null",
               stringsAsFactors = FALSE))
  d$which <- factor(d$which,
                    levels = c("Point estimate",
                               "Confidence limit nearest the null"))
  tip <- data.frame(gamma = st$tip_effect, adjusted = 1)

  lab <- if (identical(ea$confounder, "binary")) {
    sprintf("Confounder-outcome risk ratio\n(prevalence %s vs %s)",
            format(ea$exposed_prev), format(ea$unexposed_prev))
  } else {
    sprintf("Confounder-outcome risk ratio\n(standardised mean difference %s)",
            format(ea$smd))
  }

  ggplot2::ggplot(d, ggplot2::aes(x = gamma, y = adjusted)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2) +
    ggplot2::geom_line(ggplot2::aes(colour = which, linetype = which),
                       linewidth = 0.9) +
    ggplot2::geom_point(data = tip, size = 2.6) +
    ggplot2::geom_vline(xintercept = st$tip_effect, linetype = 3) +
    ggplot2::annotate("text", x = st$tip_effect, y = Inf,
                      label = sprintf("tipping point %.3f", st$tip_effect),
                      hjust = -0.05, vjust = 1.6, size = 3.2) +
    ggplot2::labs(
      x = lab,
      y = sprintf("Adjusted effect (risk-ratio scale, rare = %s)",
                  isTRUE(ea$rare)),
      colour = NULL, linetype = NULL,
      title = if (is.null(title)) sprintf("Tipping point for %s", a$treat) else title) +
    ggplot2::theme(legend.position = "bottom")
}

#' @keywords internal
#' @noRd
.sens_plt_evalue <- function(x, lim, title) {
  a  <- attr(x, "analysis")
  st <- x$stats
  rr <- .sens_hr_to_rr(st$estimate, a$evalue_args$rare)
  rr_near <- .sens_hr_to_rr(
    if (st$estimate < 1) st$conf.high else st$conf.low, a$evalue_args$rare)

  # A bias factor B is exactly explained away by the pairs
  # (RR_EU, RR_UD) with RR_UD = B (RR_EU - 1) / (RR_EU - B), RR_EU > B.
  curve_for <- function(b, label) {
    if (!is.finite(b) || b <= 1) return(NULL)
    xs <- seq(b * 1.001, max(lim[2L], b * 3), length.out = 200L)
    data.frame(rr_eu = xs, rr_ud = b * (xs - 1) / (xs - b),
               which = label, stringsAsFactors = FALSE)
  }
  b_point <- if (rr < 1) 1 / rr else rr
  b_ci    <- if (rr_near < 1) 1 / rr_near else rr_near
  d <- do.call(rbind, Filter(Negate(is.null), list(
    curve_for(b_point, "Point estimate"),
    curve_for(b_ci, "Confidence limit nearest the null"))))
  if (is.null(d) || !nrow(d))
    stop("The estimate is already at the null; there is no E-value curve to draw.",
         call. = FALSE)
  d$which <- factor(d$which,
                    levels = c("Point estimate",
                               "Confidence limit nearest the null"))
  marks <- stats::na.omit(data.frame(
    e = c(st$evalue_point, st$evalue_ci),
    which = factor(c("Point estimate", "Confidence limit nearest the null"),
                   levels = levels(d$which))))
  marks <- marks[marks$e > 1, , drop = FALSE]

  top <- max(c(lim[2L], st$evalue_point * 2), na.rm = TRUE)
  ggplot2::ggplot(d, ggplot2::aes(x = rr_eu, y = rr_ud)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2,
                         colour = "grey60") +
    ggplot2::geom_line(ggplot2::aes(colour = which, linetype = which),
                       linewidth = 0.9) +
    ggplot2::geom_point(data = marks,
                        ggplot2::aes(x = e, y = e, colour = which),
                        size = 2.6, show.legend = FALSE) +
    ggplot2::coord_cartesian(xlim = c(1, top), ylim = c(1, top)) +
    ggplot2::labs(
      x = "Confounder-exposure risk ratio",
      y = "Confounder-outcome risk ratio",
      colour = NULL, linetype = NULL,
      title = if (is.null(title)) sprintf(
        "E-value %.2f (interval %.2f) for %s",
        st$evalue_point, st$evalue_ci, a$treat) else title) +
    ggplot2::theme(legend.position = "bottom")
}


# ---- L1 public entry point -------------------------------------------------

#' Plot a sensitivity analysis
#'
#' Renders the figure that matches the backend behind a [get_sens()] result:
#' a partial \eqn{R^2} contour for the linear, DML and IV backends, and a
#' tipping-point or E-value curve for the Cox backend.
#'
#' @param x A `sens_res` object from [get_sens()].
#' @param type Which figure to draw. `"contour"` and `"extreme"` apply to
#'   `method = "lm"`; `"contour"` alone applies to `"dml"` and `"iv"`;
#'   `"tip"` and `"evalue"` apply to `"cox"`. Left unset, the default is
#'   `"tip"` for the Cox backend and `"contour"` otherwise.
#' @param estimand Length-1 character naming which estimand to draw when the
#'   backend reports several (`"ate"`/`"att"`/`"atu"` for DML, `"iv"`/`"fs"`/
#'   `"rf"` for IV). `NULL` (default) takes the first row of `x$stats`.
#' @param threshold Numeric value of the estimate the contour marks as the
#'   critical line. Default `0`, the null.
#' @param lim Length-2 numeric giving the axis window: the two partial
#'   \eqn{R^2} limits for `"contour"`, and the confounder-association range
#'   for `"tip"` and `"evalue"`. `NULL` (default) uses the window printed by
#'   `print()` on the `sens_res` object.
#' @param extreme_r2 Numeric vector of outcome-side partial \eqn{R^2} values,
#'   one curve each, for `type = "extreme"`. Default `c(1, 0.75, 0.5)`.
#' @param contour_args Named list forwarded to the upstream contour plot.
#'   \describe{
#'     \item{`sensitivity_of`}{Which quantity the contour is drawn for.
#'       `"auto"` (default) resolves to `"estimate"` for `method = "lm"` and
#'       `"lwr"` for `"dml"` and `"iv"`. Otherwise one of `"estimate"`,
#'       `"t-value"`, `"lwr"`, `"upr"`, restricted to what the backend
#'       supports: `"lm"` takes all four, `"dml"` only `"lwr"`/`"upr"`, and
#'       `"iv"` those two plus `"t-value"`.}
#'     \item{`n_levels`}{Integer, default `10L`. Number of contour levels.}
#'     \item{`grid_n`}{Integer, default `70L`. Grid resolution; DML only.}
#'     \item{`round`}{Integer, default `3L`. Digits on the contour labels.}
#'   }
#' @param title Plot title, or `NULL` (default) for a generated one. Used by
#'   `"tip"` and `"evalue"` only; the wrapped upstream plots carry their own.
#' @param save `NULL` or a list with `filename`, `width` and `height`, passed
#'   to `RegR::save_plt()` for PDF output. `list()` and `NULL` skip saving; a
#'   list naming only the file is completed with this figure's pinned size.
#'   Do not include the `plot` argument; it is supplied internally.
#'
#' @section Native figures versus wrapped ones:
#' All four backends draw with base graphics upstream, and `gridGraphics`,
#' which [ggplotify::as.ggplot()] relies on, cannot emulate labels on contour
#' lines. Wrapping therefore loses every contour level, which is the readable
#' content of the figure. Three of the four types avoid that:
#' \describe{
#'   \item{`"contour"` for `method = "lm"`}{Drawn natively from sensemakr's
#'     own [sensemakr::adjusted_estimate()], [sensemakr::adjusted_t()] and
#'     [sensemakr::adjusted_ci()], so the surface is identical to the upstream
#'     figure but the levels stay labelled and the result is a real `ggplot`
#'     that accepts further layers, scales and themes.}
#'   \item{`"tip"` and `"evalue"`}{Drawn natively; there is no upstream
#'     figure for either.}
#'   \item{`"contour"` for `"dml"`/`"iv"`, and `"extreme"`}{Wrapped with
#'     [ggplotify::as.ggplot()], because those surfaces are not exposed as
#'     grid functions. The result is a `ggplot` that saves and composes, but
#'     its contents are a fixed grob: further layers and themes do not reach
#'     it, and the DML and IV contours carry no level labels. The IV contour
#'     additionally ignores `lim`, `threshold` and every `contour_args` field,
#'     since [iv.sensemakr::ovb_contour_plot()] accepts none of them.}
#' }
#'
#' @return A `ggplot` object, carrying its pinned size in
#'   `attr(p, "plot_size")`. If `save` is non-empty, the same plot is also
#'   written to PDF through `RegR::save_plt()`.
#'
#' @seealso [get_sens()].
#'
#' @examplesIf requireNamespace("sensemakr", quietly = TRUE) && requireNamespace("ggplotify", quietly = TRUE)
#' \donttest{
#' res <- get_sens(sensemakr::darfur,
#'                 treat = "directlyharmed", outcome = "peacefactor",
#'                 adj_var = c("age", "farmer_dar", "herder_dar", "pastvoted",
#'                             "hhsize_darfur", "female", "village"),
#'                 bench_var = "female", method = "lm",
#'                 bench_args = list(k_treat = 1:3))
#' plt_sens(res, type = "contour")
#' }
#'
#' @export
plt_sens <- function(x,
                     type         = c("contour", "extreme", "tip", "evalue"),
                     estimand     = NULL,
                     threshold    = 0,
                     lim          = NULL,
                     extreme_r2   = c(1, 0.75, 0.5),
                     contour_args = list(sensitivity_of = c("auto", "estimate",
                                                            "t-value", "lwr",
                                                            "upr"),
                                         n_levels       = 10L,
                                         grid_n         = 70L,
                                         round          = 3L),
                     title        = NULL,
                     save         = list()) {

  if (!inherits(x, "sens_res"))
    stop("`x` must be a `sens_res` object from get_sens().", call. = FALSE)
  a <- attr(x, "analysis")
  spec <- .sens_plt_spec(x)

  # The default type depends on the backend, so it is resolved before, not by,
  # match.arg().
  if (missing(type)) type <- spec$type
  type <- match.arg(type)
  ok <- .SENS_TYPES[[a$method]]
  if (!type %in% ok)
    stop(sprintf("type = \"%s\" is not available for method = \"%s\"; available types are %s.",
                 type, a$method,
                 paste0("\"", ok, "\"", collapse = ", ")), call. = FALSE)

  if (!is.null(lim) && (!is.numeric(lim) || length(lim) != 2L || anyNA(lim)))
    stop("`lim` must be `NULL` or two numbers.", call. = FALSE)
  if (is.null(lim)) lim <- spec$lim
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold))
    stop("`threshold` must be a single number.", call. = FALSE)

  contour_args <- .merge_named_arg(contour_args, .SENS_CONTOUR_DEFAULTS,
                                   "contour_args")
  if (identical(type, "contour")) {
    allowed <- .SENS_SENS_OF[[a$method]]
    so <- contour_args$sensitivity_of
    if (length(so) > 1L) so <- so[1L]
    so <- match.arg(so, c("auto", allowed))
    contour_args$sensitivity_of <-
      if (identical(so, "auto")) allowed[1L] else so
  }

  if (is.null(estimand)) {
    estimand <- x$stats$estimand[1L]
  } else if (!estimand %in% x$stats$estimand) {
    stop(sprintf("`estimand` must be one of %s.",
                 paste0("\"", x$stats$estimand, "\"", collapse = ", ")),
         call. = FALSE)
  }
  x$stats <- x$stats[x$stats$estimand == estimand, , drop = FALSE]

  if (type %in% c("contour", "extreme")) {
    for (pkg in c("ggplotify", a$backend)) {
      if (!requireNamespace(pkg, quietly = TRUE))
        stop(sprintf("Package '%s' is required for plt_sens(type = \"%s\")",
                     pkg, type), call. = FALSE)
    }
  } else if (!requireNamespace("tipr", quietly = TRUE)) {
    stop(sprintf("Package 'tipr' is required for plt_sens(type = \"%s\")",
                 type), call. = FALSE)
  }

  p <- switch(
    type,
    contour = .sens_plt_contour(x, estimand, threshold, lim, contour_args),
    extreme = .sens_plt_extreme(x, threshold, extreme_r2),
    tip     = .sens_plt_tip(x, lim, title),
    evalue  = .sens_plt_evalue(x, lim, title),
    stop(sprintf("Unsupported type: '%s'", type), call. = FALSE))

  plot_size <- if (type %in% c("contour", "extreme")) c(7, 6) else c(7.5, 5.5)
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
