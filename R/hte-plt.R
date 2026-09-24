# =============================================================================
# hte-plt.R -- plots for get_hte()
# =============================================================================
#
# Architecture:
#
#   L1  plt_hte_dep()     one panel per covariate ("dep") or a two-covariate
#                         partial-dependence heat map ("heat")
#   L1  plt_hte_sub()     subgroup forest plot (forestplot, wrapped as a ggplot)
#                         of the doubly robust subgroup ATEs, recomputed with
#                         .hte_subgroup()
#   L1  plt_hte_cate()    waterfall / density of the per-patient CATE by
#                         subgroup, with the doubly robust ATE lines
#   L2  .hte_pdp()        forest CATE averaged with covariates set to grid values
#   L2  .hte_beyond()     survival patients followed past `time`
#   L2  .hte_cate_label() axis label of the CATE
#   L2  .hte_sub_var()    `sub_var` checks of plt_hte_sub() and plt_hte_cate()
#   L2  .hte_muffle_ps()  silences grf's overlap warning in recomputed estimates
#
# The doubly robust layer comes from .hte_dr_var() in hte-get.R, so the p_het
# in the strips is the one in get_hte()$importance at the default spline df.
#
# The ggplots here and in hte-rate.R take UtilsR::theme_my(base_rect_size =
# 1.5), the theme of MLR::plt_bar_per(). theme_my() reads its colours from
# ggprism, which UtilsR only suggests, so causalR imports both.
# =============================================================================


# ---- L2 partial dependence -------------------------------------------------

# Partial dependence as in StratifiedMedicine::plot_dependence(): every
# analysed row -- or an evenly spaced subset of at most `max_n`, so the result
# does not depend on the random seed -- gets the covariates in `vars` set to
# each grid combination, and the forest's CATE is averaged. A factor is set
# through its one-hot columns, which get_hte() keeps for every level, so no
# row ever carries two levels at once -- or, with factor_encoding =
# "integer", through its one column of level codes; a numeric covariate is
# one column.
#' @keywords internal
#' @noRd
.hte_pdp <- function(x, vars, grid_n, max_n) {
  fit  <- x$fit
  X    <- fit$X.orig
  d    <- x$data
  src  <- attr(x, "analysis")$covariates[attr(X, "assign")]
  rows <- if (nrow(X) > max_n) {
    unique(round(seq(1, nrow(X), length.out = max_n)))
  } else {
    seq_len(nrow(X))
  }
  X0 <- X[rows, , drop = FALSE]

  grid <- lapply(stats::setNames(vars, vars), function(v) {
    xv <- d[[v]]
    if (.hte_is_num(xv)) seq(min(xv, na.rm = TRUE), max(xv, na.rm = TRUE),
                             length.out = grid_n)
    else if (is.numeric(xv)) sort(unique(xv))
    else levels(droplevels(as.factor(xv)))
  })
  combo <- expand.grid(grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

  big <- do.call(rbind, lapply(seq_len(nrow(combo)), function(k) {
    Xk <- X0
    for (v in vars) {
      cols <- which(src == v)
      if (is.numeric(d[[v]])) {
        Xk[, cols] <- combo[[v]][k]
      } else if (identical(attr(x, "analysis")$factor_encoding, "integer")) {
        Xk[, cols] <- match(combo[[v]][k], grid[[v]])
      } else {
        Xk[, cols] <- 0
        Xk[, cols[match(combo[[v]][k], grid[[v]])]] <- 1
      }
    }
    Xk
  }))
  pred <- as.numeric(stats::predict(fit, big)$predictions)
  combo$estimate <- colMeans(matrix(pred, nrow = nrow(X0)))

  for (v in vars)
    if (!.hte_is_num(d[[v]]))
      combo[[v]] <- factor(as.character(combo[[v]]),
                           levels = as.character(grid[[v]]))
  attr(combo, "n_rows") <- nrow(X0)
  combo
}

# Survival patients followed past `time`, flagged by the rule get_hte()
# applies -- beyond `time` for S(t), up to it for RMST -- from the original
# times in `$data`, since grf keeps only the times cut at the horizon. `NULL`
# for other outcomes.
#' @keywords internal
#' @noRd
.hte_beyond <- function(x) {
  a <- attr(x, "analysis")
  if (!identical(a$outcome_type, "survival")) return(NULL)
  y <- x$data[[a$outcome[1L]]]
  if (identical(a$target, "RMST")) y >= a$time else y > a$time
}

# Axis label of the CATE, on the "diff" scale of get_hte().
#' @keywords internal
#' @noRd
.hte_cate_label <- function(x) {
  a    <- attr(x, "analysis")
  ref  <- setdiff(levels(factor(x$data[[a$cat_var]])), a$treated)[1L]
  what <- switch(a$outcome_type,
                 survival   = sprintf("%s(%s)",
                                      if (identical(a$target, "RMST")) "RMST" else "S",
                                      format(a$time)),
                 binary     = "risk",
                 continuous = "mean")
  sprintf("CATE: %s %s - %s", what, a$treated, ref)
}

# The `sub_var` a caller names, checked against `x$data`: categorical columns
# other than the exposure, with a note for any the forest does not condition
# on.
#' @keywords internal
#' @noRd
.hte_sub_var <- function(sub_var, d, a) {
  if (!is.character(sub_var) || anyNA(sub_var))
    stop("`sub_var` must be `NULL` or column names.", call. = FALSE)
  sub_var <- unique(setdiff(sub_var, a$cat_var))
  miss <- setdiff(sub_var, names(d))
  if (length(miss))
    stop(sprintf("`sub_var` names no column of `x$data`: %s.",
                 paste0("`", miss, "`", collapse = ", ")), call. = FALSE)
  num <- sub_var[vapply(d[sub_var], .hte_is_num, logical(1L))]
  if (length(num))
    stop(sprintf("%s %s continuous; cut it into groups first, for example with cut().",
                 paste0("`", num, "`", collapse = ", "),
                 if (length(num) == 1L) "is" else "are"), call. = FALSE)
  outside <- setdiff(sub_var, a$covariates)
  if (length(outside))
    cli::cli_inform(c("i" = paste(
      "{.field {outside}} {?is/are} not a forest covariate: the subgroup",
      "estimates hold for a function of the covariates or a variable that",
      "is no confounder; otherwise add it with get_hte(sub_var = ).")))
  sub_var
}

# get_hte() already reported grf's overlap warning once for this forest, so
# estimates recomputed from it leave it out.
#' @keywords internal
#' @noRd
.hte_muffle_ps <- function(expr)
  withCallingHandlers(expr, warning = function(w) {
    if (startsWith(conditionMessage(w), "Estimated treatment propensities"))
      invokeRestart("muffleWarning")
  })


# ---- L1 public entry point -------------------------------------------------

#' Dependence plots for heterogeneous treatment effects
#'
#' Draws how the conditional average treatment effect (CATE) estimated by
#' [get_hte()] varies with each covariate. The layers follow
#' `StratifiedMedicine::plot_dependence()` (partial dependence and a two-way
#' heat map) and the panel layout of `MLR::plt_shp_TM()`, with a doubly robust
#' layer that carries valid confidence intervals.
#'
#' @param x An `hte_res` object from [get_hte()].
#' @param x_var Covariates to draw: `NULL` (default) for every covariate in
#'   the order of `x$importance`, a character vector of covariate names, or
#'   `"fct"` / `"num"` for only the categorical / continuous ones. A numeric
#'   covariate with more than 5 distinct values counts as continuous.
#'   `type = "heat"` needs exactly two names. A patient missing a covariate
#'   is left out of its panel.
#' @param type `"dep"` (default) draws one panel per covariate; `"heat"` draws
#'   the partial dependence of two covariates as a heat map.
#' @param display Layers of a `"dep"` panel, any of
#'   \describe{
#'     \item{`"cate"`}{The out-of-bag CATE of every patient (grey), jittered
#'       for a categorical covariate and smoothed by loess (span
#'       `cate_smooth`) for a continuous one. Descriptive only: forest
#'       estimates are shrunk towards the overall mean.}
#'     \item{`"dr"`}{The doubly robust estimate with pointwise confidence
#'       intervals (red): the AIPW score mean per level, or a natural spline
#'       of the AIPW scores with HC3 errors. The strip adds `p_het`, the Wald
#'       test that the CATE does not vary with the covariate. Levels with
#'       fewer than two patients in either arm are left out, and for a
#'       survival outcome so are levels in which an arm has no patient
#'       followed beyond `time`; the spline of a continuous covariate covers
#'       only the values that patients followed beyond `time` reach in both
#'       arms. With nothing left the panel has no red layer and `p_het` is
#'       `NA`.}
#'     \item{`"pdp"`}{The partial dependence of the forest (blue): its CATE
#'       averaged over the patients with the covariate set to each level or
#'       grid value, other covariates as observed.}
#'   }
#'   Default `c("cate", "dr")`. Only used by `type = "dep"`.
#' @param conf_level Confidence level of the `"dr"` intervals. Default `0.95`.
#'   Only used by `type = "dep"`.
#' @param ylim `NULL` (default) or two increasing numbers giving the y range
#'   of every panel; intervals running past it are clipped.
#' @param cate_smooth Loess span of the `"cate"` line of a continuous
#'   covariate, from `0.05` to `1`: smaller follows the patients more closely,
#'   larger is smoother. Default `0.6`. `0` leaves the line out and keeps the
#'   points. The line is dark grey beside the red `"dr"` layer and red without
#'   it. Only used by `type = "dep"`.
#' @param dr_args Named list for the `"dr"` layer: `spline_df` (default `2`),
#'   the natural-spline degrees of freedom of a continuous covariate, which
#'   also sets the degrees of freedom of its `p_het`: fewer is smoother, and
#'   `1` fits a straight line. At the default the `p_het` equals the one in
#'   `get_hte()$importance`. Only used by `type = "dep"`.
#' @param pdp_args Named list for the partial dependence: `grid_n` (default
#'   `21`), the grid points of a continuous covariate, and `max_n` (default
#'   `1000`), the most patients averaged over, taken evenly spaced so the
#'   result does not depend on the random seed; `Inf` uses everyone. The
#'   forest predicts one row per grid point and patient -- for a heat map, per
#'   grid combination and patient -- so large values are slow.
#' @param axis_arg Named list with `share_y`: `"all"` (default) for one y
#'   range across panels or `"none"` for a range per panel. Without `ylim` the
#'   range covers the points and estimates but not the intervals, so a level
#'   with very few patients cannot flatten every other panel. Only used by
#'   `type = "dep"`.
#' @param title Plot title, or `NULL` (default).
#' @param save `NULL` or a list with `filename`, `width` and `height`, passed
#'   to `RegR::save_plt()` for PDF output. `list()` and `NULL` skip saving; a
#'   list naming only the file is completed with this figure's pinned size.
#'   Do not include the `plot` argument; it is supplied internally.
#'
#' @section Three views of the same CATE:
#' The layers answer different questions. The `"dr"` layer estimates the mean
#' effect among patients with a given covariate value, averaging over the
#' other covariates as they occur in that group, and is the only one with
#' confidence intervals. The `"pdp"` layer shows the shape the forest has
#' learnt when only this covariate changes and the others keep their observed
#' distribution. The `"cate"` points are the forest's own per-patient
#' estimates. When covariates are correlated the three can differ.
#'
#' @return A `ggplot` for a single covariate or a heat map, otherwise a
#'   patchwork with one panel per covariate, three per row. The y axis (the
#'   fill for a heat map) is the CATE on the `"diff"` scale of [get_hte()]:
#'   the S(t) or RMST difference for a survival outcome, the risk difference
#'   for a binary one and the mean difference otherwise. The pinned size is in
#'   `attr(p, "plot_size")`. If `save` is non-empty, the same plot is also
#'   written to PDF through `RegR::save_plt()`.
#'
#' @seealso [get_hte()]; `MLR::plt_bar_per()` for `x$importance`.
#'
#' @examplesIf requireNamespace("grf", quietly = TRUE) && requireNamespace("sandwich", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)
#' \donttest{
#' set.seed(20260923)
#' n <- 600
#' d <- data.frame(age   = round(runif(n, 20, 85)),
#'                 sex   = factor(sample(c("F", "M"), n, replace = TRUE)),
#'                 stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
#' d$z <- rbinom(n, 1, 0.5)
#' d$y <- rbinom(n, 1, plogis(-1 + 0.02 * (d$age - 50) +
#'                              d$z * (0.2 + 0.6 * (d$sex == "M"))))
#' res <- get_hte(d, cat_var = "z", adj_var = c("age", "sex", "stage"),
#'                surv = "y", grf_args = list(num.trees = 500, seed = 1))
#'
#' # Every covariate, ordered by importance
#' plt_hte_dep(res)
#'
#' # One covariate with all three layers
#' plt_hte_dep(res, x_var = "age", display = c("cate", "dr", "pdp"))
#'
#' # Two-way partial dependence
#' plt_hte_dep(res, x_var = c("age", "sex"), type = "heat")
#' }
#'
#' @export
plt_hte_dep <- function(x,
                        x_var    = NULL,
                        type     = c("dep", "heat"),
                        display  = c("cate", "dr"),
                        conf_level = 0.95,
                        ylim     = NULL,
                        cate_smooth = 0.6,
                        dr_args  = list(spline_df = 2),
                        pdp_args = list(grid_n = 21, max_n = 1000),
                        axis_arg = list(share_y = "all"),
                        title    = NULL,
                        save     = list()) {

  if (!inherits(x, "hte_res"))
    stop("`x` must be an `hte_res` object from get_hte().", call. = FALSE)
  type <- match.arg(type)
  if (type == "heat") {
    used <- c(display = !missing(display), conf_level = !missing(conf_level),
              cate_smooth = !missing(cate_smooth), dr_args = !missing(dr_args),
              axis_arg = !missing(axis_arg))
    if (any(used))
      stop(sprintf("%s only applies to type = \"dep\"; the heat map shows the partial dependence.",
                   paste0("`", names(used)[used], "`", collapse = ", ")),
           call. = FALSE)
  }
  display  <- match.arg(display, c("cate", "dr", "pdp"), several.ok = TRUE)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  if (!is.numeric(cate_smooth) || length(cate_smooth) != 1L ||
      is.na(cate_smooth) ||
      !(cate_smooth == 0 || (cate_smooth >= 0.05 && cate_smooth <= 1)))
    stop("`cate_smooth` must be 0 (no line) or a loess span from 0.05 to 1.",
         call. = FALSE)
  dr_args  <- .merge_named_arg(dr_args, list(spline_df = 2), "dr_args")
  pdp_args <- .merge_named_arg(pdp_args, list(grid_n = 21, max_n = 1000),
                               "pdp_args")
  axis_arg <- .merge_named_arg(axis_arg, list(share_y = "all"), "axis_arg")
  if (!is.character(axis_arg$share_y) || length(axis_arg$share_y) != 1L ||
      !axis_arg$share_y %in% c("all", "none"))
    stop("`axis_arg$share_y` must be \"all\" or \"none\".", call. = FALSE)
  if (!is.numeric(dr_args$spline_df) || length(dr_args$spline_df) != 1L ||
      is.na(dr_args$spline_df) || dr_args$spline_df < 1)
    stop("`dr_args$spline_df` must be a single number of at least 1.",
         call. = FALSE)
  if (!is.numeric(pdp_args$grid_n) || length(pdp_args$grid_n) != 1L ||
      is.na(pdp_args$grid_n) || pdp_args$grid_n < 2 ||
      !is.numeric(pdp_args$max_n) || length(pdp_args$max_n) != 1L ||
      is.na(pdp_args$max_n) || pdp_args$max_n < 1)
    stop("`pdp_args$grid_n` must be at least 2 and `pdp_args$max_n` at least 1.",
         call. = FALSE)
  if (!is.null(ylim) && (!is.numeric(ylim) || length(ylim) != 2L ||
                         anyNA(ylim) || ylim[1L] >= ylim[2L]))
    stop("`ylim` must be `NULL` or two increasing numbers.", call. = FALSE)
  if (!is.null(save) && !is.list(save))
    stop("`save` must be `NULL` or a list.", call. = FALSE)

  # ---- Covariates, in importance order -------------------------------------
  a        <- attr(x, "analysis")
  d        <- x$data
  vars_all <- x$importance$variable
  num      <- vapply(d[vars_all], .hte_is_num, logical(1L))
  vars <- if (is.null(x_var)) {
    vars_all
  } else if (identical(x_var, "fct")) {
    vars_all[!num]
  } else if (identical(x_var, "num")) {
    vars_all[num]
  } else {
    if (!is.character(x_var) || anyNA(x_var))
      stop("`x_var` must be `NULL`, covariate names, \"fct\" or \"num\".",
           call. = FALSE)
    miss <- setdiff(x_var, vars_all)
    if (length(miss))
      stop(sprintf("`x_var` names no covariate of `x`: %s. Covariates are %s.",
                   paste0("`", miss, "`", collapse = ", "),
                   paste0("`", vars_all, "`", collapse = ", ")), call. = FALSE)
    unique(x_var)
  }
  if (!length(vars))
    stop("No covariate of `x` matches `x_var`.", call. = FALSE)
  if (type == "heat" && length(vars) != 2L)
    stop("type = \"heat\" needs exactly two covariates in `x_var`.",
         call. = FALSE)

  if ((type == "heat" || "pdp" %in% display) &&
      !requireNamespace("grf", quietly = TRUE))
    stop("Package 'grf' is required for the partial dependence.", call. = FALSE)
  if (type == "dep" && "dr" %in% display && any(num[vars]) &&
      !requireNamespace("sandwich", quietly = TRUE))
    stop("Package 'sandwich' is required for the spline of a continuous covariate.",
         call. = FALSE)
  if (type == "dep" && length(vars) > 1L &&
      !requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required to combine several panels.",
         call. = FALSE)

  ylab <- .hte_cate_label(x)
  ate  <- x$stats$estimate[x$stats$estimand == "ATE" &
                             x$stats$measure == "diff"]

  # ---- Heat map -------------------------------------------------------------
  if (type == "heat") {
    pd    <- .hte_pdp(x, vars, pdp_args$grid_n, pdp_args$max_n)
    pd$x1 <- pd[[vars[1L]]]
    pd$x2 <- pd[[vars[2L]]]
    p <- ggplot2::ggplot(pd, ggplot2::aes(x = x1, y = x2, fill = estimate)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "navy", mid = "white",
                                    high = "firebrick", midpoint = 0,
                                    name = sub(": ", "\n", ylab, fixed = TRUE)) +
      ggplot2::labs(x = vars[1L], y = vars[2L], title = title,
                    caption = sprintf("Partial dependence of the forest CATE, averaged over %d patients",
                                      attr(pd, "n_rows"))) +
      UtilsR::theme_my(base_rect_size = 1.5) +
      # theme_my() drops legend titles, and the fill needs its own
      ggplot2::theme(legend.title = ggplot2::element_text(size = ggplot2::rel(0.8)))
    if (!is.numeric(pd$x1))
      p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,
                                                                  hjust = 1,
                                                                  vjust = 1))
    size <- c(7, 5.5)

  # ---- One panel per covariate ---------------------------------------------
  } else {
    z     <- stats::qnorm(1 - (1 - conf_level) / 2)
    w     <- x$fit$W.orig
    beyond <- .hte_beyond(x)
    fmt_p <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "< 0.001"
                         else sprintf("= %.3f", p)
    # red is the dr layer's colour, so the loess line takes it only without one
    cate_line <- if ("dr" %in% display) "grey30" else "firebrick"

    panel <- function(v) {
      is_n  <- num[[v]]
      lv    <- if (!is_n) levels(droplevels(as.factor(d[[v]])))
      xval  <- function(val) if (is_n) val else factor(as.character(val),
                                                       levels = lv)
      dr    <- if ("dr" %in% display) .hte_dr_var(d, v, w, z, dr_args$spline_df,
                                                  beyond)
      label <- if (is.null(dr)) v else sprintf("%s (p_het %s)", v, fmt_p(dr$p_het))
      yv    <- c(0, ate)

      q <- ggplot2::ggplot() +
        ggplot2::geom_hline(yintercept = 0, colour = "grey75")
      if (length(ate))
        q <- q + ggplot2::geom_hline(yintercept = ate, linetype = 2)

      if ("cate" %in% display) {
        pts <- data.frame(x = xval(d[[v]]), y = d$.cate, panel = label)
        pts <- pts[!is.na(pts$x), , drop = FALSE]
        q <- q + if (is_n) {
          list(ggplot2::geom_point(data = pts, ggplot2::aes(x = x, y = y),
                                   colour = "grey55", alpha = 0.4, size = 0.8),
               if (cate_smooth > 0)
                 ggplot2::geom_smooth(data = pts, ggplot2::aes(x = x, y = y),
                                      method = "loess", formula = y ~ x,
                                      span = cate_smooth, se = FALSE,
                                      colour = cate_line, linewidth = 0.6))
        } else {
          ggplot2::geom_point(data = pts, ggplot2::aes(x = x, y = y),
                              position = ggplot2::position_jitter(
                                width = 0.15, height = 0, seed = 1),
                              colour = "grey55", alpha = 0.4, size = 0.8)
        }
        yv <- c(yv, pts$y)
      }

      # a covariate with no level or value left to estimate gets no dr layer
      if (!is.null(dr) && (if (is_n) nrow(dr$curve) > 0L
                           else any(!is.na(dr$levels$estimate)))) {
        if (is_n) {
          cv <- data.frame(dr$curve, panel = label)
          q <- q +
            ggplot2::geom_ribbon(data = cv, ggplot2::aes(x = x, ymin = conf.low,
                                                         ymax = conf.high),
                                 fill = "firebrick", alpha = 0.15) +
            ggplot2::geom_line(data = cv, ggplot2::aes(x = x, y = estimate),
                               colour = "firebrick", linewidth = 0.9)
          yv <- c(yv, cv$estimate)
        } else {
          ok <- !is.na(dr$levels$estimate)
          ld <- data.frame(x = xval(dr$levels$level[ok]),
                           dr$levels[ok, c("estimate", "conf.low", "conf.high")],
                           panel = label)
          q <- q + ggplot2::geom_pointrange(
            data = ld, ggplot2::aes(x = x, y = estimate, ymin = conf.low,
                                    ymax = conf.high),
            colour = "firebrick")
          yv <- c(yv, ld$estimate)
        }
      }

      if ("pdp" %in% display) {
        pd <- .hte_pdp(x, v, pdp_args$grid_n, pdp_args$max_n)
        pd <- data.frame(x = xval(pd[[v]]), estimate = pd$estimate, panel = label)
        q <- q + if (is_n) {
          ggplot2::geom_line(data = pd, ggplot2::aes(x = x, y = estimate),
                             colour = "steelblue", linewidth = 0.9,
                             linetype = "longdash")
        } else {
          ggplot2::geom_point(data = pd, ggplot2::aes(x = x, y = estimate),
                              colour = "steelblue", shape = 18, size = 3)
        }
        yv <- c(yv, pd$estimate)
      }

      q <- q + ggplot2::facet_wrap(~panel) +
        ggplot2::labs(x = v, y = ylab) +
        UtilsR::theme_my(base_rect_size = 1.5)
      if (!is_n)
        q <- q + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,
                                                                    hjust = 1,
                                                                    vjust = 1))
      list(plot = q, yv = yv)
    }

    built <- lapply(vars, panel)
    rng <- if (!is.null(ylim)) {
      rep(list(ylim), length(built))
    } else if (axis_arg$share_y == "all") {
      rep(list(range(unlist(lapply(built, `[[`, "yv")), na.rm = TRUE)),
          length(built))
    } else {
      lapply(built, function(b) range(b$yv, na.rm = TRUE))
    }
    plots <- Map(function(b, r) b$plot + ggplot2::coord_cartesian(ylim = r),
                 built, rng)

    caption <- paste(c(
      if ("cate" %in% display) "grey: out-of-bag CATE per patient",
      if ("cate" %in% display && !"dr" %in% display && cate_smooth > 0 &&
          any(num[vars])) "red: its loess",
      if ("dr" %in% display)
        sprintf("red: AIPW mean or spline with %g%% CI", 100 * conf_level),
      if ("pdp" %in% display) "blue: partial dependence of the forest",
      if (length(ate)) "dashed: ATE"), collapse = "; ")
    caption <- paste0(toupper(substr(caption, 1L, 1L)), substring(caption, 2L))
    ncol <- min(3L, length(plots))
    size <- if (length(plots) == 1L) c(5, 4.2) else
      c(3.4 * ncol + 0.6, 3 * ceiling(length(plots) / ncol) + 0.8)
    # About 10 caption characters fit per inch; wrap so nothing is cut off.
    caption <- paste(strwrap(caption, width = floor(10 * size[1L])),
                     collapse = "\n")

    p <- if (length(plots) == 1L) {
      plots[[1L]] + ggplot2::labs(title = title, caption = caption)
    } else {
      patchwork::wrap_plots(plots, ncol = ncol) +
        patchwork::plot_layout(axis_titles = "collect") +
        patchwork::plot_annotation(title = title, caption = caption,
                                   theme = UtilsR::theme_my(base_rect_size = 1.5))
    }
  }

  attr(p, "plot_size") <- stats::setNames(size, c("width", "height"))
  if (!is.null(save) && length(save) > 0L) {
    if (!requireNamespace("RegR", quietly = TRUE))
      stop("Package 'RegR' is required for a non-empty `save`.", call. = FALSE)
    if (is.null(save$width))  save$width  <- size[1L]
    if (is.null(save$height)) save$height <- size[2L]
    do.call(RegR::save_plt, c(list(plot = p), save))
  }
  p
}


#' Subgroup forest plot for heterogeneous treatment effects
#'
#' Recomputes the doubly robust average treatment effect (ATE) within each
#' level of the chosen categorical variables from a [get_hte()] result and
#' draws them with \pkg{forestplot}, in the diamond style of
#' `RegR::plt_eff2()`, wrapped as a ggplot by [ggplotify::as.ggplot()] so
#' \pkg{patchwork} can combine it with other plots. The estimates come from the
#' stored forest exactly as `get_hte(sub_var = )` computes them, so nothing is
#' refitted and the subgroups can be chosen after the fact.
#'
#' @param x An `hte_res` object from [get_hte()].
#' @param sub_var Character vector of categorical columns of `x$data`, drawn
#'   in the order given; every level is one subgroup. `NULL` (default) takes
#'   every categorical covariate of the forest, in the order of `adj_var`. A
#'   numeric column with more than 5 distinct values is continuous and has to
#'   be cut into groups first. A column the forest does not condition on is
#'   accepted with a message: its subgroup estimates hold when it is a
#'   function of the covariates (as `age_55` is of `Age`) or no confounder;
#'   otherwise add it with `get_hte(sub_var = )`.
#' @param measure One of `"diff"` (default), `"ratio"` or `"OR"`, as in
#'   [get_hte()]. `"ratio"` and `"OR"` are drawn on a log axis; `"OR"` needs
#'   an outcome probability, a binary outcome or \eqn{S(t)}.
#' @param conf_level Confidence level of the intervals. Default `0.95`,
#'   whatever level `x` was computed with.
#' @param overall Logical. `TRUE` (default) adds an "All patients" row with the
#'   overall ATE.
#' @param show_n Logical. `TRUE` (default) shows the column of patients per
#'   arm, treated first; `FALSE` leaves it out.
#' @param show_pvalue Logical. `TRUE` adds a `P` column with the p-value of
#'   every row. Default `FALSE`.
#' @param show_pinter Logical. `TRUE` adds a "P for interaction" column: for
#'   every variable, the Wald test that its subgroup effects are equal (the
#'   `p_inter` of [get_hte()]). Default `FALSE`.
#' @param xlim `NULL` (default) or two increasing numbers giving the axis
#'   range, positive for `"ratio"` and `"OR"`; intervals running past it end
#'   in arrows.
#' @param ticks_at `NULL` (default) or the axis ticks, on the scale of the
#'   estimates.
#' @param title Plot title, or `NULL` (default).
#' @param fixed_size Size policy. `TRUE` (default) pins the plot to the size
#'   its text needs -- rows 0.3 in tall, the graph column a quarter of the
#'   width and at least 2 in -- so it draws the same on every device: a larger
#'   one only adds white space, a smaller one clips it. A single positive
#'   number pins the total width in inches instead, the graph column taking
#'   what the text leaves. `FALSE` stretches the plot over the device, or over
#'   its cell of a \pkg{patchwork} composition: that is what
#'   `plot_layout(widths = , heights = )` needs, since patchwork cannot give a
#'   pinned plot a proportional share.
#' @param save `NULL` or a list with `filename`, `width` and `height`, passed
#'   to `RegR::save_plt()` for PDF output. `list()` and `NULL` skip saving; a
#'   `width` or `height` left out is taken from `attr(p, "plot_size")`. Do not
#'   include the `plot` argument; it is supplied internally.
#'
#' @section Subgroup estimates:
#' Every row is the doubly robust estimate within the subgroup, as in
#' `get_hte()$subgroup`: [grf::average_treatment_effect()] with `subset` for
#' `"diff"`, and the arm scores averaged over the subgroup for `"ratio"` and
#' `"OR"`. The counts are patients per arm, treated first. A level with fewer
#' than two patients in either arm -- or, for a survival outcome, no patient
#' in an arm followed beyond `time` -- has no estimate: it is drawn with a
#' dash, with a warning.
#'
#' @return A `ggplot`: the \pkg{forestplot} drawing, captured on the PDF
#'   metrics `RegR::save_plt()` writes with and wrapped by
#'   [ggplotify::as.ggplot()]. \pkg{patchwork} combines it with other plots
#'   as one plot (`p | q`, `p / q`, [patchwork::wrap_plots()]), but it is a
#'   picture of the table: themes and scales added to it do not restyle it.
#'   `attr(p, "subgroup")` holds the subgroup estimates drawn, laid out as
#'   `get_hte()$subgroup`, and `attr(p, "plot_size")` the `c(width, height)`
#'   in inches the plot is pinned to, or a suggested size with
#'   `fixed_size = FALSE`. If `save` is non-empty, the plot is also written to
#'   PDF through `RegR::save_plt()`.
#'
#' @seealso [get_hte()]; [plt_hte_dep()] for how the CATE varies with a
#'   covariate.
#'
#' @examplesIf requireNamespace("grf", quietly = TRUE) && requireNamespace("forestplot", quietly = TRUE) && requireNamespace("ggplotify", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)
#' \donttest{
#' set.seed(20260923)
#' n <- 600
#' d <- data.frame(age   = round(runif(n, 20, 85)),
#'                 sex   = factor(sample(c("F", "M"), n, replace = TRUE)),
#'                 stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
#' d$z <- rbinom(n, 1, 0.5)
#' d$y <- rbinom(n, 1, plogis(-1 + 0.02 * (d$age - 50) +
#'                              d$z * (0.2 + 0.6 * (d$sex == "M"))))
#' res <- get_hte(d, cat_var = "z", adj_var = c("age", "sex", "stage"),
#'                surv = "y", grf_args = list(num.trees = 500, seed = 1))
#'
#' # Every categorical covariate, the overall ATE on top
#' plt_hte_sub(res)
#'
#' # Chosen variables as risk ratios, with both P columns
#' plt_hte_sub(res, sub_var = c("sex", "stage"), measure = "ratio",
#'             show_pvalue = TRUE, show_pinter = TRUE)
#'
#' # A ggplot, so patchwork stacks it over the CATE density by sex
#' patchwork::wrap_plots(plt_hte_sub(res, sub_var = "sex"),
#'                       plt_hte_cate(res, sub_var = "sex", type = "density"),
#'                       ncol = 1)
#' }
#'
#' @export
plt_hte_sub <- function(x,
                        sub_var     = NULL,
                        measure     = c("diff", "ratio", "OR"),
                        conf_level  = 0.95,
                        overall     = TRUE,
                        show_n      = TRUE,
                        show_pvalue = FALSE,
                        show_pinter = FALSE,
                        xlim        = NULL,
                        ticks_at    = NULL,
                        title       = NULL,
                        fixed_size  = TRUE,
                        save        = list()) {

  if (!inherits(x, "hte_res"))
    stop("`x` must be an `hte_res` object from get_hte().", call. = FALSE)
  measure <- match.arg(measure)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  flags <- list(overall = overall, show_n = show_n, show_pvalue = show_pvalue,
                show_pinter = show_pinter)
  for (nm in names(flags))
    if (!is.logical(flags[[nm]]) || length(flags[[nm]]) != 1L ||
        is.na(flags[[nm]]))
      stop(sprintf("`%s` must be TRUE or FALSE.", nm), call. = FALSE)
  log_x <- measure != "diff"
  if (!is.null(xlim) && (!is.numeric(xlim) || length(xlim) != 2L ||
                         anyNA(xlim) || xlim[1L] >= xlim[2L]))
    stop("`xlim` must be `NULL` or two increasing numbers.", call. = FALSE)
  if (!is.null(ticks_at) && (!is.numeric(ticks_at) || !length(ticks_at) ||
                             anyNA(ticks_at)))
    stop("`ticks_at` must be `NULL` or a numeric vector.", call. = FALSE)
  if (log_x && any(c(xlim, ticks_at) <= 0))
    stop("`xlim` and `ticks_at` must be positive on the log axis of \"ratio\" and \"OR\".",
         call. = FALSE)
  if (!(isTRUE(fixed_size) || isFALSE(fixed_size) ||
        (is.numeric(fixed_size) && length(fixed_size) == 1L &&
         is.finite(fixed_size) && fixed_size > 0)))
    stop("`fixed_size` must be TRUE, FALSE, or a single positive width in inches.",
         call. = FALSE)
  if (!is.null(save) && !is.list(save))
    stop("`save` must be `NULL` or a list.", call. = FALSE)
  # grf also has to be loaded for predict() on a forest read back from disk
  for (pkg in c("grf", "forestplot", "ggplotify"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required for plt_hte_sub().", pkg),
           call. = FALSE)

  a <- attr(x, "analysis")
  d <- x$data
  scale <- if (identical(a$outcome_type, "survival")) {
    if (identical(a$target, "RMST")) "RMST" else "S"
  } else {
    a$outcome_type
  }
  if (measure == "OR" && scale %in% c("RMST", "continuous"))
    stop("`measure = \"OR\"` needs an outcome probability: a binary outcome or S(t).",
         call. = FALSE)

  # ---- Subgroup variables ---------------------------------------------------
  sub_var <- if (is.null(sub_var)) {
    a$covariates[!vapply(d[a$covariates], .hte_is_num, logical(1L))]
  } else {
    .hte_sub_var(sub_var, d, a)
  }
  if (!length(sub_var) && !overall)
    stop("Nothing to draw: `sub_var` holds no categorical variable and `overall = FALSE`.",
         call. = FALSE)

  # ---- Estimates from the stored forest -------------------------------------
  fit  <- x$fit
  s    <- .hte_arm_scores(fit)
  z    <- stats::qnorm(1 - (1 - conf_level) / 2)
  grid <- data.frame(estimand = "ATE", measure = measure,
                     stringsAsFactors = FALSE)
  event_risk <- identical(a$target, "survival.probability")
  beyond     <- .hte_beyond(x)
  sub <- if (length(sub_var))
    .hte_muffle_ps(.hte_subgroup(fit, s, d, sub_var, grid, event_risk, z,
                                 beyond))
  ov  <- if (overall)
    .hte_muffle_ps(.hte_estimate(fit, s, rep(TRUE, nrow(d)), grid, event_risk,
                                 z, "Overall", beyond))

  # ---- Table ----------------------------------------------------------------
  ref <- setdiff(levels(factor(d[[a$cat_var]])), a$treated)[1L]
  t_x <- format(a$time)
  lab <- switch(measure,
    diff  = switch(scale, S = sprintf("S(%s) difference", t_x),
                   RMST = sprintf("RMST(%s) difference", t_x),
                   binary = "Risk difference", continuous = "Mean difference"),
    ratio = switch(scale, S = "Event risk ratio", RMST = "RMST ratio",
                   binary = "Risk ratio", continuous = "Ratio of means"),
    OR    = switch(scale, S = "Event odds ratio", binary = "Odds ratio"))
  num_f   <- paste0("%.", if (!log_x && scale %in% c("S", "binary")) 3 else 2, "f")
  fmt_est <- function(e, l, h)
    ifelse(is.na(e), "\u2014",
           sprintf(paste0(num_f, " (", num_f, ", ", num_f, ")"), e, l, h))
  fmt_p <- function(p)
    ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
  fmt_n <- function(n, nt) sprintf("%d / %d", nt, n - nt)
  row_df <- function(label, n, est, p, p_inter, mean, lower, upper, bold)
    data.frame(label = label, n = n, est = est, p = p, p_inter = p_inter,
               mean = mean, lower = lower, upper = upper, bold = bold,
               stringsAsFactors = FALSE)

  body <- if (overall)
    row_df("All patients", fmt_n(nrow(d), sum(fit$W.orig)),
           fmt_est(ov$estimate, ov$conf.low, ov$conf.high), fmt_p(ov$p.value),
           "", ov$estimate, ov$conf.low, ov$conf.high, FALSE)
  for (v in sub_var) {
    r <- sub[sub$sub_var == v, ]
    body <- rbind(body,
                  row_df(v, "", "", "", fmt_p(r$p_inter[1L]), NA, NA, NA, TRUE),
                  row_df(paste0("   ", r$level), fmt_n(r$n, r$n_treat),
                         fmt_est(r$estimate, r$conf.low, r$conf.high),
                         fmt_p(r$p.value), "", r$estimate, r$conf.low,
                         r$conf.high, FALSE))
  }
  cols <- c("label", if (show_n) "n", "est", if (show_pvalue) "p",
            if (show_pinter) "p_inter")
  text <- rbind(c("Subgroup",
                  if (show_n) sprintf("N (%s / %s)", a$treated, ref),
                  sprintf("%s (%s%% CI)", lab, format(100 * conf_level)),
                  if (show_pvalue) "P", if (show_pinter) "P for interaction"),
                as.matrix(body[cols]))
  dimnames(text) <- NULL
  mean  <- c(NA, body$mean)
  lower <- c(NA, body$lower)
  upper <- c(NA, body$upper)
  bold  <- c(TRUE, body$bold)

  # ---- Axis -----------------------------------------------------------------
  zero <- if (log_x) 1 else 0
  if (is.null(ticks_at)) {
    rng <- if (is.null(xlim)) range(c(lower, upper, zero), na.rm = TRUE) else xlim
    ticks_at <- if (log_x) {
      cand <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100)
      tk <- cand[cand >= max(c(cand[cand <= rng[1L]], cand[1L])) &
                   cand <= min(c(cand[cand >= rng[2L]], cand[length(cand)]))]
      # 1-2-5 steps crowd a wide range: thin to 1-5 steps, then to powers of
      # 10, always keeping the two outer ticks
      for (drop in list(c(0.02, 0.2, 2, 20), c(0.05, 0.5, 5, 50)))
        if (length(tk) > 7L)
          tk <- tk[!tk %in% drop | seq_along(tk) %in% c(1L, length(tk))]
      tk
    } else {
      pretty(rng, n = 4)
    }
    if (!is.null(xlim)) {
      ticks_at <- ticks_at[ticks_at >= xlim[1L] & ticks_at <= xlim[2L]]
      if (length(ticks_at) < 2L) ticks_at <- xlim
    }
  }
  # forestplot wants every log-axis value positive, so the open clip starts
  # at the smallest limit or tick rather than at 0
  clip <- if (!is.null(xlim)) {
    xlim
  } else if (log_x) {
    c(min(c(lower, ticks_at), na.rm = TRUE), Inf)
  } else {
    c(-Inf, Inf)
  }
  attr(ticks_at, "labels") <- as.character(ticks_at)   # 1, not 1.00
  xlab <- if (log_x) {
    sprintf("%s%s, %s vs %s", lab,
            if (scale == "S") sprintf(" by t = %s", t_x) else "", a$treated, ref)
  } else {
    sprintf("%s, %s - %s", lab, a$treated, ref)
  }

  # ---- forestplot, styled as RegR::plt_eff2() -------------------------------
  # forestplot() measures its text as it builds, as does the size estimate
  # below: both run on a null PDF device, so no stray window or Rplots.pdf
  # opens, and the device that was current stays current.
  on_null <- function(expr, ...) {
    prev <- grDevices::dev.cur()
    grDevices::pdf(NULL, ...)
    dev <- grDevices::dev.cur()
    on.exit({
      grDevices::dev.off(dev)
      if (prev > 1L) grDevices::dev.set(prev)
    })
    expr
  }
  n_text <- ncol(text)
  row_in <- 0.3   # row height in inches, pinned or assumed by the suggested size
  rules  <- list("1" = grid::gpar(lty = 1, lwd = 2), "2" = grid::gpar(lty = 2))
  rules[[as.character(nrow(text) + 1L)]] <-
    grid::gpar(lty = 1, lwd = 2, columns = seq_len(n_text))
  p <- on_null(forestplot::forestplot(
    labeltext  = text, mean = mean, lower = lower, upper = upper,
    is.summary = bold, zero = zero, xlog = log_x, xticks = ticks_at,
    clip = clip, graph.pos = "right", align = c("l", rep("r", n_text - 1L)),
    hrzl_lines = rules, xlab = xlab, title = title,
    fn.ci_norm = forestplot::fpDrawDiamondCI, boxsize = 0.3,
    col = forestplot::fpColors(box = "blue4", lines = "blue4",
                               zero = "black"),
    txt_gp = forestplot::fpTxtGp(label = grid::gpar(cex = 0.8),
                                 ticks = grid::gpar(cex = 0.8),
                                 xlab  = grid::gpar(cex = 0.9),
                                 title = grid::gpar(cex = 1.2)),
    lwd.zero = 1, lwd.ci = 1.5, lwd.xaxis = 2, ci.vertices = TRUE,
    ci.vertices.height = 0.2, colgap = grid::unit(6, "mm"),
    lineheight = if (isFALSE(fixed_size)) "auto" else grid::unit(row_in, "in")))

  if (isFALSE(fixed_size)) {
    # Suggested size: the widest cell of every text column (bold rows at
    # forestplot's summary size, 1.1 x the label cex), a graph column that
    # takes a quarter of the width as in RegR::plt_eff2(), and about 0.3 in a
    # row.
    text_in <- on_null(sum(vapply(seq_len(n_text), function(j)
      max(vapply(seq_len(nrow(text)), function(i) {
        gp <- if (bold[i]) grid::gpar(cex = 0.88, fontface = "bold")
              else grid::gpar(cex = 0.8)
        grid::convertWidth(grid::grobWidth(grid::textGrob(text[i, j], gp = gp)),
                           "in", valueOnly = TRUE)
      }, numeric(1L))), numeric(1L)))) + n_text * 6 / 25.4
    size <- round(c(width  = max(text_in / 0.75, text_in + 2) + 0.4,
                    height = row_in * nrow(text) + 0.9 +
                      if (is.null(title)) 0 else 0.4), 1)
  } else {
    # Pin the drawing. forestplot sizes the text columns from the rendered
    # text but leaves the graph column relative, and surrounds the table with
    # its margins, the title and the x-axis strip: one probe draw, on the PDF
    # metrics RegR::save_plt() writes with and a page wide enough for any
    # table, reads the text columns and that overhead back in inches.
    probe <- on_null({
      print(p)
      grid::seekViewport("BaseGrid")
      w <- grid::convertWidth(grid::current.viewport()$layout$widths, "in",
                              valueOnly = TRUE)
      area <- c(grid::convertWidth(grid::unit(1, "npc"), "in", valueOnly = TRUE),
                grid::convertHeight(grid::unit(1, "npc"), "in", valueOnly = TRUE))
      grid::upViewport(0)
      # the graph is the last layout column
      list(text = sum(w[-length(w)]), over = grDevices::dev.size("in") - area)
    }, width = 40, height = 40)
    graph_w <- if (is.numeric(fixed_size)) {
      fixed_size - probe$text - probe$over[1L]
    } else {
      max(probe$text / 3, 2)
    }
    if (graph_w <= 0)
      stop(sprintf("`fixed_size` leaves no room for the forest column: the labels alone need %.1f in.",
                   probe$text + probe$over[1L]), call. = FALSE)
    p$graphwidth <- grid::unit(graph_w, "in")
    # Rounded up, as a canvas short of the pinned drawing clips it; a width
    # asked for is kept as it is.
    up   <- function(v) ceiling(10 * v - 1e-8) / 10
    size <- c(width  = if (is.numeric(fixed_size)) fixed_size
                       else up(probe$text + graph_w + probe$over[1L]),
              height = up(row_in * nrow(text) + probe$over[2L]))
  }

  # ---- ggplot ---------------------------------------------------------------
  # Captured the way ggplotify captures grid plots, on a null PDF device of
  # `size`, so the text keeps the metrics measured above, then wrapped by
  # as.ggplot() for patchwork. The capture fills the panel it is drawn in: a
  # pinned plot holds that panel to `size`; otherwise it stretches as
  # forestplot does, save for the 5 mm margins, which forestplot turns into a
  # share of the capture. grid.grabExpr() hands back to the device that was
  # current, and selecting the null device opens a new one, hence on_null().
  g <- on_null(grid::grid.grabExpr(print(p), warn = 0, width = size[["width"]],
                                   height = size[["height"]]))
  # ggplotify 0.1.2 builds its ggplot with aes_(), deprecated in ggplot2
  p <- withCallingHandlers(ggplotify::as.ggplot(g),
                           lifecycle_warning_deprecated = function(w)
                             invokeRestart("muffleWarning"))
  if (!isFALSE(fixed_size))
    p <- p + ggplot2::theme(panel.widths  = grid::unit(size[["width"]], "in"),
                            panel.heights = grid::unit(size[["height"]], "in"))

  attr(p, "subgroup")  <- sub
  attr(p, "plot_size") <- size
  if (!is.null(save) && length(save) > 0L) {
    if (!requireNamespace("RegR", quietly = TRUE))
      stop("Package 'RegR' is required for a non-empty `save`.", call. = FALSE)
    if (is.null(save$width))  save$width  <- size[[1L]]
    if (is.null(save$height)) save$height <- size[[2L]]
    do.call(RegR::save_plt, c(list(plot = p), save))
  }
  p
}


#' Waterfall and density plots of the CATE by subgroup
#'
#' Draws the conditional average treatment effect (CATE) that [get_hte()]
#' estimated for every patient, after `StratifiedMedicine::plot_ple()`: a
#' waterfall of the patients sorted by their CATE, or the density of the CATE,
#' coloured by the levels of a subgroup variable. Dashed lines mark the
#' overall average treatment effect (ATE) and the subgroup ATEs, the doubly
#' robust estimates [plt_hte_sub()] draws.
#'
#' @param x An `hte_res` object from [get_hte()].
#' @param sub_var `NULL` (default) for one panel of every patient, or a
#'   character vector of categorical columns of `x$data`: one panel each,
#'   three per row, coloured by level. The columns follow the rules of
#'   [plt_hte_sub()]; a patient missing the value is left out of that panel.
#' @param type `"waterfall"` (default) draws one bar per patient, sorted by
#'   the CATE; `"density"` draws the density of the CATE per level.
#' @param show_ci Logical. `TRUE` adds grf's pointwise confidence interval of
#'   every patient's CATE (grey), which needs a forest grown with
#'   `ci.group.size` of at least 2, grf's default. Default `FALSE`. Only used
#'   by `type = "waterfall"`.
#' @param conf_level Confidence level of `show_ci`. Default `0.95`. Only used
#'   by `type = "waterfall"`.
#' @param overall Logical. `TRUE` (default) marks the overall ATE with a black
#'   dashed line.
#' @param title Plot title, or `NULL` (default).
#' @param save `NULL` or a list with `filename`, `width` and `height`, passed
#'   to `RegR::save_plt()` for PDF output. `list()` and `NULL` skip saving; a
#'   `width` or `height` left out is taken from `attr(p, "plot_size")`. Do not
#'   include the `plot` argument; it is supplied internally.
#'
#' @section Reading the plots:
#' The bars and curves are the forest's out-of-bag estimates, one per patient,
#' on the `"diff"` scale of [get_hte()]. They describe how the forest spreads
#' the effect and test nothing: forest estimates are shrunk towards the
#' overall mean, and their spread mixes real heterogeneity with estimation
#' noise. The coloured lines are the doubly robust subgroup ATEs, whose
#' intervals and `p_inter` [plt_hte_sub()] shows. A level without an estimate
#' -- fewer than two patients in an arm, or for a survival outcome no patient
#' in an arm followed beyond `time` -- keeps its bars or curve but gets no
#' line, with a warning.
#'
#' @return A `ggplot` for one panel, otherwise a patchwork. The CATE axis is
#'   the S(t) or RMST difference for a survival outcome, the risk difference
#'   for a binary one and the mean difference otherwise.
#'   `attr(p, "subgroup")` holds the subgroup ATEs drawn, laid out as
#'   `get_hte()$subgroup` (`NULL` without `sub_var`), and
#'   `attr(p, "plot_size")` a suggested `c(width, height)` in inches. If
#'   `save` is non-empty, the plot is also written to PDF through
#'   `RegR::save_plt()`.
#'
#' @seealso [plt_hte_sub()] for the subgroup ATEs with their intervals;
#'   [plt_hte_dep()] for how the CATE varies with a covariate.
#'
#' @examplesIf requireNamespace("grf", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)
#' \donttest{
#' set.seed(20260923)
#' n <- 600
#' d <- data.frame(age   = round(runif(n, 20, 85)),
#'                 sex   = factor(sample(c("F", "M"), n, replace = TRUE)),
#'                 stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
#' d$z <- rbinom(n, 1, 0.5)
#' d$y <- rbinom(n, 1, plogis(-1 + 0.02 * (d$age - 50) +
#'                              d$z * (0.2 + 0.6 * (d$sex == "M"))))
#' res <- get_hte(d, cat_var = "z", adj_var = c("age", "sex", "stage"),
#'                surv = "y", grf_args = list(num.trees = 500, seed = 1))
#'
#' # Patients sorted by their CATE, coloured by sex
#' plt_hte_cate(res, sub_var = "sex")
#'
#' # Density of the CATE by sex and by stage
#' plt_hte_cate(res, sub_var = c("sex", "stage"), type = "density")
#' }
#'
#' @export
plt_hte_cate <- function(x,
                         sub_var    = NULL,
                         type       = c("waterfall", "density"),
                         show_ci    = FALSE,
                         conf_level = 0.95,
                         overall    = TRUE,
                         title      = NULL,
                         save       = list()) {

  if (!inherits(x, "hte_res"))
    stop("`x` must be an `hte_res` object from get_hte().", call. = FALSE)
  type <- match.arg(type)
  if (type == "density") {
    used <- c(show_ci = !missing(show_ci), conf_level = !missing(conf_level))
    if (any(used))
      stop(sprintf("%s only applies to type = \"waterfall\".",
                   paste0("`", names(used)[used], "`", collapse = ", ")),
           call. = FALSE)
  }
  flags <- list(show_ci = show_ci, overall = overall)
  for (nm in names(flags))
    if (!is.logical(flags[[nm]]) || length(flags[[nm]]) != 1L ||
        is.na(flags[[nm]]))
      stop(sprintf("`%s` must be TRUE or FALSE.", nm), call. = FALSE)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  if (!is.null(save) && !is.list(save))
    stop("`save` must be `NULL` or a list.", call. = FALSE)
  # grf also has to be loaded for predict() on a forest read back from disk
  if (!requireNamespace("grf", quietly = TRUE))
    stop("Package 'grf' is required for plt_hte_cate().", call. = FALSE)

  a <- attr(x, "analysis")
  d <- x$data
  if (!is.null(sub_var)) sub_var <- .hte_sub_var(sub_var, d, a)
  if (length(sub_var) > 1L && !requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required to combine several panels.",
         call. = FALSE)

  # ---- Estimates: the lines as plt_hte_sub() computes them -------------------
  fit    <- x$fit
  s      <- .hte_arm_scores(fit)
  z      <- stats::qnorm(1 - (1 - conf_level) / 2)
  grid   <- data.frame(estimand = "ATE", measure = "diff",
                       stringsAsFactors = FALSE)
  beyond <- .hte_beyond(x)
  sub <- if (length(sub_var))
    .hte_muffle_ps(.hte_subgroup(fit, s, d, sub_var, grid, FALSE, z, beyond))
  ate <- if (overall)
    .hte_muffle_ps(.hte_estimate(fit, s, rep(TRUE, nrow(d)), grid, FALSE, z,
                                 "Overall", beyond))$estimate
  se  <- if (show_ci)
    sqrt(as.numeric(stats::predict(fit, estimate.variance = TRUE)$variance.estimates))

  # ---- One panel per variable ------------------------------------------------
  lab  <- .hte_cate_label(x)
  cate <- d$.cate
  # ggplot2's default hues, so the lines take the colour of their level
  hue  <- function(k) grDevices::hcl(h = seq(15, 375, length.out = k + 1L)[seq_len(k)],
                                     c = 100, l = 65)

  panel <- function(v) {
    g <- if (!is.null(v)) droplevels(as.factor(d[[v]]))
    o <- order(cate)
    if (!is.null(g)) o <- o[!is.na(g[o])]
    pd <- data.frame(rank = seq_along(o), cate = cate[o])
    if (show_ci) {
      pd$conf.low  <- pd$cate - z * se[o]
      pd$conf.high <- pd$cate + z * se[o]
    }
    lines <- NULL
    if (!is.null(g)) {
      pd$level <- g[o]
      pd$panel <- v
      r <- sub[sub$sub_var == v & !is.na(sub$estimate), ]
      lines <- data.frame(level = factor(r$level, levels = levels(g)),
                          estimate = r$estimate, panel = v)
    }
    ref <- if (overall) data.frame(ate = ate)

    q <- ggplot2::ggplot()
    if (type == "waterfall") {
      q <- q + ggplot2::geom_hline(yintercept = 0, colour = "grey75")
      if (show_ci)
        q <- q + ggplot2::geom_linerange(
          data = pd, ggplot2::aes(x = rank, ymin = conf.low, ymax = conf.high),
          colour = "grey70", linewidth = 0.3)
      # Outlined in their own colour: side by side, bare bars leave hairline
      # gaps once rasterised.
      q <- q + if (is.null(g)) {
        list(ggplot2::geom_col(data = pd, ggplot2::aes(x = rank, y = cate,
                                                       fill = cate,
                                                       colour = cate),
                               width = 1, linewidth = 0.3),
             ggplot2::scale_fill_gradient2(low = "navy", mid = "white",
                                           high = "firebrick", midpoint = 0,
                                           guide = "none"),
             ggplot2::scale_colour_gradient2(low = "navy", mid = "white",
                                             high = "firebrick", midpoint = 0,
                                             guide = "none"))
      } else {
        ggplot2::geom_col(data = pd, ggplot2::aes(x = rank, y = cate,
                                                  fill = level, colour = level),
                          width = 1, linewidth = 0.3)
      }
      if (overall)
        q <- q + ggplot2::geom_hline(data = ref, ggplot2::aes(yintercept = ate),
                                     linetype = 2)
      if (length(lines) && nrow(lines))
        q <- q + ggplot2::geom_hline(data = lines,
                                     ggplot2::aes(yintercept = estimate,
                                                  colour = level),
                                     linetype = 2, linewidth = 0.8)
      q <- q + ggplot2::labs(x = "Patients ranked by CATE", y = lab) +
        UtilsR::theme_my(base_rect_size = 1.5) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       panel.grid.major.x = ggplot2::element_blank(),
                       panel.grid.minor.x = ggplot2::element_blank())
    } else {
      q <- q + ggplot2::geom_vline(xintercept = 0, colour = "grey75") +
        if (is.null(g)) {
          ggplot2::geom_density(data = pd, ggplot2::aes(x = cate),
                                fill = "grey70", colour = "grey30", alpha = 0.5)
        } else {
          ggplot2::geom_density(data = pd, ggplot2::aes(x = cate, fill = level,
                                                        colour = level),
                                alpha = 0.3)
        }
      if (overall)
        q <- q + ggplot2::geom_vline(data = ref, ggplot2::aes(xintercept = ate),
                                     linetype = 2)
      if (length(lines) && nrow(lines))
        q <- q + ggplot2::geom_vline(data = lines,
                                     ggplot2::aes(xintercept = estimate,
                                                  colour = level),
                                     linetype = 2, linewidth = 0.8)
      q <- q + ggplot2::labs(x = lab, y = "Density") +
        UtilsR::theme_my(base_rect_size = 1.5)
    }
    # 1% of the axis beside the bars or curves instead of ggplot2's 5%
    q <- q + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.01))
    if (!is.null(g)) {
      cols <- stats::setNames(hue(nlevels(g)), levels(g))
      q <- q + ggplot2::scale_fill_manual(values = cols, name = NULL) +
        ggplot2::scale_colour_manual(values = cols, guide = "none") +
        ggplot2::facet_wrap(~panel) +
        ggplot2::theme(legend.position = "bottom")
    }
    q
  }
  plots <- lapply(if (length(sub_var)) sub_var else list(NULL), panel)

  caption <- paste(c(
    if (type == "waterfall") "bars: out-of-bag CATE per patient, sorted"
    else "curves: density of the out-of-bag CATE",
    if (show_ci) sprintf("grey: %g%% CI", 100 * conf_level),
    if (overall) "black dashed: ATE",
    if (length(sub_var)) "coloured dashed: subgroup ATE (doubly robust)"),
    collapse = "; ")
  caption <- paste0(toupper(substr(caption, 1L, 1L)), substring(caption, 2L))
  ncol <- min(3L, length(plots))
  size <- if (length(plots) == 1L) c(if (type == "waterfall") 7 else 6, 4.5) else
    c(3.6 * ncol + 0.6, 3.6 * ceiling(length(plots) / ncol) + 0.8)
  # About 10 caption characters fit per inch; wrap so nothing is cut off.
  caption <- paste(strwrap(caption, width = floor(10 * size[1L])),
                   collapse = "\n")

  p <- if (length(plots) == 1L) {
    plots[[1L]] + ggplot2::labs(title = title, caption = caption)
  } else {
    patchwork::wrap_plots(plots, ncol = ncol) +
      patchwork::plot_layout(axis_titles = "collect") +
      patchwork::plot_annotation(title = title, caption = caption,
                                 theme = UtilsR::theme_my(base_rect_size = 1.5))
  }

  attr(p, "subgroup")  <- sub
  attr(p, "plot_size") <- stats::setNames(size, c("width", "height"))
  if (!is.null(save) && length(save) > 0L) {
    if (!requireNamespace("RegR", quietly = TRUE))
      stop("Package 'RegR' is required for a non-empty `save`.", call. = FALSE)
    if (is.null(save$width))  save$width  <- size[1L]
    if (is.null(save$height)) save$height <- size[2L]
    do.call(RegR::save_plt, c(list(plot = p), save))
  }
  p
}
