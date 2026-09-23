# =============================================================================
# hte-plt.R -- dependence plots for get_hte()
# =============================================================================
#
# Architecture:
#
#   L1  plt_hte_dep()  one panel per covariate ("dep") or a two-covariate
#                      partial-dependence heat map ("heat")
#   L2  .hte_pdp()     forest CATE averaged with covariates set to grid values
#
# The doubly robust layer comes from .hte_dr_var() in hte-get.R, so the p_het
# in the strips is the one in get_hte()$importance at the default spline df.
# =============================================================================


# ---- L2 partial dependence -------------------------------------------------

# Partial dependence as in StratifiedMedicine::plot_dependence(): every
# analysed row -- or an evenly spaced subset of at most `max_n`, so the result
# does not depend on the random seed -- gets the covariates in `vars` set to
# each grid combination, and the forest's CATE is averaged. A factor is set
# through its one-hot columns, which get_hte() keeps for every level, so no
# row ever carries two levels at once; a numeric covariate is one column.
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
#'       for a categorical covariate and smoothed by loess for a continuous
#'       one. Descriptive only: forest estimates are shrunk towards the
#'       overall mean.}
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
#' @param ylim `NULL` (default) or two increasing numbers giving the y range
#'   of every panel; intervals running past it are clipped.
#' @param dr_args Named list for the `"dr"` layer: `spline_df` (default `3`),
#'   the natural-spline degrees of freedom of a continuous covariate, which
#'   also sets the degrees of freedom of its `p_het`. Only used by
#'   `type = "dep"`.
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
                        ylim     = NULL,
                        dr_args  = list(spline_df = 3),
                        pdp_args = list(grid_n = 21, max_n = 1000),
                        axis_arg = list(share_y = "all"),
                        title    = NULL,
                        save     = list()) {

  if (!inherits(x, "hte_res"))
    stop("`x` must be an `hte_res` object from get_hte().", call. = FALSE)
  type <- match.arg(type)
  if (type == "heat") {
    used <- c(display = !missing(display), dr_args = !missing(dr_args),
              axis_arg = !missing(axis_arg))
    if (any(used))
      stop(sprintf("%s only applies to type = \"dep\"; the heat map shows the partial dependence.",
                   paste0("`", names(used)[used], "`", collapse = ", ")),
           call. = FALSE)
  }
  display  <- match.arg(display, c("cate", "dr", "pdp"), several.ok = TRUE)
  dr_args  <- .merge_named_arg(dr_args, list(spline_df = 3), "dr_args")
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

  ref  <- setdiff(levels(factor(d[[a$cat_var]])), a$treated)[1L]
  what <- switch(a$outcome_type,
                 survival   = sprintf("%s(%s)",
                                      if (identical(a$target, "RMST")) "RMST" else "S",
                                      format(a$time)),
                 binary     = "risk",
                 continuous = "mean")
  ylab <- sprintf("CATE: %s %s - %s", what, a$treated, ref)
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
      ggplot2::theme_bw()
    if (!is.numeric(pd$x1))
      p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,
                                                                  hjust = 1))
    size <- c(7, 5.5)

  # ---- One panel per covariate ---------------------------------------------
  } else {
    z     <- stats::qnorm(1 - (1 - a$conf_level) / 2)
    w     <- x$fit$W.orig
    # patients followed past `time`, flagged by the rule get_hte() applies
    # (grf keeps only the times cut at the horizon)
    beyond <- if (identical(a$outcome_type, "survival")) {
      y <- d[[a$outcome[1L]]]
      if (identical(a$target, "RMST")) y >= a$time else y > a$time
    }
    fmt_p <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "< 0.001"
                         else sprintf("= %.3f", p)

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
               ggplot2::geom_smooth(data = pts, ggplot2::aes(x = x, y = y),
                                    method = "loess", formula = y ~ x,
                                    span = 0.6, se = FALSE, colour = "grey30",
                                    linewidth = 0.6))
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
        ggplot2::labs(x = NULL, y = ylab) +
        ggplot2::theme_bw()
      if (!is_n)
        q <- q + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,
                                                                    hjust = 1))
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
      if ("dr" %in% display)
        sprintf("red: AIPW mean or spline with %g%% CI", 100 * a$conf_level),
      if ("pdp" %in% display) "blue: partial dependence of the forest",
      if (length(ate)) "dashed: ATE"), collapse = "; ")
    caption <- paste0(toupper(substr(caption, 1L, 1L)), substring(caption, 2L))
    ncol <- min(3L, length(plots))
    size <- if (length(plots) == 1L) c(5, 4.2) else
      c(3.4 * ncol + 0.6, 3 * ceiling(length(plots) / ncol) + 0.8)
    # About 13 caption characters fit per inch; wrap so nothing is cut off.
    caption <- paste(strwrap(caption, width = floor(13 * size[1L])),
                     collapse = "\n")

    p <- if (length(plots) == 1L) {
      plots[[1L]] + ggplot2::labs(title = title, caption = caption)
    } else {
      patchwork::wrap_plots(plots, ncol = ncol) +
        patchwork::plot_layout(axis_titles = "collect") +
        patchwork::plot_annotation(title = title, caption = caption)
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
