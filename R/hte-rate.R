# =============================================================================
# hte-rate.R -- how well a ranking of the patients targets the effect (RATE)
# =============================================================================
#
# Architecture:
#
#   L1  plt_hte_rate()  rank the patients, evaluate the ranking with grf's
#                       rank-weighted average treatment effect and draw its
#                       TOC and Qini curves
#
# A ranking learnt from the outcomes -- the forest CATE -- is learnt on one
# random split and evaluated on the other, as grf's RATE vignette requires. A
# pre-specified score is not learnt from them, so it is evaluated on every
# patient with the forest get_hte() stored. With unit costs the Qini curve is
# q x TOC(q), whose area is grf's QINI, so it needs no further estimation.
# =============================================================================


#' Targeting curves for heterogeneous treatment effects
#'
#' Evaluates how well a ranking of the patients picks out those who benefit
#' most, with the rank-weighted average treatment effect (RATE) of
#' [grf::rank_average_treatment_effect()], and draws its targeting operator
#' characteristic (TOC) and Qini curves from a [get_hte()] result.
#'
#' @param x An `hte_res` object from [get_hte()].
#' @param priority One or two ranking rules; patients with larger values are
#'   treated first.
#'   \describe{
#'     \item{`"cate"` (default)}{The forest CATE. It is learnt from the
#'       outcomes, and ranking the patients it was learnt from would flatter
#'       it, so the patients are split at random: a forest refitted to a share
#'       `train_frac` of them ranks the rest, and a forest refitted to the
#'       rest evaluates that ranking. Both refits take the `grf_args` of
#'       [get_hte()] with `seed`; for a survival outcome, each half needs
#'       patients of both arms followed past `time`.}
#'     \item{A numeric or logical column of `x$data`}{A pre-specified score,
#'       such as a biomarker or a published risk score, evaluated on every
#'       patient with the doubly robust scores of the stored forest. Patients
#'       missing it are left out. As for a subgroup in [plt_hte_sub()], a
#'       column the forest does not condition on is evaluated correctly only
#'       if it is no confounder. To treat small values first, add the negated
#'       column to `x$data`.}
#'   }
#'   Two rules are drawn together and their difference is tested; with
#'   `"cate"` among them, both are evaluated on the held-out patients.
#' @param type Panels to draw, `"toc"` and / or `"qini"` (default both).
#'   \describe{
#'     \item{`"toc"`}{TOC(q): the average effect among the share q of
#'       patients ranked first minus the average effect of all of them. Its
#'       area, the AUTOC, weights the top of the ranking most, so it detects
#'       an effect held by a few patients.}
#'     \item{`"qini"`}{q x TOC(q): what treating the share q ranked first
#'       gains over treating a random share q. Its area, the QINI, weights
#'       every share alike, so it detects an effect that changes gradually.}
#'   }
#' @param smooth `0` (default) draws the TOC and its band as estimated; a
#'   number from `0.05` to `1` is a LOESS span that smooths both for display,
#'   separately for each rule: about `0.1` smooths slightly, `0.3` strongly.
#'   The Qini panel is then q times the smoothed TOC. The AUTOC, QINI and
#'   `attr(p, "rate")` do not change.
#' @param conf_level Confidence level of the bands and intervals. Default
#'   `0.95`.
#' @param train_frac Share of the patients the ranking forest is refitted to,
#'   strictly between 0 and 1. Default `0.5`. Only used when `priority`
#'   includes `"cate"`.
#' @param seed Seed of the split, the refitted forests and grf's bootstrap
#'   standard errors, so a call can be repeated exactly; `NULL` (default)
#'   takes the seed of the forest in `x`. The random number stream of the
#'   session is restored afterwards.
#' @param title Plot title, or `NULL` (default).
#' @param save `NULL` or a list with `filename`, `width` and `height`, passed
#'   to `RegR::save_plt()` for PDF output. `list()` and `NULL` skip saving; a
#'   list naming only the file is completed with this figure's pinned size.
#'   Do not include the `plot` argument; it is supplied internally.
#'
#' @section Reading the curves:
#' A ranking no better than chance keeps both curves at zero. A TOC above
#' zero on the left means the patients ranked first gain more than the
#' average patient; the Qini curve shows, in the units of the effect, what
#' treating them instead of a random share of the same size gains. The
#' curves start at q = 5%, as the top of a smaller share rests on too few
#' patients; the AUTOC and QINI use the whole ranking. They test whether the
#' ranking targets the effect at all. Whether the gain is worth acting on is
#' a clinical judgement on the size of the Qini curve.
#'
#' With `"cate"` the result depends on the split, and only the held-out
#' patients inform the test; `seed` fixes the split.
#'
#' @return A `ggplot` with one panel per `type`. The y axis is on the `"diff"`
#'   scale of [get_hte()]: the S(t) or RMST difference for a survival outcome,
#'   the risk difference for a binary one and the mean difference otherwise.
#'   `attr(p, "rate")` is a tibble with one row per rule, and for two rules
#'   their difference, and target: `rule`, `target` (`"AUTOC"` or `"QINI"`),
#'   `estimate`, `std.error` (grf's half-sample bootstrap), `conf.low`,
#'   `conf.high`, `p.value` (two-sided Wald test of zero) and `n`, the
#'   patients evaluated. The pinned size is in `attr(p, "plot_size")`. If
#'   `save` is non-empty, the plot is also written to PDF through
#'   `RegR::save_plt()`.
#'
#' @references
#' Yadlowsky S, Fleming S, Shah N, Brunskill E, Wager S (2025). Evaluating
#' treatment prioritization rules via rank-weighted average treatment
#' effects. \emph{Journal of the American Statistical Association}
#' 120(549).
#'
#' @seealso [get_hte()]; [plt_hte_dep()] and [plt_hte_sub()] for how the effect
#'   varies with one covariate.
#'
#' @examplesIf requireNamespace("grf", quietly = TRUE)
#' \donttest{
#' set.seed(20260924)
#' n <- 800
#' d <- data.frame(age    = round(runif(n, 30, 85)),
#'                 marker = rnorm(n),
#'                 stage  = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
#' d$z <- rbinom(n, 1, 0.5)
#' d$y <- rbinom(n, 1, plogis(-1 + 0.02 * (d$age - 60) +
#'                              d$z * (0.1 + 0.8 * (d$marker > 0))))
#' res <- get_hte(d, cat_var = "z", adj_var = c("age", "marker", "stage"),
#'                surv = "y", grf_args = list(num.trees = 500, seed = 1))
#'
#' # The forest CATE, learnt on half of the patients and evaluated on the rest
#' p <- plt_hte_rate(res)
#' p
#' attr(p, "rate")
#'
#' # The forest against the marker alone
#' plt_hte_rate(res, priority = c("cate", "marker"))
#' }
#'
#' @export
plt_hte_rate <- function(x,
                         priority   = "cate",
                         type       = c("toc", "qini"),
                         smooth     = 0,
                         conf_level = 0.95,
                         train_frac = 0.5,
                         seed       = NULL,
                         title      = NULL,
                         save       = list()) {

  if (!inherits(x, "hte_res"))
    stop("`x` must be an `hte_res` object from get_hte().", call. = FALSE)
  if (!requireNamespace("grf", quietly = TRUE))
    stop("Package 'grf' is required for plt_hte_rate().", call. = FALSE)
  type <- intersect(c("toc", "qini"), match.arg(type, several.ok = TRUE))
  # A LOESS fit on the 96 points of the curve needs a span of 0.05 or more.
  if (!is.numeric(smooth) || length(smooth) != 1L || is.na(smooth) ||
      !(smooth == 0 || (smooth >= 0.05 && smooth <= 1)))
    stop("`smooth` must be 0 (no smoothing) or a LOESS span from 0.05 to 1.",
         call. = FALSE)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
                         !is.finite(seed)))
    stop("`seed` must be `NULL` or a single number.", call. = FALSE)
  if (!is.null(save) && !is.list(save))
    stop("`save` must be `NULL` or a list.", call. = FALSE)

  a   <- attr(x, "analysis")
  d   <- x$data
  fit <- x$fit
  n   <- nrow(d)
  ref <- setdiff(levels(factor(d[[a$cat_var]])), a$treated)[1L]

  # ---- Ranking rules ---------------------------------------------------------
  if (!is.character(priority) || !length(priority) || anyNA(priority) ||
      length(unique(priority)) > 2L)
    stop("`priority` takes one or two ranking rules: \"cate\" or numeric columns of `x$data`.",
         call. = FALSE)
  priority <- unique(priority)
  learn    <- "cate" %in% priority
  cols     <- setdiff(priority, "cate")
  if (learn) {
    if (!is.numeric(train_frac) || length(train_frac) != 1L ||
        is.na(train_frac) || train_frac <= 0 || train_frac >= 1)
      stop("`train_frac` must be a single number strictly between 0 and 1.",
           call. = FALSE)
  } else if (!missing(train_frac)) {
    stop("`train_frac` only applies when `priority` includes \"cate\"; a pre-specified rule is evaluated on every patient.",
         call. = FALSE)
  }
  own <- intersect(cols, c(".cate", ".dr_score"))
  if (length(own))
    stop(sprintf("`priority` %s would rank the patients by the forest fitted to their own outcomes; use \"cate\", which learns the ranking on a training split.",
                 paste0("`", own, "`", collapse = ", ")), call. = FALSE)
  used <- intersect(cols, c(a$cat_var, a$outcome))
  if (length(used))
    stop(sprintf("`priority` must not use the treatment or the outcome: %s.",
                 paste0("`", used, "`", collapse = ", ")), call. = FALSE)
  miss <- setdiff(cols, names(d))
  if (length(miss))
    stop(sprintf("`priority` names no column of `x$data`: %s.",
                 paste0("`", miss, "`", collapse = ", ")), call. = FALSE)
  odd <- cols[!vapply(d[cols], function(v) is.numeric(v) || is.logical(v),
                      logical(1L))]
  if (length(odd))
    stop(sprintf("`priority` column %s must be numeric or logical, larger values ranked first; code a factor as numbers.",
                 paste0("`", odd, "`", collapse = ", ")), call. = FALSE)
  outside <- setdiff(cols, a$covariates)
  if (length(outside))
    cli::cli_inform(c("i" = paste(
      "{.field {outside}} {?is/are} not a forest covariate: the evaluation",
      "holds for a function of the covariates or a variable that is no",
      "confounder; otherwise add it with get_hte(adj_var = ).")))
  ok <- if (length(cols)) stats::complete.cases(d[cols]) else rep(TRUE, n)
  if (!all(ok)) {
    gaps <- cols[vapply(d[cols], anyNA, logical(1L))]
    cli::cli_inform(c("i" = "Left out {sum(!ok)} patient{?s} missing {.field {gaps}}."))
  }

  # One seed drives the split, the refits and grf's bootstrap; the session's
  # random number stream is put back afterwards.
  if (is.null(seed)) seed <- a$seed
  genv <- globalenv()
  old  <- if (exists(".Random.seed", envir = genv, inherits = FALSE))
    get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit({
    if (!is.null(old)) assign(".Random.seed", old, envir = genv)
    else if (exists(".Random.seed", envir = genv, inherits = FALSE))
      rm(".Random.seed", envir = genv)
  }, add = TRUE)
  set.seed(seed)

  # ---- Learn the CATE on one split, evaluate on the other --------------------
  surv <- identical(a$outcome_type, "survival")
  X    <- fit$X.orig
  W    <- fit$W.orig
  if (learn) {
    train  <- sort(sample.int(n, floor(train_frac * n)))
    rows   <- setdiff(seq_len(n), train)
    halves <- list(training = train, evaluation = rows)
    for (h in names(halves)) {
      w <- W[halves[[h]]]
      if (sum(w == 1L) < 2L || sum(w == 0L) < 2L)
        stop(sprintf("The %s half needs at least two patients in each arm of `%s`; choose `train_frac` nearer 0.5.",
                     h, a$cat_var), call. = FALSE)
    }
    if (surv) {
      # get_hte()'s follow-up rule, within each half
      beyond <- .hte_beyond(x)
      past   <- if (identical(a$target, "RMST")) "up to" else "beyond"
      arm    <- c(a$treated, ref)
      n_after <- lapply(halves, function(i)
        c(sum(beyond[i] & W[i] == 1L), sum(beyond[i] & W[i] == 0L)))
      for (h in names(halves))
        if (any(n_after[[h]] == 0L))
          stop(sprintf("No patient with `%s` = %s in the %s half is followed %s `time` = %s, so the effect at that time is not identified there; try another `seed` or `train_frac`, or a smaller `time` in get_hte().",
                       a$cat_var, arm[which(n_after[[h]] == 0L)[1L]], h, past,
                       format(a$time)), call. = FALSE)
      for (h in names(halves))
        if (any(n_after[[h]] < 10L))
          warning(sprintf("Few patients in the %s half are followed %s `time` = %s: %s.",
                          h, past, format(a$time),
                          paste(sprintf("%d with `%s` = %s", n_after[[h]],
                                        a$cat_var, arm), collapse = ", ")),
                  call. = FALSE)
    }

    ga <- a$grf_args
    if (is.null(ga))   # an hte_res from before get_hte() kept its grf_args
      ga <- c(list(num.trees = fit[["_num_trees"]]),
              if (length(fit$clusters)) list(clusters = fit$clusters),
              if (!is.null(fit$sample.weights))
                list(sample.weights = fit$sample.weights),
              if (surv) list(target = a$target))
    # grf keeps survival times cut at the horizon, so refit from `$data`
    Y <- if (surv) d[[a$outcome[1L]]] else fit$Y.orig
    D <- if (surv) as.integer(d[[a$outcome[2L]]])
    refit <- function(i) {
      args <- ga
      for (f in intersect(c("W.hat", "Y.hat", "sample.weights", "clusters"),
                          names(args)))
        if (length(args[[f]]) == n) args[[f]] <- args[[f]][i]
      args$seed <- seed
      do.call(if (surv) grf::causal_survival_forest else grf::causal_forest,
              c(list(X = X[i, , drop = FALSE], Y = Y[i], W = W[i]),
                if (surv) list(D = D[i], horizon = a$time), args))
    }
    cate   <- stats::predict(refit(train), X[rows, , drop = FALSE])$predictions
    forest <- refit(rows)
  } else {
    rows   <- seq_len(n)
    forest <- fit
  }

  # ---- RATE ------------------------------------------------------------------
  pri <- data.frame(lapply(stats::setNames(priority, priority), function(v)
    if (v == "cate") as.numeric(cate) else as.numeric(d[[v]][rows])),
    check.names = FALSE)
  sub  <- which(ok[rows])
  qs   <- (5:100) / 100
  rate <- lapply(c(AUTOC = "AUTOC", QINI = "QINI"), function(tg)
    grf::rank_average_treatment_effect(forest, pri[sub, , drop = FALSE],
                                       target = tg, q = qs, subset = sub))

  z    <- stats::qnorm(1 - (1 - conf_level) / 2)
  rule <- if (length(priority) == 2L)
    c(priority, paste(priority, collapse = " - ")) else priority
  tbl <- do.call(rbind, lapply(names(rate), function(tg) {
    est <- unname(rate[[tg]]$estimate)
    se  <- unname(rate[[tg]]$std.err)
    data.frame(rule = rule, target = tg, estimate = est, std.error = se,
               conf.low = est - z * se, conf.high = est + z * se,
               p.value = 2 * stats::pnorm(-abs(est / se)), n = length(sub),
               stringsAsFactors = FALSE)
  }))
  tbl <- tibble::as_tibble(tbl)

  # ---- Curves: the TOC, and q x TOC for the Qini panel -----------------------
  lab <- ifelse(priority == "cate", "Forest CATE", priority)
  toc <- rate$AUTOC$TOC
  toc <- toc[toc$priority %in% priority, , drop = FALSE]
  toc$rule <- factor(toc$priority, levels = priority, labels = lab)
  toc$conf.low  <- toc$estimate - z * toc$std.err
  toc$conf.high <- toc$estimate + z * toc$std.err
  # Smoothing touches only what is drawn, one rule at a time.
  if (smooth > 0)
    for (r in priority) {
      i <- toc$priority == r
      for (v in c("estimate", "conf.low", "conf.high"))
        toc[[v]][i] <- as.numeric(stats::predict(stats::loess(
          y ~ q, data = data.frame(q = toc$q[i], y = toc[[v]][i]),
          span = smooth)))
    }
  band <- function(panel, k)
    data.frame(panel = panel, q = toc$q, rule = toc$rule,
               y = k * toc$estimate, conf.low = k * toc$conf.low,
               conf.high = k * toc$conf.high)
  pd <- do.call(rbind, list(toc = band("TOC", 1), qini = band("Qini", toc$q))[type])
  pd$panel <- factor(pd$panel, levels = c(toc = "TOC", qini = "Qini")[type])
  rownames(pd) <- NULL

  what <- switch(a$outcome_type,
                 survival   = sprintf("%s(%s)",
                                      if (identical(a$target, "RMST")) "RMST" else "S",
                                      format(a$time)),
                 binary     = "risk",
                 continuous = "mean")
  fmt_p <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "< 0.001"
                       else sprintf("= %.3f", p)
  rlab <- if (length(priority) == 2L)
    c(lab, paste(lab, collapse = " - ")) else lab
  # One line per target, broken between rules rather than inside one.
  lines_of <- function(tg, width) {
    r     <- tbl[tbl$target == tg, ]
    items <- sprintf("%s %.3f (%.3f to %.3f), p %s", rlab, r$estimate,
                     r$conf.low, r$conf.high, vapply(r$p.value, fmt_p, ""))
    out <- paste0(tg, ": ", items[1L])
    for (it in items[-1L]) {
      k <- length(out)
      if (nchar(out[k]) + 2L + nchar(it) <= width) {
        out[k] <- paste0(out[k], "; ", it)
      } else {
        out <- c(out, paste0("  ", it))
      }
    }
    out
  }
  caption <- paste0(paste(c(
    if (learn) {
      sprintf("Ranking learnt on %d patients and evaluated on %d held-out patients (seed %s)",
              length(train), length(sub), format(seed))
    } else {
      sprintf("Evaluated on %s%d patients with the doubly robust scores of the stored forest",
              if (all(ok)) "all " else "", length(sub))
    },
    if (smooth > 0) sprintf("curves smoothed by LOESS span %s", format(smooth)),
    sprintf("shaded: %g%% CI", 100 * conf_level)), collapse = "; "), ".")
  size <- c(if (length(type) == 2L) 9.5 else 5.5, 4.4)
  # About 13 characters fit per inch; wrap so nothing is cut off.
  width    <- floor(13 * size[1L])
  subtitle <- unlist(lapply(c(toc = "AUTOC", qini = "QINI")[type],
                            lines_of, width = width))
  caption  <- strwrap(caption, width = width)
  size[2L] <- size[2L] + 0.2 * (length(subtitle) + length(caption) - 2L)

  pal   <- stats::setNames(c("firebrick", "steelblue")[seq_along(lab)], lab)
  strip <- c(TOC = "TOC: effect in the top q minus the ATE",
             Qini = "Qini: gain over treating a random q")
  p <- ggplot2::ggplot(pd, ggplot2::aes(x = q, y = y, colour = rule,
                                        fill = rule)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50", linetype = 2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax = conf.high),
                         alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~panel, scales = "free_y",
                        labeller = ggplot2::as_labeller(strip)) +
    ggplot2::scale_x_continuous(limits = c(0, 1),
                                labels = function(v) paste0(round(100 * v), "%"),
                                expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_colour_manual(NULL, values = pal,
                                 aesthetics = c("colour", "fill")) +
    ggplot2::labs(x = "Treated share q, highest priority first",
                  y = sprintf("%s difference, %s - %s", what, a$treated, ref),
                  title = title, subtitle = paste(subtitle, collapse = "\n"),
                  caption = paste(caption, collapse = "\n")) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = if (length(lab) > 1L) "top" else "none",
                   plot.subtitle = ggplot2::element_text(size = 9))

  attr(p, "rate")      <- tbl
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
