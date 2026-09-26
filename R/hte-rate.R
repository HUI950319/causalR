# =============================================================================
# hte-rate.R -- how well a ranking of the patients targets the effect (RATE)
# =============================================================================
#
# Architecture:
#
#   L1  plt_hte_rate()  rank the patients, evaluate the ranking with grf's
#                       rank-weighted average treatment effect and draw its
#                       TOC and Qini curves, and on request its GATES
#
# A ranking learnt from the outcomes -- the forest CATE -- is learnt on one
# random split and evaluated on the other, as grf's RATE vignette requires. A
# pre-specified score is not learnt from them, so it is evaluated on every
# patient with the forest get_hte() stored. With unit costs the Qini curve is
# q x TOC(q), whose area is grf's QINI, so it needs no further estimation.
# The GATES of a group, grf's average_treatment_effect() on its patients, is
# the ATE plus the mean slope of the Qini curve over the group's share.
# =============================================================================


#' Targeting curves for heterogeneous treatment effects
#'
#' Evaluates how well a ranking of the patients picks out those who benefit
#' most, with the rank-weighted average treatment effect (RATE) of
#' [grf::rank_average_treatment_effect()], and draws its targeting operator
#' characteristic (TOC) and Qini curves from a [get_hte()] result. On request
#' it also draws the sorted group average treatment effects (GATES): the
#' doubly robust effect within each group of the ranking.
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
#'       patients of both arms followed past `time`. With clusters, whole
#'       clusters are split and each half needs at least two clusters.}
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
#' @param type Panels to draw, any of `"toc"`, `"qini"` and `"gates"`
#'   (default `"toc"` and `"qini"`).
#'   \describe{
#'     \item{`"toc"`}{TOC(q): the average effect among the share q of
#'       patients ranked first minus the average effect of all of them. Its
#'       area, the AUTOC, weights the top of the ranking most, so it detects
#'       an effect held by a few patients.}
#'     \item{`"qini"`}{q x TOC(q): what treating the share q ranked first
#'       gains over treating a random share q. Its area, the QINI, weights
#'       every share alike, so it detects an effect that changes gradually.}
#'     \item{`"gates"`}{The average effect within each of
#'       `gates_args$n_groups` groups of equal size, highest priority first;
#'       a ranking that targets the effect gives estimates falling from left
#'       to right. For `"cate"`, a diamond marks the mean forest CATE of each
#'       group.}
#'   }
#' @param smooth `0` (default) draws the TOC and its band as estimated; a
#'   number from `0.05` to `1` is a LOESS span that smooths both for display,
#'   separately for each rule: about `0.1` smooths slightly, `0.3` strongly.
#'   The Qini panel is then q times the smoothed TOC. The AUTOC, QINI and
#'   `attr(p, "rate")` do not change.
#' @param gates_args Named list for the `"gates"` panel: `n_groups` (default
#'   `5`), the number of groups the ranking is cut into, a whole number of at
#'   least 2. Tied patients share a group, so a logical rule gives two. Only
#'   used when `type` includes `"gates"`.
#' @param conf_level Confidence level of the bands and intervals. Default
#'   `0.95`.
#' @param train_frac Share of the patients the ranking forest is refitted to,
#'   strictly between 0 and 1. Default `0.5`. Only used when `priority`
#'   includes `"cate"`. With clusters this is the share of clusters, so the
#'   share of patients can differ when cluster sizes vary.
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
#' The GATES panel draws each group's doubly robust average treatment effect,
#' [grf::average_treatment_effect()] with `subset` as in [plt_hte_sub()], over
#' the share of the ranking the group covers. A group's estimate is the ATE
#' of the evaluated patients plus the mean slope of the Qini curve over that
#' share, so the panel shows where along the ranking the effect changes. The
#' subtitle tests the top minus the bottom group, taking the two as
#' independent. For `"cate"`, diamonds spread less than the estimates mean
#' the forest shrinks the effect towards its mean.
#'
#' @return A `ggplot` with one panel per `type`. The y axis is on the `"diff"`
#'   scale of [get_hte()]: the S(t) or RMST difference for a survival outcome,
#'   the risk difference for a binary one and the mean difference otherwise.
#'   `attr(p, "rate")` is a tibble with one row per rule, and for two rules
#'   their difference, and target: `rule`, `target` (`"AUTOC"` or `"QINI"`),
#'   `estimate`, `std.error` (grf's half-sample bootstrap), `conf.low`,
#'   `conf.high`, `p.value` (two-sided Wald test of zero) and `n`, the
#'   patients evaluated. With `"gates"` in `type`, `attr(p, "gates")` is a
#'   tibble with one row per rule and group, and per rule one for the top
#'   minus the bottom group: `rule`, `group` (`"1"` ranked first, `"1 - K"`
#'   the difference), `q_from` and `q_to` (the share of the ranking the group
#'   covers, `NA` for the difference), `n`, `estimate`, `std.error`,
#'   `conf.low`, `conf.high`, `p.value` (two-sided Wald test of zero) and
#'   `cate_mean`, the mean forest CATE for `"cate"` (`NA` for a pre-specified
#'   rule). The pinned size is in `attr(p, "plot_size")`. If `save` is
#'   non-empty, the plot is also written to PDF through `RegR::save_plt()`.
#'
#' @references
#' Yadlowsky S, Fleming S, Shah N, Brunskill E, Wager S (2025). Evaluating
#' treatment prioritization rules via rank-weighted average treatment
#' effects. \emph{Journal of the American Statistical Association}
#' 120(549).
#'
#' Chernozhukov V, Demirer M, Duflo E, \enc{Fernández}{Fernandez}-Val I
#' (2025). Fisher-Schultz lecture: generic machine learning inference on
#' heterogeneous treatment effects in randomized experiments, with an
#' application to immunization in India. \emph{Econometrica} 93(4),
#' 1121-1164.
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
#'
#' # The effect within each fifth of the forest ranking
#' plt_hte_rate(res, type = "gates")
#' }
#'
#' @export
plt_hte_rate <- function(x,
                         priority   = "cate",
                         type       = c("toc", "qini"),
                         smooth     = 0,
                         gates_args = list(n_groups = 5),
                         conf_level = 0.95,
                         train_frac = 0.5,
                         seed       = NULL,
                         title      = NULL,
                         save       = list()) {

  if (!inherits(x, "hte_res"))
    stop("`x` must be an `hte_res` object from get_hte().", call. = FALSE)
  if (!requireNamespace("grf", quietly = TRUE))
    stop("Package 'grf' is required for plt_hte_rate().", call. = FALSE)
  panels <- c(toc = "TOC", qini = "Qini", gates = "GATES")
  type   <- intersect(names(panels),
                      match.arg(type, names(panels), several.ok = TRUE))
  curves <- intersect(type, c("toc", "qini"))
  gates  <- "gates" %in% type
  # A LOESS fit on the 96 points of the curve needs a span of 0.05 or more.
  if (!is.numeric(smooth) || length(smooth) != 1L || is.na(smooth) ||
      !(smooth == 0 || (smooth >= 0.05 && smooth <= 1)))
    stop("`smooth` must be 0 (no smoothing) or a LOESS span from 0.05 to 1.",
         call. = FALSE)
  if (!gates && !missing(gates_args))
    stop("`gates_args` only applies when `type` includes \"gates\".",
         call. = FALSE)
  gates_args <- .merge_named_arg(gates_args, list(n_groups = 5), "gates_args")
  n_groups <- gates_args$n_groups
  if (!is.numeric(n_groups) || length(n_groups) != 1L ||
      !is.finite(n_groups) || n_groups < 2 || n_groups != round(n_groups))
    stop("`gates_args$n_groups` must be a whole number of at least 2.",
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
    train <- if (length(fit$clusters)) {
      ids <- unique(fit$clusters)
      k <- floor(train_frac * length(ids))
      if (k < 2L || length(ids) - k < 2L)
        stop("Each split needs at least two clusters; change `train_frac` or use more clusters.",
             call. = FALSE)
      which(fit$clusters %in% ids[sample.int(length(ids), k)])
    } else {
      sort(sample.int(n, floor(train_frac * n)))
    }
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
  lab <- ifelse(priority == "cate", "Forest CATE", priority)

  # ---- GATES: the effect within each group of the ranking --------------------
  # Groups of equal size, highest priority first; tied patients share the
  # group of the middle of their tie, so a logical rule gives two. A group's
  # effect is grf's average_treatment_effect(subset = ) through .hte_estimate(),
  # as in plt_hte_sub(); the top minus the bottom group takes the two as
  # independent.
  if (gates) {
    grid <- data.frame(estimand = "ATE", measure = "diff",
                       stringsAsFactors = FALSE)
    bey  <- .hte_beyond(x)[rows]
    gt <- do.call(rbind, lapply(seq_along(priority), function(j) {
      v  <- priority[j]
      pv <- pri[[v]][sub]
      g  <- ceiling(n_groups * rank(-pv, ties.method = "average") / length(pv))
      g  <- match(g, sort(unique(g)))
      ng <- max(g)
      if (ng < 2L)
        stop(sprintf("`priority` `%s` takes one value among the evaluated patients, so it forms a single GATES group.",
                     v), call. = FALSE)
      m    <- tabulate(g, ng)
      from <- (cumsum(m) - m) / length(pv)
      one  <- do.call(rbind, lapply(seq_len(ng), function(i) {
        est <- .hte_muffle_ps(.hte_estimate(
          forest, NULL, seq_along(rows) %in% sub[g == i], grid, FALSE, z,
          sprintf("GATES of %s, group %d", lab[j], i), bey))
        data.frame(rule = v, group = as.character(i), q_from = from[i],
                   q_to = from[i] + m[i] / length(pv), n = m[i],
                   est[c("estimate", "std.error", "conf.low", "conf.high",
                         "p.value")],
                   cate_mean = if (v == "cate") mean(pv[g == i]) else NA_real_,
                   stringsAsFactors = FALSE)
      }))
      dif <- one$estimate[1L] - one$estimate[ng]
      se  <- sqrt(one$std.error[1L]^2 + one$std.error[ng]^2)
      rbind(one, data.frame(
        rule = v, group = sprintf("1 - %d", ng), q_from = NA_real_,
        q_to = NA_real_, n = m[1L] + m[ng], estimate = dif, std.error = se,
        conf.low = dif - z * se, conf.high = dif + z * se,
        p.value = 2 * stats::pnorm(-abs(dif / se)),
        cate_mean = one$cate_mean[1L] - one$cate_mean[ng],
        stringsAsFactors = FALSE))
    }))
    rownames(gt) <- NULL
    gt <- tibble::as_tibble(gt)
  }

  # ---- Curves: the TOC, and q x TOC for the Qini panel -----------------------
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
  pd <- do.call(rbind, list(toc = band("TOC", 1), qini = band("Qini", toc$q))[curves])
  if (is.null(pd)) pd <- band("TOC", 1)[0L, ]    # the GATES panel alone
  pd$panel <- factor(pd$panel, levels = panels[type])
  rownames(pd) <- NULL
  if (gates) {
    gd <- as.data.frame(gt[!is.na(gt$q_from), ])
    gd$panel <- factor("GATES", levels = panels[type])
    gd$rule  <- factor(gd$rule, levels = priority, labels = lab)
    # two rules side by side at the middle of each group
    gd$x <- (gd$q_from + gd$q_to) / 2 +
      if (length(priority) == 2L) c(-0.012, 0.012)[as.integer(gd$rule)] else 0
  }

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
  lines_of <- function(head, r, labs, width) {
    items <- sprintf("%s %.3f (%.3f to %.3f), p %s", labs, r$estimate,
                     r$conf.low, r$conf.high, vapply(r$p.value, fmt_p, ""))
    out <- paste0(head, ": ", items[1L])
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
    if (smooth > 0 && length(curves))
      sprintf("curves smoothed by LOESS span %s", format(smooth)),
    if (length(curves)) sprintf("shaded: %g%% CI", 100 * conf_level),
    if (gates) sprintf("GATES bars: %g%% CI", 100 * conf_level),
    if (gates && learn) "diamonds: mean forest CATE per group"),
    collapse = "; "), ".")
  size <- c(c(6, 10.5, 15)[length(type)], 4.4)
  # About 10 characters fit per inch; wrap so nothing is cut off.
  width    <- floor(10 * size[1L])
  subtitle <- c(
    unlist(lapply(c(toc = "AUTOC", qini = "QINI")[curves], function(tg)
      lines_of(tg, tbl[tbl$target == tg, ], rlab, width))),
    if (gates) lines_of("GATES top - bottom group", gt[is.na(gt$q_from), ],
                        lab, width))
  caption  <- strwrap(caption, width = width)
  size[2L] <- size[2L] + 0.2 * (length(subtitle) + length(caption) - 2L)

  pal   <- stats::setNames(c("firebrick", "steelblue")[seq_along(lab)], lab)
  strip <- c(TOC = "TOC: effect in the top q minus the ATE",
             Qini = "Qini: gain over treating a random q",
             GATES = "GATES: effect within each priority group")
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
    ggplot2::labs(x = "Treated fraction (q)",
                  y = sprintf("%s difference, %s - %s", what, a$treated, ref),
                  title = title, subtitle = paste(subtitle, collapse = "\n"),
                  caption = paste(caption, collapse = "\n")) +
    # A label of 100% runs about 14 pt past its panel, so the panels and the
    # right edge get room for it.
    UtilsR::theme_my(base_rect_size = 1.5, panel.spacing = 21,
                     plot.margin = c(7, 21, 7, 7)) +
    # The estimates stay left-aligned, as lines_of() indents their breaks,
    # and start at the left edge, as it never breaks inside an estimate.
    ggplot2::theme(legend.position = if (length(lab) > 1L) "top" else "none",
                   plot.subtitle = ggplot2::element_text(size = ggplot2::rel(0.8),
                                                         hjust = 0),
                   plot.title.position = "plot")
  # Each group: a line over its share and the estimate with its interval in
  # the middle; for "cate" a diamond at the group's mean forest CATE.
  if (gates) {
    ge <- gd[is.finite(gd$estimate), , drop = FALSE]
    p <- p +
      ggplot2::geom_segment(data = ge, ggplot2::aes(x = q_from, xend = q_to,
                                                    y = estimate,
                                                    yend = estimate,
                                                    colour = rule),
                            linewidth = 0.5, inherit.aes = FALSE) +
      ggplot2::geom_errorbar(data = ge, ggplot2::aes(x = x, ymin = conf.low,
                                                     ymax = conf.high,
                                                     colour = rule),
                             width = 0.02, inherit.aes = FALSE) +
      ggplot2::geom_point(data = ge, ggplot2::aes(x = x, y = estimate,
                                                  colour = rule),
                          size = 2, inherit.aes = FALSE)
    dia <- gd[!is.na(gd$cate_mean), , drop = FALSE]
    if (nrow(dia))
      p <- p + ggplot2::geom_point(data = dia, ggplot2::aes(x = x, y = cate_mean,
                                                            colour = rule),
                                   shape = 5, size = 2.5, inherit.aes = FALSE)
  }

  attr(p, "rate")      <- tbl
  if (gates) attr(p, "gates") <- gt
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
