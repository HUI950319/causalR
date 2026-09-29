# =============================================================================
# hte-icf.R -- subgroup rules from an iterative causal forest
# =============================================================================
#
# Architecture:
#
#   L1  get_hte_icf()     split the sample, cross-validate the depth, find the
#                         rules on the discovery part, estimate them on the rest
#   L2  .icf_discover()   screen, grow `n_forest` forests, keep each forest's
#                         best pruned tree per depth, vote
#   L2  .icf_leaves()     readable rule and voting key of every leaf
#   L2  .icf_assign()     leaf of each row under a partition
#   L3  print.hte_icf()
#
# The algorithm follows Wang et al. (2024). The iCF repository carries no
# licence, so nothing here is taken from its code; the differences from it
# are listed in the Details of get_hte_icf().
# =============================================================================


# ---- L2 discovery ------------------------------------------------------------

# Least drop in the pruning loss a kept split must buy, in residual variances
# of the AIPW scores: the 95% point of chi-square(1), about a 5% test that the
# two children differ.
.ICF_PENALTY <- stats::qchisq(0.95, 1)


# One discovery run on the rows `fit` was fitted to. `X` is the design matrix
# of those rows, `cols` says what each of its columns means. The forests reuse
# the nuisance estimates of `fit`. Every tree is judged on the rows it did not
# choose its splits on by the within-leaf squared error of the AIPW scores,
# inflated by m^2 / (m - 1)^2 for a leaf of m of those rows: a split survives
# only if it lowers that loss by .ICF_PENALTY residual variances. Judged on
# all rows instead, each forest's best tree kept a spurious cut in the branch
# without heterogeneity (X4 at 0.75, the same in all 5 forests). Judged
# honestly but without the penalty, as the iCF code prunes, the best of 100
# trees still kept one in 5 of 12 simulated data sets; with it, in none, and
# the weaker effects of the iCF README simulation were found as before.
#' @keywords internal
#' @noRd
.icf_discover <- function(fit, X, cols, depth, ra, seed, forest_args) {
  a   <- attr(fit, "analysis")
  f   <- fit$fit
  W   <- f$W.orig
  g   <- fit$data$.dr_score
  imp <- fit$importance

  # iCF screening: covariates at or above the mean importance, at least four
  screened <- imp$variable
  if (isTRUE(ra$screen)) {
    screened <- imp$variable[imp$importance >= mean(imp$importance)]
    k <- min(4L, nrow(imp))
    if (length(screened) < k) screened <- imp$variable[seq_len(k)]
  }
  sc <- which(cols$var %in% screened)
  Xs <- X[, sc, drop = FALSE]

  grow <- if (identical(a$outcome_type, "survival")) {
    D <- fit$data$DSS
    if (is.logical(D)) D <- as.integer(D)
    yd <- .hte_surv_yd(fit$data$time, D, a$time, a$target)
    function(s) do.call(grf::causal_survival_forest, c(
      list(X = Xs, Y = yd$Y, W = W, D = yd$D, W.hat = f$W.hat,
           horizon = a$time, target = a$target, num.trees = ra$num_trees,
           seed = s), forest_args))
  } else {
    function(s) do.call(grf::causal_forest, c(
      list(X = Xs, Y = f$Y.orig, W = W, Y.hat = f$Y.hat, W.hat = f$W.hat,
           num.trees = ra$num_trees, seed = s), forest_args))
  }

  n     <- nrow(Xs)
  dmax  <- max(depth)
  G     <- cbind(1, W, g, g^2)
  best  <- replicate(length(depth), vector("list", ra$n_forest),
                     simplify = FALSE)

  for (b in seq_len(ra$n_forest)) {
    forest <- grow(seed + b)
    top <- rep(list(list(gain = -Inf)), length(depth))
    for (t in seq_len(forest[["_num_trees"]])) {
      tr    <- grf::get_tree(forest, t)
      nodes <- tr$nodes
      leaf  <- vapply(nodes, function(x) x$is_leaf, logical(1L))
      # A tree is judged on the rows it did not choose its splits on: its
      # honest half and the rows it never drew.
      honest <- unlist(lapply(nodes[leaf], function(x) x$samples))
      e <- rep(1, n)
      e[setdiff(tr$drawn_samples, honest)] <- 0
      Ge    <- G * e
      min_n <- max(2L, ceiling(ra$min_leaf * sum(e)))
      pick  <- function(field) vapply(nodes, function(x)
        if (x$is_leaf) NA_real_ else as.numeric(x[[field]]), numeric(1L))
      var <- pick("split_variable")
      val <- pick("split_value")
      kid <- cbind(pick("left_child"), pick("right_child"))

      # Stats (patients, treated, sum and sum of squares of the scores) of the
      # node each row sits in at every depth; a row at a leaf stays there.
      cur <- rep(1L, n)
      st  <- vector("list", dmax + 1L)
      tally <- function() {
        s <- matrix(0, length(nodes), 4L)
        r <- rowsum(Ge, cur)
        s[as.integer(rownames(r)), ] <- r
        s
      }
      st[[1L]] <- tally()
      for (k in seq_len(dmax)) {
        i <- which(!leaf[cur])
        if (length(i)) {
          nd <- cur[i]
          go_left <- Xs[cbind(i, var[nd])] <= val[nd]
          cur[i] <- ifelse(go_left, kid[nd, 1L], kid[nd, 2L])
        }
        st[[k + 1L]] <- tally()
      }

      sse  <- function(s) s[4L] - s[3L]^2 / s[1L]
      root <- st[[1L]][1L, ]
      own0 <- sse(root) * root[1L]^2 / (root[1L] - 1)^2
      pen  <- .ICF_PENALTY * sse(root) / (root[1L] - 1)
      prune <- function(id, k, d) {
        s   <- st[[k + 1L]][id, ]
        own <- sse(s) * s[1L]^2 / (s[1L] - 1)^2
        if (k >= d || leaf[id]) return(list(loss = own, keep = integer()))
        cs <- st[[k + 2L]][kid[id, ], , drop = FALSE]
        if (any(cs[, 1L] < min_n | cs[, 2L] < 2 | cs[, 1L] - cs[, 2L] < 2))
          return(list(loss = own, keep = integer()))
        l <- prune(kid[id, 1L], k + 1L, d)
        r <- prune(kid[id, 2L], k + 1L, d)
        if (l$loss + r$loss + pen < own)
          list(loss = l$loss + r$loss + pen, keep = c(id, l$keep, r$keep))
        else list(loss = own, keep = integer())
      }
      for (j in seq_along(depth)) {
        p <- prune(1L, 0L, depth[j])
        gain <- (own0 - p$loss) / sum(e)
        if (gain > top[[j]]$gain)
          top[[j]] <- list(gain = gain, keep = p$keep, nodes = nodes)
      }
    }

    # Kept splits as paths from the root ("", "L", "LR", ...), columns of `X`
    for (j in seq_along(depth)) {
      nodes <- top[[j]]$nodes
      tree  <- data.frame(path = character(), col = integer(),
                          value = numeric(), stringsAsFactors = FALSE)
      todo  <- list(list(1L, ""))
      while (length(todo)) {
        id <- todo[[1L]][[1L]]
        pa <- todo[[1L]][[2L]]
        todo <- todo[-1L]
        if (id %in% top[[j]]$keep) {
          x <- nodes[[id]]
          tree[nrow(tree) + 1L, ] <- list(pa, sc[x$split_variable],
                                          x$split_value)
          todo <- c(todo, list(list(x$left_child, paste0(pa, "L")),
                               list(x$right_child, paste0(pa, "R"))))
        }
      }
      best[[j]][[b]] <- tree
    }
  }

  # Plurality vote on the partition, split values ignored; the winning
  # partition keeps its most frequent tree shape, at median split values. The
  # leaf keys are sorted, as the leaves of X1-then-X3 and X3-then-X1 trees
  # come in different orders.
  parts <- lapply(seq_along(depth), function(j) {
    trees <- best[[j]]
    keys  <- vapply(trees, function(tr)
      paste(sort(.icf_leaves(tr, cols)$key), collapse = " | "), character(1L))
    u     <- unique(keys)
    win   <- u[which.max(tabulate(match(keys, u)))]
    wins  <- which(keys == win)
    shape <- vapply(trees[wins], function(tr)
      paste(tr$path, tr$col, collapse = ";"), character(1L))
    su    <- unique(shape)
    same  <- wins[shape == su[which.max(tabulate(match(shape, su)))]]
    tree  <- trees[[same[1L]]]
    if (nrow(tree))
      tree$value <- apply(matrix(vapply(trees[same], function(tr) tr$value,
                                        numeric(nrow(tree))),
                                 nrow = nrow(tree)), 1L, stats::median)
    list(tree = tree, leaves = .icf_leaves(tree, cols),
         share = length(wins) / length(keys))
  })

  list(parts = parts, g = g, fit = fit, screened = screened)
}


# ---- L2 rules ----------------------------------------------------------------

# Leaves of a tree, left to right, with the rule a reader sees and the key the
# vote compares: the same conditions without the split values of continuous
# covariates, so the order of the splits and small moves of a cut-off do not
# split the vote.
#' @keywords internal
#' @noRd
.icf_leaves <- function(tree, cols) {
  if (!nrow(tree))
    return(data.frame(path = "", rule = "All patients", key = "",
                      stringsAsFactors = FALSE))
  lv  <- attr(cols, "levels")
  fmt <- function(x) format(signif(x, 4L), trim = TRUE, scientific = FALSE)
  kids  <- c(paste0(tree$path, "L"), paste0(tree$path, "R"))
  paths <- sort(setdiff(kids, tree$path))
  one <- function(p) {
    k     <- seq_len(nchar(p))
    j     <- match(substring(p, 1L, k - 1L), tree$path)
    dir   <- substring(p, k, k)
    col   <- tree$col[j]
    value <- tree$value[j]
    vars  <- unique(cols$var[col])
    cond  <- vapply(vars, function(v) {
      h <- cols$var[col] == v
      if (cols$type[col[h][1L]] == "num") {
        lo <- suppressWarnings(max(value[h & dir == "R"]))
        hi <- suppressWarnings(min(value[h & dir == "L"]))
        text <- if (is.finite(lo) && is.finite(hi))
          sprintf("%s < %s <= %s", fmt(lo), v, fmt(hi))
        else if (is.finite(lo)) sprintf("%s > %s", v, fmt(lo))
        else sprintf("%s <= %s", v, fmt(hi))
        return(c(text, paste0(v, if (is.finite(lo)) " >",
                              if (is.finite(hi)) " <=")))
      }
      all <- lv[[v]]
      ok  <- rep(TRUE, length(all))
      for (i in which(h)) {
        left <- switch(cols$type[col[i]],
                       bin   = c(0, 1) <= value[i],
                       code  = seq_along(all) <= value[i],
                       level = all != cols$level[col[i]])
        ok <- ok & if (dir[i] == "L") left else !left
      }
      text <- if (sum(ok) == 1L) sprintf("%s = %s", v, all[ok])
        else if (sum(!ok) == 1L) sprintf("%s != %s", v, all[!ok])
        else sprintf("%s in {%s}", v, paste(all[ok], collapse = ", "))
      c(text, paste0(v, " = ", paste(all[ok], collapse = ",")))
    }, character(2L))
    c(paste(cond[1L, ], collapse = " & "),
      paste(sort(cond[2L, ]), collapse = " & "))
  }
  out <- vapply(paths, one, character(2L))
  data.frame(path = paths, rule = make.unique(out[1L, ], sep = " #"),
             key = out[2L, ], stringsAsFactors = FALSE, row.names = NULL)
}

# Leaf path of each row of `X` under a tree of kept splits.
#' @keywords internal
#' @noRd
.icf_assign <- function(tree, X) {
  path <- rep("", nrow(X))
  repeat {
    j <- match(path, tree$path)
    i <- which(!is.na(j))
    if (!length(i)) return(path)
    j <- j[i]
    path[i] <- paste0(path[i], ifelse(X[cbind(i, tree$col[j])] <= tree$value[j],
                                      "L", "R"))
  }
}


# ---- L1 public function ------------------------------------------------------

#' Subgroup rules from an iterative causal forest
#'
#' Looks for subgroups with different treatment effects without naming them in
#' advance, following the iterative causal forest (iCF) of Wang et al. (2024),
#' and estimates their effects honestly. The data are split in two: rules such
#' as `X1 = 1 & X3 = 0` are found on the discovery part, with the tree depth
#' chosen by cross-validation, and each rule's doubly robust effect is then
#' estimated by [get_hte()] on the estimation part, which played no role in
#' finding it. Arguments follow [get_hte()] where the two overlap.
#'
#' @param data A data frame holding every column named below.
#' @param cat_var Length-1 character. The binary exposure column, coded as in
#'   [get_hte()].
#' @param adj_var Character vector of covariates the forests condition on and
#'   the rules are written in: numeric, factor, character or logical, with no
#'   missing value, since a rule cannot say where a missing value goes.
#'   Rows missing `cat_var` or the outcome are dropped, as in [get_hte()].
#' @param surv Outcome selector, as in [get_hte()]: `TRUE` (default) for the
#'   survival columns `time` and `DSS`, or a single binary (0/1) or continuous
#'   outcome column.
#' @param time Time point of a survival outcome, as in [get_hte()]. Default
#'   `120`. Only accepted with `surv = TRUE`.
#' @param depth Candidate tree depths, positive whole numbers. Default `1:3`,
#'   up to 8 subgroups. Depth 0 -- one group of every patient -- is always a
#'   candidate as well, so "no subgroups" is a possible answer.
#' @param factor_encoding How factor, character and logical covariates enter
#'   the forests, as in [get_hte()]: `"onehot"` (default) or `"integer"`. The
#'   rules are written in the original levels either way.
#' @param split_frac Share of the patients, within each arm, used for
#'   discovery. Default `0.5`; the rest are the estimation part.
#' @param rule_args Named list of discovery settings; partial overrides keep
#'   the other defaults:
#'   \describe{
#'     \item{`n_forest`}{Forests grown per discovery run, each voting once.
#'       Default `20`.}
#'     \item{`num_trees`}{Trees per forest. Default `200`.}
#'     \item{`n_folds`}{Cross-validation folds within the discovery part.
#'       Default `5`.}
#'     \item{`min_leaf`}{Smallest subgroup, as a share of the patients the
#'       rule is found on; each also needs two patients per arm. Default
#'       `0.05`.}
#'     \item{`screen`}{`TRUE` (default) grows the voting forests on the
#'       covariates whose [get_hte()] importance is at least the mean, and on
#'       at least four; `FALSE` uses all of `adj_var`.}
#'     \item{`gate`}{Largest one-sided p-value of the discovery part's
#'       calibration test (`differential.forest.prediction`, see [get_hte()])
#'       at which subgroups are reported; above it the answer is depth 0.
#'       Default `0.1`, as in Wang et al. (2024); `1` leaves the choice to the
#'       cross-validation alone.}
#'   }
#' @param grf_args Named list forwarded to every [get_hte()] fit (and so to
#'   [grf::causal_forest()] or [grf::causal_survival_forest()]) and, except
#'   `num.trees` and `seed`, to the voting forests. Per-row fields (`W.hat`,
#'   `Y.hat`, `sample.weights`, `clusters`) are not accepted, because every
#'   fit sees a different part of the rows.
#' @param seed Nonnegative whole number, default `123`. It draws the split and
#'   the folds, seeds the [get_hte()] fits unless `grf_args` sets `seed`, and
#'   seeds voting forest b with `seed + b`. The caller's random-number state is
#'   restored, also on error.
#' @param verbose Logical. `TRUE` reports the split, the screened covariates
#'   and the chosen depth. Default `FALSE`.
#'
#' @details
#' One discovery run, on some set of patients, goes as follows. [get_hte()]
#' gives the out-of-bag outcome and propensity estimates, the AIPW score of
#' every patient and the variable importance. `n_forest` causal forests are
#' grown on the screened covariates, reusing those estimates. Every tree is
#' judged on the patients it did not choose its splits on -- grf's honest half
#' of its subsample and the patients it never drew, about three quarters. It
#' is cut at each candidate depth and pruned bottom-up: a split is kept only if
#' both children hold at least `min_leaf` of those patients and two of each
#' arm, and it lowers the within-leaf squared error of their AIPW scores,
#' inflated by \eqn{m^2/(m-1)^2} for a leaf of \eqn{m} patients as in the iCF
#' implementation, by more than 3.84 (the 95% point of \eqn{\chi^2_1}) times
#' their variance -- about a 5% test that the two children differ. Without
#' that margin, the best of 100 trees kept a spurious split in 5 of 12
#' simulated data sets. Each forest's tree with the largest reduction of the
#' loss per patient is its vote. The partition with the most votes wins; its
#' share of the votes is reported as stability.
#'
#' The depth is chosen by `n_folds`-fold cross-validation within the
#' discovery part: a discovery run on the other folds, from its own
#' [get_hte()] fit, gives a partition per depth, and its leaf means of the
#' AIPW scores predict the held-out patients' scores. The depth with the
#' smallest mean squared error wins, depth 0 (the overall mean) included; a
#' depth whose final partition has one leaf is not eligible, and one whose
#' partition equals a shallower depth's is reported as that depth. The rules
#' are the discovery run on the whole discovery part at that depth. As in the
#' paper, subgroups are only reported when the calibration test of the
#' discovery part's forest finds heterogeneity at `gate`: in 12 simulated
#' data sets with a constant effect the cross-validation alone reported
#' subgroups twice, and with the gate never, while it let through all eight
#' data sets with an effect modifier.
#'
#' Compared with the iCF code, which served as a description only: the loss
#' uses AIPW scores rather than the R-loss, so survival outcomes work too;
#' depth is set by cutting trees rather than by tuning the minimum leaf size
#' for each depth, so one set of forests serves every depth; every tree is
#' judged on patients its splits were not chosen on, where the iCF code sums
#' the loss over each tree's own leaf samples, which differ in number from
#' tree to tree; the vote compares partitions rather than tree
#' shapes, so splitting X1 before X3 or after it is the same vote; every
#' cross-validation fold refits all nuisance estimates and forests on its
#' training folds only; depth 0 competes in the cross-validation besides the
#' calibration gate (the iCF README's `P_threshold = 1` turns that gate off);
#' and effects are estimated on patients not used to find the rules.
#'
#' @return An object of class `hte_icf`: a list of
#'   \describe{
#'     \item{`rules`}{Tibble with one row per subgroup of the chosen depth:
#'       `leaf`, `rule`, `n_disc` (discovery patients), and from the
#'       estimation part `n`, `n_treat`, `estimate`, `std.error`, `conf.low`,
#'       `conf.high`, `p.value` and `p_inter`, the doubly robust ATE
#'       difference with 95% Wald intervals, as in [get_hte()]'s `$subgroup`.}
#'     \item{`cv`}{Tibble with one row per depth, 0 included: `depth`,
#'       `n_leaf` of the partition found on the whole discovery part,
#'       `cv_loss`, `std.error` and `selected`.}
#'     \item{`vote`}{Tibble with one row per candidate depth: `depth`,
#'       `n_leaf`, `share` of the votes won by the partition, and `partition`,
#'       its rules.}
#'     \item{`importance`}{The discovery part's [get_hte()] importance, with
#'       `screened` marking the covariates the voting forests used.}
#'     \item{`est`}{The `hte_res` from [get_hte()] on the estimation part,
#'       whose `$data` has the rule of every patient in the factor `.rule`:
#'       `plt_hte_sub(res$est, sub_var = ".rule")` draws the subgroups.}
#'   }
#'   Analysis metadata is attached as `attr(x, "analysis")`, including the row
#'   numbers of `data` in the discovery part (`discovery`), the settings, the
#'   discovery part's calibration p-value (`calibration_p`) and whether it
#'   closed the gate (`gated`).
#'
#' @references
#' Wang T, Keil AP, Kim S, Wyss R, Htoo PT, Funk MJ, Buse JB, Kosorok MR,
#' \enc{Stürmer}{Sturmer} T (2024). Iterative causal forest: a novel algorithm
#' for subgroup identification. \emph{American Journal of Epidemiology}
#' 193(5):764-776.
#'
#' Athey S, Imbens G (2016). Recursive partitioning for heterogeneous causal
#' effects. \emph{Proceedings of the National Academy of Sciences}
#' 113(27):7353-7360.
#'
#' @seealso [get_hte()] for the forest and effect estimates; [plt_hte_sub()]
#'   to draw the subgroups; [get_hte_select()] to screen effect modifiers.
#'
#' @examplesIf requireNamespace("grf", quietly = TRUE)
#' \donttest{
#' set.seed(20260928)
#' n <- 1600
#' d <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rnorm(n),
#'                 X3 = rbinom(n, 1, 0.5), X4 = rnorm(n))
#' d$z <- rbinom(n, 1, plogis(0.4 * d$X1 - 0.3 * d$X2))
#' d$y <- d$X2 + d$z * 2 * d$X1 * d$X3 + rnorm(n)
#'
#' res <- get_hte_icf(d, cat_var = "z", adj_var = c("X1", "X2", "X3", "X4"),
#'                    surv = "y",
#'                    rule_args = list(n_forest = 5, num_trees = 100,
#'                                     n_folds = 3))
#' res
#' plt_hte_sub(res$est, sub_var = ".rule")
#' }
#'
#' @export
get_hte_icf <- function(data,
                        cat_var,
                        adj_var,
                        surv       = TRUE,
                        time       = 120,
                        depth      = 1:3,
                        factor_encoding = c("onehot", "integer"),
                        split_frac = 0.5,
                        rule_args  = list(),
                        grf_args   = list(),
                        seed       = 123,
                        verbose    = FALSE) {

  factor_encoding <- match.arg(factor_encoding)
  if (!requireNamespace("grf", quietly = TRUE))
    stop("Package 'grf' is required for get_hte_icf().", call. = FALSE)
  if (!is.data.frame(data) || !nrow(data))
    stop("`data` must be a non-empty data frame.", call. = FALSE)
  cat_var <- .sens_check_col(cat_var, data, "cat_var", n = 1L)
  adj_var <- setdiff(.sens_check_col(adj_var, data, "adj_var"), cat_var)
  if (!length(adj_var))
    stop("`adj_var` must name at least one covariate.", call. = FALSE)
  if (isFALSE(surv))
    stop("`surv = FALSE` (competing risks) is not supported: grf has no competing-risk forest.",
         call. = FALSE)
  is_surv <- isTRUE(surv)
  if (is_surv) {
    outcome <- c("time", "DSS")
    if (!all(outcome %in% names(data)))
      stop("`surv = TRUE` requires the columns `time` and `DSS`.", call. = FALSE)
  } else {
    if (!is.character(surv) || length(surv) != 1L || is.na(surv))
      stop("`surv` must be TRUE (columns `time` / `DSS`) or a single outcome column name.",
           call. = FALSE)
    outcome <- .sens_check_col(surv, data, "surv", n = 1L)
    if (!missing(time))
      stop("`time` only applies to survival outcomes (`surv = TRUE`).",
           call. = FALSE)
  }
  if (!is.numeric(depth) || !length(depth) || anyNA(depth) ||
      any(depth < 1) || any(depth != round(depth)))
    stop("`depth` must hold positive whole numbers; depth 0 is always a candidate.",
         call. = FALSE)
  depth <- sort(unique(as.integer(depth)))
  if (!is.numeric(split_frac) || length(split_frac) != 1L ||
      is.na(split_frac) || split_frac <= 0 || split_frac >= 1)
    stop("`split_frac` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  ra <- .merge_named_arg(rule_args, list(n_forest = 20L, num_trees = 200L,
                                         n_folds = 5L, min_leaf = 0.05,
                                         screen = TRUE, gate = 0.1),
                         "rule_args")
  ra$n_forest  <- .hte_select_count(ra$n_forest, "rule_args$n_forest", 1, 1e4)
  ra$num_trees <- .hte_select_count(ra$num_trees, "rule_args$num_trees", 2, 1e5)
  ra$n_folds   <- .hte_select_count(ra$n_folds, "rule_args$n_folds", 2, 100)
  if (!is.numeric(ra$min_leaf) || length(ra$min_leaf) != 1L ||
      is.na(ra$min_leaf) || ra$min_leaf <= 0 || ra$min_leaf >= 0.5)
    stop("`rule_args$min_leaf` must be a single number strictly between 0 and 0.5.",
         call. = FALSE)
  if (!isTRUE(ra$screen) && !isFALSE(ra$screen))
    stop("`rule_args$screen` must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(ra$gate) || length(ra$gate) != 1L || is.na(ra$gate) ||
      ra$gate <= 0 || ra$gate > 1)
    stop("`rule_args$gate` must be a single number in (0, 1].", call. = FALSE)
  seed <- .hte_select_count(seed, "seed", 0, .Machine$integer.max - 1e5)
  if (!is.list(grf_args))
    stop("`grf_args` must be a named list.", call. = FALSE)
  fixed <- intersect(names(grf_args), c("X", "Y", "W", "D", "horizon", "W.hat",
                                        "Y.hat", "sample.weights", "clusters"))
  if (length(fixed))
    stop(sprintf("`grf_args` cannot set %s in get_hte_icf(): the data columns set X, Y, W, D and horizon, and per-row fields cannot follow the split.",
                 paste0("`", fixed, "`", collapse = ", ")), call. = FALSE)

  # ---- Rows and design matrix -----------------------------------------------
  odd <- adj_var[!vapply(data[adj_var], function(x) is.numeric(x) ||
                           is.factor(x) || is.character(x) || is.logical(x),
                         logical(1L))]
  if (length(odd))
    stop(sprintf("Covariate column(s) %s must be numeric, factor, character or logical.",
                 paste0("`", odd, "`", collapse = ", ")), call. = FALSE)
  keep <- stats::complete.cases(data[c(cat_var, outcome)])
  data <- .sens_complete(data, c(cat_var, outcome), verbose)
  miss <- adj_var[vapply(data[adj_var], anyNA, logical(1L))]
  if (length(miss))
    stop(sprintf("get_hte_icf() needs complete covariates, since a rule cannot say where a missing value goes; %s have missing values. Impute them or drop those rows first.",
                 paste0("`", miss, "`", collapse = ", ")), call. = FALSE)
  W <- .psw_treat(data[[cat_var]], cat_var, arg = "cat_var")$z

  # The columns get_hte() builds, and what each means in a rule
  xdat <- data[adj_var]
  fac  <- adj_var[!vapply(xdat, is.numeric, logical(1L))]
  lv   <- lapply(xdat[fac], function(x) levels(droplevels(as.factor(x))))
  xdat[fac] <- lapply(xdat[fac], function(x) {
    x <- droplevels(as.factor(x))
    if (factor_encoding == "integer") as.integer(x)
    else if (nlevels(x) == 1L) rep(1, length(x)) else x
  })
  X    <- .sens_model_matrix(xdat, adj_var, one_hot = TRUE)
  src  <- adj_var[attr(X, "assign")]
  bin  <- vapply(xdat, function(x) is.numeric(x) && all(x %in% c(0, 1)),
                 logical(1L))
  type <- ifelse(src %in% fac,
                 if (factor_encoding == "integer") "code" else "level",
                 ifelse(bin[src], "bin", "num"))
  level <- rep(NA_character_, length(src))
  for (v in fac) if (factor_encoding == "onehot" && length(lv[[v]]) > 1L)
    level[src == v] <- lv[[v]]
  cols <- data.frame(var = src, type = type, level = level,
                     stringsAsFactors = FALSE)
  cols$type[cols$var %in% fac & vapply(lv[cols$var], length, 1L) == 1L] <- "num"
  attr(cols, "levels") <- c(lv, stats::setNames(rep(list(c("0", "1")),
                                                    sum(bin)), names(bin)[bin]))

  # ---- Split and folds ---------------------------------------------------------
  genv <- globalenv()
  old_seed <- if (exists(".Random.seed", envir = genv, inherits = FALSE))
    get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = genv)
    else if (exists(".Random.seed", envir = genv, inherits = FALSE))
      rm(".Random.seed", envir = genv)
  }, add = TRUE)
  set.seed(seed)
  disc <- logical(length(W))
  for (a in 0:1) {
    i <- which(W == a)
    disc[i[sample.int(length(i), round(split_frac * length(i)))]] <- TRUE
  }
  if (any(table(factor(W[disc], 0:1)) < 2L * ra$n_folds) ||
      any(table(factor(W[!disc], 0:1)) < 2L))
    stop("Too few patients per arm for this `split_frac` and `rule_args$n_folds`.",
         call. = FALSE)
  d_idx <- which(disc)
  e_idx <- which(!disc)
  fold  <- integer(length(d_idx))
  for (a in 0:1) {
    j <- which(W[d_idx] == a)
    fold[j] <- sample(rep_len(seq_len(ra$n_folds), length(j)))
  }

  # ---- Discovery runs ------------------------------------------------------------
  ga <- grf_args
  if (is.null(ga$seed)) ga$seed <- seed
  forest_args <- grf_args[setdiff(names(grf_args),
                                  c("num.trees", "seed", "target"))]
  hte <- function(rows) do.call(get_hte, c(
    list(data = rows, cat_var = cat_var, adj_var = adj_var, surv = surv,
         factor_encoding = factor_encoding, grf_args = ga),
    if (is_surv) list(time = time)))
  notes <- character()
  run <- function(idx) withCallingHandlers(
    .icf_discover(hte(data[idx, , drop = FALSE]), X[idx, , drop = FALSE],
                  cols, depth, ra, seed, forest_args),
    warning = function(w) {
      notes <<- c(notes, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))

  full <- run(d_idx)
  gD   <- full$g
  loss <- matrix(NA_real_, length(d_idx), length(depth) + 1L)
  for (k in seq_len(ra$n_folds)) {
    te   <- which(fold == k)
    tr   <- d_idx[-te]
    part <- run(tr)
    mu0  <- mean(part$g)
    loss[te, 1L] <- (gD[te] - mu0)^2
    for (j in seq_along(depth)) {
      tree <- part$parts[[j]]$tree
      mu <- if (nrow(tree)) {
        m <- tapply(part$g, .icf_assign(tree, X[tr, , drop = FALSE]), mean)
        m[.icf_assign(tree, X[d_idx[te], , drop = FALSE])]
      } else rep(mu0, length(te))
      mu[is.na(mu)] <- mu0
      loss[te, j + 1L] <- (gD[te] - mu)^2
    }
  }
  if (length(notes))
    warning(sprintf("The discovery fits raised %d warning(s); the first: %s",
                    length(notes), notes[1L]), call. = FALSE)

  n_leaf <- c(1L, vapply(full$parts, function(p) nrow(p$leaves), integer(1L)))
  cv_loss <- colMeans(loss)
  sel <- which.min(ifelse(n_leaf > 1L | seq_along(n_leaf) == 1L, cv_loss, Inf))
  # A deeper depth whose partition equals a shallower one's is reported as
  # the shallower depth.
  part_text <- c("", vapply(full$parts, function(p)
    paste(p$leaves$rule, collapse = " | "), character(1L)))
  sel <- match(part_text[sel], part_text)
  calib_p <- full$fit$calibration$p.value[2L]
  gated <- isTRUE(calib_p > ra$gate)
  if (gated) sel <- 1L
  cv <- tibble::tibble(depth = c(0L, depth), n_leaf = n_leaf,
                       cv_loss = cv_loss,
                       std.error = apply(loss, 2L, stats::sd) / sqrt(nrow(loss)),
                       selected = seq_along(n_leaf) == sel)
  vote <- tibble::tibble(
    depth = depth, n_leaf = n_leaf[-1L],
    share = vapply(full$parts, function(p) p$share, numeric(1L)),
    partition = part_text[-1L])
  chosen <- if (sel == 1L) list(tree = full$parts[[1L]]$tree[0, ],
                                leaves = .icf_leaves(full$parts[[1L]]$tree[0, ],
                                                     cols))
            else full$parts[[sel - 1L]]
  if (verbose)
    cli::cli_inform(c("i" = paste(
      "Discovery {length(d_idx)} and estimation {length(e_idx)} patients;",
      "screened {.field {full$screened}}; depth {cv$depth[sel]} selected.")))

  # ---- Estimation part -------------------------------------------------------------
  leaves <- chosen$leaves
  est_data <- data[e_idx, , drop = FALSE]
  est_data$.rule <- factor(
    leaves$rule[match(.icf_assign(chosen$tree, X[e_idx, , drop = FALSE]),
                      leaves$path)], levels = leaves$rule)
  est  <- hte(est_data)
  ea   <- attr(est, "analysis")
  sub  <- .hte_muffle_ps(.hte_subgroup(
    est$fit, .hte_arm_scores(est$fit), est$data, ".rule",
    data.frame(estimand = "ATE", measure = "diff", stringsAsFactors = FALSE),
    identical(ea$target, "survival.probability"), stats::qnorm(0.975),
    .hte_beyond(est)))
  m <- match(leaves$rule, sub$level)
  n_disc <- table(factor(.icf_assign(chosen$tree, X[d_idx, , drop = FALSE]),
                         levels = leaves$path))
  rules <- tibble::tibble(
    leaf = seq_len(nrow(leaves)), rule = leaves$rule,
    n_disc = as.integer(n_disc), n = sub$n[m], n_treat = sub$n_treat[m],
    estimate = sub$estimate[m], std.error = sub$std.error[m],
    conf.low = sub$conf.low[m], conf.high = sub$conf.high[m],
    p.value = sub$p.value[m], p_inter = sub$p_inter[m])
  rules$n[is.na(m)] <- 0L
  rules$n_treat[is.na(m)] <- 0L

  importance <- full$fit$importance
  importance$screened <- importance$variable %in% full$screened

  structure(
    list(rules = rules, cv = cv, vote = vote, importance = importance,
         est = est),
    class = c("hte_icf", "list"),
    analysis = list(
      backend_version = ea$backend_version, forest = ea$forest,
      outcome_type = ea$outcome_type, cat_var = cat_var, treated = ea$treated,
      surv = surv, outcome = outcome, time = ea$time, target = ea$target,
      adj_var = adj_var, factor_encoding = factor_encoding, depth = depth,
      depth_selected = cv$depth[sel], split_frac = split_frac,
      rule_args = ra, seed = seed, screened = full$screened,
      discovery = which(keep)[d_idx],
      calibration_p = calib_p, gated = gated,
      call = match.call()))
}


# ---- L3 print ------------------------------------------------------------------

#' @export
#' @noRd
print.hte_icf <- function(x, ...) {
  a  <- attr(x, "analysis")
  ra <- a$rule_args
  cat(sprintf("<hte_icf> iterative causal forest (%s, grf %s)\n",
              a$forest, a$backend_version))
  cat(sprintf("  discovery n = %d, estimation n = %d; %d-fold CV; %d forests x %d trees per run\n",
              sum(x$rules$n_disc), nrow(x$est$data), ra$n_folds, ra$n_forest,
              ra$num_trees))
  cat(sprintf("  screened: %s; discovery calibration p = %s\n",
              paste(a$screened, collapse = ", "),
              format(signif(a$calibration_p, 3))))
  sel <- a$depth_selected
  cat(sprintf("  selected depth = %d%s\n", sel,
              if (sel > 0L) sprintf(" (vote share %s)",
                                    format(x$vote$share[x$vote$depth == sel],
                                           digits = 2))
              else if (isTRUE(a$gated))
                sprintf(" (no subgroups: calibration p above gate = %s)",
                        format(ra$gate))
              else " (no subgroups)"))
  cat("\nCross-validated loss:\n")
  print(as.data.frame(x$cv), row.names = FALSE, digits = 4)
  cat("\nRules, estimated on the estimation part (ATE difference, 95% CI):\n")
  print(as.data.frame(x$rules[c("rule", "n_disc", "n", "estimate",
                                "conf.low", "conf.high", "p.value")]),
        row.names = FALSE, digits = 3)
  if (nrow(x$rules) > 1L)
    cat(sprintf("P for interaction = %s\n",
                format(signif(x$rules$p_inter[1L], 3))))
  invisible(x)
}
