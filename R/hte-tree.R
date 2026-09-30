# =============================================================================
# hte-tree.R -- one subgroup tree grown on doubly robust scores
# =============================================================================
#
# Architecture:
#
#   L1  get_hte_tree()   split the sample, score the discovery part with
#                        get_hte(), grow one tree by `method`, keep the splits
#                        whose children hold both arms, and estimate every
#                        node on the estimation part
#   L2  .tree_design()   the split variables as a design matrix, and what
#                        each column means in a rule (also get_hte_icf()'s)
#   L2  .tree_maxt()     the max-t test tree, one .maxt_node() per node
#   L2  .tree_engine()   ctree, mob, rpart and policy trees in the same layout
#   L2  .tree_prune()    drop the splits a child of which lacks an arm
#   L1  plt_hte_tree()   the tree drawn by plt_hte_icf()'s .hte_draw_tree()
#   L3  print.hte_tree()
#
# A tree has the layout get_hte_icf() uses: one row per split, with its path
# from the root ("", "L", "LR", ...), its column of the design matrix and its
# cut-off and whether the left side includes the cut-off. .icf_leaves(),
# .icf_assign() and .icf_party() read it. Older trees without `right` retain
# the historical `x <= value` convention.
# =============================================================================


# ---- L2 design -----------------------------------------------------------------

# The split variables as the design matrix get_hte() builds, with what each
# column means in a rule: a numeric covariate ("num"), a 0/1 covariate
# ("bin"), one level of a factor ("level", one-hot) or a factor's level codes
# ("code", integer).
#' @keywords internal
#' @noRd
.tree_design <- function(xdat, factor_encoding) {
  vars <- names(xdat)
  fac  <- vars[!vapply(xdat, is.numeric, logical(1L))]
  lv   <- lapply(xdat[fac], function(x) levels(droplevels(as.factor(x))))
  xdat[fac] <- lapply(xdat[fac], function(x) {
    x <- droplevels(as.factor(x))
    if (factor_encoding == "integer") as.integer(x)
    else if (nlevels(x) == 1L) rep(1, length(x)) else x
  })
  X    <- .sens_model_matrix(xdat, vars, one_hot = TRUE)
  src  <- vars[attr(X, "assign")]
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
  list(X = X, cols = cols)
}


# ---- L2 engines ----------------------------------------------------------------

# One node of the max-t tree. Every admissible cut of every column -- both
# sides at least `min_n` rows and `min_arm` of each arm -- gets the Welch t of
# the two sides' mean scores, and a variable's statistic is its largest |t|
# over its columns and cuts. The null distribution comes from a Rademacher
# multiplier bootstrap of the centred scores (a fixed-regressor bootstrap,
# Hansen 2000), each draw studentised as the data are; the node p-value is
# the Westfall-Young min-P over the variables. On heteroskedastic AIPW
# scores with a constant effect this held its size (5.8% at 5%) where ctree
# rejected in 12% and mob in 18% of data sets, and a Gaussian bootstrap
# divided by the observed standard errors was conservative (3.7%).
# `n_boot = 0` skips the bootstrap (p-value NA).
#' @keywords internal
#' @noRd
.maxt_node <- function(g, X, W, var, ctrl) {
  n  <- length(g)
  r  <- g - mean(g)
  tw <- sum(W)
  # |Welch t| of each cut from the left-side sums SL (cuts x draws), the
  # totals ST and the fixed sums of squares; Rademacher draws keep r^2
  tstat <- function(SL, ST, m, qm, qn) {
    mL <- SL / m
    mR <- (rep(ST, each = length(m)) - SL) / (n - m)
    vL <- (qm - m * mL^2) / (m - 1)
    vR <- (qn - qm - (n - m) * mR^2) / (n - m - 1)
    t  <- abs(mL - mR) / sqrt(pmax(vL, 0) / m + pmax(vR, 0) / (n - m))
    t[is.nan(t)] <- 0
    t
  }
  cand <- list()
  for (j in seq_len(ncol(X))) {
    o  <- order(X[, j])
    xs <- X[o, j]
    m  <- which(xs[-1L] > xs[-n])
    cw <- cumsum(W[o])[m]
    ok <- m >= ctrl$min_n & n - m >= ctrl$min_n &
      pmin(cw, m - cw, tw - cw, n - m - (tw - cw)) >= ctrl$min_arm
    m  <- m[ok]
    if (!length(m)) next
    if (length(m) > ctrl$max_cuts)
      m <- m[unique(round(seq(1, length(m), length.out = ctrl$max_cuts)))]
    ro  <- r[o]
    q   <- cumsum(ro^2)
    bin <- integer(n)
    bin[o] <- findInterval(seq_len(n), m + 1L) + 1L
    cand[[length(cand) + 1L]] <- list(
      j = j, var = var[j], m = m, value = xs[m], q = q[m], qn = q[n],
      bin = bin, t = drop(tstat(matrix(cumsum(ro)[m]), sum(ro), m, q[m], q[n])))
  }
  if (!length(cand)) return(NULL)

  cv   <- vapply(cand, `[[`, "", "var")
  vars <- unique(cv)
  tmax <- vapply(cand, function(cd) max(cd$t), numeric(1L))
  S    <- vapply(vars, function(v) max(tmax[cv == v]), numeric(1L))
  B    <- ctrl$n_boot
  p    <- NA_real_
  if (B > 0L) {
    Mb   <- matrix(0, B, length(vars))
    done <- 0L
    size <- max(1L, min(B, floor(5e6 / n)))
    while (done < B) {
      cb <- min(size, B - done)
      RE <- r * matrix(sample(c(-1, 1), n * cb, replace = TRUE), n, cb)
      ST <- colSums(RE)
      at <- done + seq_len(cb)
      for (cd in cand) {
        SL <- apply(rowsum(RE, cd$bin, reorder = TRUE), 2L, cumsum)
        tb <- tstat(SL[seq_along(cd$m), , drop = FALSE], ST, cd$m, cd$q, cd$qn)
        k  <- match(cd$var, vars)
        Mb[at, k] <- pmax(Mb[at, k], apply(tb, 2L, max))
      }
      done <- done + cb
    }
    p_var <- (1 + colSums(sweep(Mb, 2L, S, ">="))) / (B + 1)
    pb <- apply(Mb, 2L, function(b) (B + 1 - rank(b, ties.method = "min")) / B)
    p  <- (1 + sum(apply(matrix(pb, B), 1L, min) <= min(p_var))) / (B + 1)
    best <- order(p_var, -S)[1L]
  } else best <- which.max(S)
  i  <- which(cv == vars[best])
  cd <- cand[[i[which.max(tmax[i])]]]
  list(col = cd$j, value = cd$value[which.max(cd$t)], statistic = S[[best]],
       p = p)
}

# The max-t tree: a node splits at its best cut while the node test rejects
# at `alpha` (always with alpha = 1) and the depth allows. `stats` holds the
# test of every node tested, split or not.
#' @keywords internal
#' @noRd
.tree_maxt <- function(g, X, W, var, ctrl) {
  tree  <- data.frame(path = character(), col = integer(), value = numeric(),
                      right = logical(), stringsAsFactors = FALSE)
  stats <- data.frame(path = character(), statistic = numeric(),
                      split_p = numeric(), stringsAsFactors = FALSE)
  grow <- function(rows, path) {
    if (nchar(path) >= ctrl$max_depth) return(invisible())
    tst <- .maxt_node(g[rows], X[rows, , drop = FALSE], W[rows], var, ctrl)
    if (is.null(tst)) return(invisible())
    stats[nrow(stats) + 1L, ] <<- list(path, tst$statistic, tst$p)
    if (ctrl$alpha < 1 && !(tst$p <= ctrl$alpha)) return(invisible())
    tree[nrow(tree) + 1L, ] <<- list(path, tst$col, tst$value, TRUE)
    left <- X[rows, tst$col] <= tst$value
    grow(rows[left], paste0(path, "L"))
    grow(rows[!left], paste0(path, "R"))
  }
  grow(seq_along(g), "")
  list(tree = tree, stats = stats)
}

# The node model of the "_abs", "_rel" and "_aft" trees, and the R-learner
# equation of "mob_r" / "ctree_r" ("r"), as partykit::mob()
# calls a fitting function: the coefficients, the objective it minimises (the
# residual sum of squares, or minus the log-likelihood) and each row's score.
# "lm" serves a continuous outcome, a 0/1 outcome on the risk scale and
# survival pseudo-values; "logit" and "cox" work on the log-OR and log-HR
# scales. `x` holds the treatment and any adjustment columns. Only the
# scores of the columns `keep` are returned (all estimable ones by default).
# With `decorrelate`, as MOB tests a subset of parameters, the estimable
# scores are first whitened by the inverse square root of their outer
# product and the `keep` columns taken from the result; done here rather
# than by mob(parm = ), which inverts all the scores and fails once a split
# variable, also adjusted for, is constant in a child. ctree takes the raw
# scores of `keep`, as model4you's pmtree(parm = ) means to.
#' @keywords internal
#' @noRd
.tree_nodefit <- function(model, keep = NULL, decorrelate = FALSE)
  function(y, x, start = NULL, weights = NULL, offset = NULL, ...,
           estfun = FALSE, object = FALSE) {
    w <- if (length(weights)) weights else rep(1, NROW(x))
    parameter_names <- colnames(x)
    if (model == "aft") {
      # Drop columns aliased within this node, retaining their original slots
      # so `keep` still selects the treatment after an adjustment is constant.
      q <- qr(x * sqrt(w))
      active <- sort(q$pivot[seq_len(q$rank)])
      xx <- x[, active, drop = FALSE]
      m <- withCallingHandlers(
        survival::survreg(y ~ xx - 1, weights = w, dist = "weibull",
                          x = TRUE, y = TRUE),
        warning = function(w) stop(conditionMessage(w), call. = FALSE))
      if (any(!is.finite(stats::coef(m))) || !is.finite(m$scale) ||
          m$scale <= 0 || !is.finite(m$loglik[2L]))
        stop("The Weibull AFT node model did not converge.", call. = FALSE)
      parameter_names <- c(parameter_names, "Log(scale)")
      cf <- stats::setNames(rep(NA_real_, ncol(x) + 1L), parameter_names)
      cf[active] <- stats::coef(m)
      cf[length(cf)] <- log(m$scale)
      if (estfun) {
        r <- stats::residuals(m, type = "matrix")
        sc <- cbind(x * (w * r[, "dg"]), w * r[, "ds"])
      }
      obj <- -m$loglik[2L]
    } else if (model == "cox") {
      m  <- suppressWarnings(survival::coxph(y ~ x, weights = w))
      cf <- stats::coef(m)
      sc <- if (estfun) as.matrix(stats::residuals(m, type = "score"))
      obj <- -m$loglik[2L]
    } else if (model == "r") {
      # R-learner: `y` holds each row's numerator and the single column of
      # `x` its denominator, (W - e)(Y - m) and (W - e)^2 for grf's causal
      # forest; the effect solves sum(y - x * tau) = 0, the residual-on-
      # residual slope. Minus (sum y)^2 / sum x is the R-loss up to a term the
      # split does not change, so a split maximises grf's own criterion.
      m   <- NULL
      tau <- sum(w * y) / sum(w * x[, 1L])
      cf  <- stats::setNames(tau, parameter_names[1L])
      sc  <- if (estfun) cbind(w * (y - x[, 1L] * tau))
      obj <- -sum(w * y)^2 / sum(w * x[, 1L])
    } else {
      m <- suppressWarnings(if (model == "logit")
        stats::glm.fit(x, y, weights = w, family = stats::binomial())
        else stats::lm.wfit(x, y, w))
      cf <- m$coefficients
      r  <- y - m$fitted.values
      sc <- if (estfun) x * (w * r)
      obj <- if (model == "logit") m$deviance / 2 else sum(w * r^2)
    }
    ok <- which(!is.na(cf))
    if (estfun && decorrelate && !is.null(keep)) {
      s  <- sc[, ok, drop = FALSE]
      e  <- eigen(crossprod(s) / nrow(s), symmetric = TRUE)
      k  <- e$values > max(e$values) * 1e-10
      sc <- s %*% (e$vectors[, k, drop = FALSE] %*%
                     (t(e$vectors[, k, drop = FALSE]) / sqrt(e$values[k])))
      colnames(sc) <- parameter_names[ok]
      keep <- match(parameter_names[keep], colnames(sc))
    }
    list(coefficients = cf, objfun = obj,
         estfun = if (estfun) sc[, if (is.null(keep)) ok else keep, drop = FALSE],
         object = if (object) m)
  }

# The other engines on the same scores and design columns, read back into
# the tree layout. A partykit split with right = FALSE (rpart's `x < c`) is
# written as `x <= ` the largest value going left, which sends every
# discovery row the same way; ctree and mob keep their split p-values.
# "mob_model" and "ctree_model" split on the scores of the node model in
# `node` (outcome `y`, regressors `R`, `model`, and `parm`, the columns
# tested, NULL for all) instead of on `g`: MOB, and ctree with the scores as
# its transformed response, as model4you's pmtree() does. `weights` are
# rpart's case weights (the R-learner's denominators).
#' @keywords internal
#' @noRd
.tree_engine <- function(method, g, X, ctrl, node = NULL, weights = NULL) {
  tree  <- data.frame(path = character(), col = integer(), value = numeric(),
                      right = logical(), stringsAsFactors = FALSE)
  stats <- data.frame(path = character(), statistic = numeric(),
                      split_p = numeric(), stringsAsFactors = FALSE)
  if (method == "policy") {
    sgn  <- if (ctrl$better == "higher") 1 else -1
    G    <- cbind(0, sgn * g - ctrl$cost)
    step <- if (is.null(ctrl$split_step)) max(1L, length(g) %/% 1000L)
            else ctrl$split_step
    Xm   <- unname(X)
    pt   <- if (ctrl$max_depth <= 2)
      policytree::policy_tree(Xm, G, depth = ctrl$max_depth, split.step = step,
                              min.node.size = ctrl$min_n, verbose = FALSE)
    else policytree::hybrid_policy_tree(Xm, G, depth = ctrl$max_depth,
                                        search.depth = 2, split.step = step,
                                        min.node.size = ctrl$min_n,
                                        verbose = FALSE)
    walk <- function(k, path) {
      nd <- pt$nodes[[k]]
      if (isTRUE(nd$is_leaf)) return(invisible())
      tree[nrow(tree) + 1L, ] <<- list(path, as.integer(nd$split_variable),
                                       nd$split_value, TRUE)
      walk(nd$left_child, paste0(path, "L"))
      walk(nd$right_child, paste0(path, "R"))
    }
    walk(1L, "")
    return(list(tree = tree, stats = stats))
  }

  xd <- as.data.frame(unname(X))
  names(xd) <- paste0("x", seq_len(ncol(X)))
  # MOB refits the node model at every distinct value of a split variable.
  # With `max_cuts`, only that many cut-points, evenly spaced among the
  # distinct values as in "maxt", are kept, and every value is moved up to
  # the next one: `x <= cut` then sends the same rows in either scale.
  if (method %in% c("mob", "mob_model") && !is.null(ctrl$max_cuts))
    xd[] <- lapply(xd, function(x) {
      cuts <- sort(unique(x))
      cuts <- cuts[-length(cuts)]
      if (length(cuts) <= ctrl$max_cuts) return(x)
      cuts <- cuts[unique(round(seq(1, length(cuts), length.out = ctrl$max_cuts)))]
      c(cuts, max(x))[findInterval(x, cuts, left.open = TRUE) + 1L]
    })
  d  <- if (!is.null(g)) cbind(.g = g, xd)
  fit <- switch(method,
    mob_model = {
      rn <- paste0(".r", seq_len(ncol(node$R)))
      dm <- data.frame(unname(node$R), xd)
      names(dm)[seq_along(rn)] <- rn
      dm$.y <- node$y
      partykit::mob(stats::as.formula(paste(
        ".y ~ 0 +", paste(rn, collapse = " + "), "|",
        paste(names(xd), collapse = " + "))), data = dm,
        fit = .tree_nodefit(node$model, node$parm, decorrelate = TRUE),
        control = partykit::mob_control(
          alpha = ctrl$alpha, maxdepth = ctrl$max_depth + 1,
          minsize = ctrl$min_n, trim = ctrl$trim))
    },
    ctree_model = {
      fitfun <- .tree_nodefit(node$model, node$parm)
      ytrafo <- function(data, weights, control, ...)
        function(subset, weights, info = NULL, estfun = TRUE, object = FALSE) {
          f <- tryCatch(fitfun(node$y[subset], node$R[subset, , drop = FALSE],
                               estfun = TRUE), error = function(e) NULL)
          if (is.null(f) || anyNA(f$coefficients[node$parm]))
            return(list(converged = FALSE))
          ef <- matrix(0, nrow(node$R), ncol(f$estfun))
          ef[subset, ] <- f$estfun
          list(estfun = ef, converged = TRUE)
        }
      partykit::ctree(.y ~ ., data = cbind(.y = seq_len(nrow(xd)), xd),
                      ytrafo = ytrafo, control = partykit::ctree_control(
                        alpha = ctrl$alpha, testtype = ctrl$testtype,
                        maxdepth = ctrl$max_depth, minbucket = ctrl$min_n,
                        minsplit = 2L * ctrl$min_n))
    },
    ctree = partykit::ctree(.g ~ ., data = d, control = partykit::ctree_control(
      alpha = ctrl$alpha, testtype = ctrl$testtype, maxdepth = ctrl$max_depth,
      minbucket = ctrl$min_n, minsplit = 2L * ctrl$min_n)),
    # lmtree counts the root as depth 1
    mob = partykit::lmtree(.g ~ 1 | ., data = d, alpha = ctrl$alpha,
                           maxdepth = ctrl$max_depth + 1, minsize = ctrl$min_n,
                           trim = ctrl$trim),
    rpart = {
      rw <- if (is.null(weights)) rep(1, nrow(d)) else weights
      rp <- rpart::rpart(.g ~ ., data = d, weights = rw, method = "anova",
                         model = TRUE, control = rpart::rpart.control(
                           maxdepth = min(ctrl$max_depth, 30), cp = 0,
                           minbucket = ctrl$min_n, minsplit = 2L * ctrl$min_n,
                           xval = ctrl$xval, maxcompete = 0L,
                           maxsurrogate = 0L))
      cpt <- rp$cptable
      if (ctrl$xval > 0L && nrow(cpt) > 1L) {
        best <- which.min(cpt[, "xerror"])
        pick <- if (ctrl$cp_rule == "min") best else
          which(cpt[, "xerror"] <= cpt[best, "xerror"] + cpt[best, "xstd"])[1L]
        rp <- rpart::prune(rp, cp = cpt[pick, "CP"])
      }
      partykit::as.party(rp)
    })
  nm <- names(fit$data)
  walk <- function(nd, rows, path) {
    p <- partykit::info_node(nd)$p.value
    if (length(p) == 1L)
      stats[nrow(stats) + 1L, ] <<- list(path, NA_real_, as.numeric(p))
    if (partykit::is.terminal(nd)) return(invisible())
    s  <- partykit::split_node(nd)
    j  <- match(nm[partykit::varid_split(s)], names(xd))
    x  <- X[rows, j]
    b  <- partykit::breaks_split(s)
    le <- if (isFALSE(partykit::right_split(s))) x < b else x <= b
    # Keep the engine's cut-off and strictness. Replacing it with the largest
    # observed value on the left changes the partition for held-out values in
    # a gap (for example, rpart's x < 5 becomes x <= 0 when the node has only
    # x = 0 and x = 10).
    tree[nrow(tree) + 1L, ] <<- list(path, j, b,
                                      isTRUE(partykit::right_split(s)))
    kids <- partykit::kids_node(nd)
    ix   <- partykit::index_split(s)
    lo   <- if (is.null(ix) || ix[1L] == 1L) 1L else 2L
    walk(kids[[lo]], rows[le], paste0(path, "L"))
    walk(kids[[3L - lo]], rows[!le], paste0(path, "R"))
  }
  walk(partykit::node_party(fit), seq_len(nrow(X)), "")
  list(tree = tree, stats = stats)
}

# Drops every split (with its subtree) a child of which holds fewer than
# `min_arm` rows of either arm; the max-t tree never makes one, the other
# engines do not look at the arms.
#' @keywords internal
#' @noRd
.tree_prune <- function(tree, X, W, min_arm) {
  if (!nrow(tree)) return(tree)
  leaf <- .icf_assign(tree, X)
  bad  <- character()
  for (p in tree$path) {
    if (any(startsWith(p, bad))) next
    for (side in paste0(p, c("L", "R"))) {
      w <- W[startsWith(leaf, side)]
      if (sum(w == 1) < min_arm || sum(w == 0) < min_arm) {
        bad <- c(bad, p)
        break
      }
    }
  }
  tree[!vapply(tree$path, function(p) any(startsWith(p, bad)), logical(1L)), ,
       drop = FALSE]
}


# ---- L1 public function ----------------------------------------------------------

#' A subgroup tree grown on doubly robust scores
#'
#' Grows one tree of subgroups with different treatment effects, by one of
#' several recursive-partitioning methods, and estimates every node honestly.
#' The data are split in two. On the discovery part, [get_hte()] gives each
#' patient an AIPW score -- a doubly robust, confounding-adjusted effect
#' whose mean in any subgroup is that subgroup's effect -- and the tree
#' partitions the scores (the `_cate` methods: the forest's predicted
#' effects; the `_r` methods: the forest's residuals, as grf's own trees do;
#' the `_abs`, `_rel` and `_aft` methods: the outcome, through a model
#' refitted in every node). On the estimation part, which played no role in
#' finding the tree, [get_hte()] estimates the effect of every node. The
#' tree is returned as a partykit `party` object for ggparty. Arguments
#' follow [get_hte_icf()] where the two overlap.
#'
#' @param data A data frame holding every column named below.
#' @param cat_var Length-1 character. The binary exposure column, coded as in
#'   [get_hte()].
#' @param adj_var Character vector of covariates the forests condition on.
#'   Without `candidate_var` they are also the split variables. Rows missing
#'   `cat_var` or the outcome are dropped, as in [get_hte()].
#' @param candidate_var `NULL` (default), or a character vector of the only
#'   variables the tree may split on. They join `adj_var`, so a candidate is
#'   always adjusted for. Split variables must be numeric, factor, character
#'   or logical with no missing value, since a rule cannot say where a
#'   missing value goes.
#' @param surv Outcome selector, as in [get_hte()]: `TRUE` (default) for the
#'   survival columns `time` and `DSS`, or a single binary (0/1) or continuous
#'   outcome column. The scores, and so the splits (but those of the `_rel` or `_aft`
#'   methods) and every estimate, are on the difference scale: mean, risk, or
#'   \eqn{S(t)} difference (RMST difference with
#'   `grf_args = list(target = "RMST")`).
#' @param time Time point of a survival outcome, as in [get_hte()]. Default
#'   `120`. Only accepted with `surv = TRUE`.
#' @param method How the tree is grown:
#'   \describe{
#'     \item{`"maxt"`}{(default) a heteroskedasticity-robust test at every
#'       node: the Welch t of the two sides' mean scores, maximised over the
#'       cuts of each variable, with a multiplier-bootstrap null and a
#'       Westfall-Young min-P over the variables. A node splits at its best
#'       cut while the test rejects at `alpha`. See Details.}
#'     \item{`"mob_dr"`}{[partykit::lmtree()] of the scores on an intercept:
#'       parameter-instability tests with a Bonferroni correction.}
#'     \item{`"mob_cate"`}{the same on the forest's out-of-bag CATE
#'       predictions (`$data$.cate` of [get_hte()]) instead of the scores --
#'       the two-stage approach of StratifiedMedicine's `ctree_cate` and of
#'       Virtual Twins. Its p-values are optimistic; see Details.}
#'     \item{`"mob_abs"`}{model-based recursive partitioning (MOB) of a node
#'       model of the outcome, the treatment plus `adj_var`, refitted in
#'       every node, on the difference scale: [stats::lm()] of a continuous
#'       or 0/1 outcome (mean or risk difference), or for survival of the
#'       Kaplan-Meier pseudo-values of \eqn{S(t)} (RMST with
#'       `grf_args = list(target = "RMST")`) at `time`. A node splits while
#'       the treatment coefficient is unstable across a split variable, as in
#'       partykit's `lmtree()` and StratifiedMedicine's `lmtree`. No scores
#'       are used, so the discovery part needs no forest.}
#'     \item{`"mob_rel"`}{the same on the ratio scale: logistic regression
#'       (log-OR) of a binary outcome or Cox regression (log-HR) of a
#'       survival outcome, as partykit's `glmtree()` and model4you. Not for a
#'       continuous outcome.}
#'     \item{`"ctree_dr"`}{[partykit::ctree()] on the scores: permutation
#'       tests with a Bonferroni correction over the design columns.}
#'     \item{`"ctree_cate"`}{the same on the CATE predictions, as
#'       `"mob_cate"`.}
#'     \item{`"ctree_abs"`, `"ctree_rel"`}{the node models of `"mob_abs"` and
#'       `"mob_rel"`, with ctree's permutation tests on the score of the
#'       treatment coefficient, as model4you's `pmtree()`.}
#'     \item{`"rpart_dr"`}{[rpart::rpart()] on the scores, grown to `max_depth`
#'       and pruned by cross-validation; no tests.}
#'     \item{`"rpart_cate"`}{the same CART fit and pruning on the forest's
#'       out-of-bag CATE predictions; an explanatory approximation of the
#'       forest, with no split tests.}
#'     \item{`"mob_aft"`, `"ctree_aft"`}{Weibull accelerated failure time
#'       node models via [survival::survreg()], using MOB instability tests
#'       or ctree permutation tests of the treatment score. Only for
#'       `surv = TRUE`, with strictly positive times and right censoring.
#'       Splits target the log-time ratio; reported effects remain DR
#'       survival-probability or RMST differences.}
#'     \item{`"mob_r"`, `"ctree_r"`, `"rpart_r"`}{the R-learner: the tree
#'       splits on the forest's out-of-bag residuals \eqn{\tilde Y = Y -
#'       \hat m(X)} and \eqn{\tilde W = W - \hat e(X)}, minimising the R-loss
#'       \eqn{\sum (\tilde Y - \tau \tilde W)^2} as grf's own trees and the
#'       iCF do. `"mob_r"` and `"ctree_r"` test the score
#'       \eqn{\tilde W (\tilde Y - \tilde W \tau)} of the node's
#'       residual-on-residual slope; `"rpart_r"` grows CART on
#'       \eqn{\tilde Y / \tilde W} weighted by \eqn{\tilde W^2}, pruned by
#'       cross-validation. Survival outcomes use the censoring-adjusted
#'       numerator and denominator of [grf::causal_survival_forest()]. See
#'       Details.}
#'     \item{`"policy"`}{[policytree::policy_tree()] (depth up to 2) or
#'       [policytree::hybrid_policy_tree()] (deeper): the tree of exactly
#'       `max_depth` levels whose treat-or-not choice per leaf maximises the
#'       summed scores; `$rules` gains the recommended `action`.}
#'   }
#' @param max_depth Largest number of split levels, a positive whole number
#'   or `Inf`. Default `3`, up to 8 subgroups. `"policy"` needs a finite
#'   depth.
#' @param alpha Significance level of the split tests of `"maxt"` and the
#'   mob and ctree methods, in (0, 1]. Default `0.05`. `1` splits every node the
#'   depth and leaf sizes allow, so the tree grows to `max_depth`. Not
#'   accepted by the rpart methods, which prune by cross-validation, or
#'   `"policy"`.
#' @param min_leaf Smallest leaf, as a share of the discovery patients, in
#'   [0, 0.5). Default `0.05`. Every child also needs two patients of each
#'   arm; a split of another method than `"maxt"` without them is dropped.
#' @param factor_encoding How factor, character and logical split variables
#'   enter the design, as in [get_hte()]: `"integer"` (default) cuts the
#'   level codes in the order of the levels; `"onehot"` lets a split take one
#'   level off the rest. Also passed to the forests. The rules are written
#'   in the original levels either way.
#' @param split_frac Share of the patients, within each arm, used for
#'   discovery. Default `0.5`; the rest are the estimation part. `1` grows and
#'   estimates the tree on the same patients, so the intervals are no longer
#'   honest.
#' @param estimator How the estimation part estimates the effect of every
#'   node: `"aipw"` (default), grf's augmented inverse-propensity weighting,
#'   the mean of the doubly robust scores in the node; or `"tmle"`, grf's
#'   targeted maximum likelihood estimation
#'   ([grf::average_treatment_effect()] with `method = "TMLE"`), which fits
#'   the correction as a regression of each arm's residuals on its inverse
#'   propensity instead of averaging the weighted residuals. `"tmle"` needs a
#'   continuous or binary outcome, since grf's causal survival forest
#'   estimates by AIPW only. The tree itself does not change, and `$est`
#'   stays the AIPW [get_hte()] result.
#' @param tree_args Named list of settings of the chosen `method`; partial
#'   overrides keep the other defaults, and a field of another method is an
#'   error.
#'   \describe{
#'     \item{`"maxt"`}{`n_boot`, bootstrap draws per node test, default
#'       `1000` (`0` gives no p-values and needs `alpha = 1`); `max_cuts`,
#'       the most cuts tried per design column, evenly spaced among the
#'       admissible ones, default `100`.}
#'     \item{the ctree methods}{`testtype`, `"Bonferroni"` (default),
#'       `"Univariate"` or `"MonteCarlo"`, as in
#'       [partykit::ctree_control()].}
#'     \item{the `_abs`, `_rel` and `_aft` methods}{also `adjust`, `TRUE` (default)
#'       for a node model adjusted for `adj_var` or `FALSE` for the outcome
#'       on the treatment alone; and `parm`, `"treatment"` (default) to test
#'       the treatment coefficient only or `"all"` to test every coefficient
#'       (including log-scale for AFT).
#'       `adjust = FALSE, parm = "all"` is the default of partykit,
#'       StratifiedMedicine and model4you.}
#'     \item{the mob methods}{`trim`, the share of observations
#'       trimmed from the
#'       ends of a numeric split variable in the instability tests, default
#'       `0.1`; `max_cuts`, the most cut-points tried per design column,
#'       evenly spaced among its distinct values, which are moved up to the
#'       next cut-point before MOB sees them (the tests then see ties),
#'       default `100`. MOB refits the node model at every cut-point, so
#'       with every distinct value tried its time grows with the square of
#'       the rows; a value above the number of distinct values, such as
#'       `1e6`, tries them all.}
#'     \item{the rpart methods}{`xval`, cross-validation folds, default `10` (`0`
#'       keeps the tree grown to `max_depth` unpruned); `cp_rule`, `"min"`
#'       (default) prunes at the smallest cross-validated error, `"1se"` to
#'       the smallest tree within one standard error of it.}
#'     \item{`"policy"`}{`cost`, the effect a treated patient must exceed,
#'       on the score scale, default `0`; `better`, `"higher"` (default) when
#'       a larger outcome is better -- survival, RMST -- or `"lower"`, as for
#'       an adverse binary event; `split_step`, cut-points skipped between
#'       two tried (policytree's `split.step`), default `NULL`: 1 per 1000
#'       discovery patients, since the exact search grows with the square of
#'       the distinct values.}
#'   }
#' @param grf_args Named list forwarded to both [get_hte()] fits, as in
#'   [get_hte_icf()]. Per-row fields (`W.hat`, `Y.hat`, `sample.weights`,
#'   `clusters`) are not accepted, because each fit sees a different part of
#'   the rows; a single known propensity `W.hat`, as in a trial, is.
#'   Without `num.trees`, each fit grows 500 trees, or 200 above 10,000 rows
#'   as in [get_hte()]: grf's default 2000 gave the same trees and leaf
#'   effects within 0.005 in simulations, at up to four times the time.
#'   For a survival outcome with more than 100 distinct times up to `time`,
#'   `failure.times` defaults to 100 evenly spaced points from 0 to `time`
#'   and the first follow-up past it, where [get_hte()] moves the patients
#'   followed longer: in 100 simulated data sets of 4,000 patients the
#'   effects moved by at most 0.002, with the same standard errors, at a
#'   third of the time. A grid ending at `time` would put those patients on
#'   the horizon instead (effects off by up to 0.03, standard errors 19%
#'   larger).
#' @param seed Nonnegative whole number, default `123`. It draws the split,
#'   the bootstrap multipliers and rpart's folds, and seeds the [get_hte()]
#'   fits unless `grf_args` sets `seed`. The caller's random-number state is
#'   restored, also on error.
#' @param verbose Logical. `TRUE` reports the split of the sample, the
#'   scores left out and the splits dropped for a missing arm. Default
#'   `FALSE`.
#'
#' @details
#' Growing on AIPW scores rather than on predicted effects keeps the tree
#' honest about confounding: their conditional mean is the conditional effect
#' whichever of the propensity and outcome models is right, while a tree on
#' forest predictions inherits their bias and splits on confounders. Scores
#' at a propensity of 0 or 1 are undefined and left out of the growing.
#' The `_cate` methods grow on the predictions all the same, as two-stage
#' methods do, for comparison: the predictions are smooth and far less noisy
#' than the scores, but their errors follow the covariates -- where one arm
#' is rare the forest extrapolates. The `"mob_cate"` and `"ctree_cate"`
#' tests treat them as independent observations, so their p-values are
#' optimistic. `"rpart_cate"` reports no split p-values.
#'
#' The `_r` methods grow on the same out-of-bag residuals the forest was
#' fitted with. Within a node the effect is the residual-on-residual slope
#' \eqn{\sum \tilde W \tilde Y / \sum \tilde W^2} (Robinson's transformation,
#' the R-learner of Nie and Wager 2021), and a split is chosen to maximise
#' \eqn{\sum_c (\sum_{i \in c} \tilde W_i \tilde Y_i)^2 / \sum_{i \in c}
#' \tilde W_i^2}, the criterion grf's causal trees approximate. The
#' residuals involve no inverse propensity, so extreme propensities inflate
#' them far less than the AIPW scores. The slope weights each patient by
#' \eqn{\tilde W^2}, whose mean is \eqn{e(1 - e)}: the splits follow
#' overlap-weighted effects, while the reported effects stay the doubly
#' robust ATE differences of the estimation part.
#'
#' The `_abs`, `_rel` and `_aft` methods use neither: every node fits a model of the
#' outcome on the treatment and `adj_var` and tests the score of the
#' treatment coefficient (MOB first whitens the scores by their outer
#' product, as `mob(parm = )` does). A purely prognostic variable moves the
#' intercept, not the treatment coefficient, so testing that coefficient
#' and adjusting for strong prognostic and confounding covariates keeps
#' these trees on effect modifiers. With `tree_args = list(adjust = FALSE,
#' parm = "all")`, the packages' default, they split on prognostic variables
#' too. Their validity rests on the node model -- linear adjustment, an
#' effect constant within a node, proportional hazards for Cox, and for the
#' pseudo-values censoring independent of the covariates -- and a node
#' model needs about 10 patients per coefficient, which sets a floor on the
#' leaves. The `_rel` trees split on the log-OR or log-HR scale while every
#' reported effect is a difference; effect modification depends on the
#' scale, so the two can disagree.
#' The `_aft` methods assume a Weibull AFT model and conditionally independent
#' censoring. The estimated log-scale is a nuisance parameter in the MOB
#' score whitening and is also tested when `parm = "all"`. AFT splits need
#' not agree with heterogeneity on the final difference scale. No discovery
#' forest is needed unless `split_frac = 1`.
#'
#' The scores are heteroskedastic -- their variance grows where the
#' propensity nears 0 or 1 and with the outcome's variance -- which the
#' permutation tests of ctree and the instability tests of MOB assume away.
#' `"maxt"` compares, for every admissible cut, the mean scores
#' on its two sides with a Welch t, and takes the largest |t| of each
#' variable. Its null distribution is simulated by multiplying the centred
#' scores by random signs (a Rademacher multiplier bootstrap, Hansen 2000),
#' which keeps each patient's variance, and studentising each draw as the
#' data are; the node p-value is the Westfall-Young min-P over the
#' variables. The node splits at the best cut of the variable with the
#' smallest p-value. As in ctree and MOB, the tests are not corrected across
#' nodes; the root test is a test of any heterogeneity at all.
#'
#' Depth counts split levels: `max_depth = 1` is one split. Leaves hold at
#' least `min_leaf` of the discovery patients and two of each arm.
#'
#' @return An object of class `hte_tree`: a list of
#'   \describe{
#'     \item{`rules`}{Tibble with one row per leaf: `leaf`, `node` (its id in
#'       `tree`), `rule`, `n_disc` (discovery patients), and from the
#'       estimation part `n`, `n_treat`, `estimate`, `std.error`, `conf.low`,
#'       `conf.high`, `p.value` and `p_inter`, the doubly robust ATE
#'       difference with 95% Wald intervals, as in [get_hte()]'s `$subgroup`
#'       (by TMLE with `estimator = "tmle"`).
#'       `method = "policy"` adds `action`, `"Treated"` or `"Control"`.}
#'     \item{`nodes`}{Tibble with one row per node, in `tree`'s order: `node`,
#'       `parent`, `depth`, `terminal`, `rule` (leaves), `variable` and
#'       `split` (the condition sending patients left) of inner nodes,
#'       `statistic` and `split_p` of the node's split test on the discovery
#'       part (`"maxt"`; `split_p` also for the mob and ctree methods; for a leaf
#'       the test that did not reject), `n_disc`, and the estimation part's
#'       `n`, `n_treat`, `estimate`, `std.error`, `conf.low`, `conf.high` and
#'       `p.value` of the node.}
#'     \item{`tree`}{The tree as a [partykit::party()] on the estimation
#'       patients: the design columns, then the outcome, `.arm`, `.dr_score`
#'       and `.rule`. Each node's `info` is its row of `nodes`, so
#'       `ggparty::ggparty(res$tree, add_vars = list(est =
#'       "$node$info$estimate"))` maps the effects.}
#'     \item{`est`}{The `hte_res` from [get_hte()] on the estimation part,
#'       whose `$data` has the rule of every patient in the factor `.rule`:
#'       `plt_hte_sub(res$est, sub_var = ".rule")` draws the leaves.}
#'   }
#'   Analysis metadata is attached as `attr(x, "analysis")`, including the
#'   settings, the row numbers of `data` in the discovery part
#'   (`discovery`) and what each design column means (`cols`).
#'
#' @references
#' Hansen BE (2000). Testing for structural change in conditional models.
#' \emph{Journal of Econometrics} 97(1):93-115.
#'
#' Westfall PH, Young SS (1993). \emph{Resampling-Based Multiple Testing}.
#' Wiley.
#'
#' Hothorn T, Hornik K, Zeileis A (2006). Unbiased recursive partitioning: a
#' conditional inference framework. \emph{Journal of Computational and
#' Graphical Statistics} 15(3):651-674.
#'
#' Athey S, Wager S (2021). Policy learning with observational data.
#' \emph{Econometrica} 89(1):133-161.
#'
#' Nie X, Wager S (2021). Quasi-oracle estimation of heterogeneous treatment
#' effects. \emph{Biometrika} 108(2):299-319.
#'
#' Athey S, Tibshirani J, Wager S (2019). Generalized random forests.
#' \emph{Annals of Statistics} 47(2):1148-1178.
#'
#' @seealso [plt_hte_tree()] to draw the tree; [get_hte_icf()] finds rules
#'   by voting over causal-forest trees; [get_hte()] for the forest and the
#'   scores; [plt_hte_sub()] to draw the leaves as a forest plot.
#'
#' @examplesIf requireNamespace("grf", quietly = TRUE) && requireNamespace("partykit", quietly = TRUE)
#' \donttest{
#' set.seed(20260929)
#' n <- 1600
#' d <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rnorm(n), X3 = runif(n),
#'                 grp = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
#' d$z <- rbinom(n, 1, plogis(0.5 * d$X2))
#' d$y <- d$X2 + d$z * (0.5 + 1.5 * (d$X3 > 0.6)) + rnorm(n)
#'
#' res <- get_hte_tree(d, cat_var = "z", adj_var = c("X1", "X2", "X3", "grp"),
#'                     surv = "y", tree_args = list(n_boot = 500))
#' res
#' if (requireNamespace("ggparty", quietly = TRUE)) plt_hte_tree(res)
#'
#' # the same scores, other methods
#' get_hte_tree(d, "z", c("X1", "X2", "X3", "grp"), surv = "y",
#'              method = "rpart_dr")$rules
#' }
#'
#' @export
get_hte_tree <- function(data,
                         cat_var,
                         adj_var,
                         candidate_var = NULL,
                         surv       = TRUE,
                         time       = 120,
                         method     = c("maxt", "mob_dr", "mob_cate", "mob_abs",
                                        "mob_rel", "ctree_dr", "ctree_cate",
                                        "ctree_abs", "ctree_rel", "rpart_dr",
                                        "policy", "rpart_cate", "mob_aft",
                                        "ctree_aft", "mob_r", "ctree_r",
                                        "rpart_r"),
                         max_depth  = 3,
                         alpha      = 0.05,
                         min_leaf   = 0.05,
                         factor_encoding = c("integer", "onehot"),
                         split_frac = 0.5,
                         estimator  = c("aipw", "tmle"),
                         tree_args  = list(),
                         grf_args   = list(),
                         seed       = 123,
                         verbose    = FALSE) {

  method <- match.arg(method)
  # The mob and ctree methods grow on the AIPW scores (_dr), on the CATE
  # predictions (_cate), on the forest's residuals (_r, the R-learner) or on
  # the scores of a node model of the outcome on the difference (_abs),
  # ratio (_rel) or accelerated-time (_aft) scale
  engine <- sub("_(dr|cate|abs|rel|aft|r)$", "", method)
  node_scale <- if (grepl("_(abs|rel|aft)$", method)) sub("^.*_", "", method)
  model_based <- !is.null(node_scale)
  rlearner <- endsWith(method, "_r")
  factor_encoding <- match.arg(factor_encoding)
  estimator <- match.arg(estimator)
  for (pkg in c("grf", "partykit",
                switch(engine, rpart = "rpart", policy = "policytree")))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required for get_hte_tree(method = \"%s\").",
                   pkg, method), call. = FALSE)
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
    y0 <- data[[outcome]][!is.na(data[[outcome]])]
    if (identical(node_scale, "rel") && !all(y0 %in% c(0, 1)))
      stop(sprintf("`method = \"%s\"` splits on the log-OR or log-HR scale, which a continuous outcome lacks; use `method = \"%s\"`.",
                   method, sub("_rel$", "_abs", method)), call. = FALSE)
  }
  if (is_surv && estimator == "tmle")
    stop("`estimator = \"tmle\"` needs a continuous or binary outcome: grf's causal survival forest estimates by AIPW only.",
         call. = FALSE)
  if (identical(node_scale, "aft")) {
    if (!is_surv)
      stop("AFT methods require `surv = TRUE` (columns `time` / `DSS`).",
           call. = FALSE)
    tt <- data$time[!is.na(data$time)]
    if (!is.numeric(tt) || any(!is.finite(tt) | tt <= 0))
      stop("AFT methods require finite, strictly positive survival times.",
           call. = FALSE)
  }
  if (model_based && is_surv && !requireNamespace("survival", quietly = TRUE))
    stop(sprintf("Package 'survival' is required for get_hte_tree(method = \"%s\") with a survival outcome.",
                 method), call. = FALSE)
  split_var <- adj_var
  if (!is.null(candidate_var)) {
    candidate_var <- unique(.sens_check_col(candidate_var, data, "candidate_var"))
    bad <- intersect(candidate_var, c(cat_var, outcome))
    if (length(bad))
      stop(sprintf("`candidate_var` cannot hold the exposure or outcome column %s.",
                   paste0("`", bad, "`", collapse = ", ")), call. = FALSE)
    adj_var   <- union(adj_var, candidate_var)
    split_var <- candidate_var
  }
  if (!is.numeric(max_depth) || length(max_depth) != 1L || is.na(max_depth) ||
      max_depth < 1 || (is.finite(max_depth) && max_depth != round(max_depth)))
    stop("`max_depth` must be a positive whole number or Inf.", call. = FALSE)
  if (method == "policy" && !is.finite(max_depth))
    stop("`method = \"policy\"` grows exactly `max_depth` levels, so `max_depth` must be finite.",
         call. = FALSE)
  tested <- !engine %in% c("rpart", "policy")
  if (!tested && !missing(alpha))
    stop(sprintf("`alpha` does not apply to `method = \"%s\"`, which %s.", method,
                 if (engine == "rpart") "prunes by cross-validation"
                 else "grows to `max_depth`"), call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha > 1)
    stop("`alpha` must be a single number in (0, 1].", call. = FALSE)
  if (!is.numeric(min_leaf) || length(min_leaf) != 1L || is.na(min_leaf) ||
      min_leaf < 0 || min_leaf >= 0.5)
    stop("`min_leaf` must be a single number in [0, 0.5).", call. = FALSE)
  if (!is.numeric(split_frac) || length(split_frac) != 1L ||
      is.na(split_frac) || split_frac <= 0 || split_frac > 1)
    stop("`split_frac` must be a single number in (0, 1].", call. = FALSE)
  ta <- .merge_named_arg(tree_args, switch(engine,
    maxt   = list(n_boot = 1000L, max_cuts = 100L),
    ctree  = c(list(testtype = "Bonferroni"),
               if (model_based) list(adjust = TRUE, parm = "treatment")),
    mob    = c(list(trim = 0.1, max_cuts = 100L),
               if (model_based) list(adjust = TRUE, parm = "treatment")),
    rpart  = list(xval = 10L, cp_rule = "min"),
    policy = list(cost = 0, better = "higher", split_step = NULL)), "tree_args")
  one_of <- function(nm, choices)
    if (!is.character(ta[[nm]]) || length(ta[[nm]]) != 1L ||
        !ta[[nm]] %in% choices)
      stop(sprintf("`tree_args$%s` must be one of %s.", nm,
                   paste0("\"", choices, "\"", collapse = ", ")), call. = FALSE)
  switch(engine,
    maxt = {
      ta$n_boot   <- .hte_select_count(ta$n_boot, "tree_args$n_boot", 0, 1e6)
      ta$max_cuts <- .hte_select_count(ta$max_cuts, "tree_args$max_cuts", 1, 1e6)
      if (ta$n_boot == 0L && alpha < 1)
        stop("`tree_args$n_boot = 0` gives no p-values, so it needs `alpha = 1`.",
             call. = FALSE)
    },
    ctree = one_of("testtype", c("Bonferroni", "Univariate", "MonteCarlo")),
    mob = {
      if (!is.numeric(ta$trim) || length(ta$trim) != 1L || is.na(ta$trim) ||
          ta$trim < 0 || ta$trim >= 0.5)
        stop("`tree_args$trim` must be a single number in [0, 0.5).", call. = FALSE)
      ta$max_cuts <- .hte_select_count(ta$max_cuts, "tree_args$max_cuts", 1, 1e6)
    },
    rpart = {
      ta$xval <- .hte_select_count(ta$xval, "tree_args$xval", 0, 1e4)
      one_of("cp_rule", c("min", "1se"))
    },
    policy = {
      if (!is.numeric(ta$cost) || length(ta$cost) != 1L || !is.finite(ta$cost))
        stop("`tree_args$cost` must be a single finite number.", call. = FALSE)
      one_of("better", c("higher", "lower"))
      if (!is.null(ta$split_step))
        ta$split_step <- .hte_select_count(ta$split_step, "tree_args$split_step",
                                           1, 1e7)
    })
  if (model_based) {
    if (!isTRUE(ta$adjust) && !isFALSE(ta$adjust))
      stop("`tree_args$adjust` must be TRUE or FALSE.", call. = FALSE)
    one_of("parm", c("treatment", "all"))
  }
  seed <- .hte_select_count(seed, "seed", 0, .Machine$integer.max - 1e5)
  if (!is.list(grf_args))
    stop("`grf_args` must be a named list.", call. = FALSE)
  fixed <- intersect(names(grf_args), c("X", "Y", "W", "D", "horizon", "W.hat",
                                        "Y.hat", "sample.weights", "clusters"))
  if (is.numeric(grf_args[["W.hat"]]) && length(grf_args[["W.hat"]]) == 1L)
    fixed <- setdiff(fixed, "W.hat")
  if (length(fixed))
    stop(sprintf("`grf_args` cannot set %s in get_hte_tree(): the data columns set X, Y, W, D and horizon, and per-row fields cannot follow the split (a single known `W.hat` can).",
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
  miss <- split_var[vapply(data[split_var], anyNA, logical(1L))]
  if (length(miss))
    stop(sprintf("get_hte_tree() needs complete split variables, since a rule cannot say where a missing value goes; %s have missing values. Impute them or drop those rows first.",
                 paste0("`", miss, "`", collapse = ", ")), call. = FALSE)
  miss <- if (model_based && ta$adjust)
    adj_var[vapply(data[adj_var], anyNA, logical(1L))]
  if (length(miss))
    stop(sprintf("The node models of `method = \"%s\"` adjust for `adj_var`, which needs complete columns; %s have missing values. Impute them, drop those rows, or set `tree_args = list(adjust = FALSE)`.",
                 method, paste0("`", miss, "`", collapse = ", ")), call. = FALSE)
  W    <- .psw_treat(data[[cat_var]], cat_var, arg = "cat_var")$z
  des  <- .tree_design(data[split_var], factor_encoding)
  X    <- des$X
  cols <- des$cols

  # ---- Split -------------------------------------------------------------------
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
  same <- split_frac == 1
  if (any(table(factor(W[disc], 0:1)) < 2L) ||
      (!same && any(table(factor(W[!disc], 0:1)) < 2L)))
    stop("Too few patients per arm for this `split_frac`.", call. = FALSE)
  d_idx <- which(disc)
  e_idx <- if (same) d_idx else which(!disc)

  # ---- Discovery: scores and the tree -----------------------------------------
  ga <- grf_args
  if (is.null(ga$seed)) ga$seed <- seed
  # grf fits the nuisance survival and censoring curves at every distinct
  # time, which dominates a survival fit; 100 points up to `time` give the
  # same effects in a third of the time. The point after `time` keeps the
  # patients followed past it, whom get_hte() moves there, beyond the
  # horizon: grf puts every later time on the last grid point.
  if (is_surv && is.null(ga$failure.times)) {
    tm <- data$time
    if (length(unique(tm[tm <= time])) > 100L)
      ga$failure.times <- c(seq(min(0, tm), time, length.out = 100L),
                            if (any(tm > time)) min(tm[tm > time]))
  }
  # 500 trees give nearly the scores and leaf effects of grf's 2000 at a
  # quarter of the time; above 10,000 rows get_hte() grows 200
  hte <- function(rows) {
    g <- ga
    if (is.null(g$num.trees) && nrow(rows) <= 10000L) g$num.trees <- 500L
    do.call(get_hte, c(
      list(data = rows, cat_var = cat_var, adj_var = adj_var, surv = surv,
           factor_encoding = factor_encoding, grf_args = g),
      if (is_surv) list(time = time)))
  }
  notes <- character()
  # The model-based trees use no scores: the discovery forest is fitted only
  # for the other methods, or when discovery and estimation share patients
  fit_d <- if (!model_based || same)
    withCallingHandlers(hte(data[d_idx, , drop = FALSE]),
      warning = function(w) {
        notes <<- c(notes, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) invokeRestart("muffleMessage"))
  if (length(notes))
    warning(sprintf("The discovery fit raised %d warning(s); the first: %s",
                    length(notes), notes[1L]), call. = FALSE)
  node <- NULL
  if (model_based) {
    # Node model of the outcome on the treatment (and adj_var): lm on a
    # continuous or 0/1 outcome or on survival pseudo-values at `time`
    # (difference scales), logit or Cox (ratio scales), or Weibull AFT
    dd <- data[d_idx, , drop = FALSE]
    model <- if (node_scale == "aft") "aft" else if (node_scale == "abs") "lm"
             else if (is_surv) "cox" else "logit"
    y <- if (!is_surv) dd[[outcome]]
      else if (model %in% c("cox", "aft"))
        survival::Surv(dd$time, as.numeric(dd$DSS))
      # the times and events are written into the call, which pseudo()
      # evaluates again elsewhere
      else survival::pseudo(
        eval(bquote(survival::survfit(
          survival::Surv(.(dd$time), .(as.numeric(dd$DSS))) ~ 1))),
        times = time,
        type = if (identical(grf_args$target, "RMST")) "rmst" else "surv")
    adj <- NULL
    if (ta$adjust) {
      # A constant adjustment has no information in this discovery sample.
      # Drop it before model.matrix(): a one-level character/factor column
      # otherwise fails in contrasts(), and a constant numeric column only
      # creates a redundant coefficient.
      usable <- vapply(dd[adj_var], function(x) {
        if (is.factor(x)) x <- droplevels(x)
        length(unique(x)) > 1L
      }, logical(1L))
      if (any(usable)) {
        ad <- dd[adj_var[usable]]
        ad[] <- lapply(ad, function(x) if (is.factor(x)) droplevels(x) else x)
        adj <- stats::model.matrix(~ ., data = ad)[, -1L, drop = FALSE]
      }
    }
    R <- cbind(`(Intercept)` = 1, .a = W[d_idx], adj)
    if (model == "cox") R <- R[, -1L, drop = FALSE]
    node <- list(y = y, R = R, model = model,
                 parm = if (ta$parm == "treatment") match(".a", colnames(R)))
    ok <- rep(TRUE, length(d_idx))
    gd <- NULL
  } else if (rlearner) {
    # The forest's out-of-bag residuals: numerator (W - e)(Y - m) and
    # denominator (W - e)^2, or a survival forest's own censoring-adjusted
    # pair. rpart grows on num / den weighted by den (the R-loss as weighted
    # least squares); mob and ctree test the estimating equation num - den tau.
    f <- fit_d$fit
    if (inherits(f, "causal_survival_forest")) {
      num <- f[["_psi"]]$numerator
      den <- f[["_psi"]]$denominator
    } else {
      wc  <- f$W.orig - f$W.hat
      num <- wc * (f$Y.orig - f$Y.hat)
      den <- wc^2
    }
    ok   <- is.finite(num) & is.finite(den) & den > 0
    gd   <- num[ok] / den[ok]
    node <- list(y = num[ok], R = cbind(.den = den[ok]), model = "r",
                 parm = 1L)
  } else {
    g  <- fit_d$data[[if (endsWith(method, "_cate")) ".cate" else ".dr_score"]]
    ok <- is.finite(g)
    gd <- g[ok]
  }
  Xd <- X[d_idx, , drop = FALSE][ok, , drop = FALSE]
  Wd <- W[d_idx][ok]
  # A node model needs about 10 patients per coefficient
  ctrl <- c(list(max_depth = max_depth, alpha = alpha, min_arm = 2L,
                 min_n = max(2L, ceiling(min_leaf * sum(ok)),
                             if (model_based)
                               10L * (ncol(node$R) + (node$model == "aft")))), ta)
  grown <- if (method == "maxt") .tree_maxt(gd, Xd, Wd, cols$var, ctrl)
           else .tree_engine(
             if (model_based || (rlearner && engine != "rpart"))
               paste0(engine, "_model") else engine,
             gd, Xd, ctrl, node, weights = if (rlearner) node$R[, 1L])
  tree <- .tree_prune(grown$tree, Xd, Wd, ctrl$min_arm)
  if (verbose)
    cli::cli_inform(c("i" = paste(
      if (same) "Discovery and estimation on the same {length(d_idx)} patients;"
      else "Discovery {length(d_idx)} and estimation {length(e_idx)} patients;",
      "{sum(!ok)} undefined score{?s} left out;",
      "{nrow(grown$tree) - nrow(tree)} split{?s} dropped for a missing arm.")))

  # ---- Estimation part -------------------------------------------------------------
  leaves  <- .icf_leaves(tree, cols)
  pe      <- .icf_assign(tree, X[e_idx, , drop = FALSE])
  pd      <- .icf_assign(tree, X[d_idx, , drop = FALSE])
  rule_of <- factor(leaves$rule[match(pe, leaves$path)], levels = leaves$rule)
  if (same) {
    est <- fit_d
    est$data$.rule <- rule_of
  } else {
    est_data <- data[e_idx, , drop = FALSE]
    est_data$.rule <- rule_of
    est <- hte(est_data)
  }
  ea   <- attr(est, "analysis")
  s    <- .hte_arm_scores(est$fit)
  grid <- data.frame(estimand = "ATE", measure = "diff", stringsAsFactors = FALSE)
  risk <- identical(ea$target, "survival.probability")
  z    <- stats::qnorm(0.975)
  beyond <- .hte_beyond(est)
  sub <- .hte_muffle_ps(.hte_subgroup(est$fit, s, est$data, ".rule", grid,
                                      risk, z, beyond, estimator))
  m <- match(leaves$rule, sub$level)

  # Nodes depth first, as partykit numbers them
  dfs <- function(p) c(p, if (p %in% tree$path)
    c(dfs(paste0(p, "L")), dfs(paste0(p, "R"))))
  paths <- dfs("")
  inner <- paths %in% tree$path
  j     <- match(paths, tree$path)
  side  <- vapply(which(inner), function(i) .icf_leaves(
    data.frame(path = "", col = tree$col[j[i]], value = tree$value[j[i]],
               right = if ("right" %in% names(tree)) tree$right[j[i]] else TRUE),
    cols)$rule[1L], character(1L))
  est_node <- do.call(rbind, lapply(paths, function(p) {
    if (!(p %in% tree$path)) {
      k <- m[match(p, leaves$path)]
      return(if (is.na(k)) data.frame(estimate = NA_real_, std.error = NA_real_,
                                      conf.low = NA_real_, conf.high = NA_real_,
                                      p.value = NA_real_)
             else as.data.frame(sub[k, c("estimate", "std.error", "conf.low",
                                         "conf.high", "p.value")]))
    }
    .hte_muffle_ps(.hte_estimate(est$fit, s, startsWith(pe, p), grid, risk, z,
                                 if (nzchar(p)) sprintf("Node %d", match(p, paths))
                                 else "All patients", beyond, estimator))[
      c("estimate", "std.error", "conf.low", "conf.high", "p.value")]
  }))
  st <- match(paths, grown$stats$path)
  nodes <- tibble::tibble(
    node = seq_along(paths),
    parent = ifelse(nzchar(paths),
                    match(substring(paths, 1L, nchar(paths) - 1L), paths),
                    NA_integer_),
    depth = nchar(paths), terminal = !inner,
    rule = leaves$rule[match(paths, leaves$path)],
    variable = cols$var[tree$col[j]],
    split = replace(rep(NA_character_, length(paths)), inner, side),
    statistic = grown$stats$statistic[st], split_p = grown$stats$split_p[st],
    n_disc = vapply(paths, function(p) sum(startsWith(pd, p)), integer(1L),
                    USE.NAMES = FALSE),
    n = vapply(paths, function(p) sum(startsWith(pe, p)), integer(1L),
               USE.NAMES = FALSE),
    n_treat = vapply(paths, function(p) sum(W[e_idx][startsWith(pe, p)]),
                     numeric(1L), USE.NAMES = FALSE),
    est_node)
  rules <- tibble::tibble(
    leaf = seq_len(nrow(leaves)), node = match(leaves$path, paths),
    rule = leaves$rule, n_disc = nodes$n_disc[match(leaves$path, paths)],
    n = sub$n[m], n_treat = sub$n_treat[m], estimate = sub$estimate[m],
    std.error = sub$std.error[m], conf.low = sub$conf.low[m],
    conf.high = sub$conf.high[m], p.value = sub$p.value[m],
    p_inter = sub$p_inter[m])
  rules$n[is.na(m)] <- 0L
  rules$n_treat[is.na(m)] <- 0L
  if (method == "policy") {
    sgn <- if (ta$better == "higher") 1 else -1
    gain <- tapply(sgn * gd - ta$cost, factor(.icf_assign(tree, Xd),
                                              levels = leaves$path), mean)
    rules$action <- ifelse(as.vector(gain[leaves$path]) > 0, "Treated",
                           "Control")
    nodes$action <- rules$action[match(paths, leaves$path)]
  }

  # The tree as a party on the estimation patients; every node's info is its
  # row of `nodes`
  tr <- .icf_party(tree, cols, X[e_idx, , drop = FALSE], data.frame(
    est$data[outcome], .arm = factor(W[e_idx], 0:1, c("Control", "Treated")),
    .dr_score = est$data$.dr_score, .rule = est$data$.rule,
    check.names = FALSE), rules, leaves)
  info_of <- function(id) as.list(nodes[id, setdiff(names(nodes),
                                                    c("node", "parent"))])
  relabel <- function(nd) {
    id <- partykit::id_node(nd)
    if (partykit::is.terminal(nd))
      return(partykit::partynode(id, info = info_of(id)))
    partykit::partynode(id, split = partykit::split_node(nd),
                        kids = lapply(partykit::kids_node(nd), relabel),
                        info = info_of(id))
  }
  tr$node <- relabel(partykit::node_party(tr))

  structure(
    list(rules = rules, nodes = nodes, tree = tr, est = est),
    class = c("hte_tree", "list"),
    analysis = list(
      backend_version = ea$backend_version, forest = ea$forest,
      outcome_type = ea$outcome_type, cat_var = cat_var, treated = ea$treated,
      surv = surv, outcome = outcome, time = ea$time, target = ea$target,
      adj_var = adj_var, candidate_var = candidate_var, split_var = split_var,
      factor_encoding = factor_encoding, method = method,
      node_model = if (model_based) switch(node$model,
        cox = "Cox (log-HR)", logit = "logit (log-OR)",
        aft = "Weibull AFT (log-time ratio)",
        lm = if (is_surv) sprintf("lm on %s pseudo-values",
                                  if (identical(ea$target, "RMST")) "RMST" else "S(t)")
             else if (identical(ea$outcome_type, "binary")) "lm (risk difference)"
             else "lm (mean difference)"),
      max_depth = max_depth, alpha = if (tested) alpha else NA_real_,
      min_leaf = min_leaf, split_frac = split_frac, estimator = estimator,
      tree_args = ta,
      seed = seed, cols = cols, discovery = which(keep)[d_idx],
      scores_left_out = sum(!ok),
      splits_dropped = nrow(grown$tree) - nrow(tree),
      call = match.call()))
}


# ---- L1 plot -----------------------------------------------------------------------

#' Draw the subgroup tree of get_hte_tree()
#'
#' Draws `x$tree` with ggparty, as [plt_hte_icf()] draws its tree: the split
#' variable and the p-value of its split test at every inner node, the
#' condition on every edge and, above each leaf, its patients and its effect
#' with the 95% interval from `x$rules` (and, for a policy tree, the arm it
#' recommends). Beneath each leaf a panel chosen by `type` shows that
#' subgroup.
#'
#' @param x A [get_hte_tree()] result with at least one split.
#' @inheritParams plt_hte_icf
#'
#' @details The panels show the patients the effects were estimated on: the
#' estimation part, or with `split_frac = 1` every patient. The p-values are
#' those of the split tests on the discovery part, so only `"maxt"` and the
#' mob and ctree trees show them. Needs the partykit and ggparty packages.
#'
#' @return A ggplot built by ggparty, with the size it was laid out for in
#'   `attr(, "plot_size")` (inches). When `save` is non-empty, the plot is
#'   also written to PDF through `RegR::save_plt()`.
#'
#' @seealso [get_hte_tree()]; [plt_hte_sub()] draws the leaf effects as a
#'   forest plot: `plt_hte_sub(x$est, sub_var = ".rule")`.
#'
#' @examplesIf requireNamespace("grf", quietly = TRUE) && requireNamespace("ggparty", quietly = TRUE)
#' \donttest{
#' set.seed(20260929)
#' n <- 1600
#' d <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rnorm(n), X3 = runif(n))
#' d$z <- rbinom(n, 1, plogis(0.5 * d$X2))
#' d$y <- d$X2 + d$z * (0.5 + 1.5 * (d$X3 > 0.6)) + rnorm(n)
#'
#' res <- get_hte_tree(d, cat_var = "z", adj_var = c("X1", "X2", "X3"),
#'                     surv = "y", tree_args = list(n_boot = 500))
#' plt_hte_tree(res)
#' plt_hte_tree(res, type = "dr")
#' }
#'
#' @export
plt_hte_tree <- function(x, type = c("effect", "dr", "box", "bar", "km"),
                         save = list()) {
  if (!inherits(x, "hte_tree"))
    stop("`x` must be a get_hte_tree() result.", call. = FALSE)
  .hte_draw_tree(x, match.arg(type), save, "plt_hte_tree", "get_hte_tree")
}


# ---- L3 print ------------------------------------------------------------------

#' @export
#' @noRd
print.hte_tree <- function(x, ...) {
  a    <- attr(x, "analysis")
  same <- isTRUE(a$split_frac == 1)
  what <- if (a$method == "maxt") "max-t test" else sub("_.*$", "", a$method)
  cat(sprintf("<hte_tree> %s tree on %s (%s, grf %s)\n", what,
              if (!is.null(a$node_model))
                sprintf("the scores of %s node models", a$node_model)
              else if (endsWith(a$method, "_cate")) "out-of-bag CATE predictions"
              else if (endsWith(a$method, "_r")) "R-learner residuals"
              else "doubly robust scores", a$forest, a$backend_version))
  cat(sprintf("  %s\n", if (same)
    sprintf("n = %d, tree grown and estimated on the same patients",
            nrow(x$est$data))
    else sprintf("discovery n = %d, estimation n = %d", x$nodes$n_disc[1L],
                 nrow(x$est$data))))
  cat(sprintf("  max_depth = %s, %smin_leaf = %s; split on: %s\n",
              format(a$max_depth),
              if (is.na(a$alpha)) "" else sprintf("alpha = %s, ", format(a$alpha)),
              format(a$min_leaf), paste(a$split_var, collapse = ", ")))
  cat("\nNodes (split tests on the discovery part, effects on the estimation part):\n")
  shown <- c("node", "parent", "split",
             if (!all(is.na(x$nodes$split_p))) "split_p", "n", "estimate",
             "conf.low", "conf.high")
  print(as.data.frame(x$nodes[shown]), row.names = FALSE, digits = 3)
  cat(sprintf("\nRules, estimated on %s (ATE difference%s, 95%% CI):\n",
              if (same) "the patients they were found on (not honest)"
              else "the estimation part",
              if (identical(a$estimator, "tmle")) " by TMLE" else ""))
  print(as.data.frame(x$rules[c("rule", if (!same) "n_disc", "n", "estimate",
                                "conf.low", "conf.high", "p.value",
                                if (a$method == "policy") "action")]),
        row.names = FALSE, digits = 3)
  if (nrow(x$rules) > 1L && !is.na(x$rules$p_inter[1L]))
    cat(sprintf("P for interaction = %s\n",
                format(signif(x$rules$p_inter[1L], 3))))
  invisible(x)
}
