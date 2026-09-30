tree_data <- function(n = 1600L, tau = function(d) 0.5 + 1.5 * (d$X3 > 0.6),
                      seed = 7) {
  withr::local_seed(seed)
  d <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rnorm(n), X3 = runif(n),
                  grp = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  d$z <- rbinom(n, 1, plogis(0.5 * d$X2))
  d$y <- d$X2 + d$z * tau(d) + rnorm(n)
  d
}

tree_call <- function(d, ..., adj_var = c("X1", "X2", "X3", "grp"))
  get_hte_tree(d, "z", adj_var, surv = "y", grf_args = list(num.trees = 500L),
               ...)

test_that("CART methods explicitly name their discovery target", {
  methods <- eval(formals(get_hte_tree)$method)
  expect_true("rpart_dr" %in% methods)
  expect_true(all(c("mob_r", "ctree_r", "rpart_r") %in% methods))
  expect_false("rpart" %in% methods)
  expect_identical(methods[1L], "maxt")
  expect_identical(eval(formals(get_hte_tree)$factor_encoding),
                   c("integer", "onehot"))
})

test_that("rpart_cate distils CATE predictions with CART controls", {
  skip_if_not_installed("grf")
  skip_if_not_installed("rpart")
  res <- tree_call(tree_data(), method = "rpart_cate", max_depth = 1,
                   tree_args = list(cp_rule = "1se"))
  expect_s3_class(res, "hte_tree")
  expect_identical(res$nodes$variable[1L], "X3")
  expect_true(all(is.na(res$nodes$split_p)))
  expect_identical(attr(res, "analysis")$tree_args$xval, 10L)
  expect_identical(sum(res$rules$n), nrow(res$est$data))
  expect_output(print(res), "rpart tree on out-of-bag CATE predictions")
  expect_error(tree_call(tree_data(), method = "rpart_cate", alpha = 0.1),
               "alpha.*cross-validation")
})

test_that("the R-learner node model is the residual-on-residual slope", {
  withr::local_seed(1)
  wt <- rnorm(300)
  yt <- 0.7 * wt + rnorm(300)
  f <- .tree_nodefit("r")(yt * wt, cbind(.den = wt^2), estfun = TRUE)
  expect_equal(unname(f$coefficients), unname(stats::coef(stats::lm(yt ~ 0 + wt))))
  expect_equal(drop(f$estfun), wt * (yt - wt * unname(f$coefficients)))
  expect_equal(f$objfun, -sum(yt * wt)^2 / sum(wt^2))
})

test_that("_r trees split on the forest's residuals for every outcome", {
  skip_if_not_installed("grf")
  skip_if_not_installed("rpart")
  d <- tree_data()
  for (m in c("mob_r", "ctree_r", "rpart_r")) {
    res <- tree_call(d, method = m, max_depth = 1)
    expect_identical(res$nodes$variable[1L], "X3", info = m)
    expect_identical(is.na(res$nodes$split_p[1L]), m == "rpart_r", info = m)
    expect_identical(sum(res$rules$n), nrow(res$est$data))
  }
  expect_output(print(res), "rpart tree on R-learner residuals")
  expect_error(tree_call(d, method = "rpart_r", alpha = 0.1), "cross-validation")
  expect_error(tree_call(d, method = "mob_r", tree_args = list(adjust = FALSE)),
               "unknown field")

  skip_on_cran()
  withr::local_seed(3)
  ev <- rexp(nrow(d), 0.05 * exp(0.3 * d$X2 - d$z * 2 * (d$X3 > 0.6)))
  cens <- rexp(nrow(d), 0.02)
  d$time <- pmin(ev, cens)
  d$DSS <- as.integer(ev <= cens)
  rs <- get_hte_tree(d, "z", c("X1", "X2", "X3"), time = 12, method = "mob_r",
                     max_depth = 1, grf_args = list(num.trees = 500L))
  expect_s3_class(rs$est$fit, "causal_survival_forest")
  expect_identical(rs$nodes$variable[1L], "X3")
})

test_that("AFT trees find time-ratio modifiers and retain DR leaf estimates", {
  skip_if_not_installed("grf")
  skip_if_not_installed("survival")
  withr::local_seed(82)
  n <- 1600L
  d <- data.frame(x = rep(0:1, each = n / 2), z = rbinom(n, 1, 0.5),
                   prognostic = rnorm(n))
  ev <- rweibull(n, shape = 2, scale = exp(2 + 0.3 * d$prognostic +
                                           d$z * (0.1 + 1.5 * d$x)))
  cens <- rexp(n, 0.025)
  d$time <- pmin(ev, cens)
  d$DSS <- as.integer(ev <= cens)
  for (m in c("mob_aft", "ctree_aft")) {
    res <- get_hte_tree(d, "z", c("x", "prognostic"), method = m,
                        time = 8, max_depth = 2,
                        grf_args = list(num.trees = 300L, num.threads = 2L))
    expect_identical(res$nodes$variable[1L], "x", info = m)
    expect_identical(attr(res, "analysis")$node_model,
                     "Weibull AFT (log-time ratio)")
    expect_true(all(is.finite(res$rules$estimate)))
    expect_true(all(abs(res$rules$estimate) <= 1))
    expect_identical(sum(res$rules$n), nrow(res$est$data))
    expect_output(print(res), "Weibull AFT")
    expect_error(tree_call(tree_data(), method = m), "surv = TRUE")
    bad <- d
    bad$time[1L] <- 0
    expect_error(get_hte_tree(bad, "z", "x", method = m), "positive")
  }
  rmst <- get_hte_tree(d, "z", "x", method = "ctree_aft", time = 8,
                       max_depth = 1, tree_args = list(adjust = FALSE, parm = "all"),
                       grf_args = list(target = "RMST", num.trees = 300L,
                                       num.threads = 2L))
  expect_identical(attr(rmst, "analysis")$target, "RMST")
  expect_identical(rmst$nodes$variable[1L], "x")
  expect_true(all(is.finite(rmst$rules$estimate)))
  for (outcome in list("DSS", TRUE)) {
    args <- list(data = d, cat_var = "z", adj_var = c("x", "prognostic"),
                 surv = outcome, method = "rpart_cate", max_depth = 1,
                 tree_args = list(xval = 0L),
                 grf_args = list(num.trees = 300L, num.threads = 2L))
    if (isTRUE(outcome)) args$time <- 8
    cart <- do.call(get_hte_tree, args)
    expect_s3_class(cart, "hte_tree")
    expect_true(all(is.finite(cart$rules$estimate)))
  }
})

test_that("AFT node scores match the weighted Weibull likelihood", {
  skip_if_not_installed("survival")
  withr::local_seed(15)
  n <- 200L
  x <- cbind(`(Intercept)` = 1, treatment = rbinom(n, 1, 0.5),
              constant = 1, covariate = rnorm(n))
  ev <- rweibull(n, 1.8, exp(1 + 0.4 * x[, 2] + 0.2 * x[, 4]))
  cens <- rexp(n, 0.1)
  y <- survival::Surv(pmin(ev, cens), as.integer(ev <= cens))
  w <- rep(c(1, 2), length.out = n)
  f <- .tree_nodefit("aft")(y, x, weights = w, estfun = TRUE, object = TRUE)
  expect_true(is.na(f$coefficients[3L]))
  active <- c(1L, 2L, 4L)
  theta <- c(f$coefficients[active], log(f$object$scale))
  ll <- function(b) {
    shape <- exp(-b[4L])
    scale <- exp(drop(x[, active] %*% b[1:3]))
    w * ifelse(y[, 2] == 1,
               dweibull(y[, 1], shape, scale, log = TRUE),
               pweibull(y[, 1], shape, scale, lower.tail = FALSE, log.p = TRUE))
  }
  numeric_scores <- sapply(seq_along(theta), function(j) {
    step <- rep(0, length(theta)); step[j] <- 1e-5
    (ll(theta + step) - ll(theta - step)) / 2e-5
  })
  expect_equal(unname(f$estfun), unname(numeric_scores), tolerance = 1e-6)
  expect_equal(f$objfun, -sum(ll(theta)), tolerance = 1e-7)
  tr <- .tree_nodefit("aft", keep = 2L, decorrelate = TRUE)(
    y, x, weights = w, estfun = TRUE)
  expect_identical(dim(tr$estfun), c(n, 1L))
  expect_true(all(is.finite(tr$estfun)))
})

test_that("the max-t tree finds a threshold and estimates it on the other half", {
  skip_if_not_installed("grf")
  skip_if_not_installed("partykit")
  d <- tree_data()
  res <- tree_call(d, tree_args = list(n_boot = 300L))
  expect_s3_class(res, "hte_tree")
  expect_named(res, c("rules", "nodes", "tree", "est"))
  expect_identical(res$nodes$variable[1L], "X3")
  cut <- as.numeric(sub("^X3 <= ", "", res$nodes$split[1L]))
  expect_lt(abs(cut - 0.6), 0.05)
  expect_lt(res$nodes$split_p[1L], 0.01)
  expect_identical(nrow(res$rules), 2L)
  expect_lt(abs(res$rules$estimate[1L] - 0.5), 0.3)
  expect_lt(abs(res$rules$estimate[2L] - 2), 0.3)
  # Honest: the estimation half is disjoint from the discovery half
  expect_identical(res$nodes$n_disc[1L] + nrow(res$est$data), nrow(d))
  expect_identical(sum(res$rules$n), nrow(res$est$data))
  expect_true(is.factor(res$est$data$.rule))
  expect_equal(res$nodes$estimate[res$rules$node], res$rules$estimate)
  # The party carries each node's row of $nodes
  expect_s3_class(res$tree, "party")
  expect_equal(partykit::width(res$tree), nrow(res$rules))
  info <- partykit::nodeapply(res$tree, partykit::nodeids(res$tree),
                              partykit::info_node)
  expect_equal(unname(vapply(info, function(i) i$estimate, numeric(1L))),
               res$nodes$estimate)
  expect_output(print(res), "max-t test tree")
  p <- suppressMessages(plt_hte_sub(res$est, sub_var = ".rule",
                                    fixed_size = FALSE))
  expect_s3_class(p, "ggplot")
})

test_that("ggparty reads the node info of the tree", {
  skip_if_not_installed("grf")
  skip_if_not_installed("ggparty")
  res <- tree_call(tree_data(), max_depth = 1, tree_args = list(n_boot = 200L))
  g <- ggparty::ggparty(res$tree, add_vars = list(
    est = "$node$info$estimate", p = "$node$info$split_p"))
  expect_s3_class(g, "ggplot")
  expect_equal(g$data$est[order(g$data$id)], res$nodes$estimate)
  expect_equal(g$data$p[order(g$data$id)], res$nodes$split_p)
})

test_that("plt_hte_tree() draws the tree with its split p-values", {
  skip_if_not_installed("grf")
  skip_if_not_installed("ggparty")
  d <- tree_data()
  labels <- function(p) unlist(lapply(p$layers, function(l)
    if (is.data.frame(l$data)) l$data$label))
  res <- tree_call(d, max_depth = 2, tree_args = list(n_boot = 200L))
  p <- plt_hte_tree(res)
  expect_s3_class(p, "ggplot")
  expect_named(attr(p, "plot_size"), c("width", "height"))
  expect_true(any(grepl("^X3\np ", labels(p))))
  expect_s3_class(plt_hte_tree(res, type = "dr"), "ggplot")
  expect_s3_class(plt_hte_tree(res, type = "box"), "ggplot")
  expect_error(plt_hte_tree(res, type = "km"), "survival")
  expect_error(plt_hte_tree(res$rules), "get_hte_tree")
  pol <- tree_call(d, method = "policy", max_depth = 1,
                   tree_args = list(cost = 1))
  expect_true(any(grepl("\nTreated$", labels(plt_hte_tree(pol)))))
  flat <- tree_call(tree_data(tau = function(d) rep(0.5, nrow(d))),
                    tree_args = list(n_boot = 300L))
  expect_error(plt_hte_tree(flat), "no subgroups")
})

test_that("a constant effect keeps every patient together unless alpha = 1", {
  skip_if_not_installed("grf")
  d <- tree_data(tau = function(d) rep(0.5, nrow(d)))
  res <- tree_call(d, tree_args = list(n_boot = 300L))
  expect_identical(res$rules$rule, "All patients")
  expect_gt(res$nodes$split_p[1L], 0.05)
  full <- tree_call(d, alpha = 1, max_depth = 2, tree_args = list(n_boot = 0L))
  expect_identical(nrow(full$rules), 4L)
  expect_true(all(is.na(full$nodes$split_p)))
})

test_that("every method grows on the scores within max_depth", {
  skip_if_not_installed("grf")
  skip_if_not_installed("rpart")
  skip_if_not_installed("policytree")
  d <- tree_data()
  args <- list(maxt = list(n_boot = 200L), mob_dr = list(), mob_cate = list(),
               mob_abs = list(), ctree_dr = list(), ctree_cate = list(),
               ctree_abs = list(), rpart_dr = list(), policy = list(cost = 1),
               mob_r = list(), ctree_r = list(), rpart_r = list())
  for (m in names(args)) {
    one <- tree_call(d, method = m, max_depth = 1, tree_args = args[[m]])
    expect_identical(one$nodes$variable[1L], "X3", info = m)
    expect_identical(nrow(one$rules), 2L, info = m)
    two <- tree_call(d, method = m, max_depth = 2, tree_args = args[[m]])
    expect_lte(max(two$nodes$depth), 2L)
    expect_identical(attr(two, "analysis")$method, m)
  }
  # ctree and mob report their split p-values, rpart and policy none
  expect_false(anyNA(tree_call(d, method = "ctree_dr", max_depth = 1)$nodes$split_p[1L]))
  expect_true(all(is.na(tree_call(d, method = "rpart_dr", max_depth = 1)$nodes$split_p)))
  # The _cate methods grow on the forest's predictions, not on the scores
  cate <- tree_call(d, method = "mob_cate", max_depth = 1)
  expect_output(print(cate), "mob tree on out-of-bag CATE predictions")
  expect_output(print(tree_call(d, method = "mob_dr", max_depth = 1)),
                "mob tree on doubly robust scores")
})

model_data <- function(n = 2400L, seed = 7) {
  withr::local_seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rbinom(n, 1, 0.5),
                  x4 = rnorm(n))
  d$w <- rbinom(n, 1, plogis(0.8 * d$x1))
  d$y <- 1.5 * d$x2 + d$x1 + d$w * (0.5 + d$x3) + rnorm(n)
  d
}

test_that("node-model trees test the treatment coefficient, adjusted for adj_var", {
  skip_if_not_installed("grf")
  d  <- model_data()
  av <- c("x1", "x2", "x3", "x4")
  for (m in c("mob_abs", "ctree_abs")) {
    res <- get_hte_tree(d, "w", av, surv = "y", method = m, max_depth = 2,
                        grf_args = list(num.trees = 500L))
    expect_identical(unique(stats::na.omit(res$nodes$variable)), "x3", info = m)
    expect_identical(attr(res, "analysis")$node_model, "lm (mean difference)")
    # The packages' default also splits on the prognostic x2 or confounder x1
    pkg <- get_hte_tree(d, "w", av, surv = "y", method = m, max_depth = 2,
                        tree_args = list(adjust = FALSE, parm = "all"),
                        grf_args = list(num.trees = 500L))
    expect_true(any(c("x1", "x2") %in% pkg$nodes$variable), info = m)
  }
  expect_output(print(res), "ctree tree on the scores of lm")
  expect_error(get_hte_tree(d, "w", av, surv = "y", method = "mob_rel"),
               "log-OR")
  expect_error(get_hte_tree(d, "w", av, surv = "y", method = "mob_abs",
                            tree_args = list(parm = "none")), "parm")
  expect_error(get_hte_tree(d, "w", av, surv = "y", method = "mob_dr",
                            tree_args = list(adjust = FALSE)), "unknown field")
  d$x4[1] <- NA
  expect_error(get_hte_tree(d, "w", av, candidate_var = "x3", surv = "y",
                            method = "mob_abs"), "adjust = FALSE")
})

test_that("node-model trees tolerate constant adjustment columns", {
  skip_if_not_installed("grf")
  d <- model_data(n = 1600L)
  d$one_level <- "only"
  d$constant <- 1
  av <- c("x1", "x2", "x3", "x4", "one_level", "constant")
  for (m in c("mob_abs", "ctree_abs")) {
    res <- get_hte_tree(d, "w", av, surv = "y", method = m, max_depth = 1,
                        grf_args = list(num.trees = 300L))
    expect_s3_class(res, "hte_tree")
    expect_setequal(attr(res, "analysis")$adj_var, av)
  }
})

test_that("node-model trees stop splitting when a logistic fit fails", {
  withr::local_seed(21)
  d <- data.frame(x = seq(-1, 1, length.out = 800L), w = rbinom(800L, 1, 0.5))
  d$y <- as.integer(d$x > 0)
  for (method in c("ctree_rel", "mob_rel")) {
    expect_warning(res <- get_hte_tree(d, "w", "x", surv = "y", method = method,
                                       max_depth = 1,
                                       grf_args = list(num.trees = 100L,
                                                       W.hat = 0.5)),
                    "node model.*failed")
    expect_identical(nrow(res$rules), 1L)
    expect_match(attr(res, "analysis")$node_failures, "logit")
  }
  x <- cbind(`(Intercept)` = 1, .a = rep(0:1, 100L),
              x = seq(-1, 1, length.out = 200L))
  expect_error(.tree_nodefit("logit", keep = 2L)(as.integer(x[, "x"] > 0), x,
                                                  estfun = TRUE),
               "logit node model")
})

test_that("node-model trees fit logit and Cox models, or pseudo-values", {
  skip_if_not_installed("grf")
  skip_if_not_installed("survival")
  skip_on_cran()
  d  <- model_data(n = 3000L)
  av <- c("x1", "x2", "x3", "x4")
  withr::local_seed(3)
  d$yb <- rbinom(nrow(d), 1, plogis(-0.5 + 0.8 * d$x2 + 0.5 * d$x1 +
                                      d$w * (0.2 + 1.2 * d$x3)))
  for (m in c("mob_rel", "ctree_rel")) {
    rb <- get_hte_tree(d, "w", av, surv = "yb", method = m, max_depth = 1,
                       grf_args = list(num.trees = 500L))
    expect_identical(rb$nodes$variable[1L], "x3", info = m)
    expect_identical(attr(rb, "analysis")$node_model, "logit (log-OR)")
  }
  ev <- rexp(nrow(d), 0.05 * exp(0.6 * d$x2 + 0.3 * d$x1 - d$w * (0.2 + d$x3)))
  cens <- rexp(nrow(d), 0.02)
  d$time <- pmin(ev, cens)
  d$DSS <- as.integer(ev <= cens)
  rs <- get_hte_tree(d, "w", av, time = 12, method = "mob_rel", max_depth = 1,
                     grf_args = list(num.trees = 500L))
  expect_identical(rs$nodes$variable[1L], "x3")
  expect_identical(attr(rs, "analysis")$node_model, "Cox (log-HR)")
  ra <- get_hte_tree(d, "w", av, time = 12, method = "ctree_abs",
                     max_depth = 1, grf_args = list(num.trees = 500L))
  expect_identical(ra$nodes$variable[1L], "x3")
  expect_identical(attr(ra, "analysis")$node_model, "lm on S(t) pseudo-values")
})

test_that("policy trees recommend the arm by cost and direction", {
  skip_if_not_installed("grf")
  skip_if_not_installed("policytree")
  d <- tree_data()
  pol <- tree_call(d, method = "policy", max_depth = 1,
                   tree_args = list(cost = 1))
  expect_identical(pol$rules$action, c("Control", "Treated"))
  expect_identical(pol$nodes$action[pol$rules$node], pol$rules$action)
  low <- tree_call(d, method = "policy", max_depth = 1,
                   tree_args = list(cost = -1, better = "lower"))
  expect_identical(low$rules$action, c("Treated", "Control"))
  # Everyone gains more than no cost: nothing to split
  flat <- tree_call(d, method = "policy", max_depth = 1)
  expect_identical(flat$rules$rule, "All patients")
  expect_output(print(pol), "action")
})

test_that("engine trees read back the partitions their packages make", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("rpart")
  skip_if_not_installed("policytree")
  withr::local_seed(2)
  n <- 600L
  X <- cbind(a = rnorm(n), b = rbinom(n, 1, 0.5), c = round(runif(n), 1))
  g <- 2 * X[, "b"] + 1.5 * (X[, "c"] > 0.5) + rnorm(n)
  ctrl <- list(max_depth = 2, alpha = 0.05, min_n = 30L,
               testtype = "Bonferroni", trim = 0.1, xval = 0L,
               cp_rule = "min", cost = 1, better = "higher", split_step = 1L)
  same_partition <- function(a, b)
    nrow(unique(data.frame(a, b))) == length(unique(a)) &&
      length(unique(a)) == length(unique(b))
  d <- data.frame(.g = g, x1 = X[, 1L], x2 = X[, 2L], x3 = X[, 3L])
  ours <- function(m) .icf_assign(.tree_engine(m, g, X, ctrl)$tree, X)
  ct <- partykit::ctree(.g ~ ., data = d, control = partykit::ctree_control(
    maxdepth = 2, minbucket = 30L, minsplit = 60L))
  expect_true(same_partition(ours("ctree"), predict(ct, type = "node")))
  rp <- rpart::rpart(.g ~ ., data = d, control = rpart::rpart.control(
    maxdepth = 2, cp = 0, minbucket = 30L, minsplit = 60L, xval = 0L))
  expect_true(same_partition(ours("rpart"), rp$where))
  pt <- policytree::policy_tree(unname(X), cbind(0, g - 1), depth = 2,
                                min.node.size = 30L)
  expect_true(same_partition(ours("policy"),
                             predict(pt, unname(X), type = "node.id")))
})

test_that("engine cut-points keep their strictness for held-out values", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("rpart")
  X <- matrix(c(rep(0, 40), rep(10, 40)), ncol = 1L,
              dimnames = list(NULL, "x"))
  g <- c(rep(0, 40), rep(10, 40))
  ctrl <- list(max_depth = 1, alpha = 0.05, min_n = 10L,
               testtype = "Bonferroni", trim = 0.1, xval = 0L,
               cp_rule = "min", cost = 1, better = "higher", split_step = 1L)
  tree <- .tree_engine("rpart", g, X, ctrl)$tree
  held_out <- matrix(c(0, 2, 4, 5, 6, 10), ncol = 1L,
                     dimnames = list(NULL, "x"))
  expect_identical(.icf_assign(tree, held_out),
                   c("L", "L", "L", "R", "R", "R"))
  expect_false(tree$right[1L])
})

test_that("mob engines try at most max_cuts cut-points per variable", {
  skip_if_not_installed("partykit")
  withr::local_seed(3)
  n <- 800L
  X <- cbind(a = round(runif(n), 3), b = rbinom(n, 1, 0.5))
  g <- 2 * (X[, "a"] > 0.37) + rnorm(n)
  ctrl <- list(max_depth = 1, alpha = 0.05, min_n = 40L, trim = 0.1,
               max_cuts = 4L)
  cuts <- sort(unique(X[, "a"]))
  cuts <- cuts[-length(cuts)]
  kept <- cuts[unique(round(seq(1, length(cuts), length.out = 4L)))]
  node <- list(y = g, R = cbind(.den = rep(1, n)), model = "r", parm = 1L)
  for (m in c("mob", "mob_model")) {
    tree <- .tree_engine(m, g, X, ctrl, node)$tree
    expect_identical(tree$col, 1L)
    expect_true(tree$value %in% kept)
  }
  # without max_cuts every distinct value is tried
  ctrl$max_cuts <- NULL
  expect_false(.tree_engine("mob", g, X, ctrl)$tree$value %in% kept)
})

test_that("tree_args$max_cuts caps mob cut-points at 100 by default", {
  skip_if_not_installed("partykit")
  withr::local_seed(4)
  n <- 2000L
  d <- data.frame(x = round(runif(n), 4), w = rnorm(n))
  d$z <- rbinom(n, 1, 0.5)
  d$y <- d$z * (d$x > 0.5) + rnorm(n)
  fit <- function(...) get_hte_tree(d, "z", c("x", "w"), surv = "y",
                                    method = "mob_abs", max_depth = 1,
                                    grf_args = list(num.trees = 50), ...)
  kept_cuts <- function(res, k) {
    x_disc <- d$x[attr(res, "analysis")$discovery]
    cuts <- sort(unique(x_disc))
    cuts <- cuts[-length(cuts)]
    cuts[unique(round(seq(1, length(cuts), length.out = k)))]
  }
  root_cut <- function(res)
    partykit::breaks_split(partykit::split_node(partykit::node_party(res$tree)))
  res <- fit()
  expect_identical(attr(res, "analysis")$tree_args$max_cuts, 100L)
  expect_identical(res$nodes$variable[1L], "x")
  expect_true(root_cut(res) %in% kept_cuts(res, 100L))
  res3 <- fit(tree_args = list(max_cuts = 3))
  expect_identical(attr(res3, "analysis")$tree_args$max_cuts, 3L)
  expect_true(root_cut(res3) %in% kept_cuts(res3, 3L))
  # a cap above the number of distinct values tries every value
  all <- fit(tree_args = list(max_cuts = 1e6))
  expect_false(root_cut(all) %in% kept_cuts(all, 100L))
  expect_error(fit(tree_args = list(max_cuts = 0)), "max_cuts")
})

test_that("both forests grow 500 trees unless grf_args sets num.trees", {
  withr::local_seed(5)
  n <- 600L
  d <- data.frame(x = runif(n), w = rnorm(n))
  d$z <- rbinom(n, 1, 0.5)
  d$y <- d$z * (d$x > 0.5) + rnorm(n)
  fit <- function(...) get_hte_tree(d, "z", c("x", "w"), surv = "y",
                                    method = "rpart_dr", ...)
  ntree <- function(res) res$est$fit[["_num_trees"]]
  expect_identical(ntree(fit()), 500)
  expect_identical(ntree(fit(grf_args = list(num.trees = 60))), 60)
})

test_that("survival forests fit their curves on a 100-point grid by default", {
  withr::local_seed(6)
  n <- 800L
  d <- data.frame(x = runif(n), w = rnorm(n))
  d$z <- rbinom(n, 1, 0.5)
  t_ev <- rexp(n, 0.02 * exp(-0.5 * d$z * (d$x > 0.5)))
  t_c  <- runif(n, 10, 120)
  d$time <- pmin(t_ev, t_c)
  d$DSS  <- as.integer(t_ev <= t_c)
  fit <- function(...) get_hte_tree(d, "z", c("x", "w"), time = 30,
                                    method = "rpart_dr",
                                    grf_args = list(num.trees = 50, ...))
  grid <- c(seq(0, 30, length.out = 100), min(d$time[d$time > 30]))
  expect_identical(fit()$rules, fit(failure.times = grid)$rules)
  # a coarse time scale is left to grf
  d$time <- ceiling(d$time)
  expect_identical(fit()$rules,
                   fit(failure.times = sort(unique(c(0, d$time))))$rules)
})

test_that("survival grids report severe early-event compression", {
  withr::local_seed(19)
  n <- 800L
  d <- data.frame(w = rep(0:1, n / 2L), x = rep(c(0, 0, 0, 0, 1),
                                              length.out = n), z = runif(n))
  d$time <- ifelse(d$x == 1, 150, 0.01 + 0.1 * d$z + 0.85 * d$w)
  d$DSS <- 1L
  fit <- function(grid = NULL) get_hte_tree(d, "w", c("x", "z"), time = 120,
    method = "rpart_dr", max_depth = 1, tree_args = list(xval = 0L),
    grf_args = c(list(target = "RMST", num.trees = 100L, W.hat = 0.5),
                 if (!is.null(grid)) list(failure.times = grid)))
  warnings <- character()
  a <- withCallingHandlers(fit(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  expect_true(any(grepl("time grid rounds", warnings)))
  diag <- attr(a, "analysis")$time_grid
  expect_true(diag$automatic)
  expect_equal(diag$zero_fraction, 1)
  expect_length(diag$points, 101L)
  grid <- sort(unique(c(0, d$time, 120)))
  b <- suppressWarnings(fit(grid))
  expect_false(attr(b, "analysis")$time_grid$automatic)
  expect_identical(attr(b, "analysis")$time_grid$points, grid)
  expect_equal(attr(b, "analysis")$time_grid$zero_fraction, 0)
})

test_that("failed AFT candidates leave a reported unsplit tree", {
  withr::local_seed(19)
  n <- 800L
  d <- data.frame(w = rep(0:1, n / 2L), x = rep(c(0, 0, 0, 0, 1),
                                              length.out = n), z = runif(n))
  d$time <- ifelse(d$x == 1, 150, 0.01 + 0.1 * d$z + 0.85 * d$w)
  d$DSS <- 1L
  notes <- character()
  res <- withCallingHandlers(get_hte_tree(d, "w", c("x", "z"), time = 120,
    method = "mob_aft", max_depth = 1,
    grf_args = list(target = "RMST", num.trees = 100L, W.hat = 0.5)),
    warning = function(w) {
      notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  expect_identical(res$rules$rule, "All patients")
  expect_true(length(attr(res, "analysis")$node_failures) > 0L)
  expect_true(any(grepl("node model failed", notes)))
})

test_that("numeric rules retain enough cut-point precision", {
  des <- .tree_design(data.frame(x = c(10000.1, 10000.2, 10000.3)),
                      "onehot")
  tree <- data.frame(path = c("", "R"), col = c(1L, 1L),
                     value = c(10000.1, 10000.2), right = c(TRUE, TRUE))
  rules <- .icf_leaves(tree, des$cols)$rule
  expect_true(any(grepl("10000.1", rules, fixed = TRUE)))
  expect_true(any(grepl("10000.2", rules, fixed = TRUE)))
  expect_false(any(grepl("10000 < x <= 10000", rules, fixed = TRUE)))
})

test_that("splits whose children lack an arm are dropped", {
  X <- cbind(a = 1:10, b = rep(0:1, 5))
  W <- rep(c(1, 0), each = 5)
  by_a <- data.frame(path = "", col = 1L, value = 5, stringsAsFactors = FALSE)
  expect_identical(nrow(.tree_prune(by_a, X, W, 2L)), 0L)
  by_b <- data.frame(path = "", col = 2L, value = 0, stringsAsFactors = FALSE)
  expect_identical(nrow(.tree_prune(by_b, X, W, 2L)), 1L)
})

test_that("incomplete estimation leaves have diagnostics and no interaction test", {
  d <- data.frame(w = rep(0:1, 300L), x = 0)
  withr::local_seed(123)
  disc <- logical(nrow(d))
  for (a in 0:1) {
    i <- which(d$w == a)
    disc[i[sample.int(length(i), length(i) / 2L)]] <- TRUE
    j <- which(d$w == a & disc)
    d$x[j] <- rep(0:2, length.out = length(j))
  }
  e0 <- which(!disc & d$w == 0)
  e1 <- which(!disc & d$w == 1)
  d$x[e0] <- c(rep(0, 40L), rep(1, 55L), rep(2, 55L))
  d$x[e1] <- rep(1:2, each = 75L)
  d$y <- d$w * (1 + 8 * d$x) + rnorm(nrow(d), sd = 0.1)
  warnings <- character()
  res <- withCallingHandlers(get_hte_tree(d, "w", "x", surv = "y",
    method = "rpart_dr", max_depth = 2, tree_args = list(xval = 0),
    grf_args = list(num.trees = 100L, W.hat = 0.5)), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_identical(nrow(res$rules), 3L)
  expect_equal(res$rules$n_control, c(40, 55, 55))
  expect_identical(res$rules$status, c("insufficient_arm", "ok", "ok"))
  expect_true(all(is.na(res$rules$p_inter)))
  expect_true(any(grepl("Interaction test omitted", warnings)))
  expect_true(is.na(res$rules$estimate[1L]))
  expect_true(all(is.finite(res$rules$estimate[-1L])))
  expect_false(any(grepl("P for interaction", capture.output(print(res)))))
})

test_that("candidate_var names the split variables; factor rules keep levels", {
  skip_if_not_installed("grf")
  d <- tree_data(tau = function(d) 0.5 + 1.5 * (d$grp == "c"))
  for (enc in c("onehot", "integer")) {
    res <- tree_call(d, adj_var = c("X1", "X2", "X3"), candidate_var = "grp",
                     factor_encoding = enc, tree_args = list(n_boot = 200L))
    expect_setequal(res$rules$rule, c("grp != c", "grp = c"))
    a <- attr(res, "analysis")
    expect_setequal(a$adj_var, c("X1", "X2", "X3", "grp"))
    expect_identical(a$split_var, "grp")
  }
  expect_error(tree_call(d, candidate_var = "y"), "candidate_var")
})

test_that("tree prediction preserves factor encoding and strict cut-points", {
  d <- tree_data(n = 800L, tau = function(d) 5 * (d$grp == "c"))
  for (encoding in c("integer", "onehot")) {
    res <- tree_call(d, method = "rpart_dr", candidate_var = "grp",
                     max_depth = 1, factor_encoding = encoding,
                     tree_args = list(xval = 0L))
    expect_identical(predict(res, res$est$data), as.character(res$est$data$.rule))
    nd <- data.frame(grp = factor(c("a", "b", "c"), levels = c("c", "b", "a")))
    expect_identical(predict(res, nd), predict(res, transform(nd,
                                                            grp = as.character(grp))))
    expect_true(all(predict(res, nd, type = "node") %in% res$rules$node))
    expect_error(predict(res, data.frame(grp = "new")), "unseen level")
    expect_error(predict(res, data.frame(grp = NA_character_)), "missing")
  }
  res <- tree_call(tree_data(n = 800L), method = "rpart_dr", max_depth = 1,
                   candidate_var = "X3", tree_args = list(xval = 0L))
  cut <- attr(res, "analysis")$splits$value[1L]
  expect_identical(predict(res, data.frame(X3 = c(cut - 1e-9, cut))),
                   res$rules$rule)
  expect_identical(predict(res, data.frame(X3 = numeric())), character())
})

test_that("binary and survival outcomes are split on their difference scale", {
  skip_if_not_installed("grf")
  skip_on_cran()
  d <- tree_data()
  withr::local_seed(3)
  d$yb <- rbinom(nrow(d), 1, plogis(-0.5 + 0.5 * d$X2 +
                                      d$z * 2 * (d$X3 > 0.6)))
  rb <- get_hte_tree(d, "z", c("X1", "X2", "X3"), surv = "yb", max_depth = 1,
                     tree_args = list(n_boot = 200L),
                     grf_args = list(num.trees = 500L))
  expect_identical(attr(rb, "analysis")$outcome_type, "binary")
  expect_identical(rb$nodes$variable[1L], "X3")
  expect_true(all(abs(rb$rules$estimate) <= 1))

  ev <- rexp(nrow(d), 0.05 * exp(0.3 * d$X2 - d$z * 2 * (d$X3 > 0.6)))
  cens <- rexp(nrow(d), 0.02)
  d$time <- pmin(ev, cens)
  d$DSS <- as.integer(ev <= cens)
  rs <- get_hte_tree(d, "z", c("X1", "X2", "X3"), time = 12, max_depth = 1,
                     tree_args = list(n_boot = 200L),
                     grf_args = list(num.trees = 500L))
  expect_s3_class(rs$est$fit, "causal_survival_forest")
  expect_true(all(c("time", "DSS", ".arm", ".dr_score", ".rule") %in%
                    names(partykit::data_party(rs$tree))))
})

test_that("split_frac = 1 grows and estimates on every patient", {
  skip_if_not_installed("grf")
  d <- tree_data(n = 800L)
  res <- tree_call(d, split_frac = 1L, max_depth = 1,
                   tree_args = list(n_boot = 100L))
  expect_identical(nrow(res$est$data), nrow(d))
  expect_identical(res$nodes$n_disc, res$nodes$n)
  expect_output(print(res), "not honest")
})

test_that("get_hte_tree() validates its input and restores the RNG", {
  skip_if_not_installed("grf")
  d <- tree_data(n = 300L)
  expect_error(tree_call(d, tree_args = list(xval = 5)), "unknown field")
  expect_error(tree_call(d, method = "rpart_dr", alpha = 0.1), "alpha")
  expect_error(tree_call(d, method = "policy", max_depth = Inf), "finite")
  expect_error(tree_call(d, max_depth = 0), "max_depth")
  expect_error(tree_call(d, min_leaf = 0.5), "min_leaf")
  expect_error(tree_call(d, split_frac = 0), "split_frac")
  expect_error(tree_call(d, tree_args = list(n_boot = 0L)), "alpha = 1")
  expect_error(tree_call(d, method = "ctree_dr",
                         tree_args = list(testtype = "none")), "testtype")
  expect_error(tree_call(d, method = "ctree"), "should be one of")
  expect_error(get_hte_tree(d, "z", "X1", surv = "y",
                            grf_args = list(W.hat = rep(0.5, 300))), "W.hat")
  expect_error(get_hte_tree(d, "z", "X1", surv = FALSE), "competing")
  expect_error(get_hte_tree(d, "z", "X1", surv = "y", time = 5), "time")
  expect_error(tree_call(d, estimator = "iptw"), "should be one of")
  expect_error(get_hte_tree(transform(d, time = 1, DSS = 1L), "z", "X1",
                            estimator = "tmle"), "AIPW only")
  d$X2[1] <- NA
  expect_error(tree_call(d), "complete split variables")
  set.seed(5)
  before <- .Random.seed
  res <- tree_call(tree_data(n = 400L), max_depth = 1,
                   tree_args = list(n_boot = 50L))
  expect_identical(.Random.seed, before)
})

test_that("estimator = \"tmle\" estimates the same tree's nodes by grf's TMLE", {
  skip_if_not_installed("grf")
  d <- tree_data()
  aipw <- tree_call(d, method = "mob_dr", max_depth = 1)
  tmle <- tree_call(d, method = "mob_dr", max_depth = 1, estimator = "tmle")
  expect_identical(tmle$rules$rule, aipw$rules$rule)
  f <- tmle$est$fit
  for (i in seq_len(nrow(tmle$rules))) {
    a <- grf::average_treatment_effect(
      f, method = "TMLE", subset = which(tmle$est$data$.rule == tmle$rules$rule[i]))
    expect_equal(tmle$rules$estimate[i], a[["estimate"]])
    expect_equal(tmle$rules$std.error[i], a[["std.err"]])
  }
  expect_equal(tmle$nodes$estimate[1L],
               grf::average_treatment_effect(f, method = "TMLE")[["estimate"]])
  expect_false(isTRUE(all.equal(tmle$rules$estimate, aipw$rules$estimate)))
  expect_identical(attr(tmle, "analysis")$estimator, "tmle")
  expect_identical(attr(aipw, "analysis")$estimator, "aipw")
  expect_identical(attr(tmle$est, "analysis")$estimator, "tmle")
  expect_equal(tmle$est$stats$estimate, tmle$nodes$estimate[1L])
  expect_equal(tmle$est$subgroup$estimate, tmle$rules$estimate)
  expect_output(print(tmle), "ATE difference by TMLE, 95% CI", fixed = TRUE)
  expect_output(print(aipw), "(ATE difference, 95% CI)", fixed = TRUE)
})

test_that("TMLE tree effects reach subgroup and CATE plots", {
  skip_if_not_installed("forestplot")
  skip_if_not_installed("RegR")
  skip_if_not_installed("ggplotify")
  res <- tree_call(tree_data(n = 800L), method = "rpart_dr", max_depth = 1,
                   estimator = "tmle", tree_args = list(xval = 0L))
  withr::local_pdf(tempfile(fileext = ".pdf"))
  p <- suppressMessages(plt_hte_sub(res$est, sub_var = ".rule",
                                    fixed_size = FALSE))
  expect_equal(attr(p, "subgroup")$estimate, res$rules$estimate)
  q <- plt_hte_cate(res$est, sub_var = ".rule", type = "density")
  expect_equal(attr(q, "subgroup")$estimate, res$rules$estimate)
})

test_that("a single known propensity in grf_args reaches the forests", {
  skip_if_not_installed("grf")
  d <- tree_data(n = 800L)
  withr::local_seed(3)
  d$z <- rbinom(nrow(d), 1, 0.5)
  d$y <- d$X2 + d$z * (0.5 + 1.5 * (d$X3 > 0.6)) + rnorm(nrow(d))
  res <- get_hte_tree(d, "z", c("X1", "X2", "X3"), surv = "y",
                      method = "rpart_dr", max_depth = 1,
                      grf_args = list(num.trees = 200L, W.hat = 0.5))
  expect_identical(res$nodes$variable[1L], "X3")
  expect_equal(res$est$fit$W.hat, rep(0.5, nrow(res$est$data)))
  expect_error(get_hte_tree(d, "z", "X1", surv = "y",
                            grf_args = list(W.hat = c(0.4, 0.6))), "W.hat")
})

test_that("compatible HTE trees reuse both forests and reject changed inputs", {
  d <- tree_data(n = 800L)
  fit <- function(method, ...) get_hte_tree(d, "z", c("X1", "X2", "X3", "grp"),
    surv = "y", method = method, max_depth = 1,
    grf_args = list(num.trees = 100L), ...)
  first <- fit("rpart_dr", tree_args = list(xval = 0L))
  fresh <- fit("ctree_dr")
  local_mocked_bindings(get_hte = function(...) stop("unexpected forest fit"),
                        .package = "causalR")
  shared <- fit("ctree_dr", reuse = first)
  expect_equal(shared$rules, fresh$rules)
  expect_identical(shared$est$fit, first$est$fit)
  expect_identical(attr(shared, "analysis")$forests_reused,
                   c(discovery = TRUE, estimation = TRUE))
  model <- fit("mob_abs", reuse = first)
  expect_identical(model$est$fit, first$est$fit)
  expect_identical(attr(model, "analysis")$forests_reused,
                   c(discovery = FALSE, estimation = TRUE))
  expect_equal(fit("ctree_dr", reuse = model)$rules, fresh$rules)
  expect_error(fit("ctree_dr", reuse = first, seed = 124), "same data.*settings")
  d$y[1L] <- d$y[1L] + 1
  expect_error(fit("ctree_dr", reuse = first), "same data.*settings")
})

test_that("tree stability bootstraps fixed scores and full mode refits", {
  d <- tree_data(n = 600L)
  res <- get_hte_tree(d, "z", c("X1", "X2", "X3", "grp"), surv = "y",
    method = "rpart_dr", max_depth = 1, tree_args = list(xval = 0L),
    grf_args = list(num.trees = 100L, num.threads = 2L))
  withr::local_seed(55)
  before <- .Random.seed
  a <- get_hte_tree_stability(res, d, n_rep = 3L, mode = "full")
  expect_identical(.Random.seed, before)
  expect_identical(a$replicates$status, rep("ok", 3L))
  expect_identical(nrow(a$agreement), 3L)
  expect_true(all(a$agreement$jaccard >= 0 & a$agreement$jaccard <= 1))
  local_mocked_bindings(get_hte = function(...) stop("unexpected forest fit"),
                        .package = "causalR")
  b <- get_hte_tree_stability(res, d, n_rep = 3L)
  expect_identical(b$replicates$status, rep("ok", 3L))
  expect_identical(.Random.seed, before)
  expect_equal(b, get_hte_tree_stability(res, d, n_rep = 3L))
  expect_setequal(b$variable_frequency$variable, c("X1", "X2", "X3", "grp"))
  expect_true(all(b$variable_frequency$frequency >= 0 &
                   b$variable_frequency$frequency <= 1))
  expect_true(all(b$cutpoints$variable %in% b$variable_frequency$variable))
  expect_error(get_hte_tree_stability(res, d, n_rep = 1), "n_rep")
  d$y[1L] <- 0
  expect_error(get_hte_tree_stability(res, d, n_rep = 3L), "original data")
})

test_that("stability Jaccard ignores node labels and handles unsplit trees", {
  d <- tree_data(n = 600L)
  res <- get_hte_tree(d, "z", c("X1", "X2", "X3"), surv = "y",
    method = "rpart_dr", max_depth = 1, tree_args = list(xval = 0L),
    grf_args = list(num.trees = 100L, num.threads = 2L))
  permuted <- res
  permuted$rules$node <- rev(permuted$rules$node)
  local_mocked_bindings(get_hte_tree = function(...) permuted,
                        .package = "causalR")
  st <- get_hte_tree_stability(res, d, n_rep = 3L, mode = "full")
  expect_equal(st$agreement$jaccard, rep(1, 3L))
  expect_equal(st$reference_agreement$jaccard, rep(1, 3L))
  ctx <- attr(res, "analysis")$reuse
  ctx$scores[] <- 1
  attr(res, "analysis")$reuse <- ctx
  flat <- get_hte_tree_stability(res, d, n_rep = 3L)
  expect_equal(flat$agreement$jaccard, rep(1, 3L))
  expect_equal(flat$variable_frequency$frequency, rep(0, 3L))
  expect_identical(nrow(flat$cutpoints), 0L)
})

test_that("fixed-score stability supports R-learner and node-model trees", {
  d <- tree_data(n = 800L)
  fits <- lapply(c("maxt", "rpart_r", "ctree_r", "mob_abs"), function(method)
    get_hte_tree(d, "z", c("X1", "X2", "X3"), surv = "y",
      method = method, max_depth = 1,
      tree_args = if (method == "maxt") list(n_boot = 50L) else list(),
      grf_args = list(num.trees = 100L, num.threads = 2L)))
  local_mocked_bindings(get_hte = function(...) stop("unexpected forest fit"),
                        .package = "causalR")
  for (res in fits) {
    st <- get_hte_tree_stability(res, d, n_rep = 2L)
    expect_identical(st$replicates$status, rep("ok", 2L))
  }
})

test_that("stability excludes failed fits and restores the RNG on error", {
  d <- tree_data(n = 600L)
  res <- get_hte_tree(d, "z", c("X1", "X2", "X3"), surv = "y",
    method = "rpart_dr", max_depth = 1, tree_args = list(xval = 0L),
    grf_args = list(num.trees = 100L, num.threads = 2L))
  engine <- .tree_engine
  calls <- 0L
  local_mocked_bindings(.tree_engine = function(...) {
    calls <<- calls + 1L
    if (calls == 1L) stop("test fit failure")
    engine(...)
  }, .package = "causalR")
  expect_warning(st <- get_hte_tree_stability(res, d, n_rep = 3L),
                  "1 of 3 stability fits failed")
  expect_identical(st$replicates$status, c("failed", "ok", "ok"))
  expect_identical(nrow(st$agreement), 1L)
  expect_equal(st$variable_frequency$frequency,
                st$variable_frequency$n_selected / 2)
  withr::local_seed(8)
  before <- .Random.seed
  local_mocked_bindings(.tree_engine = function(...) stop("test fit failure"),
                        .package = "causalR")
  expect_error(get_hte_tree_stability(res, d, n_rep = 2L), "All stability fits failed")
  expect_identical(.Random.seed, before)
})
