# get_hte_pdp(), get_hte_ale() and get_hte_shp() explain a get_hte() forest
# without refitting it. grf, xgboost, shapviz and kernelshap are Suggests.

expl_cache <- new.env()

expl_res <- function(encoding = "integer") {
  skip_if_not_installed("grf")
  key <- paste0("res_", encoding)
  if (is.null(expl_cache[[key]])) {
    set.seed(20261006)
    n <- 400L
    d <- data.frame(age   = stats::runif(n, 20, 85),
                    stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)),
                    sex   = factor(sample(c("F", "M"), n, replace = TRUE)),
                    nodes = sample(0:3, n, replace = TRUE))
    d$z <- stats::rbinom(n, 1, 0.5)
    d$y <- 0.02 * d$age + d$z * (1 + 0.04 * (d$age - 50) +
                                   0.8 * (d$stage == "III")) +
      stats::rnorm(n, sd = 0.5)
    expl_cache[[key]] <- suppressMessages(get_hte(
      d, cat_var = "z", adj_var = c("age", "stage", "sex", "nodes"),
      surv = "y", factor_encoding = encoding,
      grf_args = list(num.trees = 200, seed = 1)))
  }
  expl_cache[[key]]
}

# Forest CATE of `rows` with design columns overwritten: a column name and one
# value, or one value per row.
pred_with <- function(fit, rows, set = list()) {
  X <- fit$X.orig[rows, , drop = FALSE]
  for (nm in names(set)) X[, nm] <- set[[nm]]
  as.numeric(stats::predict(fit, X)$predictions)
}


test_that("get_hte_pdp() matches .hte_pdp() and hand-set predictions", {
  res <- expl_res()
  all_rows <- seq_len(nrow(res$fit$X.orig))

  age <- get_hte_pdp(res, x_var = "age", grid_n = 5, max_n = Inf)
  ref <- .hte_pdp(res, "age", 5, Inf)
  # the grid follows the data: quantiles, not even steps over the range
  expect_equal(age$value, stats::quantile(res$data$age, c(0, .25, .5, .75, 1),
                                          names = FALSE))
  expect_equal(age$value, ref$age)
  expect_equal(age$estimate, ref$estimate)
  expect_true(all(is.na(age$level)))

  stage <- get_hte_pdp(res, x_var = "stage", max_n = Inf)
  expect_identical(stage$level, c("I", "II", "III"))
  expect_equal(stage$estimate[2], mean(pred_with(res$fit, all_rows,
                                                 list(stage = 2))))

  nodes <- get_hte_pdp(res, x_var = "nodes", max_n = Inf)
  expect_equal(nodes$value, 0:3)
  expect_equal(nodes$estimate[4], mean(pred_with(res$fit, all_rows,
                                                 list(nodes = 3))))

  oh <- expl_res("onehot")
  st <- get_hte_pdp(oh, x_var = "stage", max_n = Inf)
  expect_equal(st$estimate[2], mean(pred_with(oh$fit, all_rows,
    list(stageI = 0, stageII = 1, stageIII = 0))))
})

test_that("get_hte_pdp() batches covariates without changing their curves", {
  res <- expl_res()
  all <- get_hte_pdp(res, grid_n = 7, max_n = 60)
  expect_named(all, c("variable", "value", "level", "estimate"))
  expect_identical(unique(all$variable), res$importance$variable)
  for (v in c("age", "sex", "nodes"))
    expect_equal(all$estimate[all$variable == v],
                 get_hte_pdp(res, x_var = v, grid_n = 7, max_n = 60)$estimate)
  a <- attr(all, "analysis")
  expect_identical(a$method, "pdp")
  expect_identical(a$n, 60L)
  expect_identical(a$n_total, 400L)
})

test_that("get_hte_pdp() hands the forest one batch of jobs at a time", {
  res   <- expl_res()
  sizes <- integer()
  real  <- .hte_predict_set
  local_mocked_bindings(
    .hte_explain_cost = function(fit) list(a = 0, b = 0, batch = 500),
    .hte_predict_set  = function(x, row, ...) {
      sizes <<- c(sizes, length(row))
      real(x, row, ...)
    })
  get_hte_pdp(res, grid_n = 5, max_n = Inf)
  expect_lte(max(sizes), 500L)
  # age 5 grid points, stage 3 levels, sex 2, nodes 4
  expect_identical(sum(sizes), 400L * (5L + 3L + 2L + 4L))
})

test_that("time_budget sets the patient count and max_n overrides it", {
  res <- expl_res()
  expect_no_message(full <- get_hte_pdp(res, x_var = "age"))
  expect_identical(attr(full, "analysis")$n, 400L)
  expect_message(cut <- get_hte_pdp(res, x_var = "age", time_budget = 0.01),
                 "Explaining")
  expect_lt(attr(cut, "analysis")$n, 400L)
  expect_identical(attr(get_hte_pdp(res, x_var = "age", max_n = 30),
                        "analysis")$n, 30L)

  # The cost model alone: a bigger forest or more predictions per patient
  # admit fewer patients. A 10,000 x 2000-tree forest with 10 covariates
  # admits 2 patients of Kernel SHAP, as documented.
  forest <- function(n, trees, p)
    list(fit = list(`_num_trees` = trees, X.orig = matrix(0, n, p)))
  size <- function(...) suppressMessages(.hte_explain_size(...))$n
  big   <- size(forest(1e4, 2000, 10), NULL, 20, rows = 156)
  small <- size(forest(1e4, 200, 10), NULL, 20, rows = 156)
  expect_lt(big, small)
  expect_lt(size(forest(1e4, 2000, 10), NULL, 20, rows = 600), big)
  expect_identical(size(forest(1e4, 2000, 10), NULL, 20,
                        rows = 50 * (2 * (10 + 45) + 80), calls = 5), 2)
})

test_that("get_hte_ale() matches a hand computation", {
  res  <- expl_res()
  fit  <- res$fit
  rows <- .hte_explain_rows(400L, 40L)

  ale <- get_hte_ale(res, x_var = "age", n_bins = 4, max_n = 40)
  xv  <- res$data$age[rows]
  z   <- unique(stats::quantile(xv, seq(0, 1, length.out = 5), type = 1,
                                names = FALSE))
  K   <- length(z) - 1L
  j   <- pmax(1L, findInterval(xv, z, left.open = TRUE))
  dl  <- pred_with(fit, rows, list(age = z[j + 1L])) -
    pred_with(fit, rows, list(age = z[j]))
  st  <- vapply(seq_len(K), function(b) mean(dl[j == b]), numeric(1))
  nk  <- tabulate(j, K)
  A   <- c(0, cumsum(st))
  A   <- A - sum(nk * (A[-1] + A[-(K + 1)]) / 2) / sum(nk)
  expect_equal(ale$value, z)
  expect_equal(ale$ale, A)
  expect_identical(ale$n, c(NA_integer_, nk))

  cat_ale <- get_hte_ale(res, x_var = "stage", max_n = 40)
  k  <- as.integer(res$data$stage[rows])
  f0 <- pred_with(fit, rows)
  fu <- pred_with(fit, rows, list(stage = pmin(k + 1L, 3L)))
  fd <- pred_with(fit, rows, list(stage = pmax(k - 1L, 1L)))
  nk <- tabulate(k, 3L)
  st <- vapply(1:2, function(q) (sum((fu - f0)[k == q]) +
                                   sum((f0 - fd)[k == q + 1L])) /
                 (nk[q] + nk[q + 1L]), numeric(1))
  A  <- c(0, cumsum(st))
  expect_identical(cat_ale$level, c("I", "II", "III"))
  expect_equal(cat_ale$ale, A - sum(nk * A) / sum(nk))
  expect_identical(cat_ale$n, nk)

  # One-hot coding gives the same steps through the same levels.
  oh <- expl_res("onehot")
  oh_ale <- get_hte_ale(oh, x_var = c("stage", "nodes"), max_n = 40)
  expect_equal(sum(oh_ale$n[oh_ale$variable == "stage"] *
                     oh_ale$ale[oh_ale$variable == "stage"]), 0)
  expect_equal(oh_ale$value[oh_ale$variable == "nodes"], 0:3)
})

test_that("get_hte_shp() kernel SHAP is exact for up to 8 covariates", {
  skip_if_not_installed("kernelshap")
  skip_if_not_installed("shapviz")
  for (enc in c("integer", "onehot")) {
    res  <- expl_res(enc)
    sv   <- get_hte_shp(res, method = "kernel", max_n = 6, bg_n = 20)
    rows <- .hte_explain_rows(400L, 6L)
    expect_s3_class(sv, "shapviz")
    expect_identical(colnames(sv$S), c("age", "stage", "sex", "nodes"))
    expect_equal(unname(rowSums(sv$S)) + sv$baseline, pred_with(res$fit, rows),
                 tolerance = 1e-6)
    expect_true(is.factor(sv$X$stage))
    expect_identical(attr(sv, "analysis")$bg_n, 20)
  }
})

test_that("get_hte_shp() surrogate collapses one-hot columns and restores the RNG", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("shapviz")
  oh <- expl_res("onehot")
  set.seed(5)
  before <- .Random.seed
  sv <- get_hte_shp(oh, max_n = 300)
  expect_identical(.Random.seed, before)
  expect_s3_class(sv, "shapviz")
  expect_identical(dim(sv$S), c(300L, 4L))
  expect_identical(colnames(sv$S), c("age", "stage", "sex", "nodes"))
  expect_identical(levels(sv$X$stage), c("I", "II", "III"))
  a <- attr(sv, "analysis")
  expect_identical(a$method, "surrogate")
  expect_gt(a$r2, 0.8)
  # The CATE varies with age and stage only, so they carry the credit.
  imp <- colMeans(abs(sv$S))
  expect_identical(names(sort(imp, decreasing = TRUE))[1:2], c("age", "stage"))
})

test_that("explanations run on a survival forest", {
  skip_if_not_installed("grf")
  set.seed(1)
  n <- 300L
  d <- data.frame(age = stats::runif(n, 20, 85),
                  stage = factor(sample(c("I", "II"), n, replace = TRUE)))
  d$z <- stats::rbinom(n, 1, 0.5)
  ev <- stats::rexp(n, 0.02 * exp(-0.5 * d$z * (d$age > 50)))
  ce <- pmin(stats::rexp(n, 0.01), 120)
  d$time <- pmin(ev, ce)
  d$DSS  <- as.integer(ev <= ce)
  res <- suppressMessages(get_hte(d, cat_var = "z", adj_var = c("age", "stage"),
                                  surv = TRUE, time = 60,
                                  grf_args = list(num.trees = 100, seed = 1)))
  expect_true(all(is.finite(get_hte_pdp(res, max_n = 50, grid_n = 5)$estimate)))
  expect_true(all(is.finite(get_hte_ale(res, max_n = 50, n_bins = 5)$ale)))
})

test_that("invalid arguments are rejected", {
  res <- expl_res()
  expect_error(get_hte_pdp(list()), "hte_res")
  for (bad in list(0, 1.5, "a", NA))
    expect_error(get_hte_pdp(res, max_n = bad), "max_n")
  for (bad in list(0, NA, Inf, "a"))
    expect_error(get_hte_ale(res, time_budget = bad), "time_budget")
  expect_error(get_hte_pdp(res, grid_n = 1), "grid_n")
  expect_error(get_hte_ale(res, n_bins = 0), "n_bins")
  expect_error(get_hte_pdp(res, x_var = "nope"), "names no covariate")
  expect_error(get_hte_ale(res, verbose = "yes"), "verbose")
  expect_error(get_hte_shp(res, bg_n = 10), "only applies")
  expect_error(get_hte_shp(res, method = "kernel",
                           surrogate_args = list(nrounds = 5)), "only applies")
  expect_error(get_hte_shp(res, surrogate_args = list(foo = 1)), "unknown")
  expect_error(get_hte_shp(res, surrogate_args = list(learning_rate = 2)),
               "learning_rate")
})
