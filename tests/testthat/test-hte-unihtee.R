# get_hte_unihtee(): univariate projections of AIPW pseudo-outcomes. The glm
# route needs only stats + sandwich; forest tests skip without grf.

uni_args <- list(num.trees = 300, seed = 1)

# tau = w3 + 1{sex = M}, sd(w3) = 2: the TEM-VIP is 2 per SD of w3, 1 for
# M vs F and 0 for w1, w2 and w4.
uni_cont_data <- function(n = 800L, seed = 510) {
  set.seed(seed)
  d <- data.frame(w1 = stats::rnorm(n), w2 = stats::rnorm(n),
                  w3 = stats::rnorm(n, 0, 2), w4 = stats::rnorm(n),
                  sex = factor(sample(c("F", "M"), n, replace = TRUE)))
  d$a <- stats::rbinom(n, 1, stats::plogis(0.5 * d$w1))
  d$y <- d$w1 + d$w2 + d$a * (d$w3 + (d$sex == "M")) + stats::rnorm(n)
  d
}

uni_bin_data <- function(n = 1500L, seed = 511) {
  set.seed(seed)
  d <- data.frame(w1 = stats::rnorm(n), w2 = stats::rnorm(n),
                  sex = factor(sample(c("F", "M"), n, replace = TRUE)))
  d$a <- stats::rbinom(n, 1, stats::plogis(0.4 * d$w1))
  d$y <- stats::rbinom(n, 1, stats::plogis(-1 + 0.5 * d$w1 +
                                             d$a * (0.2 + 0.6 * d$w2)))
  d
}

test_that("glm route returns the documented table and recovers the truth", {
  skip_if_not_installed("sandwich")
  d   <- uni_cont_data()
  res <- get_hte_unihtee(d, cat_var = "a", adj_var = c("w1", "w2", "w3", "w4", "sex"),
                         surv = "y", method = "glm")
  expect_named(res, c("vip", "data"))
  expect_named(res$vip, c("variable", "type", "contrast", "n", "sd", "measure",
                          "estimate", "std.error", "conf.low", "conf.high",
                          "p.value", "p.adj"))
  expect_setequal(res$vip$variable, c("w1", "w2", "w3", "w4", "sex"))
  expect_equal(res$vip$p.adj, stats::p.adjust(res$vip$p.value, "BH"))
  expect_false(is.unsorted(res$vip$p.value))
  expect_setequal(res$vip$variable[1:2], c("w3", "sex"))

  v <- res$vip[match(c("w3", "sex"), res$vip$variable), ]
  expect_equal(v$type, c("continuous", "binary"))
  expect_equal(v$contrast, c("+1 SD", "M vs F"))
  expect_lt(abs(v$estimate[1] - 2), 0.3)
  expect_true(v$conf.low[2] < 1 && v$conf.high[2] > 1)

  # the estimate is the least-squares slope of the pseudo-outcome
  dd <- res$data
  expect_equal(v$estimate[1],
               unname(stats::coef(stats::lm(.score_diff ~ I(w3 / stats::sd(w3)),
                                            data = dd))[2]))
  expect_equal(v$estimate[2], mean(dd$.score_diff[dd$sex == "M"]) -
                 mean(dd$.score_diff[dd$sex == "F"]))
  g1 <- dd$.mu1 + dd$a / dd$.ps * (dd$y - dd$.mu1)
  g0 <- dd$.mu0 + (1 - dd$a) / (1 - dd$.ps) * (dd$y - dd$.mu0)
  expect_equal(dd$.score_diff, g1 - g0)
})

test_that("relative measures project log-scale pseudo-outcomes", {
  skip_if_not_installed("sandwich")
  d   <- uni_bin_data()
  res <- get_hte_unihtee(d, cat_var = "a", adj_var = c("w1", "w2", "sex"),
                         surv = "y", method = "glm",
                         measure = c("diff", "ratio", "OR"))
  expect_equal(unique(res$vip$measure), c("diff", "ratio", "OR"))
  for (m in c("diff", "ratio", "OR")) {
    r <- res$vip[res$vip$measure == m, ]
    expect_equal(r$p.adj, stats::p.adjust(r$p.value, "BH"))
    expect_equal(r$variable[1], "w2")
  }
  dd <- res$data
  g1 <- dd$.mu1 + dd$a / dd$.ps * (dd$y - dd$.mu1)
  g0 <- dd$.mu0 + (1 - dd$a) / (1 - dd$.ps) * (dd$y - dd$.mu0)
  expect_equal(dd$.score_ratio, log(dd$.mu1 / dd$.mu0) +
                 (g1 - dd$.mu1) / dd$.mu1 - (g0 - dd$.mu0) / dd$.mu0)
  r <- res$vip[res$vip$measure == "ratio" & res$vip$variable == "w2", ]
  b <- unname(stats::coef(stats::lm(.score_ratio ~ I(w2 / stats::sd(w2)), data = dd))[2])
  expect_equal(r$estimate, exp(b))
  expect_equal(c(r$conf.low, r$conf.high),
               exp(b + c(-1, 1) * stats::qnorm(0.975) * r$std.error))
})

test_that("grf route equals grf::best_linear_projection()", {
  skip_if_not_installed("grf")
  skip_if_not_installed("sandwich")
  d   <- uni_cont_data(n = 600L)
  adj <- c("w1", "w2", "w3", "w4", "sex")
  res <- get_hte_unihtee(d, cat_var = "a", adj_var = adj, surv = "y",
                         grf_args = uni_args)
  hte <- get_hte(d, cat_var = "a", adj_var = adj, surv = "y", grf_args = uni_args)
  expect_equal(res$data$.score_diff, hte$data$.dr_score)
  blp <- grf::best_linear_projection(hte$fit, A = cbind(w3 = d$w3 / stats::sd(d$w3)))
  v <- res$vip[res$vip$variable == "w3", ]
  expect_equal(v$estimate, unname(blp["w3", 1]))
  expect_equal(v$std.error, unname(blp["w3", 2]))
})

test_that("survival outcomes run through the forest on both scales", {
  skip_if_not_installed("grf")
  skip_if_not_installed("sandwich")
  set.seed(20260928)
  n <- 600
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                  grp = sample(c(TRUE, FALSE), n, replace = TRUE))
  d$z  <- stats::rbinom(n, 1, 0.5)
  ev   <- stats::rexp(n, 0.02 * exp(0.3 * d$x1 - d$z * (0.3 + 0.5 * d$x2)))
  cens <- pmin(stats::rexp(n, 0.01), 120)
  d$time <- pmin(ev, cens)
  d$DSS  <- as.integer(ev <= cens)
  res <- get_hte_unihtee(d, cat_var = "z", adj_var = c("x1", "x2", "grp"),
                         time = 60, measure = c("diff", "ratio"),
                         grf_args = uni_args)
  expect_equal(nrow(res$vip), 6L)
  expect_true(all(is.finite(res$vip$estimate)))
  expect_equal(res$vip$contrast[res$vip$variable == "grp"], rep("TRUE vs FALSE", 2))
  expect_true(all(c(".score_diff", ".score_ratio") %in% names(res$data)))
  skip_if_not_installed("patchwork")
  skip_if_not_installed("RegR")  # the default var_names
  p <- plt_hte_unihtee(res, type = "dep", x_var = c("x2", "grp"), measure = "ratio")
  expect_s3_class(p, "patchwork")
  expect_no_warning(for (q in c(p$patches$plots, list(p))) ggplot2::ggplot_build(q))
})

test_that("candidate_var narrows the table and joins the adjustment set", {
  skip_if_not_installed("sandwich")
  d   <- uni_cont_data(n = 400L)
  res <- get_hte_unihtee(d, cat_var = "a", candidate_var = "w3",
                         adj_var = c("w1", "w2"), surv = "y", method = "glm")
  expect_equal(res$vip$variable, "w3")
  expect_equal(attr(res, "analysis")$covariates, c("w1", "w2", "w3"))
})

test_that("the glm route leaves the global RNG untouched", {
  skip_if_not_installed("sandwich")
  d <- uni_cont_data(n = 300L)
  set.seed(1)
  before <- .Random.seed
  get_hte_unihtee(d, cat_var = "a", adj_var = c("w1", "w3"), surv = "y",
                  method = "glm")
  expect_identical(.Random.seed, before)
})

test_that("invalid input stops with a pointed message", {
  d <- uni_cont_data(n = 200L)
  d$stage <- factor(sample(c("I", "II", "III"), nrow(d), replace = TRUE))
  d$one   <- 1
  d$time  <- stats::rexp(nrow(d)); d$DSS <- 1L
  f <- function(...) get_hte_unihtee(d, cat_var = "a", surv = "y",
                                     method = "glm", ...)
  expect_error(f(adj_var = c("w1", "stage")), "two levels")
  expect_error(f(adj_var = c("w1", "one")), "constant")
  expect_error(f(adj_var = "w1", glm_args = list(nfold = 5)), "unknown")
  expect_error(f(adj_var = "w1", time = 60), "survival")
  expect_error(f(adj_var = "w1", candidate_var = "a"), "cat_var")
  expect_error(get_hte_unihtee(d, cat_var = "a", adj_var = "w1", method = "glm"),
               "method = \"grf\"")
})


# ---- plt_hte_unihtee() -------------------------------------------------------

uni_plot_res <- local({
  cache <- NULL
  function() {
    skip_if_not_installed("sandwich")
    skip_if_not_installed("RegR")  # the default var_names
    if (is.null(cache))
      cache <<- get_hte_unihtee(uni_cont_data(), cat_var = "a",
                                adj_var = c("w1", "w2", "w3", "w4", "sex"),
                                surv = "y", method = "glm")
    cache
  }
})

uni_layer <- function(q, geom, col) {
  for (l in q$layers)
    if (inherits(l$geom, geom) && is.data.frame(l$data) && col %in% names(l$data))
      return(l$data)
  NULL
}

test_that("bar draws every signed estimate with its interval, coloured by p.adj", {
  res <- uni_plot_res()
  v   <- res$vip
  p   <- plt_hte_unihtee(res)
  expect_s3_class(p, "ggplot")
  expect_identical(names(attr(p, "plot_size")), c("width", "height"))
  b <- uni_layer(p, "GeomCol", "estimate")
  i <- match(v$variable, b$variable)
  expect_equal(b$estimate[i], v$estimate)
  expect_equal(as.character(b$sig[i]),
               ifelse(v$p.adj < 0.05, ifelse(v$estimate > 0, "Positive", "Negative"), "n.s."))
  e <- uni_layer(p, "GeomErrorbar", "conf.low")
  expect_equal(e$conf.high[match(v$variable, e$variable)], v$conf.high)
  expect_identical(tail(names(formals(plt_hte_unihtee)), 1L), "save")
})

test_that("volcano plots -log10 p and labels only BH-significant candidates", {
  res <- uni_plot_res()
  v   <- res$vip
  p   <- plt_hte_unihtee(res, type = "volcano", sig_level = 0.01)
  pts <- uni_layer(p, "GeomPoint", "neglog10p")
  expect_equal(pts$neglog10p[match(v$variable, pts$variable)], -log10(v$p.value))
  expect_setequal(uni_layer(p, "GeomText", "variable")$variable,
                  v$variable[v$p.adj < 0.01])
})

test_that("dep draws the TEM-VIP line, a spline check and the binary contrast", {
  res <- uni_plot_res()
  v   <- res$vip
  q   <- plt_hte_unihtee(res, type = "dep", x_var = "w3")
  expect_false(inherits(q, "patchwork"))
  ln <- uni_layer(q, "GeomLine", "estimate")
  ln <- ln[order(ln$x), ]
  slope <- (ln$estimate[nrow(ln)] - ln$estimate[1]) / (ln$x[nrow(ln)] - ln$x[1])
  expect_equal(slope * v$sd[v$variable == "w3"], v$estimate[v$variable == "w3"])
  expect_false(is.null(uni_layer(q, "GeomLine", "spline")))
  pr <- uni_layer(plt_hte_unihtee(res, type = "dep", x_var = "sex"),
                  "GeomPointrange", "estimate")
  expect_equal(diff(pr$estimate), v$estimate[v$variable == "sex"])
  skip_if_not_installed("patchwork")
  expect_s3_class(plt_hte_unihtee(res, type = "dep", x_var = c("w3", "sex")),
                  "patchwork")
})

test_that("relative measures draw on the log scale", {
  skip_if_not_installed("sandwich")
  skip_if_not_installed("RegR")
  res <- get_hte_unihtee(uni_bin_data(), cat_var = "a", adj_var = c("w1", "w2", "sex"),
                         surv = "y", method = "glm", measure = c("diff", "ratio"))
  p <- plt_hte_unihtee(res, measure = "ratio")
  b <- uni_layer(p, "GeomCol", "estimate")
  r <- res$vip[res$vip$measure == "ratio", ]
  expect_equal(b$estimate[match(r$variable, b$variable)], r$estimate)
  expect_match(ggplot2::get_labs(p)$x, "ratio")
  q <- plt_hte_unihtee(res, type = "dep", x_var = "w2", measure = "ratio")
  expect_match(ggplot2::get_labs(q)$y, "log")
})

test_that("plt_hte_unihtee rejects invalid input", {
  res <- uni_plot_res()
  expect_error(plt_hte_unihtee(list(a = 1)), "get_hte_unihtee")
  expect_error(plt_hte_unihtee(res, measure = "OR"), "measure")
  expect_error(plt_hte_unihtee(res, x_var = "zz"), "x_var")
  expect_error(plt_hte_unihtee(res, sig_level = 2), "sig_level")
  expect_error(plt_hte_unihtee(res, save = "a.pdf"), "`save`")
  expect_error(plt_hte_unihtee(res, type = "dep", dr_args = list(df = 2)), "unknown")
  expect_error(plt_hte_unihtee(res, dr_args = list(spline_df = 3)), "type = \"dep\"")
})

test_that("save writes one PDF and returns the plot unchanged", {
  res <- uni_plot_res()
  expect_s3_class(plt_hte_unihtee(res, save = list()), "ggplot")
  expect_s3_class(plt_hte_unihtee(res, save = NULL), "ggplot")
  skip_if_not_installed("RegR")
  f <- tempfile(fileext = ".pdf")
  p <- plt_hte_unihtee(res, type = "volcano", save = list(filename = f))
  expect_true(file.exists(f))
  expect_s3_class(p, "ggplot")
})

test_that("dr_args$bins adds quantile-bin means of the pseudo-outcome", {
  res <- uni_plot_res()
  d   <- res$data
  q   <- plt_hte_unihtee(res, type = "dep", x_var = "w3", dr_args = list(bins = 4))
  bn  <- uni_layer(q, "GeomPointrange", "bin")
  expect_equal(bn$estimate, unname(c(tapply(
    d$.score_diff, cut(d$w3, stats::quantile(d$w3, 0:4 / 4), include.lowest = TRUE),
    mean))))
  expect_null(uni_layer(plt_hte_unihtee(res, type = "dep", x_var = "w3"),
                        "GeomPointrange", "bin"))
  expect_error(plt_hte_unihtee(res, type = "dep", dr_args = list(bins = 1)), "bins")
})

test_that("axis_arg$share_y = 'none' gives each dep panel its own range", {
  skip_if_not_installed("patchwork")
  res <- uni_plot_res()
  lim <- function(p) lapply(c(p$patches$plots, list(p)),
                            function(q) q$coordinates$limits$y)
  shared <- lim(plt_hte_unihtee(res, type = "dep", x_var = c("w3", "sex")))
  expect_equal(shared[[1]], shared[[2]])
  own <- lim(plt_hte_unihtee(res, type = "dep", x_var = c("w3", "sex"),
                             axis_arg = list(share_y = "none")))
  expect_false(isTRUE(all.equal(own[[1]], own[[2]])))
  expect_error(plt_hte_unihtee(res, axis_arg = list(share_y = "none")), "type = \"dep\"")
  expect_error(plt_hte_unihtee(res, type = "dep", axis_arg = list(share_y = "x")),
               "share_y")
  expect_identical(tail(names(formals(plt_hte_unihtee)), 4L),
                   c("axis_arg", "var_names", "title", "save"))
})

test_that("var_names relabels every type and rejects clashing labels", {
  res <- uni_plot_res()
  b <- uni_layer(plt_hte_unihtee(res, var_names = c(w3 = "Weight", sex = "Sex")),
                 "GeomCol", "estimate")
  expect_true(all(c("Weight (+1 SD)", "Sex (M vs F)", "w1 (+1 SD)") %in%
                    as.character(b$label)))
  b0 <- uni_layer(plt_hte_unihtee(res, var_names = NULL), "GeomCol", "estimate")
  expect_true("w3 (+1 SD)" %in% as.character(b0$label))
  txt <- uni_layer(plt_hte_unihtee(res, type = "volcano", var_names = c(w3 = "Weight")),
                   "GeomText", "display")
  expect_true("Weight" %in% txt$display)
  q <- plt_hte_unihtee(res, type = "dep", x_var = "w3", var_names = c(w3 = "Weight"))
  expect_match(uni_layer(q, "GeomLine", "estimate")$panel[1], "^Weight:")
  expect_identical(ggplot2::get_labs(q)$x, "Weight")
  expect_error(plt_hte_unihtee(res, var_names = c(w3 = "A", w4 = "A")),
               "same display label")
  expect_error(plt_hte_unihtee(res, var_names = list(w3 = 1)), "var_names")
})
