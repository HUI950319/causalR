# plt_hte_rate() evaluates how well a ranking of the patients targets the
# treatment effect (RATE, TOC and Qini curves) from a get_hte() result. grf is
# a Suggests; every test that needs a forest skips without it.

rate_cache <- new.env()

# Survival: only marker > 0 modifies the effect, age is prognostic.
rate_surv <- function() {
  skip_if_not_installed("grf")
  if (is.null(rate_cache$surv)) {
    set.seed(20260924)
    n <- 800L
    d <- data.frame(age    = round(stats::runif(n, 30, 85)),
                    marker = stats::rnorm(n),
                    stage  = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
    d$z <- stats::rbinom(n, 1, 0.5)
    haz <- exp(-4 + 0.03 * (d$age - 60) - d$z * (0.1 + 0.9 * (d$marker > 0)))
    tt  <- stats::rexp(n, haz)
    cc  <- stats::runif(n, 20, 150)
    d$time <- pmin(tt, cc)
    d$DSS  <- as.integer(tt <= cc)
    rate_cache$surv <- get_hte(d, "z", adj_var = c("age", "marker", "stage"),
                               time = 60,
                               grf_args = list(num.trees = 300, seed = 1))
  }
  rate_cache$surv
}

rate_cont <- function() {
  skip_if_not_installed("grf")
  if (is.null(rate_cache$cont)) {
    set.seed(3)
    n <- 600L
    d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                    risk = stats::runif(n))
    d$z <- stats::rbinom(n, 1, 0.5)
    d$y <- d$x2 + d$z * (1 + d$x1) + stats::rnorm(n)
    rate_cache$cont <- get_hte(d, "z", adj_var = c("x1", "x2"), surv = "y",
                               grf_args = list(num.trees = 300, seed = 1))
  }
  rate_cache$cont
}

rate_of <- function(p) attr(p, "rate")
gates_of <- function(p) attr(p, "gates")


test_that("get_hte() keeps the grf arguments that plt_hte_rate() refits with", {
  res <- rate_surv()
  ga  <- attr(res, "analysis")$grf_args
  expect_identical(ga$num.trees, 300)
  expect_identical(ga$seed, 1)
  expect_identical(ga$target, "survival.probability")
  expect_identical(names(formals(plt_hte_rate)),
                   c("x", "priority", "type", "smooth", "gates_args",
                     "conf_level", "train_frac", "seed", "title", "save"))
})

test_that("a pre-specified rule is evaluated on every patient by grf's RATE", {
  res <- rate_surv()
  p <- plt_hte_rate(res, priority = c("marker", "age"))
  r <- rate_of(p)
  expect_s3_class(p, "ggplot")
  expect_identical(names(r), c("rule", "target", "estimate", "std.error",
                               "conf.low", "conf.high", "p.value", "n"))
  expect_identical(r$rule, rep(c("marker", "age", "marker - age"), 2L))
  expect_identical(r$target, rep(c("AUTOC", "QINI"), each = 3L))
  expect_identical(r$n, rep(nrow(res$data), 6L))
  pr <- data.frame(marker = res$data$marker, age = res$data$age)
  for (tg in c("AUTOC", "QINI")) {
    want <- grf::rank_average_treatment_effect(res$fit, pr, target = tg)
    expect_equal(r$estimate[r$target == tg], unname(want$estimate))
  }
  z <- stats::qnorm(0.975)
  expect_equal(r$conf.low, r$estimate - z * r$std.error)
  expect_equal(r$p.value, 2 * stats::pnorm(-abs(r$estimate / r$std.error)))
  # marker drives the effect, age does not
  expect_lt(r$p.value[r$rule == "marker" & r$target == "QINI"], 0.05)
})

test_that("the forest CATE is learnt on one split and evaluated on the other", {
  res <- rate_surv()
  p <- plt_hte_rate(res, seed = 7)
  r <- rate_of(p)
  expect_identical(r$rule, c("cate", "cate"))

  # the same split and refits by hand, on the original follow-up times
  n <- nrow(res$data)
  set.seed(7)
  train <- sort(sample.int(n, floor(0.5 * n)))
  ev    <- setdiff(seq_len(n), train)
  X     <- res$fit$X.orig
  refit <- function(rows)
    grf::causal_survival_forest(X[rows, ], res$data$time[rows],
                                res$fit$W.orig[rows], res$data$DSS[rows],
                                horizon = 60, target = "survival.probability",
                                num.trees = 300, seed = 7)
  prio <- stats::predict(refit(train), X[ev, ])$predictions
  want <- grf::rank_average_treatment_effect(refit(ev), data.frame(cate = prio),
                                             target = "AUTOC")
  expect_equal(r$estimate[r$target == "AUTOC"], unname(want$estimate))
  expect_identical(r$n, rep(length(ev), 2L))

  # reproducible, and seed, train_frac move it
  expect_identical(rate_of(plt_hte_rate(res, seed = 7)), r)
  expect_false(identical(rate_of(plt_hte_rate(res, seed = 8))$estimate,
                         r$estimate))
  expect_equal(rate_of(plt_hte_rate(res, seed = 7, train_frac = 0.7))$n[1L],
               n - floor(0.7 * n))
  # seed = NULL takes the forest's seed
  expect_identical(rate_of(plt_hte_rate(res)),
                   rate_of(plt_hte_rate(res, seed = 1)))
})

test_that("the global random number stream is left untouched", {
  res <- rate_surv()
  set.seed(99)
  before <- .Random.seed
  invisible(plt_hte_rate(res, priority = "marker"))
  expect_identical(.Random.seed, before)
})

test_that("the forest CATE and a covariate are compared on the held-out half", {
  res <- rate_cont()
  expect_message(p <- plt_hte_rate(res, priority = c("cate", "risk")),
                 "not a forest covariate")
  r <- rate_of(p)
  expect_identical(r$rule, rep(c("cate", "risk", "cate - risk"), 2L))
  expect_equal(r$n, rep(nrow(res$data) - floor(0.5 * nrow(res$data)), 6L))
  expect_lt(r$p.value[r$rule == "cate - risk" & r$target == "AUTOC"], 0.05)
})

test_that("the panels hold the TOC and q times the TOC", {
  res <- rate_surv()
  p <- plt_hte_rate(res, priority = "marker")
  pd <- p$data
  expect_identical(levels(pd$panel), c("TOC", "Qini"))
  toc  <- pd[pd$panel == "TOC", ]
  qini <- pd[pd$panel == "Qini", ]
  expect_equal(qini$y, qini$q * toc$y)
  expect_equal(range(toc$q), c(0.05, 1))
  expect_match(p$labels$y, "S(60)", fixed = TRUE)
  expect_match(p$labels$caption, "all 800 patients")
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomRibbon"),
                         logical(1L))))
  # 1% of the axis beside the curves, as in plt_hte_cate()
  expect_equal(p$scales$get_scales("x")$expand, ggplot2::expansion(mult = 0.01))

  one <- plt_hte_rate(res, priority = "marker", type = "qini")
  expect_identical(levels(one$data$panel), "Qini")
  expect_false(grepl("AUTOC", one$labels$subtitle))
  expect_identical(names(attr(one, "plot_size")), c("width", "height"))

  pc <- plt_hte_rate(rate_cont(), type = "toc")
  expect_match(pc$labels$y, "mean", fixed = TRUE)
  expect_match(pc$labels$caption, "held-out")
})

test_that("smooth draws LOESS curves per rule and leaves the RATE unchanged", {
  res <- rate_surv()
  raw <- plt_hte_rate(res, priority = c("marker", "age"))
  sm  <- plt_hte_rate(res, priority = c("marker", "age"), smooth = 0.2)
  expect_identical(rate_of(sm), rate_of(raw))
  expect_identical(plt_hte_rate(res, priority = c("marker", "age"),
                                smooth = 0)$data, raw$data)

  toc_raw <- raw$data[raw$data$panel == "TOC", ]
  toc_sm  <- sm$data[sm$data$panel == "TOC", ]
  for (r in levels(toc_raw$rule)) {
    i  <- toc_raw$rule == r
    lo <- function(y) as.numeric(stats::predict(stats::loess(
      y ~ q, data = data.frame(q = toc_raw$q[i], y = y), span = 0.2)))
    expect_equal(toc_sm$y[i], lo(toc_raw$y[i]))
    expect_equal(toc_sm$conf.low[i], lo(toc_raw$conf.low[i]))
    expect_equal(toc_sm$conf.high[i], lo(toc_raw$conf.high[i]))
  }
  qini_sm <- sm$data[sm$data$panel == "Qini", ]
  expect_equal(qini_sm$y, qini_sm$q * toc_sm$y)
  expect_equal(qini_sm$conf.low, qini_sm$q * toc_sm$conf.low)
  expect_match(sm$labels$caption, "LOESS span 0.2", fixed = TRUE)
  expect_false(grepl("LOESS", raw$labels$caption))

  expect_no_warning(plt_hte_rate(res, priority = "marker", smooth = 0.05))
  for (bad in list(0.01, 1.5, -0.1, NA, "a", c(0.1, 0.2)))
    expect_error(plt_hte_rate(res, priority = "marker", smooth = bad),
                 "`smooth`")
})

test_that("GATES of a pre-specified rule are grf's subset ATEs per fifth", {
  res <- rate_surv()
  p <- plt_hte_rate(res, priority = "marker", type = "gates")
  g <- gates_of(p)
  expect_s3_class(p, "ggplot")
  expect_identical(names(g), c("rule", "group", "q_from", "q_to", "n",
                               "estimate", "std.error", "conf.low",
                               "conf.high", "p.value", "cate_mean"))
  expect_identical(g$rule, rep("marker", 6L))
  expect_identical(g$group, c(as.character(1:5), "1 - 5"))
  expect_identical(g$n, c(rep(160L, 5L), 320L))
  expect_equal(g$q_from, c((0:4) / 5, NA))
  expect_equal(g$q_to, c((1:5) / 5, NA))
  expect_true(all(is.na(g$cate_mean)))

  # the fifths by hand, highest marker first, on the stored forest
  m   <- res$data$marker
  grp <- cut(-m, stats::quantile(-m, (0:5) / 5), include.lowest = TRUE,
             labels = FALSE)
  for (k in 1:5) {
    want <- grf::average_treatment_effect(res$fit, subset = which(grp == k))
    expect_equal(g$estimate[k], unname(want[["estimate"]]))
    expect_equal(g$std.error[k], unname(want[["std.err"]]))
  }
  # the top minus the bottom fifth, the two groups taken as independent
  d  <- g$estimate[1L] - g$estimate[5L]
  se <- sqrt(g$std.error[1L]^2 + g$std.error[5L]^2)
  z  <- stats::qnorm(0.975)
  expect_equal(g$estimate[6L], d)
  expect_equal(g$std.error[6L], se)
  expect_equal(g$conf.low, g$estimate - z * g$std.error)
  expect_equal(g$p.value[6L], 2 * stats::pnorm(-abs(d / se)))
  expect_gt(d, 0)                        # marker > 0 carries the effect

  expect_match(p$labels$subtitle, "GATES top - bottom group: marker",
               fixed = TRUE)
  expect_false(grepl("AUTOC", p$labels$subtitle))
  expect_match(p$labels$caption, "GATES bars: 95% CI", fixed = TRUE)
  expect_false(grepl("diamonds|shaded", p$labels$caption))
  expect_identical(attr(p, "plot_size")[["width"]], 5.5)
  expect_null(gates_of(plt_hte_rate(res, priority = "marker")))
})

test_that("GATES of the forest CATE are estimated on the held-out half", {
  res <- rate_surv()
  g <- gates_of(plt_hte_rate(res, seed = 7, type = "gates"))
  expect_identical(g$rule, rep("cate", 6L))
  expect_identical(g$n, c(rep(80L, 5L), 160L))

  # the split and refits of the RATE test, on the original follow-up times
  n <- nrow(res$data)
  set.seed(7)
  train <- sort(sample.int(n, floor(0.5 * n)))
  ev    <- setdiff(seq_len(n), train)
  X     <- res$fit$X.orig
  refit <- function(rows)
    grf::causal_survival_forest(X[rows, ], res$data$time[rows],
                                res$fit$W.orig[rows], res$data$DSS[rows],
                                horizon = 60, target = "survival.probability",
                                num.trees = 300, seed = 7)
  prio <- stats::predict(refit(train), X[ev, ])$predictions
  fit  <- refit(ev)
  grp  <- cut(-prio, stats::quantile(-prio, (0:5) / 5), include.lowest = TRUE,
              labels = FALSE)
  for (k in 1:5) {
    want <- grf::average_treatment_effect(fit, subset = which(grp == k))
    expect_equal(g$estimate[k], unname(want[["estimate"]]))
    expect_equal(g$cate_mean[k], mean(prio[grp == k]))
  }
  expect_true(all(diff(g$cate_mean[1:5]) < 0))
  expect_equal(g$cate_mean[6L], g$cate_mean[1L] - g$cate_mean[5L])
})

test_that("adding GATES leaves the RATE alone and draws a third panel", {
  res  <- rate_cont()
  two  <- suppressMessages(plt_hte_rate(res, priority = c("cate", "risk"),
                                        seed = 7))
  all3 <- suppressMessages(plt_hte_rate(res, priority = c("cate", "risk"),
                                        seed = 7,
                                        type = c("toc", "qini", "gates")))
  expect_identical(rate_of(all3), rate_of(two))
  expect_equal(all3$data$y, two$data$y)
  g <- gates_of(all3)
  expect_identical(g$rule, rep(c("cate", "risk"), each = 6L))
  expect_true(all(is.na(g$cate_mean[g$rule == "risk"])))
  expect_false(anyNA(g$cate_mean[g$rule == "cate"]))

  b <- ggplot2::ggplot_build(all3)
  expect_identical(as.character(b$layout$layout$panel),
                   c("TOC", "Qini", "GATES"))
  expect_identical(attr(all3, "plot_size")[["width"]], 13.5)
  expect_match(all3$labels$subtitle, "GATES top - bottom group: Forest CATE",
               fixed = TRUE)
  expect_match(all3$labels$caption, "diamonds: mean forest CATE",
               fixed = TRUE)
  # the diamonds belong to the forest CATE alone
  dia <- Filter(function(l) identical(l$aes_params$shape, 5), all3$layers)
  expect_length(dia, 1L)
  expect_identical(unique(as.character(dia[[1L]]$data$rule)), "Forest CATE")
})

test_that("tied patients share a GATES group and a constant rule is refused", {
  res <- rate_cont()
  res$data$pos  <- res$data$x1 > 0
  res$data$flat <- 1
  g <- gates_of(suppressMessages(plt_hte_rate(res, priority = "pos",
                                              type = "gates")))
  expect_identical(g$group, c("1", "2", "1 - 2"))
  expect_identical(g$n, c(sum(res$data$pos), sum(!res$data$pos),
                          nrow(res$data)))
  expect_equal(g$q_to[1L], mean(res$data$pos))
  expect_error(suppressMessages(plt_hte_rate(res, priority = "flat",
                                             type = "gates")),
               "single GATES group")
})

test_that("gates_args sets the number of groups and is checked", {
  res <- rate_cont()
  g <- gates_of(plt_hte_rate(res, priority = "x1", type = "gates",
                             gates_args = list(n_groups = 4)))
  expect_identical(g$group, c(as.character(1:4), "1 - 4"))
  expect_identical(g$n, c(rep(150L, 4L), 300L))
  expect_error(plt_hte_rate(res, priority = "x1",
                            gates_args = list(n_groups = 4)),
               "only applies")
  expect_error(plt_hte_rate(res, priority = "x1", type = "gates",
                            gates_args = list(k = 4)),
               "unknown field")
  for (bad in list(1, 2.5, NA, "5", c(3, 4)))
    expect_error(plt_hte_rate(res, priority = "x1", type = "gates",
                              gates_args = list(n_groups = bad)),
                 "`gates_args$n_groups`", fixed = TRUE)
})

test_that("each half needs patients followed past `time` in both arms", {
  skip_if_not_installed("grf")
  set.seed(5)
  n <- 300L
  d <- data.frame(x = stats::rnorm(n))
  d$z <- stats::rbinom(n, 1, 0.5)
  d$time <- stats::runif(n, 1, 100)
  # only one treated patient is followed past 60
  late <- which(d$z == 1 & d$time > 60)
  d$time[late[-1L]] <- stats::runif(length(late) - 1L, 1, 59)
  d$DSS <- stats::rbinom(n, 1, 0.5)
  res <- suppressWarnings(get_hte(d, "z", adj_var = "x", time = 60,
                                  grf_args = list(num.trees = 100, seed = 1)))
  expect_error(plt_hte_rate(res), "is followed beyond `time` = 60")
  # a pre-specified rule uses every patient and needs no split
  expect_no_error(suppressWarnings(plt_hte_rate(res, priority = "x")))
})

test_that("missing priorities are left out and odd rules are refused", {
  res <- rate_surv()
  res$data$marker_na <- replace(res$data$marker, 1:10, NA)
  expect_message(p <- plt_hte_rate(res, priority = "marker_na"),
                 "10 patients")
  expect_identical(rate_of(p)$n[1L], nrow(res$data) - 10L)
  expect_message(plt_hte_rate(res, priority = "marker_na"),
                 "not a forest covariate")

  expect_error(plt_hte_rate(list()), "hte_res")
  expect_error(plt_hte_rate(res, priority = ".cate"), "training split")
  expect_error(plt_hte_rate(res, priority = "z"), "treatment or the outcome")
  expect_error(plt_hte_rate(res, priority = "time"), "treatment or the outcome")
  expect_error(plt_hte_rate(res, priority = "nope"), "nope")
  expect_error(plt_hte_rate(res, priority = "stage"), "numeric")
  expect_error(plt_hte_rate(res, priority = c("cate", "age", "marker")),
               "one or two")
  expect_error(plt_hte_rate(res, priority = "age", train_frac = 0.6),
               "train_frac")
  expect_error(plt_hte_rate(res, train_frac = 1), "train_frac")
  expect_error(plt_hte_rate(res, conf_level = 1.5), "conf_level")
  expect_error(plt_hte_rate(res, seed = "a"), "seed")
  expect_error(plt_hte_rate(res, type = "roc"), "should be one of")
  expect_error(plt_hte_rate(res, save = "a.pdf"), "`save`")
})

test_that("a result made before get_hte() kept its grf arguments still refits", {
  res <- rate_cont()
  attr(res, "analysis")$grf_args <- NULL
  r <- rate_of(plt_hte_rate(res, seed = 2))
  expect_equal(r$n[1L], nrow(res$data) - floor(0.5 * nrow(res$data)))
})

test_that("save writes one PDF and returns the plot unchanged", {
  res <- rate_cont()
  expect_s3_class(plt_hte_rate(res, priority = "x1", save = list()), "ggplot")
  expect_s3_class(plt_hte_rate(res, priority = "x1", save = NULL), "ggplot")
  skip_if_not_installed("RegR")
  f <- tempfile(fileext = ".pdf")
  p <- plt_hte_rate(res, priority = "x1", save = list(filename = f))
  expect_true(file.exists(f))
  expect_s3_class(p, "ggplot")
  unlink(f)
})
