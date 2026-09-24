# plt_hte_dep() draws from a get_hte() result. grf and sandwich are Suggests;
# every test that needs a forest skips without them.

dep_cache <- new.env()

dep_res <- function() {
  skip_if_not_installed("grf")
  skip_if_not_installed("sandwich")
  if (is.null(dep_cache$res)) {
    set.seed(20260923)
    n <- 500L
    d <- data.frame(age   = round(stats::runif(n, 20, 85)),
                    sex   = factor(sample(c("F", "M"), n, replace = TRUE)),
                    stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
    d$z <- stats::rbinom(n, 1, 0.5)
    d$y <- stats::rbinom(n, 1, stats::plogis(-1 + 0.02 * (d$age - 50) +
                                               d$z * (0.2 + 0.6 * (d$sex == "M"))))
    # sex is added to the covariates with a message
    dep_cache$res <- suppressMessages(get_hte(
      d, cat_var = "z", sub_var = "sex", adj_var = c("age", "stage"),
      surv = "y", grf_args = list(num.trees = 300, seed = 1)))
  }
  dep_cache$res
}

# A patchwork keeps its last panel as the object itself.
panels_of <- function(p)
  if (inherits(p, "patchwork")) c(p$patches$plots, list(p)) else list(p)

strip_of <- function(q) {
  for (l in q$layers)
    if (is.data.frame(l$data) && "panel" %in% names(l$data))
      return(as.character(l$data$panel[1]))
  NA_character_
}

layer_data_of <- function(q, geom, col = NULL) {
  for (l in q$layers)
    if (inherits(l$geom, geom) && (is.null(col) || col %in% names(l$data)))
      return(l$data)
  NULL
}

# Forest CATE averaged over every row with the stage columns set to one level.
manual_pdp <- function(fit, set) {
  X <- fit$X.orig
  for (nm in names(set)) X[, nm] <- set[[nm]]
  mean(stats::predict(fit, X)$predictions)
}


test_that(".hte_dr_var matches grf subsets and a hand-built spline test", {
  res <- dep_res()
  w <- res$fit$W.orig
  z <- stats::qnorm(0.975)

  st <- .hte_dr_var(res$data, "stage", w, z, 3L)
  expect_identical(st$type, "categorical")
  expect_identical(st$df, 2L)
  for (lv in c("I", "II", "III")) {
    ref <- grf::average_treatment_effect(res$fit, subset = res$data$stage == lv)
    row <- st$levels[st$levels$level == lv, ]
    expect_equal(row$estimate, unname(ref[["estimate"]]))
    expect_equal(row$std.error, unname(ref[["std.err"]]))
  }

  age <- .hte_dr_var(res$data, "age", w, z, 3L)
  fit <- stats::lm(.dr_score ~ splines::ns(age, df = 3), data = res$data)
  V   <- sandwich::vcovHC(fit, type = "HC3")
  b   <- stats::coef(fit)[-1]
  expect_identical(age$type, "continuous")
  expect_identical(age$df, 3L)
  expect_equal(age$p_het, stats::pchisq(drop(b %*% solve(V[-1, -1], b)), 3,
                                        lower.tail = FALSE))
})

test_that("get_hte() reports p_het in $importance, equal to the subgroup p_inter", {
  res <- dep_res()
  imp <- res$importance
  expect_named(imp, c("variable", "importance", "n_col", "df", "p_het"))
  expect_equal(imp$p_het[imp$variable == "sex"],
               res$subgroup$p_inter[res$subgroup$measure == "diff"][1])
  expect_identical(imp$df[imp$variable == "age"], 3L)
  expect_identical(imp$df[imp$variable == "stage"], 2L)
})

test_that("plt_hte_dep draws one panel per covariate in importance order", {
  res <- dep_res()
  p <- plt_hte_dep(res)
  expect_s3_class(p, "patchwork")
  strips <- vapply(panels_of(p), strip_of, "")
  expect_identical(sub(" \\(p_het.*", "", strips), res$importance$variable)
  expect_true(all(grepl("(p_het ", strips, fixed = TRUE)))
  expect_identical(names(attr(p, "plot_size")), c("width", "height"))

  fct <- vapply(panels_of(plt_hte_dep(res, x_var = "fct")), strip_of, "")
  expect_identical(sub(" \\(p_het.*", "", fct),
                   setdiff(res$importance$variable, "age"))
  num <- plt_hte_dep(res, x_var = "num")
  expect_false(inherits(num, "patchwork"))
  expect_match(strip_of(num), "^age ")
  expect_identical(strip_of(plt_hte_dep(res, x_var = "stage", display = "cate")),
                   "stage")
})

test_that("the dr layer is the doubly robust summary; pdp averages forest predictions", {
  res <- dep_res()
  p <- plt_hte_dep(res, x_var = "stage", display = c("dr", "pdp"))

  dr <- .hte_dr_var(res$data, "stage", res$fit$W.orig, stats::qnorm(0.975), 3L)
  expect_equal(layer_data_of(p, "GeomPointrange")$estimate, dr$levels$estimate)

  expect_identical(names(formals(plt_hte_dep))[4:6],
                   c("display", "conf_level", "ylim"))
  p90  <- plt_hte_dep(res, x_var = "stage", display = "dr", conf_level = 0.9)
  dr90 <- .hte_dr_var(res$data, "stage", res$fit$W.orig, stats::qnorm(0.95), 3L)
  expect_equal(layer_data_of(p90, "GeomPointrange")$conf.low,
               dr90$levels$conf.low)
  expect_match(p90$labels$caption, "90% CI", fixed = TRUE)

  pdp <- layer_data_of(p, "GeomPoint", "estimate")
  expect_equal(pdp$estimate[pdp$x == "II"],
               manual_pdp(res$fit, list(stageI = 0, stageII = 1, stageIII = 0)))

  cont <- plt_hte_dep(res, x_var = "age", display = "pdp",
                      pdp_args = list(grid_n = 5))
  line <- layer_data_of(cont, "GeomLine", "estimate")
  expect_identical(nrow(line), 5L)
  expect_equal(line$estimate[1], manual_pdp(res$fit, list(age = min(res$data$age))))
})

test_that("type = 'heat' tiles the two-way partial dependence", {
  res <- dep_res()
  p <- plt_hte_dep(res, x_var = c("sex", "stage"), type = "heat")
  expect_false(inherits(p, "patchwork"))
  expect_identical(nrow(p$data), 6L)
  cell <- p$data[p$data$sex == "M" & p$data$stage == "III", ]
  expect_equal(cell$estimate,
               manual_pdp(res$fit, list(sexF = 0, sexM = 1, stageI = 0,
                                        stageII = 0, stageIII = 1)))
})

test_that("pdp sets an integer-coded factor through its level code", {
  skip_if_not_installed("grf")
  set.seed(20260923)
  n <- 500L
  d <- data.frame(age   = round(stats::runif(n, 20, 85)),
                  stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
  d$z <- stats::rbinom(n, 1, 0.5)
  d$y <- stats::rbinom(n, 1, stats::plogis(-1 + d$z * (0.2 + 0.4 * (d$stage == "III"))))
  res <- get_hte(d, cat_var = "z", adj_var = c("age", "stage"), surv = "y",
                 factor_encoding = "integer",
                 grf_args = list(num.trees = 300, seed = 1))

  pdp <- layer_data_of(plt_hte_dep(res, x_var = "stage", display = "pdp"),
                       "GeomPoint", "estimate")
  expect_equal(pdp$estimate[pdp$x == "II"], manual_pdp(res$fit, list(stage = 2)))
  heat <- plt_hte_dep(res, x_var = c("age", "stage"), type = "heat",
                      pdp_args = list(grid_n = 3))
  cell <- heat$data[heat$data$stage == "III" & heat$data$age == min(d$age), ]
  expect_equal(cell$estimate,
               manual_pdp(res$fit, list(age = min(d$age), stage = 3)))
})

test_that("a covariate with missing values draws only its observed patients", {
  skip_if_not_installed("grf")
  skip_if_not_installed("sandwich")
  skip_if_not_installed("patchwork")
  set.seed(11)
  n <- 400L
  d <- data.frame(age   = stats::runif(n, 20, 85),
                  stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
  d$age[1:40]    <- NA
  d$stage[41:80] <- NA
  d$z <- stats::rbinom(n, 1, 0.5)
  d$y <- stats::rbinom(n, 1, 0.3 + 0.1 * d$z)
  res <- get_hte(d, cat_var = "z", adj_var = c("age", "stage"), surv = "y",
                 grf_args = list(num.trees = 200, seed = 1))
  expect_identical(nrow(res$data), n)

  p <- plt_hte_dep(res, display = c("cate", "dr", "pdp"))
  for (q in panels_of(p)) expect_no_warning(ggplot2::ggplot_build(q))
  pts <- layer_data_of(plt_hte_dep(res, x_var = "stage", display = "cate"),
                       "GeomPoint")
  expect_identical(nrow(pts), n - 40L)
  expect_false(anyNA(pts$x))
})

test_that("survival: levels and values no arm follows past `time` leave the dr layer", {
  skip_if_not_installed("grf")
  skip_if_not_installed("sandwich")
  skip_if_not_installed("patchwork")
  set.seed(20260924)
  n <- 800L
  d <- data.frame(age   = stats::runif(n, 20, 85),
                  stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
  d$z  <- stats::rbinom(n, 1, 0.5)
  ev   <- stats::rexp(n, 0.02)
  cens <- pmin(stats::rexp(n, 0.01), 120)
  d$time <- pmin(ev, cens)
  d$DSS  <- as.integer(ev <= cens)
  # no treated patient in stage III or over 70 is followed past 60
  cut <- d$z == 1 & (d$stage == "III" | d$age > 70) & d$time > 50
  d$time[cut] <- 50
  d$DSS[cut]  <- 0L
  expect_warning(
    res <- get_hte(d, "z", sub_var = "stage", adj_var = c("age", "stage"),
                   surv = TRUE, time = 60,
                   grf_args = list(num.trees = 300, seed = 1)),
    "stage = III: no patient in one arm")
  imp <- res$importance

  # stage III is NA in $subgroup, so it leaves p_het and the dr layer as well
  expect_equal(imp$p_het[imp$variable == "stage"],
               res$subgroup$p_inter[res$subgroup$measure == "diff"][1])
  expect_identical(imp$df[imp$variable == "stage"], 1L)
  pr <- layer_data_of(plt_hte_dep(res, x_var = "stage", display = "dr"),
                      "GeomPointrange")
  expect_identical(as.character(pr$x), c("I", "II"))

  # the age spline is fit and drawn only where both arms reach past 60
  dd   <- res$data
  past <- dd$time > 60
  lo   <- max(min(dd$age[past & dd$z == 1]), min(dd$age[past & dd$z == 0]))
  hi   <- min(max(dd$age[past & dd$z == 1]), max(dd$age[past & dd$z == 0]))
  line <- layer_data_of(plt_hte_dep(res, x_var = "age", display = "dr"),
                        "GeomLine", "estimate")
  expect_equal(range(line$x), c(lo, hi))
  fit <- stats::lm(.dr_score ~ splines::ns(age, df = 3),
                   data = dd[dd$age >= lo & dd$age <= hi, ])
  V   <- sandwich::vcovHC(fit, type = "HC3")
  b   <- stats::coef(fit)[-1]
  expect_equal(imp$p_het[imp$variable == "age"],
               stats::pchisq(drop(b %*% solve(V[-1, -1], b)), 3,
                             lower.tail = FALSE))
})

test_that("invalid requests stop with a clear message", {
  res <- dep_res()
  expect_error(plt_hte_dep(list()), "hte_res")
  expect_error(plt_hte_dep(res, x_var = "nope"), "nope")
  expect_error(plt_hte_dep(res, x_var = "age", type = "heat"), "exactly two")
  expect_error(plt_hte_dep(res, x_var = c("age", "sex"), type = "heat",
                           display = "dr"), "`display`")
  expect_error(plt_hte_dep(res, x_var = c("age", "sex"), type = "heat",
                           conf_level = 0.9), "`conf_level`")
  expect_error(plt_hte_dep(res, conf_level = 0), "conf_level")
  expect_error(plt_hte_dep(res, dr_args = list(df = 3)), "unknown field")
  expect_error(plt_hte_dep(res, axis_arg = list(share_y = "var")), "share_y")
  expect_error(plt_hte_dep(res, save = "a.pdf"), "`save`")
})

test_that("save writes one PDF and returns the plot unchanged", {
  res <- dep_res()
  expect_s3_class(plt_hte_dep(res, x_var = "stage", save = list()), "ggplot")
  expect_s3_class(plt_hte_dep(res, x_var = "stage", save = NULL), "ggplot")
  skip_if_not_installed("RegR")
  f <- tempfile(fileext = ".pdf")
  p <- plt_hte_dep(res, x_var = "stage", save = list(filename = f))
  expect_true(file.exists(f))
  expect_s3_class(p, "ggplot")
  expect_false(inherits(p, "patchwork"))
})


# ---- plt_hte_sub() ----------------------------------------------------------

# The cells of one forestplot text column, header first, and every header.
fp_col <- function(p, j)
  vapply(p$labels[[j]], function(s) paste(s, collapse = ""), "")
fp_headers <- function(p) vapply(seq_along(p$labels), function(j) fp_col(p, j)[1], "")

test_that("plt_hte_sub recomputes the get_hte() subgroup estimates", {
  skip_if_not_installed("forestplot")
  res <- dep_res()
  devs <- grDevices::dev.list()
  p <- plt_hte_sub(res, sub_var = "sex")
  expect_identical(grDevices::dev.list(), devs)   # no stray device opened
  expect_s3_class(p, "gforge_forestplot")
  expect_equal(attr(p, "subgroup"), res$subgroup[res$subgroup$measure == "diff", ])

  # a covariate get_hte() was not asked about: grf's own subset estimates
  st <- attr(plt_hte_sub(res, sub_var = "stage"), "subgroup")
  for (lv in c("I", "II", "III")) {
    ref <- grf::average_treatment_effect(res$fit, subset = res$data$stage == lv)
    expect_equal(st$estimate[st$level == lv], unname(ref[["estimate"]]))
  }
  expect_equal(unique(st$p_inter),
               res$importance$p_het[res$importance$variable == "stage"])

  # conf_level sets the intervals and the header, not the estimates
  p90 <- plt_hte_sub(res, sub_var = "stage", conf_level = 0.9)
  s90 <- attr(p90, "subgroup")
  expect_equal(s90$estimate, st$estimate)
  expect_equal(s90$conf.low, s90$estimate - stats::qnorm(0.95) * s90$std.error)
  expect_identical(fp_headers(p90)[3], "Risk difference (90% CI)")
})

test_that("measure = 'ratio' averages the arm scores on a log axis", {
  skip_if_not_installed("forestplot")
  res <- dep_res()
  p   <- plt_hte_sub(res, sub_var = "stage", measure = "ratio")
  s   <- .hte_arm_scores(res$fit)
  idx <- res$data$stage == "II"
  sg  <- attr(p, "subgroup")
  expect_equal(sg$estimate[sg$level == "II"], mean(s$g1[idx]) / mean(s$g0[idx]))
  expect_true(p$xlog)
  expect_equal(p$zero, 0)                     # log(1): forestplot logs the axis
  expect_identical(fp_headers(p)[3], "Risk ratio (95% CI)")
})

test_that("overall, show_n, show_pvalue and show_pinter set the rows and columns", {
  skip_if_not_installed("forestplot")
  res <- dep_res()
  p <- plt_hte_sub(res, sub_var = c("stage", "sex"))
  expect_identical(trimws(fp_col(p, 1)),
                   c("Subgroup", "All patients", "stage", "I", "II", "III",
                     "sex", "F", "M"))
  expect_length(p$labels, 3L)                 # no P columns by default
  expect_equal(unname(p$estimates[2, 1, 1]),
               res$stats$estimate[res$stats$measure == "diff"])

  no_all <- plt_hte_sub(res, sub_var = "stage", overall = FALSE)
  expect_false("All patients" %in% trimws(fp_col(no_all, 1)))

  expect_identical(names(formals(plt_hte_sub))[1:8],
                   c("x", "sub_var", "measure", "conf_level", "overall",
                     "show_n", "show_pvalue", "show_pinter"))
  no_n <- plt_hte_sub(res, sub_var = "stage", show_n = FALSE)
  expect_identical(fp_headers(no_n), c("Subgroup", "Risk difference (95% CI)"))

  both <- plt_hte_sub(res, sub_var = "stage", show_pvalue = TRUE,
                      show_pinter = TRUE)
  expect_identical(fp_headers(both)[4:5], c("P", "P for interaction"))
  pint <- attr(both, "subgroup")$p_inter[1]
  expect_identical(fp_col(both, 5)[trimws(fp_col(both, 1)) == "stage"],
                   if (pint < 0.001) "<0.001" else sprintf("%.3f", pint))
})

test_that("sub_var defaults to the categorical covariates and checks columns", {
  skip_if_not_installed("forestplot")
  res <- dep_res()
  cov  <- attr(res, "analysis")$covariates
  labs <- trimws(fp_col(plt_hte_sub(res, overall = FALSE), 1))
  expect_identical(labs[labs %in% cov], setdiff(cov, "age"))

  expect_error(plt_hte_sub(res, sub_var = "nope"), "nope")
  expect_error(plt_hte_sub(res, sub_var = "age"), "continuous")
  res2 <- res
  res2$data$grp <- ifelse(res2$data$age > 50, "old", "young")
  expect_message(p <- plt_hte_sub(res2, sub_var = "grp"),
                 "not a forest covariate")
  expect_identical(attr(p, "subgroup")$level, c("old", "young"))
})

test_that("survival: a subgroup no arm follows past `time` is drawn empty", {
  skip_if_not_installed("grf")
  skip_if_not_installed("forestplot")
  set.seed(20260924)
  n <- 800L
  d <- data.frame(age   = stats::runif(n, 20, 85),
                  stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)))
  d$z  <- stats::rbinom(n, 1, 0.5)
  ev   <- stats::rexp(n, 0.02)
  cens <- pmin(stats::rexp(n, 0.01), 120)
  d$time <- pmin(ev, cens)
  d$DSS  <- as.integer(ev <= cens)
  cut <- d$z == 1 & d$stage == "III" & d$time > 50
  d$time[cut] <- 50
  d$DSS[cut]  <- 0L
  res <- get_hte(d, "z", adj_var = c("age", "stage"), surv = TRUE, time = 60,
                 grf_args = list(num.trees = 300, seed = 1))

  expect_warning(p <- plt_hte_sub(res, sub_var = "stage"),
                 "stage = III: no patient in one arm")
  row <- trimws(fp_col(p, 1)) == "III"
  expect_true(is.na(attr(p, "subgroup")$estimate[3]))
  expect_identical(fp_col(p, 3)[row], "\u2014")
  expect_true(is.na(p$estimates[row, 1, 1]))
  expect_identical(fp_headers(p)[3], "S(60) difference (95% CI)")
  expect_identical(
    fp_headers(suppressWarnings(plt_hte_sub(res, sub_var = "stage",
                                            measure = "ratio")))[3],
    "Event risk ratio (95% CI)")
})

test_that("plt_hte_sub rejects invalid requests and saves a PDF", {
  skip_if_not_installed("forestplot")
  res <- dep_res()
  expect_error(plt_hte_sub(list()), "hte_res")
  expect_error(plt_hte_sub(res, overall = NA), "overall")
  expect_error(plt_hte_sub(res, conf_level = 1.5), "conf_level")
  expect_error(plt_hte_sub(res, show_pvalue = "yes"), "show_pvalue")
  expect_error(plt_hte_sub(res, xlim = c(1, 0)), "xlim")
  expect_error(plt_hte_sub(res, measure = "ratio", xlim = c(-1, 2)), "positive")
  expect_error(plt_hte_sub(res, save = "a.pdf"), "`save`")

  set.seed(2)
  dc <- data.frame(x = stats::rnorm(300),
                   g = factor(sample(c("a", "b"), 300, replace = TRUE)))
  dc$z <- stats::rbinom(300, 1, 0.5)
  dc$y <- stats::rnorm(300) + dc$z
  rc <- get_hte(dc, "z", adj_var = c("x", "g"), surv = "y",
                grf_args = list(num.trees = 100, seed = 1))
  expect_error(plt_hte_sub(rc, measure = "OR"), "OR")

  expect_identical(names(attr(plt_hte_sub(res), "plot_size")),
                   c("width", "height"))
  skip_if_not_installed("RegR")
  f <- tempfile(fileext = ".pdf")
  p <- plt_hte_sub(res, sub_var = "stage", save = list(filename = f))
  expect_true(file.exists(f))
  expect_s3_class(p, "gforge_forestplot")
})
