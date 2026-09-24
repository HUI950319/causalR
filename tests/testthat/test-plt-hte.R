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
  expect_identical(imp$df[imp$variable == "age"], 2L)
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
  # at the default spline df the strip repeats the p_het of $importance
  p_age <- res$importance$p_het[res$importance$variable == "age"]
  expect_identical(strip_of(num),
                   sprintf("age (p_het %s)", if (p_age < 0.001) "< 0.001"
                           else sprintf("= %.3f", p_age)))
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

test_that("cate_smooth sets the loess span of the cate line; 0 leaves it out", {
  res <- dep_res()
  expect_identical(names(formals(plt_hte_dep))[6:8],
                   c("ylim", "cate_smooth", "dr_args"))
  span_of <- function(q) {
    for (l in q$layers)
      if (inherits(l$geom, "GeomSmooth")) return(l$stat_params$span)
    NULL
  }
  expect_identical(span_of(plt_hte_dep(res, x_var = "age", display = "cate")), 0.6)
  expect_identical(span_of(plt_hte_dep(res, x_var = "age", display = "cate",
                                       cate_smooth = 0.3)), 0.3)
  p0 <- plt_hte_dep(res, x_var = "age", display = "cate", cate_smooth = 0)
  expect_null(span_of(p0))
  expect_identical(nrow(layer_data_of(p0, "GeomPoint")), nrow(res$data))
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
  fit <- stats::lm(.dr_score ~ splines::ns(age, df = 2),
                   data = dd[dd$age >= lo & dd$age <= hi, ])
  V   <- sandwich::vcovHC(fit, type = "HC3")
  b   <- stats::coef(fit)[-1]
  expect_equal(imp$p_het[imp$variable == "age"],
               stats::pchisq(drop(b %*% solve(V[-1, -1], b)), 2,
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
  expect_error(plt_hte_dep(res, x_var = c("age", "sex"), type = "heat",
                           cate_smooth = 0.3), "`cate_smooth` only applies")
  expect_error(plt_hte_dep(res, conf_level = 0), "conf_level")
  for (bad in list(0.01, 1.5, NA, "a", c(0.3, 0.5)))
    expect_error(plt_hte_dep(res, cate_smooth = bad), "`cate_smooth` must")
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

# Page size of a PDF in inches, from its MediaBox, and its number of pages.
pdf_size <- function(f) {
  box <- grepRaw("/MediaBox \\[[^]]*\\]", readBin(f, "raw", file.size(f)),
                 value = TRUE)
  as.numeric(strsplit(trimws(gsub("[^0-9. ]", "", rawToChar(box))), " +")[[1]][3:4]) / 72
}
pdf_pages <- function(f)
  length(grepRaw("/Type /Page[^s]", readBin(f, "raw", file.size(f)), all = TRUE))

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

  # plt_hte_cate(): that level keeps its curve but gets no subgroup line
  expect_warning(pc <- plt_hte_cate(res, sub_var = "stage", type = "density"),
                 "stage = III: no patient in one arm")
  expect_identical(as.character(layer_data_of(pc, "GeomVline", "estimate")$level),
                   c("I", "II"))
  expect_identical(nlevels(layer_data_of(pc, "GeomDensity")$level), 3L)
  expect_identical(pc$labels$x, "CATE: S(60) 1 - 0")
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
  # the page is the pinned size unless `save` names a dimension (PDF pages
  # are whole points, hence the tolerance)
  expect_equal(pdf_size(f), unname(attr(p, "plot_size")), tolerance = 0.01)
  f2 <- tempfile(fileext = ".pdf")
  plt_hte_sub(res, sub_var = "stage", save = list(filename = f2, width = 10))
  expect_equal(pdf_size(f2), c(10, attr(p, "plot_size")[["height"]]),
               tolerance = 0.01)
  unlink(c(f, f2))
})

test_that("fixed_size pins the plot to the size its text needs", {
  skip_if_not_installed("forestplot")
  res <- dep_res()
  expect_identical(tail(names(formals(plt_hte_sub)), 2L), c("fixed_size", "save"))

  # The table as drawn on a w x h page, in inches: its column widths and row
  # heights, the area forestplot gives it and the centre of that area relative
  # to the centre of the page.
  drawn <- function(p, w, h) {
    grDevices::pdf(NULL, width = w, height = h)
    on.exit(grDevices::dev.off())
    print(p)
    grid::seekViewport("BaseGrid")
    lay <- grid::current.viewport()$layout
    mid <- grid::deviceLoc(grid::unit(0.5, "npc"), grid::unit(0.5, "npc"),
                           valueOnly = TRUE)
    list(widths  = grid::convertWidth(lay$widths, "in", valueOnly = TRUE),
         heights = grid::convertHeight(lay$heights, "in", valueOnly = TRUE),
         area    = c(grid::convertWidth(grid::unit(1, "npc"), "in", TRUE),
                     grid::convertHeight(grid::unit(1, "npc"), "in", TRUE)),
         centre  = c(mid$x - w / 2, mid$y - h / 2))
  }

  p <- plt_hte_sub(res, sub_var = "stage")
  expect_s3_class(p, "hte_forestplot")
  size <- attr(p, "plot_size")
  expect_named(size, c("width", "height"))
  expect_equal(size, round(size, 1))
  small <- drawn(p, size[["width"]], size[["height"]])
  expect_equal(small, drawn(p, 14, 10), tolerance = 1e-6)
  expect_equal(small$heights, rep(0.3, 6L), tolerance = 1e-6)
  # the size is rounded up, so the table always fits its area
  expect_true(sum(small$widths) <= small$area[1] + 1e-6)
  expect_true(sum(small$widths) > small$area[1] - 0.1)

  # a width in inches: the graph column takes exactly what the text leaves
  p9 <- plt_hte_sub(res, sub_var = "stage", fixed_size = 9)
  expect_identical(attr(p9, "plot_size")[["width"]], 9)
  d9 <- drawn(p9, 12, 8)
  expect_equal(sum(d9$widths), d9$area[1], tolerance = 1e-6)
  expect_equal(d9$area[1] - small$area[1], 9 - size[["width"]], tolerance = 1e-6)
  expect_error(plt_hte_sub(res, sub_var = "stage", fixed_size = 2), "no room")
  for (bad in list("yes", c(5, 6), -1, NA))
    expect_error(plt_hte_sub(res, fixed_size = bad), "`fixed_size`")

  # FALSE stretches the table over the page, as before
  pf <- plt_hte_sub(res, sub_var = "stage", fixed_size = FALSE)
  expect_identical(pf$lineheight, "auto")
  expect_named(attr(pf, "plot_size"), c("width", "height"))
  expect_false(isTRUE(all.equal(drawn(pf, 8, 5)$heights,
                                drawn(pf, 14, 10)$heights)))

  # every print starts its own page instead of drawing over the last plot
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  plot(1)
  print(p)
  print(pf)
  grDevices::dev.off()
  expect_identical(pdf_pages(f), 3L)
  unlink(f)
})


# ---- plt_hte_cate() ---------------------------------------------------------

test_that("plt_hte_cate draws the sorted CATE with the overall and subgroup ATE", {
  res <- dep_res()
  expect_identical(names(formals(plt_hte_cate)),
                   c("x", "sub_var", "type", "show_ci", "conf_level",
                     "overall", "title", "save"))
  p <- plt_hte_cate(res, sub_var = "stage")
  expect_s3_class(p, "ggplot")
  bars <- layer_data_of(p, "GeomCol")
  expect_identical(nrow(bars), nrow(res$data))
  expect_false(is.unsorted(bars$cate))
  expect_equal(sort(bars$cate), sort(res$data$.cate))
  expect_identical(levels(bars$level), c("I", "II", "III"))
  expect_identical(strip_of(p), "stage")

  # dashed lines: the overall ATE, and the doubly robust subgroup ATEs that
  # plt_hte_sub() draws
  ate <- res$stats$estimate[res$stats$measure == "diff"]
  expect_equal(layer_data_of(p, "GeomHline", "ate")$ate, ate)
  sub <- attr(p, "subgroup")
  for (lv in c("I", "II", "III")) {
    ref <- grf::average_treatment_effect(res$fit, subset = res$data$stage == lv)
    expect_equal(sub$estimate[sub$level == lv], unname(ref[["estimate"]]))
  }
  expect_equal(layer_data_of(p, "GeomHline", "estimate")$estimate, sub$estimate)
  expect_null(layer_data_of(plt_hte_cate(res, sub_var = "stage", overall = FALSE),
                            "GeomHline", "ate"))

  # density: one curve per level and the same lines, now vertical
  pd <- plt_hte_cate(res, sub_var = "stage", type = "density")
  expect_equal(sort(layer_data_of(pd, "GeomDensity")$cate), sort(res$data$.cate))
  expect_equal(layer_data_of(pd, "GeomVline", "estimate")$estimate, sub$estimate)
  expect_equal(layer_data_of(pd, "GeomVline", "ate")$ate, ate)
  # little room left and right of the bars and curves
  expect_equal(p$scales$get_scales("x")$expand, ggplot2::expansion(mult = 0.01))
  expect_equal(pd$scales$get_scales("x")$expand, ggplot2::expansion(mult = 0.01))

  # no sub_var: one ungrouped panel; several: one panel each
  p0 <- plt_hte_cate(res)
  expect_null(attr(p0, "subgroup"))
  expect_null(layer_data_of(p0, "GeomHline", "estimate"))
  expect_false("level" %in% names(layer_data_of(p0, "GeomCol")))
  pm <- plt_hte_cate(res, sub_var = c("stage", "sex"))
  expect_identical(vapply(panels_of(pm), strip_of, ""), c("stage", "sex"))
  expect_named(attr(pm, "plot_size"), c("width", "height"))
})

test_that("plt_hte_cate adds grf's pointwise intervals and checks its input", {
  res <- dep_res()
  p  <- plt_hte_cate(res, sub_var = "stage", show_ci = TRUE, conf_level = 0.9)
  ci <- layer_data_of(p, "GeomLinerange")
  pr <- stats::predict(res$fit, estimate.variance = TRUE)
  o  <- order(as.numeric(pr$predictions))
  expect_equal(ci$cate, as.numeric(pr$predictions)[o])
  expect_equal(ci$conf.high - ci$cate,
               stats::qnorm(0.95) * sqrt(as.numeric(pr$variance.estimates))[o])
  expect_null(layer_data_of(plt_hte_cate(res, sub_var = "stage"), "GeomLinerange"))

  expect_error(plt_hte_cate(res, type = "density", show_ci = TRUE), "only applies")
  expect_error(plt_hte_cate(res, type = "density", conf_level = 0.9), "only applies")
  expect_error(plt_hte_cate(list()), "hte_res")
  expect_error(plt_hte_cate(res, sub_var = "nope"), "nope")
  expect_error(plt_hte_cate(res, sub_var = "age"), "continuous")
  expect_error(plt_hte_cate(res, show_ci = NA), "show_ci")
  expect_error(plt_hte_cate(res, conf_level = 1), "conf_level")
  expect_error(plt_hte_cate(res, save = "a.pdf"), "`save`")

  skip_if_not_installed("RegR")
  f <- tempfile(fileext = ".pdf")
  expect_s3_class(plt_hte_cate(res, sub_var = "stage", save = list(filename = f)),
                  "ggplot")
  expect_true(file.exists(f))
  unlink(f)
})
