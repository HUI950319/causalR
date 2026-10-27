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
      surv = "y", factor_encoding = "onehot",
      grf_args = list(num.trees = 300, seed = 1)))
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
  # each panel names its covariate on the x axis too
  expect_identical(plt_hte_dep(res, x_var = "stage", display = "cate")$labels$x,
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
  colour_of <- function(q) {
    for (l in q$layers)
      if (inherits(l$geom, "GeomSmooth")) return(l$aes_params$colour)
    NULL
  }
  expect_identical(colour_of(plt_hte_dep(res, x_var = "age", display = "cate")),
                   "firebrick")
  expect_identical(colour_of(plt_hte_dep(res, x_var = "age")), "grey30")
  # every curve as thick as the curves of plt_hte_rate()
  all_lines <- plt_hte_dep(res, x_var = "age", display = c("cate", "dr", "pdp"))
  widths <- unlist(lapply(all_lines$layers, function(l)
    if (inherits(l$geom, c("GeomSmooth", "GeomLine"))) l$aes_params$linewidth))
  expect_identical(unname(widths), rep(0.8, 3L))
  p0 <- plt_hte_dep(res, x_var = "age", display = "cate", cate_smooth = 0)
  expect_null(span_of(p0))
  expect_identical(nrow(layer_data_of(p0, "GeomPoint")), nrow(res$data))
})

test_that("the cate layer draws 2000 patients; its line and y range use all", {
  skip_if_not_installed("grf")
  skip_if_not_installed("sandwich")
  set.seed(20260927)
  n <- 2500L
  d <- data.frame(age = round(stats::runif(n, 20, 85)),
                  sex = factor(sample(c("F", "M"), n, replace = TRUE)))
  d$z <- stats::rbinom(n, 1, 0.5)
  d$y <- stats::rbinom(n, 1, stats::plogis(-1 + d$z * (0.2 + 0.6 * (d$sex == "M"))))
  res <- get_hte(d, cat_var = "z", adj_var = c("age", "sex"), surv = "y",
                 grf_args = list(num.trees = 50, seed = 1))
  p <- plt_hte_dep(res, display = "cate")
  for (q in panels_of(p))
    expect_identical(nrow(layer_data_of(q, "GeomPoint")), 2000L)
  expect_match(p$patches$annotation$caption, "2000 of the 2500 patients",
               fixed = TRUE)

  age <- plt_hte_dep(res, x_var = "age", display = "cate")
  smooth <- Filter(function(l) inherits(l$geom, "GeomSmooth"), age$layers)[[1L]]
  expect_identical(nrow(smooth$data), n)
  # loess's approximate trace draws the same line without an n x n pass
  expect_identical(smooth$stat_params$method.args$control$trace.hat,
                   "approximate")
  lim <- age$coordinates$limits$y
  expect_lte(lim[1L], min(res$data$.cate))
  expect_gte(lim[2L], max(res$data$.cate))
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

test_that("PDP predicts bounded batches without changing grid estimates", {
  res <- dep_res()
  cells <- numeric()
  original <- stats::predict
  local_mocked_bindings(predict = function(object, newdata, ...) {
    cells <<- c(cells, length(newdata))
    original(object, newdata, ...)
  }, .package = "stats")
  # batches of 5000 patients instead of a million design cells
  local_mocked_bindings(.hte_explain_cost = function(fit)
    list(a = 0, b = 0, batch = 5000))
  p <- plt_hte_dep(res, x_var = c("age", "stage"), type = "heat",
                   pdp_args = list(grid_n = 121, max_n = Inf))
  # whole-year ages tie, so their 121 quantiles merge into fewer points
  g <- length(unique(stats::quantile(res$data$age, seq(0, 1, length.out = 121),
                                     names = FALSE)))
  expect_lt(g, 121L)
  expect_gt(length(cells), 1L)
  expect_lte(max(cells), 5000 * ncol(res$fit$X.orig))
  expect_equal(sum(cells), g * 3 * length(res$fit$X.orig))
  for (k in c(1L, g, g + 1L, 3L * g)) {
    row <- p$data[k, ]
    value <- manual_pdp(res$fit, list(age = row$age,
      stageI = as.integer(row$stage == "I"),
      stageII = as.integer(row$stage == "II"),
      stageIII = as.integer(row$stage == "III")))
    expect_equal(row$estimate, value)
  }
})

test_that("heat-map tiles of an uneven quantile grid meet without gaps", {
  res <- dep_res()
  p   <- plt_hte_dep(res, x_var = c("age", "stage"), type = "heat",
                     pdp_args = list(grid_n = 6))
  e   <- unique(ggplot2::layer_data(p)[, c("xmin", "xmax")])
  e   <- e[order(e$xmin), ]
  expect_equal(e$xmin[-1L], e$xmax[-nrow(e)])
  g <- sort(unique(p$data$age))
  expect_true(all(g >= e$xmin & g <= e$xmax))
  # the tiles cover the covariate's range and no more
  expect_equal(range(e), range(res$data$age))
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

test_that("dependence plots reject unused limits and invalid grid sizes", {
  res <- dep_res()
  expect_error(plt_hte_dep(res, x_var = c("age", "sex"), type = "heat",
                           ylim = c(-1, 1)), "only applies to type")
  for (bad in c(2.5, Inf, -Inf, NA_real_)) {
    expect_error(plt_hte_dep(res, dr_args = list(spline_df = bad)), "spline_df")
    expect_error(plt_hte_dep(res, pdp_args = list(grid_n = bad)), "grid_n")
  }
  for (bad in c(1.5, -Inf, NA_real_))
    expect_error(plt_hte_dep(res, pdp_args = list(max_n = bad)), "max_n")
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

# The grobs of one class in the forestplot drawing a plt_hte_sub() ggplot holds
# in its annotation layer, in drawing order.
fp_grobs <- function(p, cls)
  Filter(function(k) inherits(k, cls), p$layers[[2]]$geom_params$grob$children)

# The table as drawn, one vector per text column, header first: forestplot
# draws every cell in a viewport named Label_vp_<row>_<column>. A header that
# spans several table rows, down to the one holding "Subgroup", comes back as
# one cell, its lines joined by "\n".
fp_table <- function(p) {
  cells <- Filter(function(k) !is.null(k$vp) && startsWith(k$vp$name, "Label_vp_"),
                  fp_grobs(p, "text"))
  rc  <- vapply(cells, function(k)
    as.integer(strsplit(sub("Label_vp_", "", k$vp$name, fixed = TRUE), "_")[[1]]),
    integer(2L))
  lab <- vapply(cells, function(k) paste(k$label, collapse = ""), "",
                USE.NAMES = FALSE)
  cols <- lapply(sort(unique(rc[2L, ])), function(j) {
    v <- rep("", max(rc[1L, ]))
    v[rc[1L, rc[2L, ] == j]] <- lab[rc[2L, ] == j]
    v
  })
  k <- match("Subgroup", cols[[1L]])
  lapply(cols, function(v) {
    h <- v[seq_len(k)]
    c(paste(h[nzchar(h)], collapse = "\n"), v[-seq_len(k)])
  })
}
fp_col <- function(p, j) fp_table(p)[[j]]
fp_headers <- function(p) vapply(fp_table(p), `[`, "", 1L)

# Page size of a PDF in inches, from its MediaBox, and its number of pages.
# cairo_pdf(), save_plt()'s default device, keeps the page objects in
# Flate-compressed object streams, so the MediaBox is also sought in every
# stream once inflated.
pdf_size <- function(f) {
  b <- readBin(f, "raw", file.size(f))
  at <- grepRaw("stream\r?\n", b, all = TRUE)
  at <- at + ifelse(b[at + 6L] == as.raw(10L), 7L, 8L)
  parts <- c(list(b), lapply(at, function(i)
    tryCatch(memDecompress(b[i:length(b)], "gzip"), error = function(e) raw())))
  box <- Filter(length, lapply(parts, grepRaw,
                               pattern = "/MediaBox \\[[^]]*\\]", value = TRUE))[[1]]
  as.numeric(strsplit(trimws(gsub("[^0-9. ]", "", rawToChar(box))), " +")[[1]][3:4]) / 72
}
pdf_pages <- function(f)
  length(grepRaw("/Type /Page[^s]", readBin(f, "raw", file.size(f)), all = TRUE))

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

test_that("plt_hte_sub rejects relative effects for weighted or clustered forests", {
  res <- dep_res()
  for (field in c("sample.weights", "clusters")) {
    x <- res
    x$fit[[field]] <- if (field == "clusters") rep(1:50, each = 10) else
      rep(c(1, 2), length.out = nrow(x$data))
    for (measure in c("ratio", "OR"))
      expect_error(plt_hte_sub(x, sub_var = "sex", measure = measure),
                   "do not support.*clusters.*sample.weights")
  }
})

test_that("plt_hte_sub recomputes the get_hte() subgroup estimates", {
  skip_if_not_installed("forestplot")
  skip_if_not_installed("ggplotify")
  res <- dep_res()
  devs <- grDevices::dev.list()
  p <- plt_hte_sub(res, sub_var = "sex")
  expect_identical(grDevices::dev.list(), devs)   # no stray device opened
  expect_s3_class(p, "ggplot")
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
  expect_identical(fp_headers(p90)[3], "Risk difference, %")
})

test_that("measure = 'ratio' averages the arm scores on a log axis", {
  skip_if_not_installed("forestplot")
  skip_if_not_installed("ggplotify")
  res <- dep_res()
  p   <- plt_hte_sub(res, sub_var = "stage", measure = "ratio")
  s   <- .hte_arm_scores(res$fit)
  idx <- res$data$stage == "II"
  sg  <- attr(p, "subgroup")
  expect_equal(sg$estimate[sg$level == "II"], mean(s$g1[idx]) / mean(s$g0[idx]))
  # forestplot logs the axis, labelled on the ratio scale
  ax <- fp_grobs(p, "xaxis")[[1L]]
  expect_equal(ax$at, log(as.numeric(ax$label)))
  expect_true("1" %in% ax$label)
  expect_identical(fp_headers(p)[3], "Risk ratio")
})

test_that("overall, show_n, show_pvalue and show_pinter set the rows and columns", {
  skip_if_not_installed("forestplot")
  skip_if_not_installed("ggplotify")
  res <- dep_res()
  p <- plt_hte_sub(res, sub_var = c("stage", "sex"))
  expect_identical(trimws(fp_col(p, 1)),
                   c("Subgroup", "All patients", "stage", "I", "II", "III",
                     "sex", "F", "M"))
  expect_length(fp_table(p), 3L)              # no P columns by default
  ov <- res$stats[res$stats$measure == "diff", ]
  expect_identical(fp_col(p, 3)[2], sprintf("%.1f (%.1f, %.1f)", 100 * ov$estimate,
                                            100 * ov$conf.low, 100 * ov$conf.high))
  expect_length(fp_grobs(p, "polygon"), 6L)   # a diamond for every estimate
  # "All patients" is set like a variable name; the rest of its row is not
  font <- function(lab) {
    gp <- Filter(function(k) identical(k$label, lab), fp_grobs(p, "text"))[[1L]]$gp
    c(gp$font, gp$cex)
  }
  expect_identical(font("All patients"), font("stage"))
  expect_identical(font(fp_col(p, 3)[2]), font("   I"))

  no_all <- plt_hte_sub(res, sub_var = "stage", overall = FALSE)
  expect_false("All patients" %in% trimws(fp_col(no_all, 1)))

  expect_identical(names(formals(plt_hte_sub))[1:9],
                   c("x", "sub_var", "measure", "conf_level", "effect",
                     "overall", "show_n", "show_pvalue", "show_pinter"))
  no_n <- plt_hte_sub(res, sub_var = "stage", show_n = FALSE)
  expect_identical(fp_headers(no_n), c("Subgroup", "Risk difference, %"))

  both <- plt_hte_sub(res, sub_var = "stage", show_pvalue = TRUE,
                      show_pinter = TRUE)
  expect_identical(fp_headers(both)[4:5], c("P", "P for interaction"))
  pint <- attr(both, "subgroup")$p_inter[1]
  expect_identical(fp_col(both, 5)[trimws(fp_col(both, 1)) == "stage"],
                   if (pint < 0.001) "<0.001" else sprintf("%.3f", pint))
})

test_that("effect replaces the estimate of the rows it names", {
  skip_if_not_installed("forestplot")
  skip_if_not_installed("ggplotify")
  res <- dep_res()
  st  <- attr(plt_hte_sub(res, sub_var = c("stage", "sex")), "subgroup")
  p <- plt_hte_sub(res, sub_var = c("stage", "sex"), show_pvalue = TRUE,
                   effect = list(sex = c(M = "0.120 (0.050, 0.190)"),
                                 "All patients" = "0.100 [0.070; 0.130]"))
  sg <- attr(p, "subgroup")
  m  <- sg$sub_var == "sex" & sg$level == "M"
  se <- 0.14 / (2 * stats::qnorm(0.975))
  expect_equal(c(sg$estimate[m], sg$conf.low[m], sg$conf.high[m]),
               c(0.12, 0.05, 0.19))
  expect_equal(sg$std.error[m], se)
  expect_equal(sg$p.value[m], 2 * stats::pnorm(-0.12 / se))
  expect_equal(sg[!m, ], st[!m, ])            # other rows untouched
  expect_equal(sg$p_inter, st$p_inter)        # the forest's Wald test is kept
  lab <- trimws(fp_col(p, 1))
  expect_identical(fp_col(p, 3)[lab %in% c("All patients", "M")],
                   c("10.0 (7.0, 13.0)", "12.0 (5.0, 19.0)"))
  expect_identical(fp_col(p, 4)[lab == "M"], "<0.001")

  # a ratio is read on the log scale, at the plot's conf_level, as given
  pr <- plt_hte_sub(res, sub_var = "sex", measure = "ratio", conf_level = 0.9,
                    effect = list(sex = list(F = "1.50 (1.20, 1.875)")))
  f  <- attr(pr, "subgroup")
  f  <- f[f$level == "F", ]
  expect_equal(c(f$estimate, f$conf.low, f$conf.high), c(1.5, 1.2, 1.875))
  expect_equal(f$std.error, log(1.875 / 1.2) / (2 * stats::qnorm(0.95)))
  expect_equal(f$p.value, 2 * stats::pnorm(-log(1.5) / f$std.error))

  # a bare estimate keeps the row's SE; "conf_" re-levels that row alone
  z95 <- stats::qnorm(0.975)
  z75 <- stats::qnorm(0.875)
  f0 <- st[st$sub_var == "sex" & st$level == "F", ]
  m0 <- st[st$sub_var == "sex" & st$level == "M", ]
  pe <- plt_hte_sub(res, sub_var = "sex",
                    effect = list(sex = c(F = "\u20130.050", M = "conf_0.75"),
                                  "All patients" = "conf_0.75"))
  sp <- attr(pe, "subgroup")
  f1 <- sp[sp$level == "F", ]
  m1 <- sp[sp$level == "M", ]
  expect_equal(c(f1$estimate, f1$std.error), c(-0.05, f0$std.error))
  expect_equal(c(f1$conf.low, f1$conf.high), -0.05 + c(-1, 1) * z95 * f0$std.error)
  expect_equal(f1$p.value, 2 * stats::pnorm(-0.05 / f0$std.error))
  expect_equal(m1[c("estimate", "std.error", "p.value")],
               m0[c("estimate", "std.error", "p.value")])
  expect_equal(c(m1$conf.low, m1$conf.high),
               m0$estimate + c(-1, 1) * z75 * m0$std.error)
  ov <- res$stats[res$stats$measure == "diff", ]
  expect_identical(fp_col(pe, 3)[trimws(fp_col(pe, 1)) == "All patients"],
                   sprintf("%.1f (%.1f, %.1f)", 100 * ov$estimate,
                           100 * (ov$estimate - z75 * ov$std.error),
                           100 * (ov$estimate + z75 * ov$std.error)))
  r0 <- attr(plt_hte_sub(res, sub_var = "sex", measure = "ratio"), "subgroup")
  r1 <- attr(plt_hte_sub(res, sub_var = "sex", measure = "ratio",
                         effect = list(sex = c(F = "1.20"))), "subgroup")
  expect_equal(c(r1$conf.low[1], r1$conf.high[1]),
               exp(log(1.2) + c(-1, 1) * z95 * r0$std.error[1]))

  # a level drawn under one variable only may stand alone
  nested <- plt_hte_sub(res, sub_var = c("stage", "sex"),
                        effect = list(sex = c(M = "0.120 (0.050, 0.190)"),
                                      stage = c(III = "conf_0.75")))
  bare   <- plt_hte_sub(res, sub_var = c("stage", "sex"),
                        effect = list(M = "0.120 (0.050, 0.190)", III = "conf_0.75"))
  expect_equal(attr(bare, "subgroup"), attr(nested, "subgroup"))
  res2 <- res
  res2$data$grp <- ifelse(res2$data$sex == "M", "I", "II")
  expect_error(suppressMessages(plt_hte_sub(res2, sub_var = c("stage", "grp"),
                                            effect = list(I = "0.1 (0, 0.2)"))),
               "drawn under stage, grp")

  one <- function(...) plt_hte_sub(res, sub_var = "sex", ...)
  expect_error(one(effect = list(M = "0.1 (0, 0.2)", sex = c(M = "conf_0.9"))),
               "sex = M twice")
  # a number is a bare estimate
  by_str <- attr(one(effect = list(F = "-0.05")), "subgroup")
  expect_equal(attr(one(effect = list(F = -0.05)), "subgroup"), by_str)
  expect_equal(attr(one(effect = list(sex = c(F = -0.05))), "subgroup"), by_str)
  expect_error(one(effect = list(F = Inf)), 'sex = F ("Inf")', fixed = TRUE)
  expect_error(one(effect = list(sex = c(M = "conf_1.5"))), 'sex = M ("conf_1.5")',
               fixed = TRUE)
  expect_error(one(measure = "ratio", effect = list(sex = c(M = "-0.5"))),
               "positive")
  expect_error(one(effect = "0.1 (0, 0.2)"), "`effect`")
  expect_error(one(effect = list(stage = c(I = "0.1 (0, 0.2)"))), "stage")
  expect_error(one(overall = FALSE, effect = list("All patients" = "0.1 (0, 0.2)")),
               "All patients")
  expect_error(one(effect = list(sex = "0.1 (0, 0.2)")), "level")
  expect_error(one(effect = list(sex = c(X = "0.1 (0, 0.2)"))), "X")
  expect_error(one(effect = list(sex = c(M = "abc"))), 'sex = M ("abc")',
               fixed = TRUE)
  expect_error(one(effect = list(sex = c(M = "0.3 (0, 0.2)"))), "sex = M")
  expect_error(one(measure = "ratio", effect = list(sex = c(M = "0.1 (-0.1, 0.3)"))),
               "sex = M")
})

test_that("effect moves the All patients row with the subgroup rows it replaces", {
  skip_if_not_installed("forestplot")
  skip_if_not_installed("ggplotify")
  res <- dep_res()
  N   <- nrow(res$data)
  z   <- stats::qnorm(0.975)
  s   <- .hte_arm_scores(res$fit)
  iM  <- res$data$sex == "M"
  pM  <- mean(iM)
  ov_txt <- function(p) fp_col(p, 3)[trimws(fp_col(p, 1)) == "All patients"]
  fmt <- function(e, l, h, f = "%.3f")
    sprintf(paste0(f, " (", f, ", ", f, ")"), e, l, h)
  one <- function(...) plt_hte_sub(res, sub_var = "sex", ...)
  dif <- function(...) one(axis_args = list(percent = FALSE), ...)
  ov <- res$stats[res$stats$measure == "diff", ]
  st <- attr(one(), "subgroup")
  m0 <- st[st$level == "M", ]

  # diff: the overall AIPW estimate is the patient-weighted mean of its
  # subgroups', so replacing every level gives that mean of the new values
  e_all <- (sum(!iM) * 0.10 + sum(iM) * 0.40) / N
  expect_identical(ov_txt(dif(effect = list(sex = c(F = "0.10", M = "0.40")))),
                   fmt(e_all, e_all - z * ov$std.error, e_all + z * ov$std.error))
  # a bare estimate shifts it by the row's share and keeps its SE ...
  e1 <- ov$estimate + pM * (0.40 - m0$estimate)
  expect_identical(ov_txt(dif(effect = list(M = "0.40"))),
                   fmt(e1, e1 - z * ov$std.error, e1 + z * ov$std.error))
  # ... an interval also changes its variance by the row's share squared
  se2 <- sqrt(ov$std.error^2 + pM^2 * ((0.2 / (2 * z))^2 - m0$std.error^2))
  expect_identical(ov_txt(dif(effect = list(M = "0.40 (0.30, 0.50)"))),
                   fmt(e1, e1 - z * se2, e1 + z * se2))

  # ratio and OR keep the row's control-arm mean and move its treated arm
  grid <- function(m) data.frame(estimand = "ATE", measure = m)
  A <- mean(s$g1); B <- mean(s$g0)
  AM <- mean(s$g1[iM]); BM <- mean(s$g0[iM])
  se_r <- .hte_estimate(res$fit, s, rep(TRUE, N), grid("ratio"), FALSE, z,
                        "Overall")$std.error
  th_r <- log((A + pM * (2 * BM - AM)) / B)
  expect_identical(ov_txt(one(measure = "ratio", effect = list(M = "2.00"))),
                   fmt(exp(th_r), exp(th_r - z * se_r), exp(th_r + z * se_r), "%.2f"))
  se_o <- .hte_estimate(res$fit, s, rep(TRUE, N), grid("OR"), FALSE, z,
                        "Overall")$std.error
  A1   <- A + pM * (stats::plogis(stats::qlogis(BM) + log(2)) - AM)
  th_o <- stats::qlogis(A1) - stats::qlogis(B)
  expect_identical(ov_txt(one(measure = "OR", effect = list(M = "2.00"))),
                   fmt(exp(th_o), exp(th_o - z * se_o), exp(th_o + z * se_o), "%.2f"))

  # "conf_" rows leave it; an "All patients" entry of its own wins, so
  # "conf_<conf_level>" keeps the forest's row
  forest <- fmt(ov$estimate, ov$conf.low, ov$conf.high)
  expect_identical(ov_txt(dif(effect = list(M = "conf_0.75"))), forest)
  expect_identical(ov_txt(dif(effect = list(M = "0.40", "All patients" = "conf_0.95"))),
                   forest)
  # rows of two variables would move it two ways
  expect_error(plt_hte_sub(res, sub_var = c("stage", "sex"),
                           effect = list(M = "0.40", I = "0.10")),
               "All patients")
  expect_identical(
    ov_txt(plt_hte_sub(res, sub_var = c("stage", "sex"),
                       axis_args = list(percent = FALSE),
                       effect = list(M = "0.40", I = "conf_0.75"))),
    fmt(e1, e1 - z * ov$std.error, e1 + z * ov$std.error))
})

test_that("label_column_list renames rows; the SEER dictionary names variables only", {
  skip_if_not_installed("forestplot")
  skip_if_not_installed("ggplotify")
  skip_if_not_installed("RegR")
  res <- dep_res()
  res$data$age_55 <- ifelse(res$data$age < 55, "Below 55", "55 and above")
  p <- suppressMessages(plt_hte_sub(
    res, sub_var = c("age_55", "sex"),
    label_column_list = c(sex = "Sex", F = "Female", "All patients" = "Overall")))
  # age_55 from the SEER dictionary; the level M is not its "M stage"
  expect_identical(trimws(fp_col(p, 1)),
                   c("Subgroup", "Overall", "Age category (year)", "55 and above",
                     "Below 55", "Sex", "Female", "M"))
  expect_identical(unique(attr(p, "subgroup")$sub_var), c("age_55", "sex"))
  # a user key wins over the dictionary
  q <- suppressMessages(plt_hte_sub(res, sub_var = "age_55", overall = FALSE,
                                    label_column_list = list(age_55 = "Age")))
  expect_identical(trimws(fp_col(q, 1))[2], "Age")
  # a name given twice keeps its first label, as RegR's colname_map does
  q2 <- suppressMessages(plt_hte_sub(res, sub_var = "age_55", overall = FALSE,
                                     label_column_list = list(age_55 = "Age",
                                                              age_55 = "Age group")))
  expect_identical(trimws(fp_col(q2, 1))[2], "Age")

  expect_identical(names(formals(plt_hte_sub))[9:11],
                   c("show_pinter", "label_column_list", "header_list"))
  expect_error(plt_hte_sub(res, label_column_list = "Sex"), "label_column_list")
  expect_error(plt_hte_sub(res, label_column_list = list(sex = 1)), "label_column_list")
})

test_that("sub_var defaults to the categorical covariates and checks columns", {
  skip_if_not_installed("forestplot")
  skip_if_not_installed("ggplotify")
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
  skip_if_not_installed("ggplotify")
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
  expect_length(fp_grobs(p, "polygon"), 3L)   # no diamond for III
  expect_identical(fp_headers(p)[3], "5-year DSS difference, %")   # months
  # effect fills the empty level in, diamond and all
  pf <- suppressWarnings(plt_hte_sub(res, sub_var = "stage",
                                     effect = list(stage = c(III = "-0.050 (-0.150, 0.050)"))))
  expect_identical(fp_col(pf, 3)[row], "-5.0 (-15.0, 5.0)")
  expect_length(fp_grobs(pf, "polygon"), 4L)
  # ... and "All patients" follows it from the level's raw arm scores
  s3 <- .hte_arm_scores(res$fit)
  i3 <- res$data$stage == "III"
  e3 <- res$stats$estimate[res$stats$measure == "diff"] +
    mean(i3) * (-0.05 - mean((s3$g1 - s3$g0)[i3]))
  expect_true(startsWith(fp_col(pf, 3)[trimws(fp_col(pf, 1)) == "All patients"],
                         sprintf("%.1f (", 100 * e3)))
  # a bare estimate or a level needs the row's own estimate
  expect_error(suppressWarnings(plt_hte_sub(res, sub_var = "stage",
                                            effect = list(stage = c(III = "0.05")))),
               "no estimate")
  expect_identical(
    fp_headers(suppressWarnings(plt_hte_sub(res, sub_var = "stage",
                                            measure = "ratio")))[3],
    "Event risk ratio")

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
  skip_if_not_installed("ggplotify")
  res <- dep_res()
  expect_error(plt_hte_sub(list()), "hte_res")
  expect_error(plt_hte_sub(res, overall = NA), "overall")
  expect_error(plt_hte_sub(res, conf_level = 1.5), "conf_level")
  expect_error(plt_hte_sub(res, show_pvalue = "yes"), "show_pvalue")
  expect_error(plt_hte_sub(res, axis_args = list(xlim = c(1, 0))), "xlim")
  expect_error(plt_hte_sub(res, measure = "ratio",
                           axis_args = list(xlim = c(-1, 2))), "positive")
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
  expect_s3_class(p, "ggplot")
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
  skip_if_not_installed("ggplotify")
  res <- dep_res()
  expect_identical(tail(names(formals(plt_hte_sub)), 6L),
                   c("fixed_size", "row_height", "ci_shape", "box_col", "...",
                     "save"))

  p <- plt_hte_sub(res, sub_var = "stage")
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

  # row_height sets the pinned row pitch, and the height follows it
  p21 <- plt_hte_sub(res, sub_var = "stage", row_height = 0.21)
  size21 <- attr(p21, "plot_size")
  expect_equal(drawn(p21, 14, 10)$heights, rep(0.21, 6L), tolerance = 1e-6)
  expect_identical(size21[["width"]], size[["width"]])
  expect_equal(size[["height"]] - size21[["height"]], 0.09 * 6, tolerance = 0.1)
  for (bad in list("yes", c(0.2, 0.3), 0, NA))
    expect_error(plt_hte_sub(res, row_height = bad), "`row_height`")

  # FALSE stretches the table over the page, as before
  pf <- plt_hte_sub(res, sub_var = "stage", fixed_size = FALSE)
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

# Every text the drawing writes, header cells, axis title and labels included.
fp_texts <- function(p)
  unname(vapply(fp_grobs(p, "text"), function(k) paste(k$label, collapse = ""), ""))
fp_fills <- function(p, cls)
  unname(vapply(fp_grobs(p, cls), function(k) paste(k$gp$fill, collapse = ","), ""))

test_that("plt_hte_sub shows percentages, plain titles and direction labels", {
  skip_if_not_installed("forestplot")
  skip_if_not_installed("ggplotify")
  res <- dep_res()
  ov  <- res$stats[res$stats$measure == "diff", ]
  # a difference of risks reads in percentage points by default
  p <- plt_hte_sub(res, sub_var = "stage")
  expect_identical(fp_headers(p)[3], "Risk difference, %")
  expect_identical(fp_col(p, 3)[2],
                   sprintf("%.1f (%.1f, %.1f)", 100 * ov$estimate,
                           100 * ov$conf.low, 100 * ov$conf.high))
  ax <- fp_grobs(p, "xaxis")[[1L]]
  expect_identical(ax$label, as.character(round(100 * ax$at, 8)))
  # z is 0/1, so the arms are named after it
  expect_true("Risk difference (%), z - no z" %in% fp_texts(p))
  expect_true(all(c("Higher with no z", "Higher with z") %in% fp_texts(p)))
  # the estimates themselves keep their scale
  expect_equal(attr(p, "subgroup"),
               attr(plt_hte_sub(res, sub_var = "stage",
                                axis_args = list(percent = FALSE)), "subgroup"))

  p0 <- plt_hte_sub(res, sub_var = "stage", conf_level = 0.9,
                    axis_args = list(percent = FALSE, direction_labels = FALSE))
  expect_identical(fp_headers(p0)[3], "Risk difference")
  expect_identical(fp_col(p0, 3)[2],
                   sprintf("%.3f (%.3f, %.3f)", ov$estimate,
                           ov$estimate - stats::qnorm(0.95) * ov$std.error,
                           ov$estimate + stats::qnorm(0.95) * ov$std.error))
  expect_false(any(c("Higher with no z", "Higher with z") %in% fp_texts(p0)))
  expect_true("Risk difference, z - no z" %in% fp_texts(p0))

  p1 <- plt_hte_sub(res, sub_var = "stage",
                    header_list = list(effect = "ARD (95% CI)"),
                    xlab = "Absolute risk difference",
                    axis_args = list(direction_labels = c("Favours control",
                                                          "Favours z")))
  expect_identical(fp_headers(p1)[3], "ARD (95% CI)")
  expect_true(all(c("Absolute risk difference", "Favours control",
                    "Favours z") %in% fp_texts(p1)))
  # each line of a header takes a table row of its own; the direction labels
  # sit in the last header row
  q <- plt_hte_sub(res, sub_var = "stage", conf_level = 0.9,
                   header_list = list(effect = "Risk difference, %\n(90% CI)"))
  expect_identical(fp_headers(q)[3], "Risk difference, %\n(90% CI)")
  expect_length(drawn(p, 14, 10)$heights, 6L)
  expect_length(drawn(q, 14, 10)$heights, 7L)
  vp_of <- function(q, lab) Filter(function(k) identical(k$label, lab),
                                   fp_grobs(q, "text"))[[1L]]$vp$name
  expect_match(vp_of(p, "Higher with z"), "^Line_1_")
  expect_match(vp_of(q, "Higher with z"), "^Line_2_")
  p3 <- plt_hte_sub(res, sub_var = "stage",
                    header_list = list(effect = "Risk\ndifference\n(95% CI)"))
  expect_identical(fp_headers(p3)[3], "Risk\ndifference\n(95% CI)")
  expect_length(drawn(p3, 14, 10)$heights, 8L)
  # ratios are no percentages
  pr <- plt_hte_sub(res, sub_var = "stage", measure = "ratio")
  expect_identical(fp_headers(pr)[3], "Risk ratio")

  expect_error(plt_hte_sub(res, axis_args = list(percent = NA)),
               "`axis_args$percent`", fixed = TRUE)
  expect_error(plt_hte_sub(res, header_list = list(effect = 1)),
               "`header_list$effect`", fixed = TRUE)
  expect_error(plt_hte_sub(res, xlab = c("a", "b")), "`xlab`")
  expect_error(plt_hte_sub(res, axis_args = list(direction_labels = "left")),
               "`axis_args$direction_labels`", fixed = TRUE)
})

test_that("plt_hte_sub takes header_list and deprecates est_label", {
  skip_if_not_installed("forestplot")
  skip_if_not_installed("ggplotify")
  res <- dep_res()
  p <- plt_hte_sub(res, sub_var = "stage", show_pvalue = TRUE, show_pinter = TRUE,
                   header_list = list(count = "N", p.value = "P value",
                                      p.interaction = "P inter"))
  # the headers left out keep their defaults
  expect_identical(fp_headers(p),
                   c("Subgroup", "N", "Risk difference, %", "P value", "P inter"))
  # fp_table() finds the header by "Subgroup", so a renamed one is looked up
  lab <- plt_hte_sub(res, sub_var = "stage", header_list = list(label = "Stage"))
  expect_true("Stage" %in% fp_texts(lab))
  expect_false("Subgroup" %in% fp_texts(lab))
  expect_error(plt_hte_sub(res, header_list = list(estimate = "x")),
               "unknown field")
  expect_error(plt_hte_sub(res, est_lab = "x"), "unused argument")

  expect_false("est_label" %in% names(formals(plt_hte_sub)))
  expect_warning(q <- plt_hte_sub(res, sub_var = "stage", est_label = "ARD"),
                 "use `header_list$effect` instead", fixed = TRUE)
  expect_identical(fp_headers(q)[3], "ARD")
  expect_warning(q <- plt_hte_sub(res, sub_var = "stage", est_label = "ARD",
                                  header_list = list(effect = "RD")),
                 "`header_list$effect` wins", fixed = TRUE)
  expect_identical(fp_headers(q)[3], "RD")
})

test_that("axis_args gathers the axis settings of plt_hte_sub", {
  skip_if_not_installed("forestplot")
  skip_if_not_installed("ggplotify")
  res <- dep_res()
  expect_identical(eval(formals(plt_hte_sub)$axis_args),
                   list(xlim = NULL, xticks = NULL, percent = TRUE,
                        direction_labels = TRUE))
  # xticks sits on the scale of the estimates; percent labels it
  p <- plt_hte_sub(res, sub_var = "stage",
                   axis_args = list(xticks = c(-0.2, 0, 0.2)))
  ax <- fp_grobs(p, "xaxis")[[1L]]
  expect_equal(ax$at, c(-0.2, 0, 0.2))
  expect_identical(ax$label, c("-20", "0", "20"))
  # a partial list keeps the other defaults: the direction labels stay on
  expect_true("Higher with z" %in% fp_texts(p))
  # xlim and xticks together, as the notes draw them
  q <- plt_hte_sub(res, sub_var = "stage",
                   axis_args = list(xlim = c(-0.4, 0.4), xticks = c(-0.4, 0, 0.4),
                                    percent = FALSE))
  expect_identical(fp_grobs(q, "xaxis")[[1L]]$label, c("-0.4", "0", "0.4"))

  expect_error(plt_hte_sub(res, axis_args = "x"), "`axis_args`")
  expect_error(plt_hte_sub(res, axis_args = list(ticks_at = 0)), "ticks_at")
  expect_error(plt_hte_sub(res, axis_args = list(xlim = c(0, 1), xlim = c(0, 2))),
               "duplicated")
  expect_error(plt_hte_sub(res, axis_args = list(xticks = "a")), "xticks")
  # the arguments it replaces are gone
  for (old in c("xlim", "ticks_at", "percent", "direction_labels", "overall_shape"))
    expect_false(old %in% names(formals(plt_hte_sub)))
})

test_that("ci_shape and box_col draw a meta-analysis style", {
  skip_if_not_installed("forestplot")
  skip_if_not_installed("ggplotify")
  res <- dep_res()
  # default: a blue4 diamond for every estimate
  p <- plt_hte_sub(res, sub_var = "stage")
  expect_identical(fp_fills(p, "polygon"), rep("blue4", 4L))

  expect_identical(eval(formals(plt_hte_sub)$ci_shape), "diamond")
  m <- plt_hte_sub(res, sub_var = "stage", ci_shape = c("normal", "summary"),
                   box_col = c("blue4", "#B2182B"))
  # one summary diamond for "All patients", squares for the levels
  expect_identical(fp_fills(m, "polygon"), "#B2182B")
  sq <- Filter(function(k) identical(k$gp$fill, "blue4"), fp_grobs(m, "rect"))
  expect_length(sq, 3L)
  # the most precise level gets the largest square
  se <- attr(m, "subgroup")$std.error
  grDevices::pdf(NULL)   # converting units needs a device, not the default one
  side <- vapply(sq, function(k) grid::convertWidth(k$width, "in", TRUE), 0)
  grDevices::dev.off()
  expect_identical(unname(which.max(side)), which.min(se))
  # the summary row is bold across its columns
  bold <- function(lab) Filter(function(k) identical(k$label, lab),
                               fp_grobs(m, "text"))[[1L]]$gp$font
  expect_identical(unname(bold(fp_col(m, 3)[2])), 2L)

  # one ci_shape serves the overall row too; a second colour marks it
  n <- plt_hte_sub(res, sub_var = "stage", ci_shape = "normal")
  expect_length(fp_grobs(n, "polygon"), 0L)
  expect_length(Filter(function(k) identical(k$gp$fill, "blue4"),
                       fp_grobs(n, "rect")), 4L)
  d2 <- plt_hte_sub(res, sub_var = "stage", box_col = c("blue4", "red"))
  expect_identical(sort(fp_fills(d2, "polygon")), c(rep("blue4", 3L), "red"))

  expect_error(plt_hte_sub(res, ci_shape = "square"), "`ci_shape`")
  expect_error(plt_hte_sub(res, ci_shape = c("normal", "box")), "`ci_shape`")
  expect_error(plt_hte_sub(res, ci_shape = "summary"), "`ci_shape`")
  expect_error(plt_hte_sub(res, ci_shape = c("normal", "normal", "summary")),
               "`ci_shape`")
  expect_error(plt_hte_sub(res, box_col = c("red", "blue", "green")), "`box_col`")
  expect_error(plt_hte_sub(res, box_col = "notacolour"), "`box_col`")
  expect_identical(names(formals(plt_hte_sub))[10:14],
                   c("label_column_list", "header_list", "xlab", "axis_args",
                     "title"))
})

test_that("plt_hte_sub returns a ggplot patchwork can combine", {
  skip_if_not_installed("forestplot")
  skip_if_not_installed("ggplotify")
  skip_if_not_installed("patchwork")
  res <- dep_res()
  # the wrapping raises no warning of its own
  expect_no_warning(p <- plt_hte_sub(res, sub_var = "stage"))
  expect_s3_class(p, "ggplot")
  q  <- plt_hte_cate(res, sub_var = "stage", type = "density")
  pw <- patchwork::wrap_plots(p, q, ncol = 1)
  expect_s3_class(pw, "patchwork")
  # the pinned table draws in the composition as it does alone
  keep <- c("widths", "heights", "area")
  expect_equal(drawn(pw, 9, 8)[keep], drawn(p, 9, 8)[keep], tolerance = 1e-6)

  # unpinned, it takes a proportional share
  pf <- plt_hte_sub(res, sub_var = "stage", fixed_size = FALSE)
  even <- drawn(patchwork::wrap_plots(pf, q, ncol = 1), 9, 8)
  more <- drawn(patchwork::wrap_plots(pf, q, ncol = 1, heights = c(2, 1)), 9, 8)
  expect_gt(more$area[2], even$area[2])
})


# ---- plt_hte_cate() ---------------------------------------------------------

test_that("CATE-only plots do not compute unused arm scores", {
  res <- dep_res()
  local_mocked_bindings(.hte_arm_scores = function(...) stop("unused scores"))
  for (type in c("waterfall", "density")) {
    p <- plt_hte_cate(res, type = type, overall = FALSE)
    expect_s3_class(p, "ggplot")
    expect_null(attr(p, "subgroup"))
    expect_no_error(ggplot2::ggplot_build(p))
  }
  p <- plt_hte_cate(res, show_ci = TRUE, overall = FALSE)
  expect_true(all(is.finite(layer_data_of(p, "GeomLinerange")$conf.low)))
})

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
  # Intervals must be drawn over the opaque bars, including without groups.
  for (q in list(p, plt_hte_cate(res, show_ci = TRUE, overall = FALSE))) {
    geoms <- vapply(q$layers, function(l) class(l$geom)[1L], "")
    expect_gt(which(geoms == "GeomLinerange"), which(geoms == "GeomCol"))
    built <- ggplot2::ggplot_build(q)$data
    expect_equal(built[[which(geoms == "GeomLinerange")]]$ymin,
                 layer_data_of(q, "GeomLinerange")$conf.low)
  }
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


# ---- Theme ------------------------------------------------------------------

test_that("CATE plots exclude SPSS user-missing subgroup values", {
  res <- dep_res()
  values <- as.numeric(res$data$stage)
  values[1:8] <- c(999, 95:99, NA, 999)
  res$data$aux <- structure(values,
                            class = c("haven_labelled_spss", "haven_labelled"),
                            na_values = 999, na_range = c(95, 99))
  original <- res$data$aux
  ref <- res
  values[1:8] <- NA_real_
  ref$data$aux <- values
  for (type in c("waterfall", "density")) {
    p <- suppressMessages(plt_hte_cate(res, sub_var = "aux", type = type))
    want <- suppressMessages(plt_hte_cate(ref, sub_var = "aux", type = type))
    geom <- if (type == "waterfall") "GeomCol" else "GeomDensity"
    expect_identical(nrow(layer_data_of(p, geom)), nrow(res$data) - 8L)
    expect_equal(ggplot2::ggplot_build(p)$data, ggplot2::ggplot_build(want)$data)
    expect_equal(attr(p, "subgroup"), attr(want, "subgroup"))
    expect_identical(res$data$aux, original)
  }
})

test_that("CATE plots diagnose empty subgroups and retain usable panels", {
  skip_if_not_installed("patchwork")
  res <- dep_res()
  res$data$empty <- NA
  for (type in c("waterfall", "density")) {
    expect_error(suppressMessages(plt_hte_cate(res, sub_var = "empty", type = type)),
                 "No non-missing observations.*empty")
    expect_warning(p <- suppressMessages(plt_hte_cate(
      res, sub_var = c("stage", "empty", "sex"), type = type)),
      "Omitted.*empty.*no non-missing observations")
    ref <- plt_hte_cate(res, sub_var = c("stage", "sex"), type = type)
    expect_s3_class(p, "patchwork")
    expect_identical(vapply(panels_of(p), strip_of, ""), c("stage", "sex"))
    expect_equal(attr(p, "subgroup"), attr(ref, "subgroup"))
    expect_equal(attr(p, "plot_size"), attr(ref, "plot_size"))
    expect_no_error(patchwork::patchworkGrob(p))
    expect_warning(single <- suppressMessages(plt_hte_cate(
      res, sub_var = c("empty", "stage"), type = type)), "Omitted.*empty")
    expect_false(inherits(single, "patchwork"))
    expect_identical(strip_of(single), "stage")
  }
  for (value in list(rep(NA_character_, nrow(res$data)),
                     factor(rep(NA_character_, nrow(res$data)), levels = "unused"),
                     structure(rep(99, nrow(res$data)),
                               class = c("haven_labelled_spss", "haven_labelled"),
                               na_values = 99))) {
    res$data$empty <- value
    expect_error(suppressMessages(plt_hte_cate(res, sub_var = "empty")),
                 "No non-missing observations.*empty")
  }
  # An explicit NA factor level is a group, not a missing factor code.
  res$data$empty <- factor(rep(NA_character_, nrow(res$data)), exclude = NULL)
  p <- suppressMessages(plt_hte_cate(res, sub_var = "empty"))
  expect_identical(nrow(layer_data_of(p, "GeomCol")), nrow(res$data))
})

# The fonts and line widths a theme draws with.
look_of <- function(th) {
  el <- function(name) ggplot2::calc_element(name, th)
  list(family = el("text")$family, face = el("text")$face,
       size = el("text")$size, axis_text = el("axis.text.y")$size,
       ticks = el("axis.ticks")$linewidth,
       border = el("panel.border")$linewidth,
       grid = el("panel.grid.major")$linetype)
}

test_that("the ggplots draw with the theme of MLR::plt_bar_per()", {
  res  <- dep_res()
  want <- look_of(UtilsR::theme_my(base_rect_size = 1.5))
  dep  <- plt_hte_dep(res)
  cate <- plt_hte_cate(res, sub_var = c("stage", "sex"))
  heat <- plt_hte_dep(res, x_var = c("age", "sex"), type = "heat")
  for (q in c(panels_of(dep), panels_of(cate), list(heat),
              list(plt_hte_cate(res, type = "density")),
              list(plt_hte_rate(res, priority = "age", type = c("toc", "gates")))))
    expect_identical(look_of(q$theme), want)
  # a patchwork draws its title and caption with it too
  for (pw in list(dep, cate))
    expect_identical(look_of(pw$patches$annotation$theme), want)
  # theme_my() drops legend titles; the heat map keeps the title of its fill
  expect_false(inherits(heat$theme$legend.title, "element_blank"))
})
