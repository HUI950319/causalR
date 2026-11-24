# save_hte() stacks a top panel over plt_hte_sub(). grf, sandwich, forestplot,
# ggplotify and patchwork are Suggests; every test that needs a forest skips
# without them, and the importance panel without MLR.

save_cache <- new.env()

save_res <- function() {
  for (p in c("grf", "sandwich", "forestplot", "ggplotify", "patchwork"))
    skip_if_not_installed(p)
  if (is.null(save_cache$res)) {
    set.seed(20261115)
    n <- 400L
    d <- data.frame(age = round(stats::runif(n, 20, 85)),
                    sex = factor(sample(c("F", "M"), n, replace = TRUE)),
                    grp = factor(sample(c("a", "b"), n, replace = TRUE)),
                    x1  = stats::rnorm(n), x2 = stats::rnorm(n))
    d$PORT <- factor(sample(c("No", "Yes"), n, replace = TRUE),
                     levels = c("No", "Yes"))
    ev   <- stats::rexp(n, 0.02 * exp(-0.8 * (d$PORT == "Yes") * (d$sex == "M")))
    cens <- pmin(stats::rexp(n, 0.01), 120)
    d$time <- pmin(ev, cens)
    d$DSS  <- as.integer(ev <= cens)
    save_cache$res <- suppressMessages(get_hte(
      d, cat_var = "PORT", sub_var = "sex",
      adj_var = c("age", "grp", "x1", "x2"), surv = TRUE, time = 60,
      grf_args = list(num.trees = 200, seed = 1)))
  }
  save_cache$res
}

# The arguments save_hte() hands to one of its plot functions, and its plot.
# The mock keeps the formals of the function, which save_hte() checks its
# `*_args` by.
capture_call <- function(fun, pkg, ...) {
  got  <- NULL
  orig <- get(fun, envir = asNamespace(pkg))
  mock <- orig
  body(mock) <- quote({
    cl  <- match.call()
    got <<- lapply(as.list(cl)[-1L], eval, parent.frame())
    cl[[1L]] <- orig
    eval(cl, parent.frame())
  })
  environment(mock) <- environment()
  do.call(local_mocked_bindings,
          c(stats::setNames(list(mock), fun), list(.package = pkg, .env = environment())))
  p <- save_hte(...)
  list(args = got, plot = p)
}
capture_sub <- function(...) capture_call("plt_hte_sub", "causalR", ...)

test_that("save_hte() has the usage of its contract", {
  expect_identical(names(formals(save_hte)),
                   c("x", "sub_var", "panel", "name_map", "sub_args",
                     "imp_args", "cate_args", "top_height", "color", "save"))
})

test_that("color colours the three plots from one palette, and a colour in *_args wins", {
  res <- save_res()
  pal <- UtilsR::pal_get("PiYG")
  expect_length(pal, 11L)
  args_of <- function(fun, pkg, panel, ...)
    suppressWarnings(capture_call(fun, pkg, res, sub_var = "sex", panel = panel, ...))$args

  # the forest plot: the subgroup rows, then the "All patients" row
  expect_identical(args_of("plt_hte_sub", "causalR", "cate", color = pal)$box_col,
                   pal[c(3, 9)])
  # a palette name of UtilsR::pal_get() is that palette
  expect_identical(args_of("plt_hte_sub", "causalR", "cate", color = "PiYG")$box_col,
                   pal[c(3, 9)])
  # without it the plot keeps its own colours, and its own colour wins
  expect_false("box_col" %in% names(args_of("plt_hte_sub", "causalR", "cate")))
  expect_identical(args_of("plt_hte_sub", "causalR", "cate", color = pal,
                           sub_args = list(box_col = "red"))$box_col, "red")

  # the CATE waterfall: from the most negative to the most positive
  expect_identical(args_of("plt_hte_cate", "causalR", "cate", color = pal)$color,
                   pal[c(9, 3)])
  expect_false("color" %in% names(args_of("plt_hte_cate", "causalR", "cate")))
  expect_identical(args_of("plt_hte_cate", "causalR", "cate", color = pal,
                           cate_args = list(color = c("blue", "red")))$color,
                   c("blue", "red"))

  expect_error(save_hte(res, sub_var = "sex", color = pal[1:8]), "at least 9")
  expect_error(save_hte(res, sub_var = "sex", color = "no such palette"),
               "`color`", fixed = TRUE)
  expect_error(save_hte(res, sub_var = "sex", color = 1:11), "`color`", fixed = TRUE)

  # the importance bars and their rings
  skip_if_not_installed("MLR")
  expect_identical(args_of("plt_bar_per", "MLR", "importance", color = pal)$bar_per_palette,
                   list(pal[5:3], pal[7:9]))
  expect_false("bar_per_palette" %in% names(args_of("plt_bar_per", "MLR", "importance")))
  own <- list(pal[1:3], pal[9:11])
  expect_identical(args_of("plt_bar_per", "MLR", "importance", color = pal,
                           imp_args = list(bar_per_palette = own))$bar_per_palette, own)
})

test_that("the forest plot gets two-line headers named after the arms", {
  res <- save_res()
  out <- suppressWarnings(capture_sub(res, sub_var = "sex", panel = "cate"))
  a   <- out$args
  expect_identical(a$header_list$count, "N\n(PORT / no PORT)")
  expect_identical(a$header_list$effect, "5-year DSS\ndifference, %")
  expect_identical(a$xlab, "Difference in 5-year DSS (%),\nPORT - no PORT")
  expect_identical(a[c("show_n", "fixed_size", "row_height", "ci_shape")],
                   list(show_n = TRUE, fixed_size = TRUE, row_height = 0.5,
                        ci_shape = c("normal", "summary")))
  expect_s3_class(out$plot, "patchwork")
})

test_that("a partial sub_args keeps the other defaults and its own headers win", {
  res <- save_res()
  out <- suppressWarnings(capture_sub(
    res, sub_var = "sex", panel = "cate",
    sub_args = list(box_col = "red", header_list = list(count = "N"),
                    axis_args = list(xlim = c(-0.3, 0.3)))))
  a <- out$args
  expect_identical(a$box_col, "red")
  expect_identical(a$row_height, 0.5)
  expect_identical(a$header_list$count, "N")
  expect_identical(a$header_list$effect, "5-year DSS\ndifference, %")
  expect_identical(a$axis_args, list(xlim = c(-0.3, 0.3)))
})

test_that("other measures keep the one-line defaults of plt_hte_sub()", {
  res <- save_res()
  out <- suppressWarnings(capture_sub(res, sub_var = "sex", panel = "cate",
                                      sub_args = list(measure = "ratio")))
  expect_null(out$args$header_list)
  expect_null(out$args$xlab)
})

test_that("the page is the forest plot's width plus the margin of the top panel", {
  res <- save_res()
  sub <- suppressWarnings(capture_sub(res, sub_var = "sex", panel = "cate"))
  sz  <- attr(suppressWarnings(do.call(plt_hte_sub, sub$args)), "plot_size")
  expect_equal(unname(attr(sub$plot, "plot_size")),
               unname(c(sz[[1]] + 1.2, 3.6 + sz[[2]] + 0.6)))
  expect_named(attr(sub$plot, "plot_size"), c("width", "height"))

  tall <- suppressWarnings(save_hte(res, sub_var = "sex", panel = "cate",
                                    top_height = 5))
  expect_equal(unname(attr(tall, "plot_size")[2]), 5 + sz[[2]] + 0.6)

  skip_if_not_installed("MLR")
  imp <- suppressWarnings(save_hte(res, sub_var = "sex", panel = "importance"))
  expect_equal(unname(attr(imp, "plot_size")),
               unname(c(sz[[1]] + 0.5, 4.1 + sz[[2]] + 0.6)))
})

test_that("panel = 'auto' draws the importance bars only for a top-three variable", {
  skip_if_not_installed("MLR")
  res  <- save_res()
  imp  <- res$importance
  top3 <- utils::head(imp$variable[order(-imp$importance)], 3L)
  for (v in c("sex", "grp")) {
    p   <- suppressWarnings(save_hte(res, sub_var = v))
    pad <- if (v %in% top3) 0.5 else 1.2
    sz  <- attr(suppressWarnings(capture_sub(res, sub_var = v, panel = "cate"))$plot,
                "plot_size")
    expect_equal(unname(attr(p, "plot_size")[1]), unname(sz[1] - 1.2 + pad))
  }
})

test_that("bad arguments are rejected before any plot is drawn", {
  res <- save_res()
  expect_error(save_hte(res, sub_var = c("sex", "grp")), "single column name")
  expect_error(save_hte(res, sub_var = "sex", top_height = 0), "top_height")
  expect_error(save_hte(res, sub_var = "sex", sub_args = list(foo = 1)),
               "unknown field")
  expect_error(save_hte(res, sub_var = "sex", sub_args = list(name_map = list())),
               "unknown field")
  expect_error(save_hte(res, sub_var = "sex", panel = "cate",
                        cate_args = list(save = list())), "unknown field")
  expect_error(save_hte(list(), sub_var = "sex"), "hte_res")
})

test_that("save completes a bare label from RegR::set_paper()", {
  skip_if_not_installed("RegR")
  res <- save_res()
  old <- options(RegR.paper = NULL)
  on.exit(options(old), add = TRUE)
  suppressMessages(RegR::set_paper(treat = "postoperative radiotherapy",
                                   modifier = "sex"))
  seen <- NULL
  local_mocked_bindings(save_plt = function(plot, filename, ...) {
    seen <<- filename
    invisible(plot)
  }, .package = "RegR")
  draw <- function(x) suppressWarnings(suppressMessages(save_hte(
    x, sub_var = "sex", panel = "cate",
    save = list(filename = "Supplementary Figure 13"))))
  draw(res)
  expect_identical(seen, paste(
    "Supplementary Figure 13. Heterogeneity analysis of the effect of",
    "postoperative radiotherapy on disease-specific survival based on causal",
    "survival forest"))
  # a result that does not say which forest it is gets no name guessed
  bare <- res
  attr(bare, "analysis")$forest <- NULL
  expect_error(draw(bare), "does not record its forest")
})


test_that("save writes a PDF at the page size", {
  skip_if_not_installed("RegR")
  res <- save_res()
  withr::local_dir(withr::local_tempdir())
  p <- suppressWarnings(save_hte(res, sub_var = "sex", panel = "cate",
                                 save = list(filename = "sf5")))
  expect_true(file.exists("sf5.pdf"))
  expect_s3_class(p, "patchwork")
})
