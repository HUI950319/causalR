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

# The arguments save_hte() hands to plt_hte_sub(), and its plot. The mock
# keeps the formals of plt_hte_sub(), which save_hte() checks `sub_args` by.
capture_sub <- function(...) {
  got  <- NULL
  orig <- plt_hte_sub
  mock <- orig
  body(mock) <- quote({
    cl  <- match.call()
    got <<- lapply(as.list(cl)[-1L], eval, parent.frame())
    cl[[1L]] <- orig
    eval(cl, parent.frame())
  })
  environment(mock) <- environment()
  local_mocked_bindings(plt_hte_sub = mock, .package = "causalR")
  p <- save_hte(...)
  list(args = got, plot = p)
}

test_that("save_hte() has the usage of its contract", {
  expect_identical(names(formals(save_hte)),
                   c("x", "sub_var", "panel", "name_map", "sub_args",
                     "imp_args", "cate_args", "top_height", "save"))
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

test_that("save writes a PDF at the page size", {
  skip_if_not_installed("RegR")
  res <- save_res()
  withr::local_dir(withr::local_tempdir())
  p <- suppressWarnings(save_hte(res, sub_var = "sex", panel = "cate",
                                 save = list(filename = "sf5")))
  expect_true(file.exists("sf5.pdf"))
  expect_s3_class(p, "patchwork")
})
