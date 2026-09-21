# MatchIt is an Imports, so the core paths never skip. optmatch (optimal /
# full) and halfmoon (balance) are Suggests and skip when absent.

psm_data <- function(n = 400L) {
  set.seed(20260921)
  d <- data.frame(x1 = stats::rnorm(n),
                  x2 = stats::rbinom(n, 1, 0.4),
                  x3 = stats::runif(n))
  # controls must outnumber treated, or MatchIt warns that not every treated
  # unit can be matched
  lp <- -0.9 - 1.1 * d$x1 + 0.8 * d$x2 - 1.2 * d$x3
  d$z <- stats::rbinom(n, 1, stats::plogis(lp))
  d
}

psm_adj <- c("x1", "x2", "x3")

psm_stats_cols <- c("method", "estimand", "n", "n_treat", "n_ctrl",
                    "n_unmatched", "n_discarded", "pct_retained", "n_pairs",
                    "ess", "ess_treat", "ess_ctrl", "ess_pct", "w_max",
                    "w_cv", "smd_max", "smd_over")

psm_ps <- function(d) {
  unname(stats::fitted(stats::glm(z ~ x1 + x2 + x3, stats::binomial(),
                                  data = d)))
}


test_that("get_PSM returns the documented structure", {
  d   <- psm_data()
  res <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE)

  expect_s3_class(res, "psm_res")
  expect_named(res, c("data", "stats", "balance", "fit"))
  expect_null(res$balance)
  expect_named(res$fit, "nearest")
  expect_s3_class(res$fit$nearest, "matchit")

  expect_named(res$stats, psm_stats_cols)
  expect_identical(nrow(res$stats), 1L)
  expect_identical(res$stats$method, "nearest")
  expect_true(all(c("ps", "w_nearest", "s_nearest") %in% names(res$data)))
  expect_identical(nrow(res$data), nrow(d))      # every input row kept

  a <- attr(res, "analysis")
  expect_identical(a$treat, "z")
  expect_identical(a$wcols, "w_nearest")
  expect_identical(a$estimand, "ATT")
})


test_that("the weights and subclasses are MatchIt's own", {
  d   <- psm_data()
  res <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE)

  obj <- MatchIt::matchit(z ~ x1 + x2 + x3, data = d, method = "nearest",
                          estimand = "ATT", distance = psm_ps(d),
                          normalize = FALSE, ratio = 1)
  expect_equal(res$data$w_nearest, unname(obj$weights))
  expect_equal(res$data$s_nearest,
               as.integer(as.character(obj$subclass)))
  expect_equal(res$data$ps, psm_ps(d), tolerance = 1e-10)
})


test_that("one score is shared by every matching scheme", {
  skip_if_not_installed("optmatch")

  d   <- psm_data()
  res <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                 method = c("nearest", "full"))

  expect_identical(res$stats$method, c("nearest", "full"))
  expect_identical(nrow(res$stats), 2L)
  expect_named(res$fit, c("nearest", "full"))
  expect_true(all(c("w_nearest", "w_full", "s_nearest", "s_full") %in%
                    names(res$data)))
  # the two matchit objects were handed the same distance
  expect_equal(unname(res$fit$nearest$distance),
               unname(res$fit$full$distance))
  expect_equal(res$data$ps, unname(res$fit$full$distance))

  # full matching keeps everyone, nearest does not
  expect_equal(res$stats$pct_retained[res$stats$method == "full"], 100)
  expect_lt(res$stats$pct_retained[res$stats$method == "nearest"], 100)
})


test_that("unmatched units keep their row with weight 0 and no subclass", {
  d   <- psm_data()
  res <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                 caliper = 0.1)
  w <- res$data$w_nearest

  expect_identical(nrow(res$data), nrow(d))
  expect_true(all(is.na(res$data$s_nearest[w == 0])))
  expect_false(anyNA(res$data$s_nearest[w > 0]))
  expect_identical(res$stats$n, sum(w > 0))
  expect_identical(res$stats$n_treat, sum(w > 0 & d$z == 1))

  # match.data() drops exactly the rows we keep at weight 0
  md <- MatchIt::match.data(res$fit$nearest)
  expect_identical(nrow(md), sum(w > 0))
})


test_that("matchit objects carry a readable call, not the data", {
  # do.call() would inline the matchit() body and a data snapshot into each
  # $call, so summary() printed the whole data frame and the object was mostly
  # that. match.data() must still find the data on its own.
  d   <- psm_data()
  res <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                 caliper = 0.2)
  cl  <- res$fit$nearest$call

  expect_identical(cl[[1L]], quote(MatchIt::matchit))
  expect_identical(cl$data, quote(data))
  expect_identical(cl$distance, quote(ps))
  expect_lt(sum(nchar(deparse(cl))), 300L)
  expect_lt(as.numeric(utils::object.size(res$fit$nearest)),
            30 * as.numeric(utils::object.size(d)))

  out <- utils::capture.output(summary(res$fit$nearest))
  expect_lt(length(out), 60L)
  expect_identical(nrow(MatchIt::match.data(res$fit$nearest)),
                   sum(res$data$w_nearest > 0))
})


test_that("unmatched and discarded are counted separately", {
  d   <- psm_data()
  res <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                 caliper = 0.1, match_args = list(discard = "both"))
  disc <- res$fit$nearest$discarded

  expect_gt(sum(disc), 0L)
  expect_identical(res$stats$n_discarded, sum(disc))
  expect_identical(res$stats$n_unmatched,
                   sum(res$data$w_nearest == 0 & !disc))
  # every discarded unit is also unmatched, but not the reverse
  expect_true(all(res$data$w_nearest[disc] == 0))
  expect_gt(res$stats$n_unmatched, 0L)
})


test_that("weights are unnormalised, so an ATT control column sums to n_treated", {
  d   <- psm_data()
  nt  <- sum(d$z == 1)

  for (r in c(1, 2)) {
    res <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                   ratio = r)
    w <- res$data$w_nearest
    expect_equal(sum(w[d$z == 0]), nt, tolerance = 1e-8,
                 info = paste("ratio", r))
    expect_equal(sum(w[d$z == 1]), nt, info = paste("ratio", r))
  }

  # MatchIt's own default would break that at ratio = 2
  norm <- MatchIt::matchit(z ~ x1 + x2 + x3, data = d, method = "nearest",
                           estimand = "ATT", distance = psm_ps(d), ratio = 2)
  expect_gt(sum(norm$weights[d$z == 0]), nt)
})


test_that("weights stop being 0/1 with replacement and with full matching", {
  d <- psm_data()

  plain <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE)
  expect_setequal(unique(plain$data$w_nearest), c(0, 1))

  repl <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                  replace = TRUE)
  expect_gt(length(unique(repl$data$w_nearest)), 2L)
  expect_false(all(repl$data$w_nearest %in% c(0, 1)))

  skip_if_not_installed("optmatch")
  full <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                  method = "full")
  expect_gt(length(unique(full$data$w_full)), 10L)
})


test_that("ESS agrees with MatchIt's own and is unaffected by normalisation", {
  d   <- psm_data()
  res <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                 replace = TRUE)
  nn  <- summary(res$fit$nearest)$nn

  expect_equal(res$stats$ess_ctrl, unname(nn["Matched (ESS)", "Control"]),
               tolerance = 1e-8)
  expect_equal(res$stats$ess_treat, unname(nn["Matched (ESS)", "Treated"]),
               tolerance = 1e-8)
})


test_that("matching improves balance and the table is computed on the cohort", {
  skip_if_not_installed("halfmoon")

  d   <- psm_data()
  res <- get_PSM(d, treat = "z", adj_var = psm_adj, caliper = 0.2)

  expect_s3_class(res$balance, "data.frame")
  expect_setequal(unique(res$balance$method), c("observed", "w_nearest"))
  expect_false(anyNA(res$balance$estimate))

  obs <- max(abs(res$balance$estimate[res$balance$method == "observed"]))
  expect_lt(res$stats$smd_max, obs)

  # computed on the full cohort, not on the matched subset -- the subset's
  # smaller SD would inflate every SMD
  full <- suppressMessages(do.call(halfmoon::check_balance,
    list(.data = res$data, .vars = psm_adj, .exposure = "z",
         .weights = "w_nearest", .metrics = "smd", na.rm = TRUE)))
  expect_equal(res$balance$estimate, full$estimate)

  sub <- res$data[res$data$w_nearest > 0, ]
  onsub <- suppressMessages(do.call(halfmoon::check_balance,
    list(.data = sub, .vars = psm_adj, .exposure = "z",
         .weights = "w_nearest", .metrics = "smd", na.rm = TRUE)))
  expect_false(isTRUE(all.equal(
    res$balance$estimate[res$balance$method == "w_nearest"],
    onsub$estimate[onsub$method == "w_nearest"])))
})


test_that("exact matching needs discrete covariates; a failing method is named", {
  d <- psm_data()

  ex <- get_PSM(d, treat = "z", adj_var = "x2", method = "exact",
                balance = FALSE)
  expect_identical(ex$stats$n_pairs, 2L)
  expect_identical(ex$stats$n, nrow(d))

  # a continuous covariate has no exact matches, and the error says which
  # method it was that failed, since the whole call stops there
  expect_error(
    get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
            method = c("nearest", "exact")),
    "matchit(method = \"exact\") failed", fixed = TRUE)
})


test_that("ATC matches each control, and the set counts mean what they say", {
  d <- psm_data()

  atc <- suppressWarnings(
    get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
            estimand = "ATC"))
  w <- atc$data$w_nearest
  expect_identical(atc$stats$estimand, "ATC")
  expect_true(all(w[d$z == 0] %in% c(0, 1)))
  expect_equal(sum(w[d$z == 1]), atc$stats$n_ctrl)

  att <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE)
  expect_identical(att$stats$n_pairs, att$stats$n_treat)   # 1:1 pairs
  expect_equal(att$stats$ess_pct, 1)                        # equal weights
  expect_equal(att$stats$w_cv, 0)

  sub <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                 method = "subclass")
  expect_identical(sub$stats$n_pairs, 6L)                   # MatchIt default
})


test_that("estimand is validated per method, naming the offender", {
  d <- psm_data()

  expect_error(get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                       estimand = "ATE"),
               "not available for method \"nearest\"")
  expect_error(get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                       method = c("nearest", "full"), estimand = "ATE"),
               "\"nearest\"")

  skip_if_not_installed("optmatch")
  ok <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                method = "full", estimand = "ATE")
  expect_identical(ok$stats$estimand, "ATE")
})


test_that("the caliper narrows the matched set and is recorded", {
  d    <- psm_data()
  none <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE)
  tight <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                   caliper = 0.05)

  expect_lt(tight$stats$n, none$stats$n)
  expect_null(attr(none, "analysis")$caliper)
  expect_identical(attr(tight, "analysis")$caliper, 0.05)
  expect_identical(attr(tight, "analysis")$specs$nearest$caliper, 0.05)
})


test_that("knobs reach only the methods that accept them", {
  skip_if_not_installed("optmatch")
  d   <- psm_data()
  # full matching ignores `ratio`; asking for both methods must still work
  res <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                 method = c("nearest", "full"), ratio = 2)
  sp <- attr(res, "analysis")$specs
  expect_identical(sp$nearest$ratio, 2)
  expect_null(sp$full$ratio)

  # optimal ignores caliper and replace, subclass ignores all three
  sp2 <- attr(get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                      method = c("optimal", "subclass"),
                      ratio = 2, caliper = 0.2, replace = TRUE),
              "analysis")$specs
  expect_identical(sp2$optimal$ratio, 2)
  expect_null(sp2$optimal$caliper)
  expect_null(sp2$optimal$replace)
  expect_identical(sp2$subclass, list())
})


test_that("no method is handed an argument MatchIt would ignore", {
  # MatchIt warns "the argument `x` is not used with method `y`" once per
  # ignored argument; a multi-method call would otherwise be noisy, and the
  # warning is the only signal that a knob silently did nothing.
  skip_if_not_installed("optmatch")
  d <- psm_data()

  expect_no_warning(
    get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
            method = c("nearest", "optimal", "full", "subclass"),
            ratio = 2, caliper = 0.2, replace = TRUE))

  # cem and exact ignore the score entirely, so it is not passed to them
  expect_no_warning(
    get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
            method = "cem"))
  expect_false("distance" %in%
                 names(attr(get_PSM(d, treat = "z", adj_var = psm_adj,
                                    balance = FALSE,
                                    method = "cem"), "analysis")$specs$cem))
})


test_that("factor and character exposures take the second level as treated", {
  d    <- psm_data()
  d$zf <- factor(d$z, levels = c(0, 1), labels = c("no", "yes"))
  d$zc <- ifelse(d$z == 1, "yes", "no")

  num <- get_PSM(d, treat = "z",  adj_var = psm_adj, balance = FALSE)
  fct <- get_PSM(d, treat = "zf", adj_var = psm_adj, balance = FALSE)
  chr <- get_PSM(d, treat = "zc", adj_var = psm_adj, balance = FALSE)

  expect_equal(fct$data$w_nearest, num$data$w_nearest)
  expect_equal(chr$data$w_nearest, num$data$w_nearest)

  # the column comes back as supplied; the 0/1 coding stays internal, and
  # match.data() sees that internal coding
  expect_s3_class(fct$data$zf, "factor")
  expect_type(chr$data$zc, "character")
  expect_setequal(unique(MatchIt::match.data(fct$fit$nearest)$zf), c(0L, 1L))

  # the arm taken as treated is recorded and printed
  expect_identical(attr(num, "analysis")$treated, "1")
  expect_identical(attr(fct, "analysis")$treated, "yes")
  expect_match(paste(utils::capture.output(print(chr)), collapse = "\n"),
               "treat = zc (treated = yes)", fixed = TRUE)

  # alphabetical order puts "control" second; the label makes that visible
  d$zr <- ifelse(d$z == 1, "active", "control")
  rev  <- suppressWarnings(
    get_PSM(d, treat = "zr", adj_var = psm_adj, balance = FALSE))
  expect_identical(attr(rev, "analysis")$treated, "control")
  expect_identical(.psw_treat(d$zr, "zr")$z, 1L - d$z)

  # MatchIt on its own cannot take a character exposure at all
  expect_error(MatchIt::matchit(zc ~ x1 + x2 + x3, data = d,
                                method = "nearest"),
               "y values must be")
})


test_that("a supplied score bypasses the model", {
  d      <- psm_data()
  d$myps <- psm_ps(d)

  res <- get_PSM(d, treat = "z", adj_var = psm_adj, ps = "myps",
                 balance = FALSE)
  expect_equal(res$data$ps, d$myps)
  expect_true(is.na(attr(res, "analysis")$ps_method))
  expect_error(get_PSM(d, treat = "z", adj_var = psm_adj, ps = "myps",
                       balance = FALSE, ps_method = "gbm"),
               "does not apply when `ps` is supplied")
})


test_that("invalid input is rejected rather than absorbed", {
  d <- psm_data()

  expect_error(get_PSM(d, treat = "x1", adj_var = c("x2", "x3"),
                       balance = FALSE), "must be 0/1, logical")
  expect_error(get_PSM(d, treat = "z", balance = FALSE), "`adj_var` is required")
  expect_error(get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                       method = "nope"), "should be one of")
  expect_error(get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                       match_args = list(normalize = TRUE)), "may not set")
  expect_error(get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                       caliper = -1), "single positive number")
  expect_error(get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE,
                       ratio = 0), "at least 1")

  clash <- d
  clash$w_nearest <- 1
  expect_error(get_PSM(clash, treat = "z", adj_var = psm_adj,
                       balance = FALSE), "already has column")
})


test_that("print reports the table and echoes a plt_PSM call", {
  d   <- psm_data()
  res <- get_PSM(d, treat = "z", adj_var = psm_adj, balance = FALSE)
  out <- paste(utils::capture.output(print(res)), collapse = "\n")

  expect_match(out, "<psm_res>")
  expect_match(out, "plt_PSM: type = \"ess\"")
  expect_match(out, "nearest")
})
