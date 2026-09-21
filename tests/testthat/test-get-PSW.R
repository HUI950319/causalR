# Weights are computed in-package; propensity / halfmoon are Suggests and the
# tests that need them skip when absent.

psw_data <- function(n = 400L) {
  set.seed(20260921)
  d <- data.frame(x1 = stats::rnorm(n),
                  x2 = stats::rbinom(n, 1, 0.4),
                  x3 = stats::runif(n))
  lp <- -0.3 + 0.8 * d$x1 - 0.6 * d$x2 + 1.1 * d$x3
  d$z <- stats::rbinom(n, 1, stats::plogis(lp))
  d
}

psw_adj <- c("x1", "x2", "x3")

stats_cols <- c("estimand", "n", "n_treat", "n_ctrl", "ess", "ess_treat",
                "ess_ctrl", "ess_pct", "w_min", "w_max", "w_mean", "w_sd",
                "w_cv", "smd_max", "smd_over")

all_wcols <- c("w_ate", "w_att", "w_atc", "w_ato", "w_atm", "w_ew")

# Weighted mean of `v` over the units selected by `g`.
wmean <- function(v, w, g) sum(v[g] * w[g]) / sum(w[g])


test_that("get_PSW returns the documented structure", {
  d   <- psw_data()
  res <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE)

  expect_s3_class(res, "psw_res")
  expect_named(res, c("data", "stats", "balance", "fit"))
  expect_null(res$balance)
  expect_s3_class(res$fit, "weightit")

  expect_named(res$stats, stats_cols)
  expect_identical(res$stats$estimand,
                   c("ATE", "ATT", "ATC", "ATO", "ATM", "EW"))
  expect_true(all(all_wcols %in% names(res$data)))
  expect_true(all(c("ps", ".trimmed") %in% names(res$data)))
  expect_false(any(res$data$.trimmed))
  expect_identical(nrow(res$data), nrow(d))

  a <- attr(res, "analysis")
  expect_identical(a$treat, "z")
  expect_identical(a$wcols, all_wcols)
  expect_identical(a$n_trimmed, 0L)
})


test_that("the six weights match propensity::wt_*()", {
  skip_if_not_installed("propensity")

  d   <- psw_data()
  res <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE)
  e   <- res$data$ps
  z   <- d$z

  up <- list(
    w_ate = propensity::wt_ate(e, z),
    w_att = propensity::wt_att(e, z),
    w_atc = propensity::wt_atc(e, z),
    w_ato = propensity::wt_ato(e, z),
    w_atm = propensity::wt_atm(e, z),
    w_ew  = propensity::wt_entropy(e, z))

  for (nm in names(up))
    expect_equal(res$data[[nm]], as.numeric(up[[nm]]),
                 tolerance = 1e-12, info = nm)
})


test_that("the five shared estimands equal WeightIt's own weights", {
  # With method = "glm" the score does not depend on the estimand, so our
  # tilting functions must reproduce weightit()'s weights exactly -- same
  # values, not merely proportional. EW has no WeightIt counterpart.
  d    <- psw_data()
  res  <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE)
  form <- stats::as.formula(paste("z ~", paste(psw_adj, collapse = " + ")))

  for (e in c("ATE", "ATT", "ATC", "ATO", "ATM")) {
    o <- suppressMessages(
      WeightIt::weightit(form, data = d, method = "glm", estimand = e))
    expect_equal(res$data$ps, unname(as.numeric(o$ps)),
                 tolerance = 1e-12, info = e)
    expect_equal(res$data[[paste0("w_", tolower(e))]], unname(o$weights),
                 tolerance = 1e-10, info = e)
  }

  expect_error(
    suppressMessages(
      WeightIt::weightit(form, data = d, method = "glm", estimand = "EW")),
    "allowable estimand")
})


test_that("score-adaptive backends deliberately share one ATE-fitted score", {
  # cbps refits the score per estimand, so get_PSW() and weightit() are not
  # the same estimator there. This pins that the divergence comes from the
  # score, not from the weight formula. (No skip_on_cran(): under test_dir()
  # NOT_CRAN is unset, so it would skip locally and never run at all.)
  d    <- psw_data()
  form <- stats::as.formula(paste("z ~", paste(psw_adj, collapse = " + ")))

  ps_ate <- suppressMessages(
    WeightIt::weightit(form, data = d, method = "cbps", estimand = "ATE"))$ps
  ps_ato <- suppressMessages(
    WeightIt::weightit(form, data = d, method = "cbps", estimand = "ATO"))$ps
  expect_false(isTRUE(all.equal(unname(as.numeric(ps_ate)),
                                unname(as.numeric(ps_ato)))))

  res <- suppressMessages(
    get_PSW(d, treat = "z", adj_var = psw_adj, method = "cbps",
            estimand = c("ATE", "ATO"), balance = FALSE))
  expect_equal(res$data$ps, unname(as.numeric(ps_ate)), tolerance = 1e-10)
  # and the weights are still our closed form applied to that one score
  expect_equal(res$data$w_ato,
               ifelse(d$z == 1, 1 - res$data$ps, res$data$ps),
               tolerance = 1e-12)
})


test_that("overlap weights balance the covariate means exactly", {
  # Li, Morgan & Zaslavsky (2018): with a logistic score fitted on exactly
  # these covariates, ATO weights equalise their two-arm means to machine
  # precision. No other weight in the family does.
  d   <- psw_data()
  res <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE)
  tr  <- d$z == 1
  ct  <- d$z == 0

  for (v in psw_adj)
    expect_equal(wmean(d[[v]], res$data$w_ato, tr),
                 wmean(d[[v]], res$data$w_ato, ct),
                 tolerance = 1e-8, info = v)

  worst <- max(vapply(psw_adj, function(v)
    abs(wmean(d[[v]], res$data$w_ate, tr) -
        wmean(d[[v]], res$data$w_ate, ct)), numeric(1)))
  expect_gt(worst, 1e-3)
})


test_that("each weight has its defining property", {
  d   <- psw_data()
  res <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE)
  w   <- res$data
  tr  <- d$z == 1
  ct  <- d$z == 0
  e   <- w$ps

  expect_equal(w$w_att[tr], rep(1, sum(tr)))          # ATT leaves treated as is
  expect_equal(w$w_atc[ct], rep(1, sum(ct)))          # ATC leaves controls as is
  expect_true(all(w$w_atm <= 1 + 1e-12))              # matching weights cap at 1
  expect_equal(w$w_ato[tr], 1 - e[tr])                # overlap weight is 1 - e
  expect_equal(w$w_ato[ct], e[ct])
  expect_true(all(w$w_ate >= 1 - 1e-12))              # IPTW is at least 1

  # The entropy tilt peaks at e = 0.5 with h = log 2
  expect_equal(max(.psw_tilt("EW", seq(0.001, 0.999, by = 0.001))), log(2),
               tolerance = 1e-6)
})


test_that("stabilize recentres the ATE weight and leaves the others alone", {
  d  <- psw_data()
  a  <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE)
  b  <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE,
                stabilize = TRUE)
  p1 <- mean(d$z)

  expect_equal(mean(b$data$w_ate), 1, tolerance = 0.05)
  expect_equal(b$data$w_ate, a$data$w_ate * ifelse(d$z == 1, p1, 1 - p1))
  for (nm in setdiff(all_wcols, "w_ate"))
    expect_equal(b$data[[nm]], a$data[[nm]], info = nm)

  # No selected estimand can use it -> hard error, not a silent no-op
  expect_error(
    get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE,
            estimand = c("ATO", "ATM"), stabilize = TRUE),
    "defined for")
})


test_that("trimming keeps rows, blanks the weights and refits the score", {
  d <- psw_data()

  keep_all <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE)
  res <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE,
                 trim_args = list(method = "ps", lower = 0.3, upper = 0.7))

  expect_identical(nrow(res$data), nrow(d))
  out <- keep_all$data$ps < 0.3 | keep_all$data$ps > 0.7
  expect_identical(res$data$.trimmed, unname(out))
  expect_true(all(is.na(res$data$ps[out])))
  for (nm in all_wcols) expect_true(all(is.na(res$data[[nm]][out])))
  expect_identical(attr(res, "analysis")$n_trimmed, sum(out))
  expect_identical(res$stats$n[[1L]], sum(!out))

  # refit really re-estimates: the retained scores move
  no_refit <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE,
                      trim_args = list(method = "ps", lower = 0.3,
                                       upper = 0.7, refit = FALSE))
  expect_equal(no_refit$data$ps[!out], keep_all$data$ps[!out])
  expect_false(isTRUE(all.equal(res$data$ps[!out], no_refit$data$ps[!out])))

  # the Crump rule returns a symmetric window
  cr <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE,
                trim_args = list(method = "cr"))
  b <- attr(cr, "analysis")$trim_bounds
  expect_equal(b[[1L]], 1 - b[[2L]])
})


test_that("the Crump cut-off equals its alpha-by-alpha definition", {
  # The one-pass implementation is checked against the literal rule: the
  # smallest observed alpha whose retained set satisfies the criterion.
  crump_ref <- function(ps) {
    v <- 1 / (ps * (1 - ps))
    for (a in sort(unique(pmin(ps, 1 - ps)))) {
      if (a >= 0.5) break
      keep <- ps >= a & ps <= 1 - a
      if (1 / (a * (1 - a)) <= 2 * mean(v[keep])) return(a)
    }
    0
  }

  for (s in 1:5) {
    set.seed(s)
    e <- stats::plogis(stats::rnorm(2000, 0, 2.5))     # poor overlap
    expect_identical(.psw_crump(e), crump_ref(e), info = paste("seed", s))
    expect_gt(.psw_crump(e), 0)
  }
  set.seed(6)
  e <- stats::plogis(stats::rnorm(2000, 0, 0.5))       # good overlap
  expect_identical(.psw_crump(e), crump_ref(e))
  expect_identical(.psw_crump(rep(0.5, 10)), 0)         # nothing to cut
})


test_that("truncation clamps the score without dropping anyone", {
  d   <- psw_data()
  raw <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE)
  res <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE,
                 trunc_args = list(method = "ps", lower = 0.2, upper = 0.8))

  expect_false(any(res$data$.trimmed))
  expect_gte(min(res$data$ps), 0.2)
  expect_lte(max(res$data$ps), 0.8)
  expect_lt(max(res$data$w_ate), max(raw$data$w_ate))

  q <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE,
               trunc_args = list(method = "pctl", lower = 0.05, upper = 0.95))
  expect_equal(range(q$data$ps),
               unname(stats::quantile(raw$data$ps, c(0.05, 0.95))))
})


test_that("estimand selects columns and rows", {
  d   <- psw_data()
  res <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE,
                 estimand = c("EW", "ATO"))

  expect_identical(res$stats$estimand, c("ATO", "EW"))   # canonical order
  expect_true(all(c("w_ato", "w_ew") %in% names(res$data)))
  expect_false(any(c("w_ate", "w_att", "w_atc", "w_atm") %in% names(res$data)))
  expect_identical(attr(res, "analysis")$wcols, c("w_ato", "w_ew"))
})


test_that("the weightit object carries a readable call, not the data", {
  # do.call() would inline the function body and the data frame into $call;
  # the object is then mostly that, and printing the call dumps the data.
  d   <- psw_data()
  res <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE)
  cl  <- res$fit$call

  expect_identical(cl[[1L]], quote(WeightIt::weightit))
  expect_identical(cl$data, quote(data))
  expect_lt(sum(nchar(deparse(cl))), 200L)
  expect_lt(as.numeric(utils::object.size(res$fit)),
            20 * as.numeric(utils::object.size(d)))
  expect_no_error(summary(res$fit))
})


test_that("a supplied score bypasses the model", {
  d      <- psw_data()
  d$myps <- stats::fitted(stats::glm(z ~ x1 + x2 + x3, binomial, data = d))

  res <- get_PSW(d, treat = "z", adj_var = psw_adj, ps = "myps",
                 balance = FALSE)
  expect_null(res$fit)
  expect_equal(res$data$ps, unname(d$myps))
  expect_true(is.na(attr(res, "analysis")$method))

  expect_error(get_PSW(d, treat = "z", adj_var = psw_adj, ps = "myps",
                       balance = FALSE, method = "gbm"),
               "does not apply when `ps` is supplied")
  expect_error(get_PSW(d, treat = "z", adj_var = psw_adj, ps = "myps",
                       balance = FALSE,
                       trim_args = list(method = "pctl")),
               "needs a propensity model to refit")
})


test_that("factor and character exposures take the second level as treated", {
  skip_if_not_installed("halfmoon")
  d    <- psw_data()
  d$zf <- factor(d$z, levels = c(0, 1), labels = c("no", "yes"))
  d$zc <- ifelse(d$z == 1, "yes", "no")

  num <- get_PSW(d, treat = "z",  adj_var = psw_adj, estimand = "ATO")
  fct <- get_PSW(d, treat = "zf", adj_var = psw_adj, estimand = "ATO")
  chr <- get_PSW(d, treat = "zc", adj_var = psw_adj, estimand = "ATO")

  expect_equal(fct$data$w_ato, num$data$w_ato)
  expect_equal(chr$data$w_ato, num$data$w_ato)

  # the balance table is computed on the 0/1 coding, so its signs do not
  # follow the labelling, and the column comes back exactly as supplied
  expect_equal(fct$balance$estimate, num$balance$estimate)
  expect_equal(chr$balance$estimate, num$balance$estimate)
  expect_s3_class(fct$data$zf, "factor")
  expect_type(chr$data$zc, "character")

  # the arm taken as treated is recorded and printed
  expect_identical(attr(num, "analysis")$treated, "1")
  expect_identical(attr(fct, "analysis")$treated, "yes")
  expect_match(paste(utils::capture.output(print(chr)), collapse = "\n"),
               "treat = zc (treated = yes)", fixed = TRUE)

  # alphabetical order puts "control" second; the label makes that visible
  d$zr <- ifelse(d$z == 1, "active", "control")
  rev  <- get_PSW(d, treat = "zr", adj_var = psw_adj, balance = FALSE,
                  estimand = "ATT")
  expect_identical(attr(rev, "analysis")$treated, "control")
  expect_identical(rev$stats$n_treat, sum(d$z == 0))
})


test_that("the default score still equals stats::glm's", {
  # method = "glm" routes through WeightIt like every other backend; this pins
  # it to the logistic fit users expect, and `ps_args` to WeightIt's own knobs.
  d   <- psw_data()
  res <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE)
  expect_equal(
    res$data$ps,
    unname(stats::fitted(stats::glm(z ~ x1 + x2 + x3, stats::binomial(),
                                    data = d))),
    tolerance = 1e-10)

  probit <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE,
                    ps_args = list(link = "probit"))
  expect_equal(
    probit$data$ps,
    unname(stats::fitted(stats::glm(z ~ x1 + x2 + x3,
                                    stats::binomial("probit"), data = d))),
    tolerance = 1e-10)
  expect_false(isTRUE(all.equal(probit$data$ps, res$data$ps)))
})


test_that("invalid input is rejected rather than absorbed", {
  d <- psw_data()

  expect_error(get_PSW(d, treat = "x1", adj_var = c("x2", "x3"),
                       balance = FALSE),
               "must be 0/1, logical")
  expect_error(get_PSW(d, treat = "z", adj_var = "nope", balance = FALSE),
               "not found in `data`")
  expect_error(get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE,
                       estimand = "ATZ"), "should be one of")
  expect_error(get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE,
                       method = "ebal"), "`method` must be one of")
  expect_error(get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE,
                       trim_args = list(nope = 1)), "unknown field")
  expect_error(get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE,
                       ps_args = list(formula = z ~ x1)), "may not set")
  expect_error(get_PSW(d, treat = "z", balance = FALSE),
               "Supply `adj_var`")

  supplied <- d
  supplied$myps <- stats::fitted(
    stats::glm(z ~ x1 + x2 + x3, stats::binomial(), data = d))
  expect_error(get_PSW(supplied, treat = "z", ps = "myps", balance = TRUE),
               "`adj_var` is required for the balance table")

  bad <- d
  bad$bad_ps <- 0
  expect_error(get_PSW(bad, treat = "z", adj_var = psw_adj, ps = "bad_ps",
                       balance = FALSE), "strictly between 0 and 1")

  clash <- d
  clash$w_ate <- 1
  expect_error(get_PSW(clash, treat = "z", adj_var = psw_adj,
                       balance = FALSE), "already has column")
})


test_that("balance diagnostics fill smd_max and smd_over", {
  skip_if_not_installed("halfmoon")

  d   <- psw_data()
  res <- get_PSW(d, treat = "z", adj_var = psw_adj)

  expect_s3_class(res$balance, "data.frame")
  expect_true(all(c("variable", "method", "metric", "estimate") %in%
                    names(res$balance)))
  expect_setequal(unique(res$balance$method), c("observed", all_wcols))

  expect_false(anyNA(res$stats$smd_max))
  expect_false(anyNA(res$stats$smd_over))

  # ATO balances exactly, so it must have the smallest worst-case SMD
  expect_identical(res$stats$estimand[which.min(res$stats$smd_max)], "ATO")
  expect_lt(res$stats$smd_max[res$stats$estimand == "ATO"], 1e-6)
})


test_that("trimming does not turn the balance table into NA", {
  skip_if_not_installed("halfmoon")

  # check_balance() reports NA for any weight column holding an NA, even with
  # na.rm = TRUE, so the trimmed rows have to be excluded from it. Without
  # that, every smd_max below is silently NA.
  d   <- psw_data()
  res <- get_PSW(d, treat = "z", adj_var = psw_adj,
                 trim_args = list(method = "ps", lower = 0.3, upper = 0.7))

  expect_gt(attr(res, "analysis")$n_trimmed, 0L)
  expect_false(anyNA(res$balance$estimate))
  expect_false(anyNA(res$stats$smd_max))
  expect_identical(nrow(res$data), nrow(d))          # rows still aligned
  expect_lt(res$stats$smd_max[res$stats$estimand == "ATO"], 1e-6)
})


test_that("the overlap weight keeps the most information", {
  d   <- psw_data()
  res <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE)
  ess <- stats::setNames(res$stats$ess, res$stats$estimand)

  expect_identical(names(which.max(ess)), "ATO")
  expect_lt(ess[["ATE"]], ess[["ATO"]])
  expect_true(all(res$stats$ess <= res$stats$n + 1e-8))
})


test_that("print reports the table and echoes a plt_PSW call", {
  d   <- psw_data()
  res <- get_PSW(d, treat = "z", adj_var = psw_adj, balance = FALSE)
  out <- paste(utils::capture.output(print(res)), collapse = "\n")

  expect_match(out, "<psw_res>")
  expect_match(out, "plt_PSW: type = \"ess\"")
  expect_match(out, "ATO")
})
