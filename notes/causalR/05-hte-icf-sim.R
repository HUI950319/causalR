# =============================================================================
# 05-hte-icf-sim.R -- 手册第 5 页（get_hte_icf()）的模拟与计时
# =============================================================================
#
# 生成两个文件，页面读取它们作表：
#   data/hte_icf_sim.csv     四个场景下的识别结果（每行一次模拟）
#   data/hte_icf_timing.csv  不同样本量下一次默认调用的耗时
#
# 全部使用 get_hte_icf() 的默认设置（20 个森林 x 200 棵树 x 5 折），
# 2026-09-28 在 Windows R 4.4.3、grf 2.6.1、Intel Core Ultra 9 275HX
# （24 核）、128 GB 内存上运行。计时部分须单独运行，不与其他任务并行。
#
# 用法（在项目根目录）：
#   Rscript notes/causalR/05-hte-icf-sim.R sim      # 识别结果，约 1 小时
#   Rscript notes/causalR/05-hte-icf-sim.R timing   # 耗时，约半小时
# 也可以只跑某个场景：Rscript ... sim inter null
# =============================================================================

library(causalR)

# ---- 数据生成 ----------------------------------------------------------------
# 每个数据集都带真值列 tau：连续结局是个体效应，生存结局是 S1(12) - S0(12)。
# tau 不进 adj_var，只用来核对估计集上的亚组效应。

gen <- function(scenario, n, seed) {
  set.seed(seed)
  if (scenario == "readme") {
    # iCF README 的模拟（Setoguchi 等的设定加上 W:X1:X3 交互）
    cor1 <- function(x, rho)
      rho * (x - mean(x)) / sd(x) + sqrt(1 - rho^2) * rnorm(length(x))
    bin <- function(x) as.numeric(x > mean(x))
    X2 <- rnorm(n); X4 <- rnorm(n); X7 <- rnorm(n); X10 <- rnorm(n)
    X1 <- bin(rnorm(n)); X3 <- bin(rnorm(n))
    X5 <- bin(cor1(X1, 0.2)); X6 <- bin(cor1(X2, 0.9))
    X8 <- bin(cor1(X3, 0.2)); X9 <- bin(cor1(X4, 0.9))
    z <- rbinom(n, 1, plogis(0.8 * X1 - 0.25 * X2 + 0.6 * X3 - 0.4 * X4 -
                               0.8 * X5 - 0.5 * X6 + 0.7 * X7))
    tau <- -0.4 + 0.3 * X1 + 0.4 * X3 + 0.4 * X1 * X3
    y <- -3.85 + 0.3 * X1 - 0.36 * X2 - 0.73 * X3 - 0.2 * X4 + 0.71 * X8 -
      0.19 * X9 + 0.26 * X10 + z * tau + 0.2 * X1 * X3 + rnorm(n)
    return(data.frame(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, z, y, tau))
  }
  # 其余三个场景共用 6 个协变量：X1、X3、X6 二分类，X2、X4、X5 正态
  d <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rnorm(n),
                  X3 = rbinom(n, 1, 0.5), X4 = rnorm(n),
                  X5 = rnorm(n), X6 = rbinom(n, 1, 0.3))
  d$z <- rbinom(n, 1, plogis(0.4 * d$X1 - 0.3 * d$X2))
  if (scenario == "surv") {
    # 治疗把 X1 = X3 = 1 者的风险乘以 exp(-1.2)，其余人无效应
    rate <- 0.05 * exp(-1.2 * d$z * d$X1 * d$X3)
    ev   <- rexp(n, rate)
    cens <- rexp(n, 0.02)
    d$time <- pmin(ev, cens)
    d$DSS  <- as.integer(ev <= cens)
    d$tau  <- exp(-0.6 * exp(-1.2 * d$X1 * d$X3)) - exp(-0.6)
    return(d)
  }
  d$tau <- if (scenario == "inter") 2 * d$X1 * d$X3 else rep(0.5, n)
  d$y   <- d$X2 + d$z * d$tau + rnorm(n)
  d
}

# ---- 一次调用，整理成一行 ----------------------------------------------------

one <- function(scenario, n, seed, rule_args = list()) {
  d  <- gen(scenario, n, seed)
  xv <- grep("^X", names(d), value = TRUE)
  t0 <- Sys.time()
  r <- if (scenario == "surv")
    get_hte_icf(d, "z", xv, time = 12, seed = seed, rule_args = rule_args)
  else get_hte_icf(d, "z", xv, surv = "y", seed = seed, rule_args = rule_args)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  a <- attr(r, "analysis")

  # 关掉校准闸门、只看交叉验证时会选哪个深度（与 get_hte_icf() 内部同一规则）
  txt <- c("", r$vote$partition)
  ok  <- r$cv$n_leaf > 1 | r$cv$depth == 0
  cv_sel   <- match(txt[which.min(ifelse(ok, r$cv$cv_loss, Inf))], txt)
  vars_of  <- function(s) paste(sort(unique(unlist(
    regmatches(s, gregexpr("X[0-9]+", s))))), collapse = "+")

  # 估计集上每条规则的真值：该组患者 tau 的平均
  est   <- r$est$data
  truth <- tapply(est$tau, est$.rule, mean)[r$rules$rule]
  cover <- r$rules$conf.low <= truth & truth <= r$rules$conf.high

  data.frame(
    scenario = scenario, n = n, seed = seed,
    screen = a$rule_args$screen,
    depth = a$depth_selected, n_leaf = nrow(r$rules),
    vars = vars_of(r$rules$rule),
    depth_cv = r$cv$depth[cv_sel], vars_cv = vars_of(txt[cv_sel]),
    calib_p = a$calibration_p, gated = a$gated,
    share = if (a$depth_selected > 0)
      r$vote$share[r$vote$depth == a$depth_selected] else NA_real_,
    p_inter = r$rules$p_inter[1L],
    n_cover = sum(cover, na.rm = TRUE), max_err = max(abs(r$rules$estimate - truth)),
    rules = paste(r$rules$rule, collapse = " | "),
    secs = round(secs, 1), stringsAsFactors = FALSE)
}

# 每次调用后追加写入，跑到一半也能看到进度
run <- function(jobs, file) {
  for (i in seq_len(nrow(jobs))) {
    j <- jobs[i, ]
    row <- one(j$scenario, j$n, j$seed,
               rule_args = if (isFALSE(j$screen)) list(screen = FALSE) else list())
    utils::write.table(row, file, sep = ",", row.names = FALSE,
                       col.names = !file.exists(file), append = file.exists(file))
    message(sprintf("%s n=%d seed=%d screen=%s: depth %d, %s (%.0f s)",
                    j$scenario, j$n, j$seed, j$screen, row$depth, row$rules,
                    row$secs))
  }
}

args <- commandArgs(trailingOnly = TRUE)
what <- if (length(args)) args[1L] else "sim"
pick <- args[-1L]

if (what == "sim") {
  jobs <- rbind(
    data.frame(scenario = "inter",  n = 1600,  seed = 1:20, screen = TRUE),
    data.frame(scenario = "null",   n = 1600,  seed = 1:20, screen = TRUE),
    data.frame(scenario = "readme", n = 3000,  seed = 1:10, screen = TRUE),
    data.frame(scenario = "readme", n = 10000, seed = 1:10, screen = TRUE),
    data.frame(scenario = "surv",   n = 3000,  seed = 1:10, screen = TRUE),
    data.frame(scenario = "surv",   n = 3000,  seed = 1:10, screen = FALSE))
  if (length(pick)) jobs <- jobs[jobs$scenario %in% pick, ]
  run(jobs, "data/hte_icf_sim.csv")
}

if (what == "timing") {
  jobs <- rbind(
    data.frame(scenario = "inter", n = c(1000, 2000, 5000, 10000, 20000, 50000),
               seed = 1, screen = TRUE),
    data.frame(scenario = "surv",  n = c(1000, 2000, 5000, 10000, 20000),
               seed = 1, screen = TRUE))
  if (length(pick)) jobs <- jobs[jobs$scenario %in% pick, ]
  run(jobs, "data/hte_icf_timing.csv")
}
