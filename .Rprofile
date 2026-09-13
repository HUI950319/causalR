# causalR · 项目级 .Rprofile
#
# 只有当工作目录 = 项目根时才会被加载（用 causalR.Rproj 打开即可）。
#
# 【重要】项目级 .Rprofile 会屏蔽用户级的 ~/.Rprofile（Windows 下是
# Documents/.Rprofile），所以那边的 btw / mcptools 配置必须在这里重复一遍，
# 否则在本项目里 BTW MCP 连不上 RStudio 会话。

local({

  ## ---- CRAN 镜像 -------------------------------------------------------
  ## 仅在未配置时生效。clone 本仓库的人可自行改成就近镜像。
  repos <- getOption("repos")
  if (is.null(repos[["CRAN"]]) || identical(unname(repos[["CRAN"]]), "@CRAN@")) {
    repos["CRAN"] <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
    options(repos = repos)
  }

  if (!interactive()) return(invisible(NULL))

  ## ---- btw / mcptools --------------------------------------------------
  ## 这里【只设连接参数，不自动注册会话】。
  ## 自动注册会在 .Rprofile 被重复执行时留下多个 nanonext listener，
  ## 表现为「能列出会话、但工具调用永久挂起」。注册走下面的 btw_connect()。
  has_btw <- requireNamespace("btw", quietly = TRUE) &&
    requireNamespace("mcptools", quietly = TRUE)

  if (has_btw) {
    state <- get("the", envir = asNamespace("mcptools"))
    state$socket_url <- "tcp://127.0.0.1:4777"
    options(btw.run_r.enabled = TRUE)
  }

  ## ---- 会话工具（挂搜索路径，不污染 globalenv）-------------------------
  tools <- new.env()

  ## 向 Claude 注册本 R 会话。重复调用是安全的：已注册时只提示，不再监听。
  tools$btw_connect <- function() {
    if (!has_btw) {
      message("btw / mcptools 未安装，无法注册会话。")
      return(invisible(NULL))
    }
    state <- get("the", envir = asNamespace("mcptools"))
    if (!is.null(state$session)) {
      message("本会话已注册（slot ", state$session, "）。",
              "要重连请先 Session > Restart R，再跑 btw_connect()。")
      return(invisible(state$session))
    }
    btw::btw_mcp_session()
    message("已监听 ", sprintf("%s%d", state$socket_url, state$session),
            "，Claude 现在可以看到这个会话了。")
    invisible(state$session)
  }

  ## 本会话的 BTW 状态，连不上时先看这个。
  tools$btw_status <- function() {
    if (!has_btw) return(c(btw = "未安装"))
    state <- get("the", envir = asNamespace("mcptools"))
    c(socket_url = state$socket_url,
      slot       = if (is.null(state$session)) "未注册" else as.character(state$session),
      run_r      = as.character(isTRUE(getOption("btw.run_r.enabled"))),
      wd         = getwd())
  }

  ## 手动 source(".Rprofile") 调试时不叠加挂载
  if ("causalR-tools" %in% search()) detach("causalR-tools", character.only = TRUE)
  attach(tools, name = "causalR-tools", warn.conflicts = FALSE)

  ## ---- 自动注册到 Claude ----------------------------------------------
  ## 注：这一步偏离了 .claude/docs/rstudio-btw-workflow.md 里「不自动注册」的
  ## 约定。那条约定防的是重复调用留下多个 nanonext listener，而 btw_connect()
  ## 内部已有幂等保护，所以这里可以安全地自动跑。
  ## 不想开机自动连的话删掉下面这一行即可，btw_connect() 仍可手动调用。
  slot <- if (has_btw) suppressMessages(try(tools$btw_connect(), silent = TRUE)) else NULL

  conn <- if (is.numeric(slot)) {
    paste0("已连 Claude（slot ", slot, "）")
  } else if (has_btw) {
    "未连上 Claude —— 手动试 btw_connect()"
  } else {
    "btw / mcptools 未安装"
  }

  ## ---- 启动提示 --------------------------------------------------------
  cat("causalR · 24 周因果推断学习\n",
      "  计划 learn/00-roadmap.md   笔记 learn/notes/   本周代码 weeks/wNN-*/\n",
      "  ", conn, "   自检：btw_status()\n",
      sep = "")
})
