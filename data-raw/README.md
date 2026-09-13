# 数据

> ## 本仓库是公开仓库，这个目录里除本文件和 `synthetic/` 外，全部不进版本库。
>
> · **MIMIC-IV** 受 PhysioNet DUA 约束，明文禁止再分发
> · **院内 HIS / 自有队列** 含患者信息，任何形式都不得入库
> · 误提交后即使 revert，历史里仍然可以取回

`.gitignore` 对本目录采用「默认拒绝 + 白名单」，`.githooks/pre-commit`
再拦一道。两道都绕过去需要 `--no-verify`，那时候请确认你知道自己在做什么。

## 目录

| 目录 | 内容 | 入库 |
|---|---|---|
| `raw/` | 原始数据：NHEFS、MIMIC-IV 抽数结果、HIS 导出 | 否 |
| `derived/` | 清洗后、克隆扩展后、插补后的分析数据集 | 否 |
| `synthetic/` | 小体量合成 / 完全脱敏的演示数据 | **是** |

## 怎么拿到数据

| 数据 | 来源 | 说明 |
|---|---|---|
| NHEFS | https://miguelhernan.org/whatifbook | 随《What If》免费下载，放 `raw/nhefs/` |
| causaldata | `install.packages("causaldata")` | R 包直接带，不用落盘 |
| MIMIC-IV | https://physionet.org/content/mimiciv/ | 需 CITI 培训 + DUA 签署，审批 2–4 周 |
| MIMIC 抽数代码 | https://github.com/MIT-LCP/mimic-code | 官方 SQL |
| 自有队列 | 院内 | 伦理批件 + 脱敏后放 `raw/`，**脱敏脚本本身可以入库** |

## 复现约定

`raw/` 一旦落盘就**锁定不再改动**。任何清洗、派生都走 `weeks/` 或 `R/` 下的脚本
写进 `derived/`，这样别人拿到代码 + 自己申请的数据能跑出同样结果。
