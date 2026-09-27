# causalR (development version)

* For a survival probability, `get_hte()` passes patients followed beyond
  `time` to grf as events just after it, as grf itself does for RMST. The
  estimand is unchanged, but with rare events grf's nuisance survival forest
  could not split at all -- one leaf holding every patient -- which was slow
  and ignored the covariates: at 130,923 patients with 1,043 events, 41 s
  instead of 236 s, with a CATE correlated 0.82 across seeds instead of
  0.44. In 200 simulated data sets each with 1.5% and 24% events, bias,
  standard errors and 95% coverage did not change. Survival estimates
  change numerically. `plt_hte_rate()` refits the same way, and its new
  `max_n` (default `10000`, after `train_frac`) sets the most patients the
  forest CATE is learnt and evaluated on; `Inf` uses every patient, which
  rare events need.

* Large samples run faster. `get_hte()` grows 200 trees instead of grf's
  2000 when more than 10,000 rows are analysed and `grf_args` sets no
  `num.trees`, with a message: at 100,000 survival rows 13 s instead of
  83 s, the ATE unchanged. Above 10,000 patients `plt_hte_rate()` learns and
  evaluates the forest CATE on a random 10,000 of them (whole clusters) with
  500 trees, with a message: 2 s instead of a minute at 100,000, at the
  cost of power, as the test rests on about 5,000 held-out patients.
  `plt_hte_dep()` draws at most 2000 evenly spaced patients in the `"cate"`
  layer; its loess line, still fitted to every patient, skips loess's exact
  O(n^2) trace, which `se = FALSE` never used, so the line is unchanged:
  2 s instead of 122 s at 100,000 patients.

* New `get_bal()` draws the covariate balance of several propensity score
  schemes on one love plot, styled like `RegR::get_ps()`. `methods` takes
  shorthands (`"PSM"`, `"ATE"`, `"ATT"`, `"ATC"`, `"ATO"`, `"ATM"`, `"EW"`)
  or a named list of up to 14 schemes, each `design = "matching"` or
  `"weighting"` plus `get_PSM()` / `get_PSW()` arguments under their own
  names. Every scheme is standardised by the unadjusted pooled SD
  (`cobalt::bal.tab()`), the same denominator halfmoon gives `plt_PSM()` /
  `plt_PSW()` except for its variance over n rather than n - 1. `cat_smd`
  sets how a factor with three or more levels is summarised: `"overall"`
  (default) gives one unsigned Yang & Dalton (2012) row, the number
  `gtsummary::add_difference()` reports; `"level"` one row per level, as
  cobalt does. `tbl = TRUE` adds `$tbl`, a gtsummary `tbl_merge` with one
  spanner per sample (both arms, gtsummary's own SMD, p-value; schemes
  through `tbl_svysummary()`); gtsummary and survey join Suggests. Two save
  lists, in place of the usual `save`: `save_plt` (`filename`, `width`,
  `height`) writes the plot through `RegR::save_plt()`, and `save_tbl`
  (`path`, `title`, `note`, ...) writes the table through `RegR::save_tb()`,
  building it even when `tbl = FALSE`. Returns
  `list(plt, balance, data)`, plus `tbl`; cobalt, already a
  WeightIt dependency, joins Imports. Plot settings live in `love_args`
  (`threshold`, `colors`, `shapes`, `size`, `line`, `var_order`,
  `base_size`, `ref_color`, `legend_position`, `legend_justification`).
  `var_names`, default `RegR::name_map_seer` with your labels merged in
  first, relabels the plot and `$tbl` but not `$balance`; two covariates
  under one label are an error, as they would share a plot row, and so is
  a scheme labelled `"Un"`, cobalt's name for the unadjusted column.
  Weighting schemes that differ only in `estimand` are fitted in one
  `get_PSW()` call, the default six in one instead of six. `cores`
  (default `NULL`, automatic) builds the `$tbl` tables in parallel, forked
  on Linux and WSL and on a PSOCK cluster on Windows: at 20,000 rows under
  WSL `get_bal(tbl = TRUE)` takes 15 s instead of 49 s, the table
  unchanged. `parallel` joins Imports.

* `get_hte_select()` with `sel_metric` now selects only a step whose training
  scores are nonconstant and which strictly beats a constant-effect baseline
  (0 for AUTOC/QINI/score SD/IQR; the best constant effect on the evaluation
  split for R-loss/DR-loss). Without such a step `selected` is `character(0)`,
  `analysis$best_step` is 0 and the combined plot omits the red optimum.
  Previously constant-effect data still selected the first candidate, so the
  result depended on candidate order. `analysis$sel_baseline` reports the
  reference; default `sel_metric = NULL` and manual `n_select` are unchanged.

* `get_hte_select()` reuses the top single-variable fit as forward step 1
  instead of refitting it, reducing model fits from 2p to 2p - 1. Results are
  unchanged; that fit's warnings are still recorded for both stages.

* `get_hte_select()` deals PS-route LASSO folds within treatment-arm by
  event (survival) or outcome-class (binary) strata, and stops before fitting
  when any training fold would lack an arm, events or an outcome class.
  Sparse events previously could share one fold, failing binary fits and
  silently returning empty Cox models. Fold IDs therefore differ from earlier
  versions for the same seed.

* `get_hte_select()` renames `rank_metric` to `imp_metric` (variable importance)
  and `select_metric` to `sel_metric` (variable selection), including the
  corresponding `analysis` fields. Defaults and calculations are unchanged;
  callers using the old argument names must update them.
* `get_hte_select()` adds independent `imp_metric` and `sel_metric`
  controls, benefit-score IQR, held-out AUTOC/QINI, and continuous-outcome
  R-loss/DR-loss. `eval_args` controls a fixed training/evaluation split and
  a shared GRF evaluator; matched pairs stay together. Selection minimizes
  losses and maximizes other metrics, with manual `n_select` taking precedence.
  All plots follow the chosen metrics; defaults retain manual SD screening.
* `get_hte_select()` adds signed `score_mean` and `plots$combined`, aligning
  ranked SD bars with cumulative mean scores on a separate top axis. Red marks
  the first mean-score maximum without automatically selecting variables;
  a different `n_select` is marked separately. The example has 40 candidates.
* `get_hte_select()` ranks candidates by personalized benefit-score SD and
  returns the full fixed-order accumulation path for survival, continuous and
  binary outcomes, with existing 1:1 matching or fixed propensity scores.
* `plt_hte_dep()` rejects unused heat-map `ylim` and invalid non-integer or
  non-finite grid sizes and spline degrees; `max_n = Inf` remains supported.
* `plt_hte_cate(overall = FALSE)` skips unused arm-score calculations when no
  subgroup is requested, including when individual CATE intervals are drawn.
* PDP and heat-map predictions use bounded matrix batches instead of expanding
  the entire grid-by-patient matrix in memory, retaining the same estimates.
* GATES groups, population shares and mean CATE use the forest's analysis
  weights, matching TOC/Qini; zero-weight observations are excluded.
* GATES top-bottom contrasts include covariance from shared clusters in
  their standard errors, confidence intervals and p-values.
* `plt_hte_rate()` separates whole clusters when learning a CATE ranking,
  requiring at least two clusters in both training and evaluation samples.
* `plt_hte_sub()` rejects relative effects for weighted or clustered forests,
  matching the restrictions of `get_hte()`.
* `get_hte()` skips unused spline prediction grids when assembling importance
  tests; `plt_hte_dep()` still constructs the full curves on demand.
* `get_hte()` accepts single-level factor, character and logical covariates
  after complete-case filtering, retaining their missing values in the design.
* `get_hte()` subgroup `cate_mean` combines target-population weights with
  sample weights or equal cluster weights, matching the analysis population.
* `get_hte(estimand = "ATO")` accepts boundary propensities without losing
  valid overlap estimates; unavailable ordinary AIPW diagnostics are marked NA.
* `get_hte()` retains average effects when spline or calibration diagnostics
  are degenerate, returning unavailable diagnostics with a warning.
* `get_hte()` and `plt_hte_sub()` account for shared clusters when testing
  differences between subgroup effects, retaining grf's marginal variances.
* `get_hte()` and `plt_hte_dep()` now use the forest's observation weights
  and clusters for covariate heterogeneity tests and doubly robust curves.
