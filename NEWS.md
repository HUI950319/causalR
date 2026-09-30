# causalR (development version)

* `predict.hte_tree()` assigns new patients to the original rules or
  terminal node ids using the training factor encoding and exact split
  boundaries, without refitting. Missing values and unseen levels are rejected.

* `get_hte_tree()` records its survival time grid and event compression
  in `analysis$time_grid`. It warns when at least ten positive events,
  comprising at least 10% of events up to the horizon, round to zero;
  the existing default and user-supplied grids are retained.

* `get_hte_tree()` reports control counts, follow-up support and an
  estimation status for each leaf. Its overall interaction test is now
  omitted when any leaf cannot be estimated, including in downstream plots.

* `get_hte_tree()` rejects failed logistic, Cox and non-finite node fits
  before using their scores. A failed root remains an unsplit tree;
  node-model failures are warned once and recorded in analysis metadata.

* `get_hte_tree(estimator = "tmle")` now keeps its stored overall and
  subgroup results, `plt_hte_sub()` and `plt_hte_cate()` on the same TMLE
  estimator as its node and rule tables.

* `plt_sens()` draws its native figures (the `"lm"` contour, `"tip"` and
  `"evalue"`) with `UtilsR::theme_my(base_rect_size = 1.5)`, the theme of
  `MLR::plt_bar_per()`, and gains `colors = c(main, highlight)`: the contour
  lines and critical line, or the point-estimate and confidence-limit
  curves. The DML contour takes it as the upstream `col.contour` /
  `col.thr.line`; the IV contour and `"extreme"` refuse it, as upstream fixes
  their colours. The Cox curves now default to firebrick and steelblue.

* `plt_sens()` gains `legend_position`. The tip and E-value legends now sit
  inside the panel by default -- top right for `"evalue"`; for `"tip"` the
  first corner the curves, null line and tipping-point line leave free, or
  below the panel when none is. The E-value points and the tipping point
  carry `ggrepel` labels of their values; the tipping point's replaces the
  text atop its line, which now rises from the axis only to the point, and
  the point takes the confidence-limit colour. The tip plot's curves,
  reference lines, point and label are drawn larger.
  causalR now imports ggrepel (already required through UtilsR) and needs
  ggplot2 >= 3.5.0 for inside legends.

* A `failure.times` grid in `grf_args` no longer biases survival estimates
  of `get_hte()`, `get_hte_icf()` and `plt_hte_rate()`. Patients followed
  past `time` reach grf as events at the first such follow-up, and grf
  moves every time down to the grid point at or before it: a grid with no
  point between `time` and that follow-up, such as one ending at `time`,
  counted them as deaths at `time` (effects off by up to 0.03 in
  simulations, standard errors 19% larger). That point is now added to the
  grid.

* For a survival outcome with more than 100 distinct times up to `time`,
  `get_hte_tree()` fits grf's survival and censoring curves on 100 evenly
  spaced points up to `time` plus the first follow-up past it, unless
  `grf_args` sets `failure.times`. In 100 simulated data sets the effects
  moved by at most 0.002 with the same standard errors, at a third of the
  time. With the 500 trees below, a survival tree at 20,000 patients falls
  from 37 to 10 s (`"maxt"`) and from about 30 to 5 s (ctree, rpart).

* `get_hte_tree()` grows 500 trees in each of its two forests unless
  `grf_args` sets `num.trees` (above 10,000 rows `get_hte()` still grows
  200). The forests take most of the time; at 20,000 patients a continuous
  `"maxt"` tree falls from 11 to 4 s with the same rules and leaf effects
  within 0.005. Pass `grf_args = list(num.trees = 2000)` for grf's default.

* The mob methods of `get_hte_tree()` accept `tree_args$max_cuts`, the most
  cut-points tried per split variable, as `"maxt"` does. MOB refits the node
  model at every distinct value, so its split search grew with the square
  of the rows (a Weibull AFT tree took 8 minutes at 20,000 patients). The
  default `100` makes it about linear; MOB cut-points now fall on one of 100
  evenly spaced values, so earlier mob trees can move slightly. Pass
  `tree_args = list(max_cuts = 1e6)` to try every value, as before.

* `get_hte()` and `get_hte_tree()` now default to
  `factor_encoding = "integer"`, as `get_hte_icf()` already did: a factor
  enters the forest (and the tree) as one column of level codes. Pass
  `factor_encoding = "onehot"` for the previous one-column-per-level design;
  results with multi-level factors change under the new default.

* `get_hte_tree()` gains `estimator`: `"aipw"` (default, unchanged) or
  `"tmle"`, grf's targeted maximum likelihood estimate of every node and leaf
  (`grf::average_treatment_effect(method = "TMLE")`) for continuous and binary
  outcomes. The tree itself does not change; survival outcomes stay AIPW,
  the only estimator grf's causal survival forest has.

* `get_hte_tree()` and `get_hte_icf()` accept a single known propensity in
  `grf_args` (`W.hat = 0.5`, as in a trial), which every forest then uses.
  A per-row `W.hat` is still refused, since it cannot follow the split.

* `get_hte_tree()` adds the R-learner methods `"mob_r"`, `"ctree_r"` and
  `"rpart_r"`. They split on the forest's out-of-bag residuals
  (`Y - Y.hat` and `W - W.hat`), minimising the R-loss as grf's causal
  trees and the iCF do: MOB and ctree test the score of the node's
  residual-on-residual slope, and CART grows on
  `(Y - Y.hat) / (W - W.hat)` weighted by `(W - W.hat)^2` with
  cross-validation pruning. Survival outcomes use the censoring-adjusted
  numerator and denominator of the causal survival forest. With no inverse
  propensity in the residuals, extreme propensities inflate them far less
  than the AIPW scores; leaf effects stay the doubly robust ATE differences.

* `get_hte_tree()` renames `method = "rpart"` to `"rpart_dr"` to explicitly
  distinguish CART on DR scores from `"rpart_cate"`. Update existing calls
  to the new name; the fitting and pruning behavior is unchanged.

* `get_hte_tree()` adds `"mob_aft"` and `"ctree_aft"`: Weibull AFT node
  models for positive, right-censored survival times. They test treatment
  effects on the log-time-ratio scale while retaining DR survival or RMST
  differences for leaf estimates. Node fits handle aliased adjustments and
  include the estimated log-scale in the nuisance scores.

* `get_hte_tree(method = "rpart_cate")` now fits and cross-validates a CART
  approximation of the forest's out-of-bag CATE predictions. It retains
  independent leaf effect estimation and uses the existing CART controls.

* New `get_hte_tree()` grows one subgroup tree on the AIPW scores of
  `get_hte()` (continuous, binary or survival outcomes, as `surv` selects)
  and estimates every node on the other half of the patients. `method`
  chooses the partitioning: `"maxt"` (default), a heteroskedasticity-robust
  test at every node -- the Welch t of the two sides' mean scores maximised
  over the cuts, a Rademacher multiplier bootstrap and a Westfall-Young
  min-P over the variables; `"mob_dr"`, `"ctree_dr"` (partykit on the
  scores); `"mob_cate"`, `"ctree_cate"` (the same on the forest's CATE
  predictions, the two-stage approach, for comparison); `"mob_abs"`,
  `"mob_rel"`, `"ctree_abs"`, `"ctree_rel"` (MOB and model4you-style ctree
  on the scores of a node model of the outcome -- lm, or Kaplan-Meier
  pseudo-values for survival, on the difference scale; logistic or Cox
  regression on the ratio scale -- adjusted for `adj_var` and testing the
  treatment coefficient only, both switchable in `tree_args`); `"rpart_dr"`
  (pruned by cross-validation) or `"policy"` (policytree, with the
  recommended arm per leaf). `max_depth`, `alpha` and `min_leaf` stop the
  recursion. `$tree` is a partykit `party` whose node `info` holds the split
  test and the honest effect of each node, for ggparty. rpart and
  policytree join Suggests.

* New `plt_hte_tree()` draws that tree as `plt_hte_icf()` does, with the
  p-value of each split test at the inner nodes and a policy tree's
  recommended arm at the leaves; the two share one drawing.

* `get_hte_icf()` returns the selected tree as `$tree`, a partykit `party`
  on the estimation patients (design columns, outcome, `.arm`, `.dr_score`,
  `.rule`; each leaf's `info` holds its row of `$rules`), which
  `ggparty::ggparty()` and `plot()` also draw.

* New `plt_hte_icf()` draws that tree with ggparty: the split variable at
  every inner node, the conditions of the rules on the edges (`= 1`,
  `<= 0.4647`, `!= c`) and, above each leaf, its patients and its effect
  with the 95% interval from `$rules`. `type` sets the panel beneath each
  leaf: `"effect"` (default, the effect and interval on one axis, red above
  zero and blue below), `"dr"` (box plot of the patients' AIPW scores),
  `"box"` (the outcome by arm, continuous outcomes), `"bar"` (mean or event
  share by arm with 95% intervals) or `"km"` (Kaplan-Meier curves by arm
  up to `time`, survival outcomes). partykit and ggparty join Suggests.

* New `get_hte_icf(candidate_var = )` names the only covariates the rules
  may split on: the voting forests are grown on them in place of the
  importance screening, while the outcome and propensity estimates and the
  calibration test use the union of `adj_var` and `candidate_var`. The
  default `NULL` keeps the screening.

* `get_hte_icf(depth = )` with a single value now fixes the depth: that
  depth's voted partition is reported without cross-validation or
  calibration gate, and, unless `split_frac` is given, the rules are found
  and estimated on every patient. Previously a single depth was
  cross-validated against depth 0 on half the patients. At 1,600 patients
  `depth = 2` takes 6 s instead of 23 s for `depth = 1:3`.

* `get_hte_icf()` defaults to `factor_encoding = "integer"`, as the iCF code
  codes its categorical covariates: the levels of a factor, character or
  logical covariate are numbered in level order and a split keeps
  neighbouring levels together. `"onehot"`, `get_hte()`'s default and
  previously this one's, lets any set of levels form a subgroup; pass it for
  a nominal covariate whose levels have no order. Binary covariates give the
  same rules either way.

* `get_hte_icf()` gains `style = c("causalR", "icf")`. `"icf"` sets every
  step where it departs from the iCF code back to that code: the R-loss
  judged on each tree's leaf samples, no pruning margin, the unpruned best
  tree, forests grown per depth with `min.node.size = n / 25, 45, 65, 85`
  for D2-D5, a vote on tree shape with averaged split values, the
  cross-validated MSE of `lm(Y* ~ W + G + W:G + X)` on the IPW-transformed
  outcome with depth 0 left to the gate, each depth's rules from the fold
  majority, and IPTW effects (a lasso propensity score within each subgroup)
  on the same patients (`split_frac = 1`, now allowed). Each step is also a
  `rule_args` field (`loss`, `eval`, `penalty`, `prune`, `vote`, `grow`,
  `cv_loss`, `cv_zero`, `cv_rules`, `estimate`, and `screen = "icf"`), so
  one can be switched on its own. The folds still grow their forests on
  their training folds only, where the iCF code grows them on every patient.
  Needs a binary or continuous outcome; glmnet joins Suggests. The default
  results are unchanged.

* `get_hte_icf()` counts a partition split on X1 first and the same
  partition split on X3 first as one vote, as its help page said; they were
  counted apart.

* New `plt_hte_unihtee()` draws a `get_hte_unihtee()` screen for one
  measure: `type = "bar"` (default) the signed estimates with their
  intervals, `"volcano"` estimate against -log10(p) with the BH-significant
  candidates labelled, and `"dep"` one panel per candidate with the
  projection behind its row -- the line whose slope times `sd` is the
  estimate, or the two level means of a binary candidate -- beside a natural
  spline of the same pseudo-outcome that shows what the slope leaves out
  (a U-shaped modifier has a flat line and a curved spline). Works for both
  `method`s, every measure (ratio and OR on a log axis) and survival
  outcomes, with no forest refitted. `get_hte_unihtee()` and the plot share
  one projection, so every drawn line matches its table row.
  `dr_args$bins` (default `0`, off) adds quantile-bin means of the
  pseudo-outcome as a second shape check, `axis_arg = list(share_y =
  "none")` gives each panel its own y range, as in `plt_hte_dep()`, and
  `var_names` relabels the candidates as in `get_bal()` (default
  `RegR::name_map_seer`).

* New `get_hte_unihtee()` screens candidate effect modifiers with the
  treatment effect modifier variable importance parameter (TEM-VIP) of
  Boileau et al. (2025), the estimand of the `unihtee` package, written
  without depending on it: the least-squares slope of cross-fitted AIPW
  pseudo-outcomes on each candidate alone, per SD for a continuous
  candidate and as a level difference for a binary one, with HC3 intervals
  and Benjamini-Hochberg `p.adj`. The scores come from the `get_hte()`
  causal forest (`method = "grf"`, also for survival outcomes; a row equals
  `grf::best_linear_projection()`) or from a logistic propensity model and
  per-arm outcome GLMs with 5-fold cross-fitting (`method = "glm"`, no
  forest). `measure = c("diff", "ratio", "OR")` projects the difference,
  log ratio or log odds ratio. In simulations the glm route covered at
  0.94-0.97 on all three scales and rejected 4.5% of null candidates at
  p < 0.05. Returns `list(vip, data)`.

* New `get_hte_icf()` finds subgroup rules such as `X1 = 1 & X3 = 0`
  without naming them in advance, after the iterative causal forest (iCF)
  of Wang et al. (2024), written from the paper since the iCF code carries
  no licence. Half the patients (`split_frac`) discover the rules: forests
  grown on the covariates `get_hte()` ranks at or above mean importance,
  each tree cut at every candidate `depth` (default `1:3`), pruned and
  judged on the patients it did not choose its splits on by the squared
  error of the AIPW scores, each forest voting with its best tree, and the
  depth chosen by cross-validation with depth 0 ("no subgroups") among the
  candidates. As in the paper, subgroups are only reported when the
  forest's calibration test gives p <= 0.1 (`rule_args$gate`). The other
  half estimates every rule's doubly robust effect through `get_hte()`, so
  `plt_hte_sub(res$est, sub_var = ".rule")` draws them. Continuous, binary
  and survival outcomes work. In simulations the X1 x X3 partition was
  found in 12 of 12 data sets with no spurious split and no subgroups were
  reported in 12 of 12 without heterogeneity; the default 20 forests x 200
  trees x 5 folds took 49 s at 5,000 patients. Returns
  `list(rules, cv, vote, importance, est)`.

* `get_bal(methods = )` accepts a list that mixes named scheme
  specifications with unnamed shorthands, e.g.
  ``list(`PSM 1:2` = list(design = "matching", ratio = 2), "ATE", "ATO")``;
  each shorthand keeps its usual legend label.

* `plt_hte_sub()` heads the effect column with the effect alone, e.g.
  `Risk difference` instead of `Risk difference (95% CI)`, so the plot no
  longer states `conf_level`.

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
