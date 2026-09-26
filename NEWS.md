# causalR (development version)

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
