# causalR (development version)

* `get_hte(estimand = "ATO")` accepts boundary propensities without losing
  valid overlap estimates; unavailable ordinary AIPW diagnostics are marked NA.
* `get_hte()` retains average effects when spline or calibration diagnostics
  are degenerate, returning unavailable diagnostics with a warning.
* `get_hte()` and `plt_hte_sub()` account for shared clusters when testing
  differences between subgroup effects, retaining grf's marginal variances.
* `get_hte()` and `plt_hte_dep()` now use the forest's observation weights
  and clusters for covariate heterogeneity tests and doubly robust curves.
