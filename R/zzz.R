# Columns of the local data frames built inside the sens-plt.R, psw-plt.R,
# psm-plt.R and hte-plt.R renderers and referenced through ggplot2::aes().
utils::globalVariables(c(
  # sens-plt.R
  "gamma",     # .sens_plt_tip:     confounder-outcome effect grid
  "adjusted",  # .sens_plt_tip:     adjusted effect on that grid
  "which",     # .sens_plt_tip / .sens_plt_evalue: curve identity
  "rr_eu",     # .sens_plt_evalue:  confounder-exposure risk ratio
  "rr_ud",     # .sens_plt_evalue:  confounder-outcome risk ratio
  "e",         # .sens_plt_evalue:  E-value marker on the diagonal
  "id",        # .sens_contour_lm:  contour path identifier
  "label",     # .sens_contour_lm:  contour and bound labels
  "x",         # .sens_contour_lm:  treatment-side partial R2
  "y",         # .sens_contour_lm:  outcome-side partial R2
  # psw-plt.R
  "smd",       # .psw_plt_love:     absolute standardised mean difference
  "variable",  # .psw_plt_love:     covariate name from halfmoon
  "Weighting", # .psw_plt_love:     weighting scheme, shown in the legend
  "ess",       # .psw_plt_ess:      effective sample size
  "ess_pct",   # .psw_plt_ess:      the same as a share of n
  "estimand",  # .psw_plt_ess:      bar identity
  "weight",    # .psw_plt_weight:   the weight itself
  "arm",       # .psw_plt_weight:   treated / control
  "ps",        # .psw_plt_ps:       propensity score
  # psm-plt.R (reuses smd / variable / ess / weight / arm / ps above)
  "Matching",  # .psm_plt_love:     matching scheme, shown in the legend
  "method",    # .psm_plt_ess:      bar identity
  "n",         # .psm_plt_ess:      matched count in the bar label
  # hte-plt.R (reuses x / y above)
  "estimate",  # plt_hte_dep:       AIPW mean, spline or partial dependence
  "conf.low",  # plt_hte_dep:       lower interval bound
  "conf.high", # plt_hte_dep:       upper interval bound
  "panel",     # plt_hte_dep:       strip label of a panel
  "x1",        # plt_hte_dep:       first covariate of the heat map
  "x2"         # plt_hte_dep:       second covariate of the heat map
))
