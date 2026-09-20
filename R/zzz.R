# Columns of the local data frames built inside sens-plt.R renderers and
# referenced through ggplot2::aes().
utils::globalVariables(c(
  "gamma",     # .sens_plt_tip:     confounder-outcome effect grid
  "adjusted",  # .sens_plt_tip:     adjusted effect on that grid
  "which",     # .sens_plt_tip / .sens_plt_evalue: curve identity
  "rr_eu",     # .sens_plt_evalue:  confounder-exposure risk ratio
  "rr_ud",     # .sens_plt_evalue:  confounder-outcome risk ratio
  "e",         # .sens_plt_evalue:  E-value marker on the diagonal
  "id",        # .sens_contour_lm:  contour path identifier
  "label",     # .sens_contour_lm:  contour and bound labels
  "x",         # .sens_contour_lm:  treatment-side partial R2
  "y"          # .sens_contour_lm:  outcome-side partial R2
))
