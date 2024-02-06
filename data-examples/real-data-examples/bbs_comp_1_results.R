################################################################################
### BBS DATA ANALYSIS - Compare three copula NBNNMP models
################################################################################

rm(list = ls())

inpath <- "/path/to/mcmc output/"
outpath <- "/path/to/output/"

library(readr)

load(paste0(inpath, "comp_res_20.RData"))


#-------------------------------------------------------------------------------
# Prediction at the hold-out set locations
#-------------------------------------------------------------------------------
# Prediction
set.seed(42)
gaus_nbnnmp_pred <-  predict(nnmp_gaus_out, 
                             nonref_coords = valid_dat[,1:2], 
                             nonref_covars = valid_XX,
                             probs = c(0.025, 0.5, 0.975),
                             predict_sam = TRUE)
  
set.seed(42)
gum_nbnnmp_pred <- predict(nnmp_gum_out, 
                           nonref_coords = valid_dat[,1:2], 
                           nonref_covars = valid_XX,
                           probs = c(0.025, 0.5, 0.975),
                           predict_sam = TRUE)

set.seed(42)
clay_nbnnmp_pred <-predict(nnmp_clay_out, 
                           nonref_coords = valid_dat[,1:2], 
                           nonref_covars = valid_XX,
                           probs = c(0.025, 0.5, 0.975),
                           predict_sam = TRUE)

gaus_nbnnmp_qq <- gaus_nbnnmp_pred[['obs_qq']]
gum_nbnnmp_qq <- gum_nbnnmp_pred[['obs_qq']]
clay_nbnnmp_qq <- clay_nbnnmp_pred[['obs_qq']]


#-------------------------------------------------------------------------------
# Table 3 of the supplementary material
#-------------------------------------------------------------------------------

# RMSPE
gaus_nbnnmp_rmspe <- sqrt(mean((gaus_nbnnmp_qq[2,] - valid_dat[,3])^2))
gum_nbnnmp_rmspe <- sqrt(mean((gum_nbnnmp_qq[2,] - valid_dat[,3])^2))
clay_nbnnmp_rmspe <- sqrt(mean((clay_nbnnmp_qq[2,] - valid_dat[,3])^2))

# 95% CI width
gaus_nbnnmp_width <- mean(gaus_nbnnmp_qq[3,] - gaus_nbnnmp_qq[1,])
gum_nbnnmp_width <- mean(gum_nbnnmp_qq[3,] - gum_nbnnmp_qq[1,])
clay_nbnnmp_width <- mean(clay_nbnnmp_qq[3,] - clay_nbnnmp_qq[1,])

# 95% CI Coverage rate
valid_num <- nrow(valid_dat)
gaus_cover_boolen <- sapply(1:(valid_num), function(x) {
  valid_dat[x,3] >= gaus_nbnnmp_qq[1,x] & valid_dat[x,3] <= gaus_nbnnmp_qq[3,x]
})
gaus_nbnnmp_cover <- mean(gaus_cover_boolen)

gum_cover_boolen <- sapply(1:(valid_num), function(x) {
  valid_dat[x,3] >= gum_nbnnmp_qq[1,x] & valid_dat[x,3] <= gum_nbnnmp_qq[3,x]
})
gum_nbnnmp_cover <- mean(gum_cover_boolen)

clay_cover_boolen <- sapply(1:(valid_num), function(x) {
  valid_dat[x,3] >= clay_nbnnmp_qq[1,x] & valid_dat[x,3] <= clay_nbnnmp_qq[3,x]
})
clay_nbnnmp_cover <- mean(clay_cover_boolen)

# CRPS
gaus_nbnnmp_crps <- mean(scoringRules::crps_sample(valid_dat[, 3], gaus_nbnnmp_pred$sam))
gum_nbnnmp_crps <- mean(scoringRules::crps_sample(valid_dat[, 3], gum_nbnnmp_pred$sam))
clay_nbnnmp_crps <- mean(scoringRules::crps_sample(valid_dat[, 3], clay_nbnnmp_pred$sam))

# ES
gaus_nbnnmp_es <- scoringRules::es_sample(valid_dat[, 3], gaus_nbnnmp_pred$sam)
gum_nbnnmp_es <- scoringRules::es_sample(valid_dat[, 3], gum_nbnnmp_pred$sam)
clay_nbnnmp_es <- scoringRules::es_sample(valid_dat[, 3], clay_nbnnmp_pred$sam)

# VS, p = 1
gaus_nbnnmp_vs <- scoringRules::vs_sample(valid_dat[, 3], gaus_nbnnmp_pred$sam, p = 1)
gum_nbnnmp_vs <- scoringRules::vs_sample(valid_dat[, 3], gum_nbnnmp_pred$sam, p = 1)
clay_nbnnmp_vs <- scoringRules::vs_sample(valid_dat[, 3], clay_nbnnmp_pred$sam, p = 1)


tbl <- array(NA, dim = c(3, 6))
tbl[1,] <- c(gaus_nbnnmp_rmspe, gaus_nbnnmp_cover, gaus_nbnnmp_width, gaus_nbnnmp_crps, gaus_nbnnmp_es, gaus_nbnnmp_vs)
tbl[2,] <- c(gum_nbnnmp_rmspe, gum_nbnnmp_cover, gum_nbnnmp_width, gum_nbnnmp_crps, gum_nbnnmp_es, gum_nbnnmp_vs)
tbl[3,] <- c(clay_nbnnmp_rmspe, clay_nbnnmp_cover, clay_nbnnmp_width, clay_nbnnmp_crps, clay_nbnnmp_es, clay_nbnnmp_vs)

rownames(tbl) <- c("Gaussian", "Gumbel", "Clayton")
colnames(tbl) <- c("RMSPE", "95\\% CI cover", "95\\% CI width", "CRPS", "ES", "VS")

tbl_latex <- knitr::kable(tbl, digits = 2, align = "c", format = "latex", escape = FALSE)
readr::write_file(tbl_latex, paste0(outpath, "nbnnmp_cop_comp_tbl3.txt"))
