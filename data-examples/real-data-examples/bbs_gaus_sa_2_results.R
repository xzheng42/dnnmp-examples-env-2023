################################################################################
### BBS DATA ANALYSIS - Sensitivity analysis of L; predictive performances
### Section D.2.1 of the supplementary material
################################################################################

rm(list = ls())

inpath <- "/path/to/mcmc output/"
outpath <- "/path/to/output/"

library(dnnmp)
library(readr)

load(paste0(inpath, "sa2_gaus_nbnnmp_out_L_5to20.RData"))


#-------------------------------------------------------------------------------
# Prediction at the hold-out set locations
#-------------------------------------------------------------------------------
ne_sizes <- sapply(nnmp_out_list, function(x) x$mod_spec$neighbor_size)
pred_list <- vector("list", length = length(ne_sizes))

for (i in seq_along(ne_sizes)) {
  
  ne_size <- ne_sizes[i]
  print(paste0("Start with L = ", ne_size))
  
  nnmp_out <- nnmp_out_list[[i]]
  set.seed(42)
  gaus_nbnnmp_pred <- predict(nnmp_out,
                             nonref_coords = valid_dat[,1:2], 
                             nonref_covars = valid_XX,
                             probs = c(0.025, 0.5, 0.975),
                             predict_sam = TRUE)
  
  pred_list[[i]] <- gaus_nbnnmp_pred
  print(paste0("Done with L = ", ne_size))
  
}

#-------------------------------------------------------------------------------
# Table 2 of the supplementary material 
#-------------------------------------------------------------------------------
tbl <- array(NA, dim = c(length(ne_sizes), 6))
rownames(tbl) <- paste0("L = ", ne_sizes)
colnames(tbl) <- c("RMPSE", "95 \\% CI", "95 \\%CI width", "CRPS", "ES", "VS")

for (i in seq_along(ne_sizes)) {
  
  gaus_nbnnmp_pred <- pred_list[[i]]
  
  gaus_nbnnmp_qq <- gaus_nbnnmp_pred[['obs_qq']]
  
  # RMSPE
  gaus_nbnnmp_rmspe <- sqrt(mean((gaus_nbnnmp_qq[2,] - valid_dat[,3])^2))
  
  # 95% CI width
  gaus_nbnnmp_width <- mean(gaus_nbnnmp_qq[3,] - gaus_nbnnmp_qq[1,])
  
  # 95% CI Coverage rate
  valid_num <- nrow(valid_dat)
  gaus_cover_boolen <- sapply(1:(valid_num), function(x) {
    valid_dat[x,3] >= gaus_nbnnmp_qq[1,x] & valid_dat[x,3] <= gaus_nbnnmp_qq[3,x]
  })
  gaus_nbnnmp_cover <- sum(gaus_cover_boolen) / length(gaus_cover_boolen)
  
  # CRPS
  gaus_nbnnmp_crps <- mean(scoringRules::crps_sample(valid_dat[, 3], gaus_nbnnmp_pred$sam))
  
  # ES
  gaus_nbnnmp_es <- scoringRules::es_sample(valid_dat[, 3], gaus_nbnnmp_pred$sam)
  
  # VS, p = 1
  gaus_nbnnmp_vs <- scoringRules::vs_sample(valid_dat[, 3], gaus_nbnnmp_pred$sam, p = 1)
  
  tbl[i,] <- c(gaus_nbnnmp_rmspe, gaus_nbnnmp_cover, gaus_nbnnmp_width, gaus_nbnnmp_crps, gaus_nbnnmp_es, gaus_nbnnmp_vs)
  
  print(paste0("Done with L = ", ne_size))
  
}

tbl_latex <- knitr::kable(tbl, digits = 2, align = "c", format = "latex", escape = FALSE)

readr::write_file(tbl_latex, paste0(outpath, "nbnnmp_gaus_sa_tbl2.txt"))

