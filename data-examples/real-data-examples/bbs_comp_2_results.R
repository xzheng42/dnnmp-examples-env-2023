################################################################################
### SIMULATION DATA EXAMPLE 2 
### Model Comparison
################################################################################

rm(list = ls())

inpath <- "/path/to/mcmc output/"
outpath <- "/path/to/output/"

load(paste0(inpath, "nbnnmp_out_20.RData"))
load(paste0(inpath, "sglmgp_res_500.RData"))

library(scoringRules)
library(spBayes)
library(dnnmp)


#-------------------------------------------------------------------------------
# Parameter estimates
#-------------------------------------------------------------------------------
# NBNNMP
nbnnmp_beta_est <- summarize(nnmp_out$post_samples$regcoef)

# SGLMM-GP 
n.samples <- nrow(sim.s1$p.beta.theta.samples)
idx <- seq(30001, n.samples, by = 5)
sglmgp_beta_est <- summarize(as.matrix(sim.s1$p.beta.theta.samples)[idx,1:2])

# Table 
tbl_beta_est <- cbind(nbnnmp_beta_est, sglmgp_beta_est)
colnames(tbl_beta_est) <- c("NBNNMP", "SGLMM-GP")
rownames(tbl_beta_est) <- c("$\\beta_0$", "$\\beta_1$")

#-------------------------------------------------------------------------------
# Out-of-sample prediction
#-------------------------------------------------------------------------------
valid_num <- nrow(valid_dat)

# NBNNMP ---------------------------------------------------------------------
set.seed(42)
nbnnmp_pred <- predict(nnmp_out,
                       nonref_coords = valid_dat[,1:2], 
                       nonref_covars = valid_XX,
                       probs = c(0.025, 0.5, 0.975),
                       predict_sam = TRUE)

nbnnmp_qq <- nbnnmp_pred[['obs_qq']]

# RMSPE
nbnnmp_rmspe <- sqrt(mean((nbnnmp_qq[2,] - valid_dat[,3])^2))

# 95% CI width
nbnnmp_width <- mean(nbnnmp_qq[3,] - nbnnmp_qq[1,])

# 95% CI Coverage rate
valid_num <- nrow(valid_dat)
cover_boolen <- sapply(1:(valid_num), function(x) {
  valid_dat[x,3] >= nbnnmp_qq[1,x] & valid_dat[x,3] <= nbnnmp_qq[3,x]
})
nbnnmp_cover <- mean(cover_boolen)

# CRPS
nbnnmp_crps <- mean(scoringRules::crps_sample(valid_dat[, 3], nbnnmp_pred$sam))

# ES
nbnnmp_es <- scoringRules::es_sample(valid_dat[, 3], nbnnmp_pred$sam)

# VS, p = 1
nbnnmp_vs <- scoringRules::vs_sample(valid_dat[, 3], nbnnmp_pred$sam, p = 1)

# Table the metrics
runtime <- nnmp_out$runtime
nbnnmp_metric <- as.matrix(c(nbnnmp_rmspe, nbnnmp_cover, nbnnmp_width, nbnnmp_crps, nbnnmp_es, nbnnmp_vs, runtime[3]/60))

# SGLMM-GP  ------------------------------------------------------------------
load(paste0(inpath, "sglmgp_res_500.RData"))
set.seed(42)
sglm_gp_pred <- spPredict(sim.s1, 
                          pred.coords = valid_dat[,1:2], 
                          pred.covars = as.matrix(cbind(1, valid_dat[,2])),
                          start = 30001,
                          end = n.samples,
                          thin = 5)

save(sglm_gp_pred, file = paste0(outpath, "slgm_gp_500_pred.RData"))

sglm_gp_sam <- apply(sglm_gp_pred$p.y.predictive.samples, 2, function(x) rpois(valid_num, x))

sglm_gp_q50 <- apply(sglm_gp_sam, 1, median)
sglm_gp_q025 <- apply(sglm_gp_sam, 1, function(x) quantile(x, probs = 0.025))
sglm_gp_q975 <- apply(sglm_gp_sam, 1, function(x) quantile(x, probs = 0.975))

sglm_gp_width <- mean(sglm_gp_q975 - sglm_gp_q025)

sglm_gp_rmspe <- sqrt(mean((sglm_gp_q50 - valid_dat[,3])^2))

cover_boolen_1 <- sapply(1:valid_num, function(x) {
  valid_dat[x,3] >= sglm_gp_q025[x] & valid_dat[x,3] <= sglm_gp_q975[x]
})
sglm_gp_cover <- mean(cover_boolen_1)

sglm_gp_crps <- mean(scoringRules::crps_sample(valid_dat[, 3], sglm_gp_sam))

sglm_gp_es <- scoringRules::es_sample(valid_dat[, 3], sglm_gp_sam)

sglm_gp_vs <- scoringRules::vs_sample(valid_dat[, 3], sglm_gp_sam, p = 1)


sglm_gp_metric <- as.matrix(c(sglm_gp_rmspe, sglm_gp_cover, sglm_gp_width, sglm_gp_crps,
                              sglm_gp_es, sglm_gp_vs, (sim.s1$run.time[1]+sim.s1$run.time[2])/60))


#-------------------------------------------------------------------------------
# Table 4 of the supplementary material
#-------------------------------------------------------------------------------
tbl_pred <- round(cbind(nbnnmp_metric, sglm_gp_metric), 2)
rownames(tbl_pred) <- c("RMSPE", "95% CI", "95% CI width", "CRPS", "ES", "VS", "Time (mins)")
tbl <- rbind(tbl_beta_est, tbl_pred)
tbl_latex <- knitr::kable(tbl, align = "c", digits = 3, format = "latex", escape = FALSE)
write_file(tbl_latex, paste0(outpath, "nbnnmp_sglmm_comp_tbl4.txt"))



