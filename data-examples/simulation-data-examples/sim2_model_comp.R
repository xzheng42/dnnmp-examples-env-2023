################################################################################
### SIMULATION DATA EXAMPLE 2 
### Model Comparison (NBNNMP, SGLMM-GP, and SGLMM-GPP)
################################################################################

rm(list = ls())

inpath <- "/path/to/mcmc output/"
outpath <- "/path/to/output/"

library(scoringRules)
library(spBayes)
library(readr)
library(dnnmp)

load(paste0(inpath, "nbnnmp_res.RData"))
load(paste0(inpath, "sglmgp_res.RData"))
load(paste0(inpath, "sglmgpp_res.RData"))


#-------------------------------------------------------------------------------
# Parameter estimates
#-------------------------------------------------------------------------------
n.samples <- nrow(sim.s1$p.beta.theta.samples)
idx <- seq(20001, n.samples, by = 5)
nbnnmp_beta_est <- dnnmp::summarize(nnmp_out$post_samples$regcoef)
sglmgp_beta_est <- dnnmp::summarize(t(as.matrix(sim.s1$p.beta.theta.samples)[idx,1:3]))
sglmgpp_beta_est <- dnnmp::summarize(t(as.matrix(sim.s2$p.beta.theta.samples)[idx,1:3]))
tbl_beta_est <- cbind(nbnnmp_beta_est, sglmgp_beta_est, sglmgpp_beta_est)
colnames(tbl_beta_est) <- c("NBNNMP", "SGLMM-GP", "SGLMM-GPP")
rownames(tbl_beta_est) <- c("$\beta_0$", "$\beta_1$", "$\beta_2$")

#-------------------------------------------------------------------------------
# Out-of-sample prediction
#-------------------------------------------------------------------------------
test_num <- nrow(test_dat)


# NBNNMP  ----------------------------------------------------------------------
set.seed(42)
nbnnmp_pred <- predict(nnmp_out, nonref_covars = cbind(1, test_dat[, 1:2]),
                       nonref_coords = test_dat[,1:2],
                            probs = c(0.025, 0.5, 0.975), predict_sam = TRUE)

nbnnmp_qq <- nbnnmp_pred[['obs_qq']]

### RMSPE
nbnnmp_rmspe <- sqrt(mean((nbnnmp_qq[2,] - test_dat[,3])^2))

### 95% CI width
nbnnmp_width <- mean(nbnnmp_qq[3,] - nbnnmp_qq[1,])

### 95% CI Coverage rate
cover_boolen <- sapply(1:(test_num), function(x) {
  test_dat[x,3] >= nbnnmp_qq[1,x] & test_dat[x,3] <= nbnnmp_qq[3,x]
})
nbnnmp_cover <- sum(cover_boolen) / length(cover_boolen)

### CRPS
nbnnmp_crps <- mean(scoringRules::crps_sample(test_dat[, 3], nbnnmp_pred$sam))

### ES
nbnnmp_es <- scoringRules::es_sample(test_dat[, 3], nbnnmp_pred$sam)

### VS, p = 1
nbnnmp_vs <- scoringRules::vs_sample(test_dat[, 3], nbnnmp_pred$sam, p = 1)

### Table the metrics
nbnnmp_metric <- as.matrix(c(nbnnmp_rmspe, nbnnmp_cover, nbnnmp_width, 
                             nbnnmp_crps, nbnnmp_es, nbnnmp_vs, nnmp_out$runtime[3]/60))

# SGLMM-GP  --------------------------------------------------------------------
set.seed(42)
sglm_gp_pred <- spPredict(sim.s1, 
                          pred.coords = test_dat[,1:2], 
                          pred.covars = as.matrix(cbind(1, test_dat[,1:2])),
                          start = 20001,
                          end = n.samples,
                          thin = 5)

sglm_gp_sam <- apply(sglm_gp_pred$p.y.predictive.samples, 2, function(x) rpois(test_num, x))

sglm_gp_q50 <- apply(sglm_gp_sam, 1, median)
sglm_gp_q025 <- apply(sglm_gp_sam, 1, function(x) quantile(x, probs = 0.025))
sglm_gp_q975 <- apply(sglm_gp_sam, 1, function(x) quantile(x, probs = 0.975))

sglm_gp_width <- mean(sglm_gp_q975 - sglm_gp_q025)

sglm_gp_rmspe <- sqrt(mean((sglm_gp_q50 - test_dat[,3])^2))

cover_boolen_1 <- sapply(1:test_num, function(x) {
  test_dat[x,3] >= sglm_gp_q025[x] & test_dat[x,3] <= sglm_gp_q975[x]
})
sglm_gp_cover <- sum(cover_boolen_1) / length(cover_boolen_1)

sglm_gp_crps <- mean(scoringRules::crps_sample(test_dat[, 3], sglm_gp_sam))

sglm_gp_es <- scoringRules::es_sample(test_dat[, 3], sglm_gp_sam)

sglm_gp_vs <- scoringRules::vs_sample(test_dat[, 3], sglm_gp_sam, p = 1)

# Table the metrics
sglm_gp_metric <- as.matrix(c(sglm_gp_rmspe, sglm_gp_cover, sglm_gp_width, sglm_gp_crps,
                              sglm_gp_es, sglm_gp_vs, (sim.s1$run.time[1]+sim.s1$run.time[2])/60))


# SGLMM-GPP  -------------------------------------------------------------------
set.seed(42)
sglm_gpp_pred <- spPredict(sim.s2, 
                             pred.coords = test_dat[,1:2], 
                             pred.covars = as.matrix(cbind(1, test_dat[,1:2])),
                             start = 20001,
                             end = n.samples,
                             thin = 5)

sglm_gpp_sam <- apply(sglm_gpp_pred$p.y.predictive.samples, 2, function(x) rpois(test_num, x))

sglm_gpp_q50 <- apply(sglm_gpp_sam, 1, median)
sglm_gpp_q025 <- apply(sglm_gpp_sam, 1, function(x) quantile(x, probs = 0.025))
sglm_gpp_q975 <- apply(sglm_gpp_sam, 1, function(x) quantile(x, probs = 0.975))

sglm_gpp_width <- mean(sglm_gpp_q975 - sglm_gpp_q025)

sglm_gpp_rmspe <- sqrt(mean((sglm_gpp_q50 - test_dat[,3])^2))

cover_boolen_1 <- sapply(1:test_num, function(x) {
  test_dat[x,3] >= sglm_gpp_q025[x] & test_dat[x,3] <= sglm_gpp_q975[x]
})
sglm_gpp_cover <- sum(cover_boolen_1) / length(cover_boolen_1)

sglm_gpp_crps <- mean(scoringRules::crps_sample(test_dat[, 3], sglm_gpp_sam))

sglm_gpp_es <- scoringRules::es_sample(test_dat[, 3], sglm_gpp_sam)

sglm_gpp_vs <- scoringRules::vs_sample(test_dat[, 3], sglm_gpp_sam, p = 1)

sglm_gpp_metric <- as.matrix(c(sglm_gpp_rmspe, sglm_gpp_cover, sglm_gpp_width, sglm_gpp_crps, 
                               sglm_gpp_es, sglm_gpp_vs, sim.s2$run.time[3]/60))

#-------------------------------------------------------------------------------
# Table 3 of the main paper
#-------------------------------------------------------------------------------
tbl_pred <- round(cbind(nbnnmp_metric, sglm_gp_metric, sglm_gpp_metric), 2)
rownames(tbl_pred) <- c("RMSPE", "95% CI cover", "95% CI width", "CRPS", "ES", "VS", "Time (mins)")
tbl <- rbind(tbl_beta_est, tbl_pred)
tbl_latex <- knitr::kable(tbl, align = "c", digits = 2, format = "latex", escape = TRUE)
write_file(tbl_latex, paste0(outpath, "nbnnmp_comp_est_pred.txt"))
