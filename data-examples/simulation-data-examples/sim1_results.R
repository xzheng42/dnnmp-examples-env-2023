################################################################################
### SIMULATION DATA EXAMPLE 1
################################################################################

rm(list = ls())

inpath <- "/path/to/mcmc output/"
outpath <- "/path/to/output/"

library(ggplot2)
library(readr)
library(dnnmp)

load(paste0(inpath, "sim1_res.RData"))

myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))


#-------------------------------------------------------------------------------
# Plot true fields; first row in Figure 1 of the supplementary material 
#-------------------------------------------------------------------------------
plot_save_list <- vector("list", length = length(skewparams))

for (j in seq_along(skewparams)) {
  
  dat <- dat_save[[j]]
  
  plot_save <- ggplot(as.data.frame(dat), aes(x = x, y = y, fill = z)) +
    geom_raster(interpolate = T) + theme_bw() + 
    scale_fill_gradientn(colours = myPalette(100), limits = range(dat[,3])) + 
    labs(fill = "", x = "Easting", y = "Northing") +  
    theme(plot.margin = unit(c(0.5,0,0,0), "cm")) +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.key.height= unit(2.2, 'cm'),
          axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
          axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent"))
  
  plot_save_list[[j]] <- plot_save
  
}

png(paste0(outpath,  "pois_skew_", skewparams[1], ".png"), units = "in", width = 7, height = 6, res = 300, bg = "transparent")
plot_save_list[[1]]
dev.off()

png(paste0(outpath,  "pois_skew_", skewparams[2], ".png"), units = "in", width = 7, height = 6, res = 300, bg = "transparent")
plot_save_list[[2]]
dev.off()

png(paste0(outpath,  "pois_skew_", skewparams[3], ".png"), units = "in", width = 7, height = 6, res = 300, bg = "transparent")
plot_save_list[[3]]
dev.off()

#-------------------------------------------------------------------------------
# Posterior means and 95% credible interval estimates; Table 2 of the main paper
#-------------------------------------------------------------------------------
tbl_th <- matrix(NA, nrow = 3, ncol = 3)

for (j in seq_along(skewparams)) {
  tbl_th[1, j] <- summarize(nnmp_out_list[[j]]$gaus$post_samples$rate)
  tbl_th[2, j] <- summarize(nnmp_out_list[[j]]$gum$post_samples$rate)
  tbl_th[3, j] <- summarize(nnmp_out_list[[j]]$clay$post_samples$rate)
}

colnames(tbl_th) <- c("$\\sigma_1 = 1$", "$\\sigma_1 = 3$", "$\\sigma_1 = 10$")
rownames(tbl_th) <- c("Gaussian", "Gumbel", "Clayton")
tbl_latex <- knitr::kable(tbl_th, align = "c", format = "latex", escape = FALSE)
write_file(tbl_latex, paste0(outpath, "ponnmp_comp_rate_est.txt"))

#-------------------------------------------------------------------------------
# Predictive performances; Table 2 of the main paper
#-------------------------------------------------------------------------------
tbl_rpmse <- tbl_ci_width <- tbl_ci_cover <- tbl_crps <- matrix(NA, nrow = 3, ncol = 3)
tbl_es <- tbl_vs_05 <- tbl_vs_1 <- tbl_vs_2 <- matrix(NA, nrow = 3, ncol = 3)

for (j in seq_along(skewparams)) {
  
  test_dat <- test_dat_save[[j]]
  
  nnmp_out_gaus <- nnmp_out_list[[j]]$gaus
  nnmp_out_gum <- nnmp_out_list[[j]]$gum
  nnmp_out_clay <- nnmp_out_list[[j]]$clay
  
  set.seed(42)
  gaus_nnmp_test <- predict(nnmp_out_gaus, nonref_coords = test_dat[,1:2],
                            probs = c(0.025, 0.5, 0.975), predict_sam = TRUE,
                            compute_control = 
                              list(lb = -1, ub = max(nnmp_out_gaus$orig_data$response) + 100,
                                   tol = .Machine$double.eps^0.25, maxit = 1000))
  
  set.seed(42)
  gum_nnmp_test <- predict(nnmp_out_gum, nonref_coords = test_dat[,1:2],
                            probs = c(0.025, 0.5, 0.975), predict_sam = TRUE,
                            compute_control = 
                              list(lb = -1, ub = max(nnmp_out_gaus$orig_data$response) + 100,
                                   tol = .Machine$double.eps^0.25, maxit = 1000))

  set.seed(42)
  clay_nnmp_test <- predict(nnmp_out_clay, nonref_coords = test_dat[,1:2],
                            probs = c(0.025, 0.5, 0.975), predict_sam = TRUE,
                            compute_control = 
                              list(lb = -1, ub = max(nnmp_out_gaus$orig_data$response) + 100,
                                   tol = .Machine$double.eps^0.25, maxit = 1000))

  gaus_nnmp_test_qq <- gaus_nnmp_test[['obs_qq']]
  gum_nnmp_test_qq <- gum_nnmp_test[['obs_qq']]
  clay_nnmp_test_qq <- clay_nnmp_test[['obs_qq']]
  
  # CRPS
  gaus_nnmp_crps <- mean(scoringRules::crps_sample(test_dat[, 3],gaus_nnmp_test$sam))
  gum_nnmp_crps <- mean(scoringRules::crps_sample(test_dat[, 3], gum_nnmp_test$sam))
  clay_nnmp_crps <- mean(scoringRules::crps_sample(test_dat[, 3], clay_nnmp_test$sam))
  
  # ES
  gaus_nnmp_es <- scoringRules::es_sample(test_dat[, 3],gaus_nnmp_test$sam)
  gum_nnmp_es <- scoringRules::es_sample(test_dat[, 3], gum_nnmp_test$sam)
  clay_nnmp_es <- scoringRules::es_sample(test_dat[, 3], clay_nnmp_test$sam)
  
  # VS, p = 1
  gaus_nnmp_vs_1 <- scoringRules::vs_sample(test_dat[, 3],gaus_nnmp_test$sam, p = 1)
  gum_nnmp_vs_1 <- scoringRules::vs_sample(test_dat[, 3], gum_nnmp_test$sam, p = 1)
  clay_nnmp_vs_1 <- scoringRules::vs_sample(test_dat[, 3], clay_nnmp_test$sam, p = 1)
  
  tbl_crps[1, j] <- gaus_nnmp_crps
  tbl_crps[2, j] <- gum_nnmp_crps
  tbl_crps[3, j] <- clay_nnmp_crps
  
  tbl_es[1, j] <- gaus_nnmp_es
  tbl_es[2, j] <- gum_nnmp_es
  tbl_es[3, j] <- clay_nnmp_es
  
  tbl_vs_1[1, j] <- gaus_nnmp_vs_1
  tbl_vs_1[2, j] <- gum_nnmp_vs_1
  tbl_vs_1[3, j] <- clay_nnmp_vs_1
  
}

tbl <- cbind(tbl_crps[,1], tbl_es[,1], round(tbl_vs_1[,1]),
             tbl_crps[,2], tbl_es[,2], round(tbl_vs_1[,2]),
             tbl_crps[,3], tbl_es[,3], round(tbl_vs_1[,3]))
tbl_latex <- knitr::kable(tbl, align = "c", digits = 2, format = "latex", escape = TRUE)
write_file(tbl_latex, paste0(outpath, "ponnmp_comp_pred.txt"))


#-------------------------------------------------------------------------------
# Prediction
#-------------------------------------------------------------------------------
nonref_grid <- as.matrix(expand.grid(seq(0.01, 1, by = 1 / 120), 
                                     seq(0.01, 1, by = 1 / 120)))
plot_gaus_nnmp_list <- vector("list", length(skewparams))
plot_gum_nnmp_list <- vector("list", length(skewparams))
plot_clay_nnmp_list <- vector("list", length(skewparams))

for (j in seq_along(skewparams)) {
  
  nnmp_out_gaus <- nnmp_out_list[[j]]$gaus
  nnmp_out_gum <- nnmp_out_list[[j]]$gum
  nnmp_out_clay <- nnmp_out_list[[j]]$clay
  
  set.seed(42)
  gaus_nnmp_pred <- predict(nnmp_out_gaus, nonref_coords = nonref_grid,
                            probs = c(0.025, 0.5, 0.975), predict_sam = TRUE,
                            compute_control = 
                              list(lb = -1, ub = max(nnmp_out_gaus$orig_data$response) + 100,
                                   tol = .Machine$double.eps^0.25, maxit = 1000))
  
  set.seed(42)
  gum_nnmp_pred <- predict(nnmp_out_gum, nonref_coords = nonref_grid,
                           probs = c(0.025, 0.5, 0.975), predict_sam = TRUE,
                           compute_control = 
                             list(lb = -1, ub = max(nnmp_out_gum$orig_data$response) + 100,
                                  tol = .Machine$double.eps^0.25, maxit = 1000))
  
  set.seed(42)
  clay_nnmp_pred <- predict(nnmp_out_clay, nonref_coords = nonref_grid,
                            probs = c(0.025, 0.5, 0.975), predict_sam = TRUE,
                            compute_control = 
                              list(lb = -1, ub = max(nnmp_out_clay$orig_data$response) + 100,
                                   tol = .Machine$double.eps^0.25, maxit = 1000))
  
  gaus_nnmp_surf <- as.data.frame(cbind(nonref_grid, gaus_nnmp_pred[[1]], t(gaus_nnmp_pred[[2]])))
  gum_nnmp_surf <- as.data.frame(cbind(nonref_grid, gum_nnmp_pred[[1]], t(gum_nnmp_pred[[2]])))
  clay_nnmp_surf <- as.data.frame(cbind(nonref_grid, clay_nnmp_pred[[1]], t(clay_nnmp_pred[[2]])))
  colnames(gaus_nnmp_surf) <- c("x", "y", "mu", "q025", "q50", "q975")
  colnames(gum_nnmp_surf) <- c("x", "y", "mu", "q025", "q50", "q975")
  colnames(clay_nnmp_surf) <- c("x", "y", "mu", "q025", "q50", "q975")
  
  #-----------------------------------------------------------------------------
  # Gaussian copula
  #-----------------------------------------------------------------------------
  plot_gaus_nnmp <- ggplot() + theme_bw() +
    geom_raster(aes(x = x, y = y, fill = q50), gaus_nnmp_surf, interpolate = T) +
    scale_fill_gradientn(colours = myPalette(100), limits = c(0, max(dat[,3]))) + 
    labs(fill = "", x = "Easting", y = "Northing") +  
    theme(plot.margin = unit(c(0.5,0,0,0), "cm")) +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.key.height= unit(2.2, 'cm'),
          axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
          axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent"))
  
  #-----------------------------------------------------------------------------
  # Gumbel copula
  #-----------------------------------------------------------------------------
  plot_gum_nnmp <- ggplot() + theme_bw() +
    geom_raster(aes(x = x, y = y, fill = q50), gum_nnmp_surf, interpolate = T) +
    scale_fill_gradientn(colours = myPalette(100), limits = c(0, max(dat[,3]))) +
    labs(fill = "", x = "Easting", y = "Northing") +
    theme(plot.margin = unit(c(0.5,0,0,0), "cm")) +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.key.height= unit(2.2, 'cm'),
          axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
          axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent"))
  
  #-----------------------------------------------------------------------------
  # Clayton copula
  #-----------------------------------------------------------------------------
  plot_clay_nnmp <- ggplot() + theme_bw() +
    geom_raster(aes(x = x, y = y, fill = q50), clay_nnmp_surf, interpolate = T) +
    scale_fill_gradientn(colours = myPalette(100), limits = c(0, max(dat[,3]))) +
    labs(fill = "", x = "Easting", y = "Northing") +
    theme(plot.margin = unit(c(0.5,0,0,0), "cm")) +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.key.height= unit(2.2, 'cm'),
          axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
          axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent"))
  
  plot_gaus_nnmp_list[[j]] <- plot_gaus_nnmp
  plot_gum_nnmp_list[[j]] <- plot_gum_nnmp
  plot_clay_nnmp_list[[j]] <- plot_clay_nnmp
  
}



#-------------------------------------------------------------------------------
# Second to Fourth rows in Figure 1 of the supplementary material 
#-------------------------------------------------------------------------------
nne <- 10
for (j in seq_along(skewparams)) {
  
  png(paste0(outpath, "skew_", skewparams[j], "_", "PoNNMP_gaus_median_", nne, ".png"), 
      units = "in", width = 7, height = 6, res = 300, bg = "transparent")
  print(plot_gaus_nnmp_list[[j]])
  dev.off()
  
  png(paste0(outpath, "skew_", skewparams[j], "_", "PoNNMP_gum_median_", nne, ".png"), 
      units = "in", width = 7, height = 6, res = 300, bg = "transparent")
  print(plot_gum_nnmp_list[[j]])
  dev.off()
  
  png(paste0(outpath, "skew_", skewparams[j], "_", "PoNNMP_clay_median_", nne, ".png"), 
      units = "in", width = 7, height = 6, res = 300, bg = "transparent")
  print(plot_clay_nnmp_list[[j]])
  dev.off()
  
}

