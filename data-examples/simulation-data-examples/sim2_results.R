################################################################################
### SIMULATION DATA EXAMPLE 2 
### Prediction (NBNNMP, SGLMM-GP, and SGLMM-GPP)
################################################################################

rm(list = ls())

inpath <- "/path/to/mcmc output/"
outpath <- "/path/to/output/"

library(spBayes)
library(ggplot2)
library(dnnmp)

load(paste0(inpath, "nbnnmp_res.RData"))
load(paste0(inpath, "sglmgp_res.RData"))
load(paste0(inpath, "sglmgpp_res.RData"))

myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))


#-------------------------------------------------------------------------------
# NBNNMP prediction
#-------------------------------------------------------------------------------
set.seed(42)

new_grid <- as.matrix(expand.grid(seq(0.01, 1, by = 1/120), seq(0.01, 1, by = 1/120)))

nbnnmp_pred <- predict(nnmp_out, nonref_covars = cbind(1, new_grid[, 1:2]),
                       nonref_coords = new_grid[,1:2],
                       probs = c(0.025, 0.5, 0.975), predict_sam = TRUE)

nbnnmp_surf <- as.data.frame(cbind(new_grid, nbnnmp_pred[[1]], t(nbnnmp_pred[[2]])))
colnames(nbnnmp_surf) <- c("x", "y", "mu", "q025", "q50", "q975")

plot_save_nbnnmp <- ggplot(as.data.frame(nbnnmp_surf), aes(x = x, y = y, fill = q50)) +
  geom_raster(interpolate = T) +
  scale_fill_gradientn(colours = myPalette(100), limits = range(dat[,3])) + 
  theme_bw() +
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

# setEPS()
# postscript(paste0(outpath, "NBNNMP_median_", nnmp_out$neighbor_size, ".eps"), width = 7, height = 6, bg = "transparent")
# plot_save_nbnnmp
# dev.off()
png(paste0(outpath, "NBNNMP_median_", nnmp_out$neighbor_size, ".png"), units = "in", width = 7, height = 6, res = 300,
    bg = "transparent")
plot_save_nbnnmp
dev.off()



#-------------------------------------------------------------------------------
# SGLMM-GP
#-------------------------------------------------------------------------------
set.seed(42)
n.samples <- n.batch * batch.length
sglm_gp_pred <- spPredict(sim.s1,
                          pred.coords = new_grid,
                          pred.covars = as.matrix(cbind(1, new_grid)),
                          start = 20001,
                          end = n.samples,
                          thin = 5)

save(sglm_gp_pred, file = paste0(outpath, "sglmgp_pred_120x120.RData"))
# load(paste0(outpath, "sglmgp_pred_120x120.RData"))

sglm_gp_sam <- apply(sglm_gp_pred$p.y.predictive.samples, 2, function(x) rpois(nrow(new_grid), x))
sglm_gp_q50 <- apply(sglm_gp_sam, 1, median)

sglm_gp_surf <- as.data.frame(cbind(new_grid, sglm_gp_q50))
colnames(sglm_gp_surf) <- c("x", "y", "z")

plot_save_sglm_gp <- ggplot(as.data.frame(sglm_gp_surf), aes(x = x, y = y, fill = z)) +
  geom_raster(interpolate = T) +
  scale_fill_gradientn(colours = myPalette(100), limits = range(dat[,3])) + 
  theme_bw() +
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

# setEPS()
# postscript(paste0(outpath, "sglm_gp_median", ".eps"), width = 7, height = 6, bg = "transparent")
# plot_save_sglm_gp
# dev.off()
png(paste0(outpath, "sglm_gp_median",".png"), units = "in", width = 7, height = 6, res = 300, bg = "transparent")
plot_save_sglm_gp
dev.off()

#-------------------------------------------------------------------------------
# SGLMM-GPP
#-------------------------------------------------------------------------------
# Prediction
set.seed(42)
sglm_gpp_pred <- spPredict(sim.s2, 
                           pred.coords = new_grid,
                           pred.covars = as.matrix(cbind(1, new_grid)),
                           start = 20001,
                           end = n.samples,
                           thin = 5)

sglm_gpp_sam <- apply(sglm_gpp_pred$p.y.predictive.samples, 2, function(x) rpois(nrow(new_grid), x))
sglm_gpp_q50 <- apply(sglm_gpp_sam, 1, median)

sglm_gpp_surf <- as.data.frame(cbind(new_grid, sglm_gpp_q50))
colnames(sglm_gpp_surf) <- c("x", "y", "z")

plot_save_sglm_gpp <- ggplot(as.data.frame(sglm_gpp_surf), aes(x = x, y = y, fill = z)) +
  geom_raster(interpolate = T) +
  scale_fill_gradientn(colours = myPalette(100), limits = range(dat[,3])) + 
  theme_bw() +
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

# setEPS()
# postscript(paste0(outpath, "sglm_gpp_median", ".eps"), width = 7, height = 6, bg = "transparent")
# plot_save_sglm_gpp
# dev.off()
png(paste0(outpath, "sglm_gpp_median",".png"), units = "in", width = 7, height = 6, res = 300, bg = "transparent")
plot_save_sglm_gpp
dev.off()

