################################################################################
### BBS DATA ANALYSIS 
### Section 5.3
################################################################################

rm(list = ls())

inpath <- "/path/to/mcmc output/"
outpath <- "/path/to/output/"

library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggplot2)
library(dnnmp)
library(readr)

myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))

load(paste0(inpath, "best_res_20.RData"))
ref_obs <- nnmp_out$orig_data$response


#-------------------------------------------------------------------------------
# Figure 2 in Section 5.3
#-------------------------------------------------------------------------------
usa <- rnaturalearth::ne_states(returnclass = "sf", country = "United States of America")
bbox <- c(minlon = -102, maxlon = -66, minlat = 24, maxlat = 50)
lon_grid <- seq(bbox["minlon"], bbox["maxlon"], length = 120)
lat_grid <- seq(bbox["minlat"], bbox["maxlat"], length = 120)
whole_grid <- expand.grid(lon_grid, lat_grid)
colnames(whole_grid) <- c("lon", "lat")
grid_in_box <- whole_grid[which(whole_grid$lon >= bbox["minlon"] & whole_grid$lon <= bbox["maxlon"] & 
                                whole_grid$lat >= bbox["minlat"] & whole_grid$lat <= bbox["maxlat"]), ]
grid_in_land_idx <- !is.na(maps::map.where("usa", grid_in_box[,1, drop = TRUE], grid_in_box[,2, drop = TRUE]))
grid_loc <- as.matrix(grid_in_box[grid_in_land_idx, ])


set.seed(42)
nbnnmp_pred <- predict(nnmp_out, nonref_covars = cbind(1, grid_loc[, 2]),
                       nonref_coords = grid_loc[, 1:2],
                       probs = c(0.025, 0.5, 0.975), 
                       predict_sam = TRUE)

nbnnmp_surf <- as.data.frame(cbind(grid_loc, nbnnmp_pred[[1]], t(nbnnmp_pred[[2]])))
colnames(nbnnmp_surf) <- c("x", "y", "mu", "q025", "q50", "q975")


base <- ggplot(data = usa) + geom_sf() +
  coord_sf(xlim = c(bbox["minlon"], bbox["maxlon"]),
           ylim = c(bbox["minlat"], bbox["maxlat"]),
           expand = FALSE)


plot_save_nbnnmp <- base + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = q50), nbnnmp_surf, interpolate = T) +
  scale_fill_gradientn(colours = myPalette(100), limits = c(0, max(ref_obs))) + 
  labs(fill = "", x = "Longitude", y = "Latitude") + 
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
setEPS()
postscript(paste0(outpath, "NBNNMP_median_", ".eps"), width = 7, height = 6, bg = "transparent")
plot_save_nbnnmp
dev.off()


plot_save_nbnnmp <- base + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = q975 - q025), nbnnmp_surf, interpolate = T) +
  scale_fill_gradientn(colours = myPalette(100)) + 
  labs(fill = "", x = "Longitude", y = "Latitude") + 
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
setEPS()
postscript(paste0(outpath, "NBNNMP_95CI_width_", ".eps"), width = 7, height = 6, bg = "transparent")
plot_save_nbnnmp
dev.off()


bb_sam <- nnmp_out$post_samples$regcoef
mu_surf_sam <- sapply(1:nrow(bb_sam), function(x) 
                      exp(cbind(1, grid_loc[,2]) %*% bb_sam[,x]))
mu_surf <- data.frame(x = grid_loc[,1], y = grid_loc[,2], z = rowMeans(mu_surf_sam))
plot_save_nbnnmp <- base + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = z), mu_surf, interpolate = T) +
  scale_fill_gradientn(colours = myPalette(100), limits = c(0, max(ref_obs))) + 
  labs(fill = "", x = "Longitude", y = "Latitude") + 
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
setEPS()
postscript(paste0(outpath, "NBNNMP_marg_mu_", ".eps"), width = 7, height = 6, bg = "transparent")
plot_save_nbnnmp
dev.off()