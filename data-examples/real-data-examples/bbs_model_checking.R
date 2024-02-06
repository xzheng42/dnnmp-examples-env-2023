################################################################################
### BBS DATA ANALYSIS: Model checking
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

myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))

load(paste0(inpath, "best_res_20.RData"))


#-------------------------------------------------------------------------------
# Figure 3 in Section 5.3
#-------------------------------------------------------------------------------
nbnnmp_rqr <- rqr(nnmp_out)
coords_ord <- nnmp_out$ord_data$coords_ord

q_mu <- rowMeans(nbnnmp_rqr)
q_025 <- apply(nbnnmp_rqr, 1, function(x) quantile(x, 0.025))
q_975 <- apply(nbnnmp_rqr, 1, function(x) quantile(x, 0.975))
qq <- qqnorm(q_mu[order(q_mu)], plot.it = FALSE)
resid <- data.frame(lon = coords_ord[,1], lat = coords_ord[,2],
                    z = qq$x, mu = qq$y, q025 = q_025, q975 = q_975)

# qq-plot
plot_gg_save <- ggplot(data = resid) + theme_bw() + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_line(aes(x = z, y = mu), linetype = "dotted", size = 0.8) +
  geom_line(aes(x = z, y = q025[order(q025)]), linetype = "dashed", size = 0.8) + 
  geom_line(aes(x = z, y = q975[order(q975)]), linetype = "dashed", size = 0.8) +
  scale_x_continuous(breaks = seq(-3.5, 3.5, by = 1)) + 
  scale_y_continuous(breaks = seq(-3.5, 3.5, by = 1)) + 
  labs(x = "Theoretical quantiles", y = "Sample quantiles") + 
  theme(plot.margin = unit(c(1, 0.5, 1, 0.5), "cm")) + 
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.key.height= unit(2.2, 'cm'),
        axis.text.x = element_text(margin = unit(c(0.3, 0, 0.3, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.3, 0, 0.3), "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent")) 
setEPS()
postscript(paste0(outpath, "NBNNMP_qqplot_", ".eps"), bg = "transparent", width = 7, height = 6)
plot_gg_save
dev.off()

# histogram
xx <- seq(-3.5, 3.5, length = 512)
zz <- dnorm(xx)
kde_sam_y <- apply(nbnnmp_rqr, 2, function(x) density(x, from = -3.5, to = 3.5, bw = 0.5)$y)
kde_mu <- rowMeans(kde_sam_y)
kde <- data.frame(xx = xx, zz = zz, mu = kde_mu)

plot_hist_save <- ggplot(data = resid, aes(x = mu)) + theme_bw() + 
  geom_histogram(aes(y=..density..), binwidth = 0.7, fill="white", color="darkgrey", size = 0.8) + 
  labs(x = "Posterior mean residuals", y = "Density") + 
  scale_x_continuous(breaks = seq(-3.5, 3.5, by = 1)) + ylim(0, 0.45) + 
  geom_line(data = kde, aes(x = xx, y = zz), size = 0.8) +
  geom_line(data = kde, aes(x = xx, y = mu), linetype = "dashed", size = 0.8) +
  theme(plot.margin = unit(c(1, 0.5, 1, 0.5), "cm")) + 
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.key.height= unit(2.2, 'cm'),
        axis.text.x = element_text(margin = unit(c(0.3, 0, 0.3, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.3, 0, 0.3), "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent")) 
setEPS()
postscript(paste0(outpath, "NBNNMP_kde_", ".eps"), bg = "transparent", width = 7, height = 6)
plot_hist_save
dev.off()


# residual map
usa <- rnaturalearth::ne_states(returnclass = "sf", country = "United States of America")
bbox <- c(minlon = -102, maxlon = -66, minlat = 24, maxlat = 50)
base <- ggplot(data = usa) + geom_sf() +
  coord_sf(xlim = c(bbox["minlon"], bbox["maxlon"]),
           ylim = c(bbox["minlat"], bbox["maxlat"]),
           expand = FALSE)


plot_resid_save <- base + 
  geom_point(aes(lon, lat, color = mu), data = resid) + theme_bw() + 
  scale_colour_gradientn(colours = myPalette(100), limits = c(-3.5, 3.5)) + 
  labs(color = "", x = "Longitude", y = "Latitude") + 
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.key.height= unit(2.2, 'cm'),
        axis.text.x = element_text(margin = unit(c(0.3, 0, 0.3, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.3, 0, 0.3), "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent")) 
setEPS()
postscript(paste0(outpath, "NBNNMP_resid_map_", ".eps"), bg = "transparent", width = 7, height = 6)
plot_resid_save
dev.off()

