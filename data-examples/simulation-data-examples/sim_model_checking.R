###################################################################################
### SIMULATION STUDY
### Data set generated from a Poisson spatial generalized linear mixed model (SGLMM)
###################################################################################

rm(list = ls())

inpath <- "/path/to/mcmc output/"
outpath <- "/path/to/output/"

library(ggplot2)
library(gridExtra)
library(dnnmp)


#-------------------------------------------------------------------------------
# Function to plot randomized quantile residuals
#-------------------------------------------------------------------------------
plotRQR <- function(rqr, coords){
  
  q_mu <- rowMeans(rqr)
  q_025 <- apply(rqr, 1, function(x) quantile(x, 0.025))
  q_975 <- apply(rqr, 1, function(x) quantile(x, 0.975))
  qq <- qqnorm(q_mu[order(q_mu)], plot.it = FALSE)
  resid <- data.frame(lon = coords[,1], lat = coords[, 2],
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
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
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
  
  # histogram
  xx <- seq(-3.5, 3.5, length = 512)
  zz <- dnorm(xx)
  kde_sam_y <- apply(rqr, 2, function(x) density(x, from = -3.5, to = 3.5, bw = 0.5)$y)
  kde_mu <- rowMeans(kde_sam_y)
  kde <- data.frame(xx = xx, zz = zz, mu = kde_mu)
  
  plot_hist_save <- ggplot(data = resid, aes(x = mu)) + theme_bw() + 
    geom_histogram(aes(y=..density..), binwidth = 0.7, fill="white", color="darkgrey", size = 0.8) + 
    labs(x = "Posterior mean residuals", y = "Density") + 
    scale_x_continuous(breaks = seq(-3.5, 3.5, by = 1)) + ylim(0, 0.5) + 
    geom_line(data = kde, aes(x = xx, y = zz), size = 0.8) +
    geom_line(data = kde, aes(x = xx, y = mu), linetype = "dashed", size = 0.8) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
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
  
  # residual map
  myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))
  plot_resid_save <- ggplot() + theme_bw() +
    geom_point(aes(lon, lat, color = mu), data = resid) + theme_bw() + 
    scale_colour_gradientn(colours = myPalette(100), limits = c(-3.5, 3.5)) + 
    labs(color = "", x = "Easting", y = "Northing") + xlim(0, 1) + ylim(0, 1) + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
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
  
  list(plot_gg_save, plot_hist_save, plot_resid_save)
  
}

#-------------------------------------------------------------------------------
# Figures 4, 5, and 6 of the supplementary material
#-------------------------------------------------------------------------------

load(paste0(inpath1, "sim1_res.RData"))

for (j in seq_along(nnmp_out_list)) {
  
  nnmp_out_gaus <- nnmp_out_list[[j]][[1]]
  nnmp_out_gum <- nnmp_out_list[[j]][[2]]
  nnmp_out_clay <- nnmp_out_list[[j]][[3]]
  
  gaus_rqr <- rqr(nnmp_out_gaus)
  gum_rqr <- rqr(nnmp_out_gum)
  clay_rqr <- rqr(nnmp_out_clay)
  
  gaus_rar_plot <- plotRQR(gaus_rqr, nnmp_out_gaus$ord_data$coords_ord)
  gum_rar_plot <- plotRQR(gum_rqr, nnmp_out_gum$ord_data$coords_ord)
  clay_rar_plot <- plotRQR(clay_rqr, nnmp_out_clay$ord_data$coords_ord)
  
  setEPS()
  postscript(paste0(outpath, "sim1_s", j, "_gaus_rqr", ".eps"), bg = "transparent", width = 21, height = 6)
  grid.arrange(gaus_rar_plot[[1]], gaus_rar_plot[[2]], gaus_rar_plot[[3]], nrow = 1)
  dev.off()
  
  setEPS()
  postscript(paste0(outpath, "sim1_s", j, "_gum_rqr", ".eps"), bg = "transparent", width = 21, height = 6)
  grid.arrange(gum_rar_plot[[1]], gum_rar_plot[[2]], gum_rar_plot[[3]], nrow = 1)
  dev.off()
  
  setEPS()
  postscript(paste0(outpath, "sim1_s", j, "_clay_rqr", ".eps"), bg = "transparent", width = 21, height = 6)
  grid.arrange(clay_rar_plot[[1]], clay_rar_plot[[2]], clay_rar_plot[[3]], nrow = 1)
  dev.off()
 
}


#-------------------------------------------------------------------------------
# Figure 7 of the supplementary material
#-------------------------------------------------------------------------------

load(paste0(inpath2, "nbnnmp_res.RData"))

gaus_rqr <- rqr(nnmp_out)
gaus_rar_plot <- plotRQR(gaus_rqr, nnmp_out$ord_data$coords_ord)

setEPS()
postscript(paste0(outpath, "sim2_nbnnmp_rqr", ".eps"), bg = "transparent", width = 21, height = 6)
grid.arrange(gaus_rar_plot[[1]], gaus_rar_plot[[2]], gaus_rar_plot[[3]], nrow = 1)
dev.off()





















