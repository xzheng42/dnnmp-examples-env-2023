################################################################################
### BBS DATA ANALYSIS - Analysis of L
### Section D.2.1 of the supplementary material
################################################################################

rm(list = ls())

inpath <- "/path/to/mcmc output/"
outpath <- "/path/to/output/"

library(ggplot2)
library(dnnmp)
library(readr)

load(paste0(inpath, "sa_1_gaus_nbnnmp_out.RData"))


#-------------------------------------------------------------------------------
# Table 1 of the supplementary material 
#-------------------------------------------------------------------------------
ne_sizes <- sapply(nnmp_out_list, function(x) x$mod_spec$neighbor_size)
tbl_est <- array(NA, dim = c(10, length(ne_sizes)))
rownames(tbl_est) <- c("$\\beta_0$", "$\\beta_1$", "$r$", "$\\phi$", "$\\zeta$", 
                       "$\\gamma_0$", "$\\gamma_1$", "$\\gamma_2$", "$\\kappa^2$",
                       "Time")
colnames(tbl_est) <- paste0("L = ", ne_sizes)

for (i in seq_along(ne_sizes)) {
  
  nnmp_out <- nnmp_out_list[[i]]
  runtime <- nnmp_out$runtime
  
  tbl_est[,i] <- c(summary(nnmp_out), round(runtime[3]/60, 2))

}

tbl_est_latex <- knitr::kable(tbl_est, digits = 2, align = "c", format = "latex", escape = FALSE)
readr::write_file(tbl_est_latex, paste0(outpath, "nbnnmp_gaus_sa_tbl1.txt"))

#-------------------------------------------------------------------------------
# Figures 2 and 3 of the supplementary material
#-------------------------------------------------------------------------------
set.seed(42)
obs_idx1 <- sample(21:200, 5)
obs_idx2 <- sample(1312:1512, 5)
obs_idx <- c(obs_idx1, obs_idx2)

ord <- nnmp_out_list[[1]]$ord_data$ord
ref_dat_ord <- ref_dat[ord, ]
coord <- ref_dat_ord[,1:2]

dd <- 5000
zeta_prior <- nnmp_out_list[[1]]$mod_spec$priors$zeta_invgamma
ga_prior <- nnmp_out_list[[1]]$mod_spec$priors$ga_gaus
kasq_prior <- nnmp_out_list[[1]]$mod_spec$priors$kasq_invgamma

# Calculate prior means and 95% CIs using prior samples
set.seed(42)
zeta_prior_sam <- rgamma(dd, zeta_prior[1], zeta_prior[2])
ga_prior_sam <- t(mvtnorm::rmvnorm(dd, ga_prior[[1]], ga_prior[[2]]))
g0_mu_prior_sam <- ga_prior_sam[1,]
kasq_prior_sam <- 1 / rgamma(dd, kasq_prior[1], kasq_prior[2])
xx_grid <- seq(0, 1, length = 100)
yy_prior_sam <- sapply(1:dd, function(x) logitnorm::dlogitnorm(xx_grid, g0_mu_prior_sam[x], sqrt(kasq_prior_sam[x])))
yy_prior_mu <- rowMeans(yy_prior_sam)
yy_prior_qq <- apply(yy_prior_sam, 1, function(x) quantile(x, probs = c(0.025, 0.975)))

# Plot prior means and 95% CIs and posterior means and 95% CIs of the weights
corrFunc <- function(d, xx) return(exp(-d / xx)) 
for (i in seq_along(ne_sizes)) {
  
  nne <- ne_sizes[i]
  post_sams <- nnmp_out_list[[i]]$post_samples
  ne_info <- neighbor(nne, ref_dat_ord[,1:2], ref_dat_ord[,3])
  dist_mat <- ne_info$ne_dist
  
  for (j in obs_idx) {
    
    g0_mu_sam <- apply(post_sams$ga, 2, function(x) sum(c(1, ref_dat_ord[j,1:2]) * x))
    prior_un_cutoff_sam <- sapply(zeta_prior_sam, function(x) corrFunc(dist_mat[j, ], x))
    prior_cutoff_sam <- apply(prior_un_cutoff_sam, 2, function(x) c(0, cumsum(x / sum(x))))
    prior_cutoff_sam[(nne + 1), ] <- 1
    un_cutoff <- corrFunc(dist_mat[j, ], Inf)
    cutoff <- c(0, cumsum(un_cutoff) / sum(un_cutoff))
    
    prior_w_j_sam <- sapply(1:dd, function(x) 
      diff(pnorm(qlogis(prior_cutoff_sam[,x]), g0_mu_prior_sam[x], sqrt(kasq_prior_sam[x]))))
    un_cutoff_sam <- sapply(post_sams$zeta, function(x) corrFunc(dist_mat[j, ], x))
    cutoff_sam <- apply(un_cutoff_sam, 2, function(x) c(0, cumsum(x / sum(x))))
    cutoff_sam[nne + 1, ] <- 1
    w_j_sam <- post_sams$weight[j,,] 
    
    prior_mean <- rowMeans(prior_w_j_sam, na.rm = TRUE)
    prior_qq <- apply(prior_w_j_sam, 1, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
    w_j_mu <- rowMeans(w_j_sam)
    w_j_qq <- apply(w_j_sam, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
    weight_j <- as.data.frame(cbind(prior_mean, t(prior_qq), w_j_mu, t(w_j_qq)))
    weight_j <- rbind(weight_j, tail(weight_j, 1))
    weight_j <- cbind(0:nne, weight_j)
    colnames(weight_j) <- c("idx", "prior_mu", "prior_q_025", "prior_q_975", "post_mu", "q_025", "q_975")
    
    plot_weight_save <- ggplot(weight_j, aes(x = idx)) + 
      pammtools::geom_stepribbon(aes(ymin = prior_q_025, ymax = prior_q_975, fill = "Prior 95% CI"), alpha = 0.3) + 
      pammtools::geom_stepribbon(aes(ymin = q_025, ymax = q_975, fill = "Posterior 95% CI"), alpha = 0.3) + 
      geom_step(aes(y = prior_mu, color = "Prior mean"), linetype = "dashed", linewidth = 1.5) + 
      geom_step(aes(y = post_mu, color = "Posterior mean"), size = 1.5) + 
      scale_fill_manual("", values = c("blue", "brown1")) + 
      scale_colour_manual("", values = c("blue", "brown1")) + 
      scale_x_continuous(breaks = 0:nne, labels = 0:nne) + 
      theme_bw() + ylim(0, 1) + 
      labs(fill = "", x = "Neighbors", y = "Probability", 
           title = substitute(paste("Weights at site ", s_j, " (L = ", nne, ")"), list(s_j = j, nne = nne))) + 
      theme(legend.position = "top",
            plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 25),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 12.5))

    png(paste0(outpath, "L_", nne, "_weights_", j, ".png"), width=600, height=500, pointsize=20)
    print(plot_weight_save)
    dev.off()
    
  }
  
}
