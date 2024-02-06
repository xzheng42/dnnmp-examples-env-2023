################################################################################
### SIMULATION DATA EXAMPLE 1
################################################################################
rm(list = ls())

outpath <- "/path/to/output/"

library(SpatialExtremes)
library(sn)
library(dnnmp)


#-------------------------------------------------------------------------------
# Generate synthetic data
#-------------------------------------------------------------------------------
set.seed(42)

# Simulate latent processes on a grid
x_seq <- seq(0, 1, length = 120)
y_seq <- seq(0, 1, length = 120)
grid_loc <- as.matrix(expand.grid(x_seq, y_seq))
nn <- nrow(grid_loc)

zz2 <- t(rgp(1, grid_loc, cov.mod = "powexp", nugget =  0, sill = 1, range = 0.1, smooth = 1))
zz1 <- abs(t(rgp(1, grid_loc, cov.mod = "powexp", nugget =  0, sill = 1, range = 0.1, smooth = 1)))

# Simulate data
skewparams <- c(1, 3, 10)
dat_save <- vector("list", length = length(skewparams))

for (j in seq_along(skewparams)) {
  
  # Generate latent field
  skewparam <- skewparams[j]
  zz <- skewparam * zz1  +  zz2
  
  # Generate data
  pp <- psn(zz, 0, sqrt(skewparam^2 + 1), skewparam)
  yy <- qpois(pp, 5)
  dat <- cbind(grid_loc, yy)
  colnames(dat) <- c("x", "y", "z")

  # Save data
  dat_save[[j]] <- dat
  
}

#-------------------------------------------------------------------------------
# Fit Poisson NNMP with different copulas
#-------------------------------------------------------------------------------
mcmc_param <- list(n_iter = 20000, n_burn = 4000, n_thin = 4, n_report = 2000)

priors <- list("rate_gamma" = c(1, 1), 
               "phi_invgamma" = c(3, 1),
               "zeta_invgamma" = c(3, 1))

starting <- list("phi" = 0.01, "zeta" = 0.1, 
                 "ga" = rep(0, 3), "kasq" = 1, "rate" = 1)

se_rate_mat <- matrix(c(.1, .08, .07, 
                        .1, .08, .06, 
                        .08, .07, .06), 
                      nrow = 3, ncol = 3, byrow = TRUE)
se_phi_mat <- matrix(c(.15, .1, .1, 
                       .15, .1, .1, 
                       .15, .1, .1), 
                      nrow = 3, ncol = 3, byrow = TRUE)
se_zeta_mat <- matrix(c(.6, .5, .6, 
                        .6, .5, .7, 
                        .6, .5, .8), 
                      nrow = 3, ncol = 3, byrow = TRUE)

nne <- 10

nnmp_out_list <- vector("list", length = length(skewparams))
test_dat_save <- vector("list", length = length(skewparam))

for (j in seq_along(skewparams)) {
 
  dat <- dat_save[[j]]
  
  set.seed(42)
  # Randomly select part of the gridded data as observations
  obs_num <- 1000
  obs_idx <- sample(1:nrow(grid_loc), obs_num, FALSE)
  obs_dat <- dat[obs_idx,]
  
  # Partition the observations into training and testing sets
  # Use the whole training set as reference set for NNMP
  ref_num <- 800
  test_num <- obs_num - ref_num
  ref_idx <- sample(1:obs_num, ref_num, FALSE)
  ref_dat <- obs_dat[ref_idx,]
  test_dat_save[[j]] <- obs_dat[-ref_idx,]
  
  # Random ordering
  ord <- sample(1:ref_num, replace = FALSE)
  
  # MCMC
  set.seed(42)
  tuning <- list("rate" = se_rate_mat[j, 1], 
                 "phi" = se_phi_mat[j, 1], 
                 "zeta" = se_zeta_mat[j, 1])
  nnmp_out_gaus <- dnnmp(response = ref_dat[,3],
                         coords = ref_dat[,1:2],
                         neighbor_size = nne,
                         marg_family = "poisson",
                         cop_family = "gaussian",
                         priors = priors,
                         starting = starting,
                         tuning = tuning,
                         ord = ord,
                         mcmc_settings = mcmc_param, 
                         weights_samples = TRUE,
                         verbose = TRUE)
  
  cat(paste0("Done the ", j, "th run (", nnmp_out_gaus$mod_spec$cop_family, " copula).\n"))
  
  set.seed(42)
  tuning <- list("rate" = se_rate_mat[j, 2], 
                 "phi" = se_phi_mat[j, 2], 
                 "zeta" = se_zeta_mat[j, 2])
  nnmp_out_gum <- dnnmp(response = ref_dat[,3],
                        coords = ref_dat[,1:2],
                        neighbor_size = nne,
                        marg_family = "poisson",
                        cop_family = "gumbel",
                        priors = priors,
                        starting = starting,
                        tuning = tuning,
                        ord = ord,
                        mcmc_settings = mcmc_param, 
                        weights_samples = TRUE,
                        verbose = TRUE)
  
  cat(paste0("Done the ", j, "th run (", nnmp_out_gum$mod_spec$cop_family, " copula).\n"))
  
  set.seed(42)
  tuning <- list("rate" = se_rate_mat[j, 3], 
                 "phi" = se_phi_mat[j, 3], 
                 "zeta" = se_zeta_mat[j, 3])
  nnmp_out_clay <- dnnmp(response = ref_dat[,3],
                         coords = ref_dat[,1:2],
                         neighbor_size = nne,
                         marg_family = "poisson",
                         cop_family = "clayton",
                         priors = priors,
                         starting = starting,
                         tuning = tuning,
                         ord = ord,
                         mcmc_settings = mcmc_param, 
                         weights_samples = TRUE,
                         verbose = TRUE)
  
  nnmp_out_list[[j]] <- list(gaus = nnmp_out_gaus, 
                             gum = nnmp_out_gum, 
                             clay = nnmp_out_clay)
  
  cat(paste0("Done the ", j, "th run (", nnmp_out_clay$mod_spec$cop_family, " copula).\n"))
  
}

save(nnmp_out_list, dat_save, test_dat_save, skewparams, 
     file = paste0(outpath, "sim1_res.RData"))

