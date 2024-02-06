################################################################################
### BBS DATA ANALYSIS - Analysis of L; model estimation
### Section D.2.1 of the supplementary material
################################################################################

rm(list = ls())

outpath <- "/path/to/output/"

library(ggplot2)
library(readr)
library(dnnmp)


#-------------------------------------------------------------------------------
# Input data
#-------------------------------------------------------------------------------
dat <- read_csv(file = "/path/to/load/bbs_dat.csv")

#-------------------------------------------------------------------------------
# Fit the model with different values of L
#-------------------------------------------------------------------------------

ref_dat <- as.matrix(dat)
ref_XX <- cbind(1, ref_dat[, 2])

# Ordering
set.seed(42)
ord <- sample(1:nrow(ref_dat), replace = FALSE)

# Estimation
prior <- list("scale_gamma" = c(1, 1))
starting <- list(phi = 1, zeta  = 0.1, ga = rep(0, 3), kasq = 1, scale = 1,
                 regcoef = c(log(mean(ref_dat[,3])), rep(0, ncol(ref_XX) - 1)))
tuning <- list(regcoef = c(0.1, 0.003), scale = 0.12, phi = 0.12, zeta = 0.15)
mcmc_param <- list(n_iter = 30000, n_burn = 10000, n_thin = 5)

se_zeta <- c(0.25, 0.25, 0.2, 0.15)
ne_sizes <- c(5, 10, 15, 20)
nnmp_out_list <- vector("list", length = length(ne_sizes))

for (i in seq_along(ne_sizes)) {
  
  nne <- ne_sizes[i]
  print(paste0("Start with L = ", nne))
  
  set.seed(42)
  tuning$zeta <- se_zeta[i]
  nnmp_out <- dnnmp(response = ref_dat[,3],
                    covars = ref_XX,
                    coords = ref_dat[, 1:2],
                    neighbor_size = nne,
                    marg_family = "negative binomial",
                    cop_family = "gaussian",
                    priors = prior,
                    starting = starting,
                    tuning = tuning,
                    ord = ord,
                    mcmc_settings = mcmc_param,
                    weights_samples = TRUE,
                    verbose = TRUE)
  
  nnmp_out_list[[i]] <- nnmp_out
  
  print(paste0("Done with L = ", nne))
  
}

save(nnmp_out_list, ref_dat, file = paste0(outpath, "sa_1_gaus_nbnnmp_out.RData"))
