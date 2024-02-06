################################################################################
### BBS DATA ANALYSIS - model comparison
### Section D.2.2 of the supplementary material
################################################################################

rm(list = ls())

outpath <- "/path/to/output/"

library(dplyr)
library(readr)
library(dnnmp)


#-------------------------------------------------------------------------------
# Input data
#-------------------------------------------------------------------------------
dat <- read_csv(file = "/path/to/bbs_dat.csv")
  
#-------------------------------------------------------------------------------
# Fit three copula NBNNMP models
#-------------------------------------------------------------------------------
ref_num <- nrow(dat)
valid_num <- 300
selec_num <- ref_num - valid_num

set.seed(42)
selec_idx <- sample(1:ref_num, selec_num, FALSE)

ref_dat <- as.matrix(dat)
ref_XX <- cbind(1, ref_dat[, 2])

selec_dat <- ref_dat[selec_idx,]
selec_XX <- ref_XX[selec_idx, ]

valid_dat <- ref_dat[-selec_idx,]
valid_XX <- ref_XX[-selec_idx, ]

# Ordering
ord <- sample(1:selec_num, replace = FALSE)

# MCMC
prior <- list("scale_gamma" = c(1, 1))
starting <- list(phi = 1, zeta  = 0.1, ga = rep(0, 3), kasq = 1, scale = 1,
                 bb = c(log(mean(selec_dat[,3])), rep(0, ncol(selec_XX) - 1)))
tuning <- list(regcoef = c(0.1, 0.003), scale = 0.12,  phi = 0.15, zeta = 0.15)
mcmc_param <- list(n_iter = 30000, n_burn = 10000, n_thin = 5)

nne <- 20

set.seed(42)
nnmp_gaus_out <- dnnmp(response = selec_dat[,3],
                       covars = selec_XX,
                       coords = selec_dat[, 1:2],
                       neighbor_size = nne,
                       marg_family = "negative binomial",
                       cop_family = "gaussian",
                       priors = prior,
                       starting = starting,
                       tuning = tuning,
                       ord = ord,
                       mcmc_settings = mcmc_param,
                       verbose = TRUE)

set.seed(42)
tuning$phi <- 0.1
tuning$regcoef <- c(0.1, 0.002)
nnmp_gum_out <- dnnmp(response = selec_dat[,3],
                      covars = selec_XX,
                      coords = selec_dat[, 1:2],
                      neighbor_size = nne,
                      marg_family = "negative binomial",
                      cop_family = "gumbel",
                      priors = prior,
                      starting = starting,
                      tuning = tuning,
                      ord = ord,
                      mcmc_settings = mcmc_param,
                      verbose = TRUE)

set.seed(42)
tuning$regcoef <- c(0.08, 0.002)
nnmp_clay_out <- dnnmp(response = selec_dat[,3],
                       covars = selec_XX,
                       coords = selec_dat[, 1:2],
                       neighbor_size = nne,
                       marg_family = "negative binomial",
                       cop_family = "clayton",
                       priors = prior,
                       starting = starting,
                       tuning = tuning,
                       ord = ord,
                       mcmc_settings = mcmc_param,
                       verbose = TRUE)

save(nnmp_gaus_out, nnmp_gum_out, nnmp_clay_out, valid_dat, valid_XX, 
     file = paste0(outpath, "comp_res_20.RData"))


#-------------------------------------------------------------------------------
# Fit the best model
#-------------------------------------------------------------------------------

ref_dat <- as.matrix(dat)
ref_XX <- cbind(1, ref_dat[, 2])

# Ordering
ord <- sample(1:nrow(ref_dat), replace = FALSE)

# MCMC
prior <- list("scale_gamma" = c(1, 1))
starting <- list(phi = 1, zeta  = 0.1, ga = rep(0, 3), kasq = 1, scale = 1,
                 bb = c(log(mean(ref_dat[,3])), rep(0, ncol(ref_XX) - 1)))
tuning <- list(regcoef = c(0.1, 0.003), scale = 0.12,  phi = 0.12, zeta = 0.15)
mcmc_param <- list(n_iter = 30000, n_burn = 10000, n_thin = 5)
nne <- 20

set.seed(42)
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
 
save(nnmp_out, file = paste0(outpath, "best_res_20.RData"))
