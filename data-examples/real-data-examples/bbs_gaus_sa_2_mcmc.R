################################################################################
### BBS DATA ANALYSIS - Sensitivity analysis of L; predictive performance
### Section D.2.1 of the supplementary material
################################################################################

rm(list = ls())

outpath <- "/path/to/output/"

library(readr)
library(dnnmp)


#-------------------------------------------------------------------------------
# Input data
#-------------------------------------------------------------------------------
dat <- read_csv(file = "/path/to/load/bbs_dat.csv")

#-------------------------------------------------------------------------------
# Fit the model with different values of L
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

ne_sizes <- 5:20
nnmp_out_list <- vector("list", length = length(ne_sizes))

for (i in seq_along(ne_sizes)) {
  
  nne <- ne_sizes[i]

  print(paste0("Start with L = ", nne))
  
  set.seed(42)
  nnmp_out <- dnnmp(response = selec_dat[,3],
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
  
  nnmp_out_list[[i]] <- nnmp_out
  
  print(paste0("Done with L = ", nne))
  
}

save(nnmp_out_list, valid_dat, valid_XX, 
     file = paste0(outpath, "sa2_gaus_nbnnmp_out_L_5to20.RData"))

