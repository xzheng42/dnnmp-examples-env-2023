################################################################################
### SIMULATION DATA EXAMPLE 2 - NBNNMP
################################################################################

rm(list = ls())

outpath <- "/path/to/output/"

library(readr)
library(spBayes)
library(dnnmp)


#-------------------------------------------------------------------------------
# Input data
#-------------------------------------------------------------------------------
dat <- read_csv(file = "/path/to/load/bbs_dat.csv")


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

#-------------------------------------------------------------------------------
# Gaussian copula NBNNMP
#-------------------------------------------------------------------------------
# Ordering
ord <- sample(1:selec_num, replace = FALSE)

# Estimation
prior <- list("scale_gamma" = c(1, 1),
              "phi_invgamma" = c(3, 1),
              "zeta_invgamma" = c(3, 1))
starting <- list(phi = 0.01, zeta  = 0.01, ga = rep(0, 3), kasq = 1, scale = 1,
                 regcoef = as.matrix(c(0.1, 0.1)))
tuning <- list(regcoef = c(0.12, 0.003), scale = 0.12, phi = 0.14, zeta = 0.2)
mcmc_param <- list(n_iter = 30000, n_burn = 10000, n_thin = 5)

nne <- 20

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

save(nnmp_out, valid_dat, valid_XX, file = paste0(outpath, "nbnnmp_out_20.RData"))

#-------------------------------------------------------------------------------
# SGLMM-GP
#-------------------------------------------------------------------------------
beta_starting <-  coefficients(glm(selec_dat[,3] ~ selec_XX[,-1], family="poisson"))
n.batch <- 500
batch.length <- 100
n.samples <- n.batch * batch.length
starting <- list("beta" = beta_starting, "phi" = 3 / 0.5, "sigma.sq" = 1, "w" = 0)
tuning <- list("beta" = c(0.01, 0.001), "phi" = 1,  "sigma.sq" = 0.1, "w" = 0.5)
priors <- list("beta.Flat", "phi.Unif" = c(3/1, 3/0.1), "sigma.sq.IG" = c(2,1))

set.seed(42)
sim.s1 <- spGLM(formula = selec_dat[,3] ~ selec_XX[, -1],
                family = "poisson",
                coords = selec_dat[,1:2],
                starting = starting,
                tuning = tuning,
                priors = priors,
                cov.model = "exponential",
                n.samples = n.samples,
                amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
                n.report = 10)

save(sim.s1, n.batch, batch.length, starting, tuning, priors, dat, selec_dat, selec_XX,
     file = paste0(outpath, "sglmgp_res_", n.batch, ".RData"))


