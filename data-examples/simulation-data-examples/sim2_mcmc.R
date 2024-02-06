################################################################################
### SIMULATION DATA EXAMPLE 2 - NBNNMP
################################################################################

rm(list = ls())

outpath <- "/path/to/output/"

library(SpatialExtremes)
library(ggplot2)
library(spBayes)
library(dnnmp)


#-------------------------------------------------------------------------------
# Simulate data
#-------------------------------------------------------------------------------

# Generate latent Gaussian random field on a grid
x_seq <- seq(0, 1, length = 120)
y_seq <- seq(0, 1, length = 120)
grid_loc <- as.matrix(expand.grid(x_seq, y_seq))
nn <- nrow(grid_loc)

set.seed(42)
gp <- t(rgp(1, grid_loc, cov.mod = "powexp", nugget =  0, sill = 0.2, range = 1/12, smooth = 1))
latent <- cbind(grid_loc, gp)
colnames(latent) <- c("x", "y", "z")
zz <- latent[,3]

# Set up the mean function (log link for the poisson mean)
set.seed(42)
XX <- cbind(1, grid_loc)
bb <- as.matrix((c(1.5, 1, 2)))

# Generate data
yy <- rpois(nn, exp(XX %*% bb + zz))
dat <- cbind(grid_loc, yy)
colnames(dat) <- c("x", "y", "z")

# Plot the true field
myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))

plot_save <- ggplot(as.data.frame(dat), aes(x = x, y = y, fill = z)) +
  geom_raster(interpolate = T) + theme_bw() +
  scale_fill_gradientn(colours = myPalette(100), limits = range(dat[,3])) +
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

setEPS()
postscript(paste0(outpath, "pois_glmm", ".eps"), width = 7, height = 6, bg = "transparent")
print(plot_save)
dev.off()

#-------------------------------------------------------------------------------
# Fit negative binomial NNMP
#-------------------------------------------------------------------------------
set.seed(42)

# Randomly select part of the gridded data as observations
obs_num <- 1000
obs_idx <- sample(1:nrow(grid_loc), obs_num, FALSE)
obs_dat <- dat[obs_idx,]
obs_XX <- XX[obs_idx,]

# Partition the observations into training and testing sets
# Use the whole training set as reference set for NNMP
ref_num <- 800
test_num <- obs_num - ref_num

ref_idx <- sample(1:obs_num, ref_num, FALSE)
ref_dat <- obs_dat[ref_idx,]
ref_XX <- obs_XX[ref_idx,]

test_dat <- obs_dat[-ref_idx,]
test_XX <- XX[-ref_idx, ]

# Random ordering
ord <- sample(1:ref_num, replace = FALSE)

# MCMC
priors <- list("scale_gamma" = c(1, 1),
               "phi_invgamma" = c(3, 1),
               "zeta_invgamma" = c(3, 1))

starting <- list("phi" = 0.01, "zeta" = 0.01, 
                 "ga" = rep(0, 3), "kasq" = 1, "scale" = 1, 
                 "regcoef" = as.matrix(c(0.1, 0.1, 0.1)))

tuning <- list("regcoef" = c(0.15, 0.2, 0.2), "scale" = 0.15, "phi" = 0.15, "zeta" = 0.8)

mcmc_param <- list(n_iter = 20000, n_burn = 4000, n_thin = 4, n_report = 2000)

nne <- 10

set.seed(42)
nnmp_out <- dnnmp(response = ref_dat[,3],
                  covars = ref_XX,
                  coords = ref_dat[,1:2],
                  neighbor_size = nne,
                  marg_family = "negative binomial",
                  cop_family = "gaussian",
                  priors = priors,
                  starting = starting,
                  tuning = tuning,
                  ord = ord,
                  mcmc_settings = mcmc_param, 
                  weights_samples = TRUE,
                  verbose = TRUE)

save(nnmp_out, dat, test_dat, test_XX, file = paste0(outpath, "nbnnmp_res.RData"))

#-------------------------------------------------------------------------------
# SGLMM
#-------------------------------------------------------------------------------
beta_starting <-  coefficients(glm(ref_dat[,3] ~ ref_XX[,-1], family="poisson"))
n.batch <- 400
batch.length <- 100
n.samples <- n.batch * batch.length
starting <- list("beta"=beta_starting, "phi"=3/0.5, "sigma.sq"=1, "w"=0)
tuning <- list("beta"=c(0.05,0.1,0.1), "phi"=0.2,  "sigma.sq"=0.5, "w"=0.5)
priors <- list("beta.Flat", "phi.Unif"=c(3/1, 3/0.1), "sigma.sq.IG"=c(2,1))

# SGLMM-GP
set.seed(42)
sim.s1 <- spGLM(formula = ref_dat[,3] ~ ref_XX[, -1],
                family = "poisson",
                coords = ref_dat[,1:2],
                starting = starting,
                tuning = tuning,
                priors = priors,
                cov.model = "exponential",
                n.samples = n.samples,
                amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
                n.report = 100)

save(sim.s1, n.batch, batch.length, starting, tuning, priors, dat, file = paste0(outpath, "sglmgp_res.RData"))


# SGLMM-GPP
set.seed(42)
sim.s2 <- spGLM(formula = ref_dat[,3] ~ ref_XX[, -1], 
                family = "poisson",
                coords = ref_dat[,1:2],
                knots = c(10, 10),
                starting = starting, 
                tuning = tuning, 
                priors = priors, 
                cov.model = "exponential",
                n.samples = n.samples, 
                amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
                n.report = 10)

save(sim.s2, n.batch, batch.length, starting, tuning, priors, dat, file = paste0(outpath, "sglmgpp_res.RData"))


