copNNMP <- function(yy, XX = NULL, coords, ne_info, marg_family, cop_family, 
                    sp_func, priors, starting, tuning, mcmc_settings, 
                    weights_samples, verbose) {
  
  #--------------------------------------------------------------
  # Convert cop_family to integer for faster comparison in loops
  #--------------------------------------------------------------
  if (cop_family == "gaussian") {
    
    cop_fam <- 1
    
  } else if (cop_family == "gumbel") {
    
    cop_fam <- 2
    
  } else if (cop_family == "clayton") {
    
    cop_fam <- 3
    
  } 
  
  #----------------------------------------------------
  # Fit models
  #----------------------------------------------------
  if (marg_family == "poisson"){
    
    marg_fam <- 3
    copPoisNNMP(yy, XX, coords, ne_info, marg_fam, cop_fam,
                sp_func, priors, starting, tuning, mcmc_settings, 
                weights_samples, verbose)
    
  } else if (marg_family == "negative binomial"){
    
    marg_fam <- 4
    copNegbinNNMP(yy, XX, coords, ne_info, marg_fam, cop_fam,
                  sp_func, priors, starting, tuning, mcmc_settings, 
                  weights_samples, verbose)
    
  }

}

# Update phi of the copula NNMP
updateFullCopPhi <- function(phi, se_phi, u_phi, v_phi, 
                           ce_dat, ce_dat_ne, dat_rho, data_label, dist_mat, rho_mat,
                           mar_param, marg_family, cop_family, sp_func) {
  
  nn <- length(ce_dat)
  prop_log_phi <- rnorm(1, log(phi), se_phi)
  prop_phi <- exp(prop_log_phi)
  prop_rho_mat <- sp_func(dist_mat, prop_phi)
  prop_dat_rho <- c(0, as.numeric(t(sapply(2:nn, function(x) prop_rho_mat[x, data_label[x]]))))
  
  prop_loglik <-
    dgamma(1 / prop_phi, u_phi, v_phi, log = TRUE) + 2 * log(1 / prop_phi) +
    sum(dBiCopMar_cpp2(ce_dat[-1], ce_dat_ne[-1], cop_family, prop_dat_rho[-1], marg_family, mar_param, TRUE, TRUE)) + prop_log_phi
  cur_loglik <-
    dgamma(1 / phi, u_phi, v_phi, log = TRUE) + 2 * log(1 / phi) +
    sum(dBiCopMar_cpp2(ce_dat[-1], ce_dat_ne[-1], cop_family, dat_rho[-1], marg_family, mar_param, TRUE, TRUE)) + log(phi)
  
  diff_loglik <- prop_loglik - cur_loglik
  if (diff_loglik > log(runif(1))) {
    return(list(phi = prop_phi, rho_mat = prop_rho_mat, dat_rho = prop_dat_rho, accept = TRUE))
  } else {
    return(list(phi = phi, rho_mat = rho_mat, dat_rho = dat_rho, accept = FALSE))
  }
  
}

# Update zeta of the copula NNMP 
updateFullCopZeta <- function(zeta, se_zeta, nne, dist_mat, mu_t, kasq, data_label, u_zeta, v_zeta){
  
  nn <- nrow(dist_mat)
  prop_log_zeta <- rnorm(1, log(zeta), se_zeta)
  prop_zeta <- exp(prop_log_zeta)
  prop_weight_res <- logitGausWeight(nne, dist_mat, prop_zeta, mu_t, rep(sqrt(kasq), nn - 2), trun = FALSE)
  prop_weight_mat <- prop_weight_res$weights
  prop_cutoff <- prop_weight_res$cutoff
  weight_res <- logitGausWeight(nne, dist_mat, zeta, mu_t, rep(sqrt(kasq), nn - 2), trun = FALSE)
  weight_mat <- weight_res$weights
  cutoff <- weight_res$cutoff
  valid_idx <- !is.na(prop_weight_mat)
  if (any(prop_weight_mat[valid_idx]==0)) prop_weight_mat[prop_weight_mat==0] <- .Machine$double.eps
  valid_idx <- !is.na(weight_mat)
  if (any(weight_mat[valid_idx]==0)) weight_mat[weight_mat==0] <- .Machine$double.eps
  
  prop_loglik <-
    dgamma(1 / prop_zeta, u_zeta, v_zeta, log = TRUE) + 2 * log(1 / prop_zeta) +
    sum(log(sapply(3:nn, function(x) prop_weight_mat[x, data_label[x]]))) + prop_log_zeta
  cur_loglik <-
    dgamma(1 / zeta, u_zeta, v_zeta, log = TRUE) + 2 * log(1 / zeta) +
    sum(log(sapply(3:nn, function(x) weight_mat[x, data_label[x]]))) + log(zeta)
  
  diff_loglik <- prop_loglik - cur_loglik
  if (diff_loglik > log(runif(1))) {
    return(list(zeta = prop_zeta, weight_mat = prop_weight_mat, cutoff = prop_cutoff, accept = TRUE))
  } else {
    return(list(zeta = zeta, weight_mat = weight_mat, cutoff = cutoff, accept = FALSE))
  }
  
}

# Update gamma of the copula NNMP 
updateCopGa <- function(latent_t, DD, kasq, mu_ga, V_ga) {
  V_ga_inv <- chol2inv(chol(V_ga))
  V_1 <- chol2inv(chol(t(DD) %*% DD / kasq + V_ga_inv))
  mu_ga_1 <- V_1 %*% (t(DD) %*% latent_t / kasq + V_ga_inv %*% mu_ga)
  ga <- t(mvtnorm::rmvnorm(1, mu_ga_1, V_1))
  return(ga)
}

# Update kappa2 of the copula NNMP 
updateCopKasq <- function(latent_t, mu_t, u_kasq, v_kasq) {
  kk <- length(mu_t)
  uu <- u_kasq + kk / 2
  vv <- v_kasq + .5 * sum((latent_t - mu_t)^2)
  kasq <- 1 / rgamma(1, uu, vv)
  return(kasq)
}


