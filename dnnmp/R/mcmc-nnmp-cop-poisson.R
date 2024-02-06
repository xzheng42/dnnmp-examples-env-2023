# Gibbs sampler for Poisson NNMP model 
copPoisNNMP <- function(yy, XX, coords, ne_info, marg_family, cop_family, 
                        sp_func, priors, starting, tuning, mcmc_settings, 
                        weights_samples, verbose) {
  
  #-----------------------------------
  # Check priors except for marginals
  #-----------------------------------
  u_phi <- priors$phi_invgamma[1]
  v_phi <- priors$phi_invgamma[2]
  u_zeta <- priors$zeta_invgamma[1]
  v_zeta <- priors$zeta_invgamma[2]
  
  mu_ga <- priors$ga_gaus$mean_vec
  V_ga <- priors$ga_gaus$var_mat
  u_kasq <- priors$kasq_invgamma[1]
  v_kasq <- priors$kasq_invgamma[2]
  
  #----------------------------------------------
  # Check tuning parameters except for marginals
  #----------------------------------------------
  se_phi <- tuning$phi
  se_zeta <- tuning$zeta
  
  #----------------------------------------------
  # Check starting values except for marginals
  #----------------------------------------------
  phi <- starting$phi
  zeta <- starting$zeta
  ga <- starting$ga
  kasq <- starting$kasq
  
  #------------------------------------------------------------
  # Check priors, tuning and starting for marginals, and MCMC
  #------------------------------------------------------------
  
  if (is.null(XX)) {
    
    # Check prior
    if (is.null(priors$rate_gamma)) {
      stop("error: prior for the rate parameter must be specified")
    }    
    u_th <- priors$rate_gamma[1]
    v_th <- priors$rate_gamma[2]
    
    if (u_th <= 0) stop("error: gamma prior for the rate parameter has positive hyperparameters")
    if (v_th <= 0) stop("error: gamma prior for the rate parameter has positive hyperparameters")
        
    # Check tuning parameters
    se_th <- tuning$rate
    if (is.null(tuning$rate)) {
      stop("error: random walk proposal distribution variance for log rate must be specified")
    }    

    # Check starting
    if (is.null(starting$rate)) {
      th <- 1
    } else {
      th <- starting$rate
    }
    
    
    runtime <- system.time(
      mcmc_out <- mcmcCopPoisNNMP_simple(yy, coords, ne_info, 
                                         marg_family, cop_family, sp_func, 
                                         mcmc_settings, weights_samples, verbose,
                                          
                                         u_th, v_th, 
                                         u_phi, v_phi, u_zeta, v_zeta,
                                         mu_ga, V_ga, u_kasq, v_kasq,
                                         
                                         se_th, se_phi, se_zeta,
                                         
                                         th, phi, zeta, ga, kasq
      )
    )
  }
  else 
  {
    stop("error: for Poisson marginals, the package only supports models without covariates.")
  }
  
  mcmc_out$runtime <- runtime
  
  mcmc_out
  
}


#-------------------------------------------------------------------
# MCMC for stationary Poisson NNMP
#-------------------------------------------------------------------
mcmcCopPoisNNMP_simple <- function(yy, coords, ne_info, 
                                   marg_family, cop_family, sp_func, 
                                   mcmc_settings, weights_samples, verbose,
                                   u_th, v_th,
                                   u_phi, v_phi, u_zeta, v_zeta,
                                   mu_ga, V_ga, u_kasq,v_kasq,
                                   se_th, se_phi, se_zeta,
                                   th, phi, zeta, ga, kasq
                                   ) {
  
  
  #--------------------------
  # MCMC settings
  #--------------------------
  if (verbose) {
    
    if (is.null(mcmc_settings$n_report)) {
      nreport <- 1000
    } else {
      nreport <- mcmc_settings$n_report
    }
    
  }
  niter <- mcmc_settings$n_iter
  nburn <- mcmc_settings$n_burn
  nthin <- mcmc_settings$n_thin
  
  #--------------------------
  # Initialization
  #--------------------------
  nne <- ne_info$ne_size
  nn <- length(yy)
  dist_mat <- ne_info$ne_dist
  ne_idx_mat <- ne_info$ne_index
  yy_ne_mat <- ne_info$ne_obs
  DD <- as.matrix(cbind(1, coords))[-(1:2),]
  mu_t <- DD %*% ga

  weight_res <- logitGausWeight(nne, dist_mat, zeta, mu_t, rep(sqrt(kasq), nn - 2), FALSE)
  weight_mat <- weight_res$weights
  cutoff <- weight_res$cutoff
  valid_idx <- !is.na(weight_mat)
  if (any(weight_mat[valid_idx]==0)){
    weight_mat[weight_mat==0] <- .Machine$double.eps
  }

  yy_label <- array(NA, dim = nn)
  yy_label[1] <- 0
  yy_label[2] <- 1
  for (i in 3:nne) {
    yy_label[i] <- sample(1:(i-1), size = 1, prob = weight_mat[i, 1:(i-1)])
  }
  for (i in (nne+1):nn) {
    yy_label[i] <- sample(1:nne, size = 1, prob = weight_mat[i, ])
  }
  rho_mat <- sp_func(dist_mat, phi)
  yy_rho <- c(0, as.numeric(t(sapply(2:nn, function(x) rho_mat[x, yy_label[x]]))))
  ne_idx <- as.numeric(t(sapply(1:nn, function(x) ne_idx_mat[x, yy_label[x]])))
  yy_ne <- c(0, yy[ne_idx[-1]])
  
  oo <- runif(nn)
  oo_ne <- c(0, oo[ne_idx[-1]])
  ce_yy <- yy - oo
  ce_yy_ne <- yy_ne - oo_ne
  ce_yy_ne_mat <- t(sapply(1:nn, function(x) ce_yy[ne_idx_mat[x, ]]))
  
  # empty stuff to save samples
  oo_save <- array(NA, dim = c(nn, niter - nburn))
  th_save <- array(NA, dim = niter - nburn)
  phi_save <- array(NA, dim = niter - nburn)
  zeta_save <- array(NA, dim = niter - nburn)
  ga_save <- array(NA, dim = c(ncol(DD), niter - nburn))
  kasq_save <- array(NA, dim = niter - nburn)
  if (weights_samples) {
    weight_save <- array(NA, dim = c(nn, nne, niter - nburn))
  }
  
  th_acct_save <- array(0, dim = niter)
  phi_acct_save <- array(0, dim = niter)
  zeta_acct_save <- array(0, dim = niter)
  oo_acct_save <- array(0, dim = c(nn, niter))

  
  if (verbose) {
    cat("----------------------------------------\n")
    cat("\t  Running MCMC\n");
    cat("----------------------------------------\n")
  }

  #--------------------------
  # MCMC
  #--------------------------
  block_runtime <- 0
  
  for (iter in 1:niter) {
    
    start_time <- Sys.time()
    #--------------------------
    # Update th
    #--------------------------
    th_res <- updateCopPoisTh(th, se_th, u_th, v_th, 
                                yy, ce_yy, ce_yy_ne, yy_rho, 
                                marg_family, cop_family)
    th <- th_res$th
    if (th_res$accept) th_acct_save[iter] <- 1
    
    #--------------------------
    # Update phi
    #--------------------------
    mar_param <- cbind(rep(th, nn - 1), rep(th, nn - 1))
    phi_res <- updateFullCopPhi(phi, se_phi, u_phi, v_phi, 
                                ce_yy, ce_yy_ne, yy_rho, yy_label, dist_mat, rho_mat,
                                mar_param, marg_family, cop_family, sp_func)
    phi <- phi_res$phi
    rho_mat <- phi_res$rho_mat
    yy_rho <- phi_res$dat_rho
    if (phi_res$accept) phi_acct_save[iter] <- 1
    
    #--------------------------
    # Update oo
    #--------------------------
    oo_res <- updateCeCopAux(yy, oo, yy_rho, ne_idx, cop_family, marg_family, rep(th, nn), rep(th, nn))
    oo <- oo_res$oo
    oo_acct_save[, iter] <- oo_res$oo_acct_save
    oo <- as.vector(oo)
    oo_ne <- c(0, oo[ne_idx[-1]])
    ce_yy <- yy - oo
    ce_yy_ne_mat <- t(sapply(1:nn, function(x) ce_yy[ne_idx_mat[x,]]))
    
    #--------------------------
    # Update labels
    #--------------------------
    yy_param1 <- rep(th, nn)
    yy_param2 <- rep(th, nn)
    yy_ne_param1 <- matrix(rep(yy_param1, nne), ncol = nne)
    yy_ne_param2 <- matrix(rep(yy_param2, nne), ncol = nne)
    labels <- updateFullCopLabel(ce_yy, ce_yy_ne_mat, nne, weight_mat, rho_mat,
                                 cop_family, marg_family, 
                                 yy_param1, yy_param2, yy_ne_param1, yy_ne_param2,
                                 mu_t, rep(sqrt(kasq), length(mu_t)), cutoff)
    yy_label <- as.numeric(labels$data_label)
    latent_t <- labels$latent_t[-(1:2)]
    
    yy_rho <- c(0, as.numeric(t(sapply(2:nn, function(x) rho_mat[x, yy_label[x]]))))
    ne_idx <- as.numeric(t(sapply(1:nn, function(x) ne_idx_mat[x, yy_label[x]])))
    yy_ne <- c(0, yy[ne_idx[-1]])
    oo_ne <- c(0, oo[ne_idx[-1]])
    ce_yy_ne <- yy_ne - oo_ne
    
    #--------------------------
    # Update ga
    #--------------------------
    ga <- updateCopGa(latent_t, DD, kasq, mu_ga, V_ga)
    mu_t <- as.numeric(DD %*% ga)
    
    #--------------------------
    # Update kasq
    #--------------------------
    kasq <- updateCopKasq(latent_t, mu_t, u_kasq, v_kasq)
    
    #--------------------------
    # Update zeta and weights
    #--------------------------
    zeta_res <- updateFullCopZeta(zeta, se_zeta, nne, dist_mat, mu_t, kasq, yy_label, u_zeta, v_zeta)
    zeta <- zeta_res$zeta
    weight_mat <- zeta_res$weight_mat
    cutoff <- zeta_res$cutoff
    if (zeta_res$accept) zeta_acct_save[iter] <- 1
    
    end_time <- Sys.time()
    block_runtime <- block_runtime + as.numeric(end_time - start_time)

    #--------------------------
    # Print MCMC progress
    #--------------------------
    if (verbose) {
      if (iter %% nreport == 0) {

        cat(paste0("Iterations: ", iter, "/", niter, 
                   "  Percentage: ", specRound(iter / niter * 100, 2), "%\n"))

        ert <- (block_runtime / nreport) * (niter - iter)
        cat(paste0("Estimated remaining time: ", 
                   specRound(ert / 60, 2), " minutes \n"))
        block_runtime <- 0

        cat(paste0("Metropolis Hastings acceptance rates: \n"))
        cat(paste0("  Model parameter           Acceptance rate\n"))
        
        cat(paste0("  Poisson rate              ", 
                   specRound(sum(th_acct_save[1:iter]) / iter * 100, 2), "%\n"))
        cat(paste0("  phi                       ", 
                   specRound(sum(phi_acct_save[1:iter]) / iter * 100, 2), "%\n"))
        cat(paste0("  zeta                      ", 
                   specRound(sum(zeta_acct_save[1:iter]) / iter * 100, 2), "%\n"))

        cat("--------------------------------------------\n")

      }
    }  
    
    #--------------------------
    # Save samples
    #--------------------------
    if (iter > nburn){
      oo_save[, iter - nburn] <- oo
      th_save[iter - nburn] <- th
      phi_save[iter - nburn] <- phi
      zeta_save[iter - nburn] <- zeta
      ga_save[, iter - nburn] <- ga
      kasq_save[iter - nburn] <- kasq
      if (weights_samples) {
        weight_save[, , iter - nburn] <- weight_mat
      }
    }
    
  }
  
  #--------------------------
  # Thinning
  #--------------------------
  selc_index <- seq(1, niter - nburn, by = nthin)
  
  post_sams <- list(
    oo = oo_save[, selc_index],
    rate = th_save[selc_index],
    phi = phi_save[selc_index],
    zeta = zeta_save[selc_index],
    ga = ga_save[,selc_index],
    kasq = kasq_save[selc_index])
  
  if (weights_samples) {
    post_sams$weight <- weight_save[, , selc_index]
  }
  
  mh_acct_save <- list(rate_mh = th_acct_save,
                       oo_mh = oo_acct_save,
                       phi_mh = phi_acct_save,
                       zeta_mh = zeta_acct_save)
  
  list(post_sams = post_sams, mh_acct_save = mh_acct_save)
  
}


# Update shape parameter 1 of the copula NNMP with beta marginals
updateCopPoisTh <- function(th, se_th, u_th, v_th, 
                            yy, ce_yy, ce_yy_ne, dat_rho, 
                            marg_family, cop_family) {
  
  nn <- length(yy)
  prop_log_th <- rnorm(1, log(th), se_th)
  prop_th <- exp(prop_log_th)
  prop_mar_param <- cbind(rep(prop_th, nn - 1), rep(prop_th, nn - 1))
  cur_mar_param <- cbind(rep(th, nn - 1), rep(th, nn - 1))
  
  prop_loglik <- 
    dgamma(prop_th, u_th, v_th, log = TRUE) +
    sum(dBiCopMar_cpp2(ce_yy[-1], ce_yy_ne[-1], cop_family, dat_rho[-1], 
                       marg_family, prop_mar_param, TRUE, TRUE)) +
    sum(dpois(yy, prop_th, log = TRUE)) + prop_log_th
  cur_loglik <-  
    dgamma(th, u_th, v_th, log = TRUE) +
    sum(dBiCopMar_cpp2(ce_yy[-1], ce_yy_ne[-1], cop_family, dat_rho[-1], 
                       marg_family, cur_mar_param, TRUE, TRUE)) +
    sum(dpois(yy, th, log = TRUE)) + log(th)
  
  diff_loglik <- prop_loglik - cur_loglik
  if (diff_loglik > log(runif(1))) {
    return(list(th = prop_th, accept = TRUE))
  } else {
    return(list(th = th, accept = FALSE))
  }
  
}