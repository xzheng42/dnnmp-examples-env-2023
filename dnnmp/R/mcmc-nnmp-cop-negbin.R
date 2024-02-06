# Gibbs sampler for negative binomial NNMP model 
copNegbinNNMP <- function(yy, XX, coords, ne_info, marg_family, cop_family, 
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
  
  if (!is.null(XX)) {
    
    # Check priors
    if (is.null(priors$scale_gamma)) {
      stop("error: prior for the scale parameter must be specified")
    }
    if (is.null(priors$scale_gamma)) {
      stop("error: prior for the scale parameter must be specified")
    }

    u_sz <- priors$scale_gamma[1]
    v_sz <- priors$scale_gamma[2]
    
    if (u_sz <= 0) stop("error: gamma prior for the scale parameter has positive hyperparameters")
    if (v_sz <= 0) stop("error: gamma prior for the scale parameter has positive hyperparameters")
    
    # Check tuning parameters
    if (is.null(tuning$regcoef)) {
      stop("error: random walk proposal distribution variance for regression coefficients must be specified")
    }
    if (is.null(tuning$scale)) {
      stop("error: random walk proposal distribution variance for log scale must be specified")
    }
    
    se_bb <- tuning$regcoef
    se_sz <- tuning$scale
    
    # Check starting
    if (is.null(starting$regcoef)) {
      if (identical(as.numeric(XX[, 1]), rep(1, length(yy)))) {
        bb <- c(log(mean(yy)), rep(0, ncol(XX) - 1))
      } else {
        bb <- rep(0, ncol(XX))
      }
    } else {
      bb <- starting$regcoef  
    }
    
    if (is.null(starting$scale)) {
      sz <- 1
    } else {
      sz <- starting$scale  
    }
    


    runtime <- system.time(
      mcmc_out <- mcmcCopNegbinNNMP_covars(yy, XX, coords, ne_info, 
                                           marg_family, cop_family, sp_func, 
                                           mcmc_settings, weights_samples, verbose,
                                          
                                           u_sz, v_sz, 
                                           u_phi, v_phi, u_zeta, v_zeta,
                                           mu_ga, V_ga, u_kasq, v_kasq,
                                         
                                           se_bb, se_sz, se_phi, se_zeta,
                                         
                                           bb, sz, phi, zeta, ga, kasq
      )
    )
  }
  else 
  {
    stop("error: for negative binomial marginals, the package only supports models with covariates.")
  }
  
  mcmc_out$runtime <- runtime
  
  mcmc_out
  
}


#-------------------------------------------------------------------
# MCMC for stationary Poisson NNMP
#-------------------------------------------------------------------
mcmcCopNegbinNNMP_covars <- function(yy, XX, coords, ne_info, 
                                     marg_family, cop_family, sp_func, 
                                     mcmc_settings, weights_samples, verbose,
                                     u_sz, v_sz, 
                                     u_phi, v_phi, u_zeta, v_zeta,
                                     mu_ga, V_ga, u_kasq,v_kasq,
                                     se_bb, se_sz, se_phi, se_zeta,
                                     bb, sz, phi, zeta, ga, kasq
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

  th <- exp(as.vector(XX %*% bb))
  th_ne <- c(0, th[ne_idx[-1]])
  
  # empty stuff to save samples
  oo_save <- array(NA, dim = c(nn, niter - nburn))
  bb_save <- array(NA, dim = c(ncol(XX), niter - nburn))
  sz_save <- array(NA, dim = niter - nburn)
  phi_save <- array(NA, dim = niter - nburn)
  zeta_save <- array(NA, dim = niter - nburn)
  ga_save <- array(NA, dim = c(ncol(DD), niter - nburn))
  kasq_save <- array(NA, dim = niter - nburn)
  if (weights_samples) {
    weight_save <- array(NA, dim = c(nn, nne, niter - nburn))
  }
  
  sz_acct_save <- array(0, dim = niter)
  phi_acct_save <- array(0, dim = niter)
  zeta_acct_save <- array(0, dim = niter)
  oo_acct_save <- array(0, dim = c(nn, niter))
  bb_acct_save <- array(0, dim = c(ncol(XX), niter))


  #--------------------------
  # MCMC
  #--------------------------
  if (verbose) {
    cat("----------------------------------------\n")
    cat("\t  Running MCMC\n");
    cat("----------------------------------------\n")
  }
  
  block_runtime <- 0
  
  for (iter in 1:niter) {
    
    start_time <- Sys.time()

    #--------------------------------------
    # Update regression coefficient beta
    #--------------------------------------
    bb_res <- updateCopNegbinCoef(bb, se_bb,
                                  yy, ce_yy, ce_yy_ne, yy_rho, XX, sz, th, th_ne,
                                  ne_idx, marg_family, cop_family)
    bb <- bb_res$bb
    bb_acct_save[, iter] <- bb_res$accept
    th <- bb_res$th
    th_ne <- bb_res$th_ne
    
    #--------------------------
    # Update scale
    #--------------------------
    sz_res <- updateCopNegbinSz(sz, se_sz, u_sz, v_sz,
                                yy, ce_yy, ce_yy_ne, yy_rho,
                                th, th_ne, marg_family, cop_family)
    sz <- sz_res$sz
    if (sz_res$accept) sz_acct_save[iter] <- 1
    
    #--------------------------
    # Update phi
    #--------------------------
    mar_param <- cbind(rep(sz, nn - 1), th[-1], rep(sz, nn - 1), th_ne[-1])
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
    oo_res <- updateCeCopAux(yy, oo, yy_rho, ne_idx, cop_family, marg_family, rep(sz, nn), th)
    oo <- oo_res$oo
    oo_acct_save[, iter] <- oo_res$oo_acct_save
    oo <- as.vector(oo)
    oo_ne <- c(0, oo[ne_idx[-1]])
    ce_yy <- yy - oo
    ce_yy_ne_mat <- t(sapply(1:nn, function(x) ce_yy[ne_idx_mat[x,]]))
    
    #--------------------------
    # Update labels
    #--------------------------
    yy_param1 <- rep(sz, nn)
    yy_param2 <- th
    yy_ne_param1 <- t(sapply(1:nn, function(x) yy_param1[ne_idx_mat[x,]]))
    yy_ne_param2 <- t(sapply(1:nn, function(x) yy_param2[ne_idx_mat[x,]]))
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
    th_ne <- c(0, th[ne_idx[-1]])
    
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
                   ". Percentage: ", specRound(iter / niter * 100, 2), "%\n"))
        ert <- (block_runtime / nreport) * (niter - iter)
        cat(paste0("Estimated remaining time: ", 
                   specRound(ert / 60, 2), " minutes \n"))
        block_runtime <- 0
        
        cat(paste0("Metropolis Hastings acceptance rates: \n"))
        cat(paste0("  Model parameter           Acceptance rate\n"))

        for (j in 1:ncol(XX)) {
          cat(paste0("  NegBin beta_", j, "             ", 
                     specRound(sum(bb_acct_save[j,1:iter]) / iter * 100, 2), "%\n"))
        }                   
        cat(paste0("  NegBin scale              ", 
                   specRound(sum(sz_acct_save[1:iter]) / iter * 100, 2), "%\n"))
        cat(paste0("  phi                       ", 
                   specRound(sum(phi_acct_save[1:iter]) / iter * 100, 2), "%\n"))
        cat(paste0("  zeta                      ", 
                   specRound(sum(zeta_acct_save[1:iter]) / iter * 100, 2), "%\n"))
      }
    }  
    
    #--------------------------
    # Save samples
    #--------------------------
    if (iter > nburn){
      oo_save[, iter - nburn] <- oo
      bb_save[, iter - nburn] <- bb
      sz_save[iter - nburn] <- sz
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
    regcoef = bb_save[, selc_index],
    scale = sz_save[selc_index],
    phi = phi_save[selc_index],
    zeta = zeta_save[selc_index],
    ga = ga_save[,selc_index],
    kasq = kasq_save[selc_index])
  
  if (weights_samples) {
    post_sams$weight <- weight_save[, , selc_index]
  }
  
  mh_acct_save <- list(regcoef_mh = bb_acct_save,
                       scale_mh = sz_acct_save,
                       oo_mh = oo_acct_save,
                       phi_mh = phi_acct_save,
                       zeta_mh = zeta_acct_save)
  
  list(post_sams = post_sams, mh_acct_save = mh_acct_save)
  
}


# Update regression parameters of the negative binomial marginal
updateCopNegbinCoef <- function(bb, se_bb, 
                                dat, ce_dat, ce_dat_ne, dat_rho, XX, sz, th, th_ne,
                                ne_idx, marg_family, cop_family) {
  
  nn <- length(ce_dat)
  accept <- array(0, dim = length(bb))

  for (j in seq_along(bb)) {

    prop_bb <- bb
    prop_bb[j] <- rnorm(1, bb[j], se_bb[j])
    prop_th <- as.vector(exp(XX %*% prop_bb))
    prop_th_ne <- c(0, prop_th[ne_idx[-1]])
    prop_mar_param <- cbind(rep(sz, nn - 1), prop_th[-1], rep(sz, nn - 1), prop_th_ne[-1])
    cur_mar_param <- cbind(rep(sz, nn - 1), th[-1], rep(sz, nn - 1), th_ne[-1])
    
    prop_loglik <-
      sum(dBiCopMar_cpp2(ce_dat[-1], ce_dat_ne[-1], cop_family, dat_rho[-1], marg_family, prop_mar_param, TRUE, TRUE)) +
      sum(dnbinom(dat, size = sz, mu = prop_th, log = TRUE))
    cur_loglik <-
      sum(dBiCopMar_cpp2(ce_dat[-1], ce_dat_ne[-1], cop_family, dat_rho[-1], marg_family, cur_mar_param, TRUE, TRUE)) +
      sum(dnbinom(dat, size = sz, mu = th, log = TRUE))
    
    diff_loglik <- prop_loglik - cur_loglik
    if (diff_loglik > log(runif(1))) {
      bb <- prop_bb
      th <- prop_th
      th_ne <- prop_th_ne
      accept[j] <- 1
    }

  }
  
  list(bb = bb, th = th, th_ne = th_ne, accept = accept)

}


# Update dispersion of the negative binomial marginals
updateCopNegbinSz <- function(sz, se_sz, u_sz, v_sz,
                              dat, ce_dat, ce_dat_ne, dat_rho, 
                              th, th_ne, marg_family, cop_family) {
  
  nn <- length(ce_dat)
  prop_log_sz <- rnorm(1, log(sz), se_sz)
  prop_sz <- exp(prop_log_sz)
  prop_mar_param <- cbind(rep(prop_sz, nn - 1), th[-1], rep(prop_sz, nn - 1), th_ne[-1])
  cur_mar_param <- cbind(rep(sz, nn - 1), th[-1], rep(sz, nn - 1), th_ne[-1])
  
  prop_loglik <-
    dgamma(prop_sz, u_sz, v_sz, log = TRUE) +
    sum(dBiCopMar_cpp2(ce_dat[-1], ce_dat_ne[-1], cop_family, dat_rho[-1], marg_family, prop_mar_param, TRUE, TRUE)) +
    sum(dnbinom(dat, size = prop_sz, mu = th, log = TRUE)) + prop_log_sz
  cur_loglik <-
    dgamma(sz, u_sz, v_sz, log = TRUE) +
    sum(dBiCopMar_cpp2(ce_dat[-1], ce_dat_ne[-1], cop_family, dat_rho[-1], marg_family, cur_mar_param, TRUE, TRUE)) +
    sum(dnbinom(dat, size = sz, mu = th, log = TRUE)) + log(sz)
  
  diff_loglik <- prop_loglik - cur_loglik
  if (diff_loglik > log(runif(1))) {
    return(list(sz = prop_sz, accept = TRUE))
  } else {
    return(list(sz = sz, accept = FALSE))
  }
  
}

