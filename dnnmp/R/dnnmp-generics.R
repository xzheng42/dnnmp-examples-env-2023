#' @title Function for prediction at new locations
#'
#' @description This function produces posterior predictive samples for a set of new locations using a DNNMP.
#'
#' @param object an object of class dnnmp.
#' @param nonref_coords an \eqn{m \times 2} matrix of the corresponding spatial coordinates.
#' @param nonref_covars an \eqn{m \times p} matrix of covariates for \eqn{m} new locations.
#' @param probs a vector of probabilities with values in \eqn{[0,1]} to generate the corresponding quantiles.
#' @param predict_sam logical; if true, posterior predictive samples for each new location are returned.
#' @param nonref_par_labels a vector of numeric values of length \eqn{m}, indicating which partition
#'                          each new location belongs to.
#' @param compute_control This argument is used to specify parameters for
#'                        numerically inverting cdf of the continuous-extended marginals for prediction.
#'                        Currently, the argument is a list that comprises: lb: lower bound of the support;
#'                        ub: upper bound of the support; maxit: maximum number of iterations;
#'                        tol: the desired accuracy (convergence tolerance). The default parameters 
#'                        are lb = -1, ub = maximum of the data + 100, maxit = 1000, and tol = 1e-5.
#'                        Parameters (cop_maxit and cop_tol) for generating samples from a conditional 
#'                        copula (e.g., Gumbel) using numerical methods are fixed now. 
#'                        Currently, cop_maxit = 1000 and cop_tol = 1e-5.
#' @param verbose logical; if true, progress of the prediction is printed to the screen.
#' @param nreport If \code{verbose = TRUE}, \code{n_report} defines the interval to report
#'                the progress of prediction.
#' @param ... additional arguments. No additional arguments are supported currently.
#' 
#' @return
#' The return object is a list that comprises:
#' 
#' @author Xiaotian Zheng \email{xiaotianzheng42@gmail.com}
#' 
#' @references 
#' Zheng, X., Kottas, A., and Sansó, B. (2023),
#' "Bayesian geostatistical modeling for discrete-valued processes".
#' Environmetrics, 34(7), e2805.
#' DOI: \href{https://doi.org/10.1002/env.2805}{https://doi.org/10.1002/env.2805}.
#' 
#' Zheng, X., Kottas, A., and Sansó, B (2023). 
#' "Nearest-neighbor mixture models for non-Gaussian spatial processes".
#' Bayesian Anal. 18 (4) 1191 - 1222, December 2023.
#' DOI: \href{https://doi.org/10.1214/23-BA1405}{https://doi.org/10.1214/23-BA1405}.
#' 
#' @exportS3Method 
predict.dnnmp <- function(object, 
                          nonref_coords, 
                          nonref_covars = NULL, 
                          probs = c(0.025, 0.5, 0.975), 
                          predict_sam = FALSE,
                          nonref_par_labels = NULL, 
                          compute_control = NULL, 
                          verbose = TRUE, 
                          nreport = NULL,
                          ...) {
  
  if (missing(object)) {
    stop("error: need dnnmp object for prediction\n")
    }
  if (!inherits(object, "dnnmp")) {
    stop("error: require dnnmp object for prediction\n")
  }
  
  
  nne <- object$mod_spec$neighbor_size
  ref_coords <- object$ord_data$coords_ord
  ref_yy <- object$ord_data$yy_ord
  
  
  nonref_ne_info <- neighbor(nne, ref_coords, ref_yy, ref = FALSE, nonref_coords)
  nonref_ne_idx <- nonref_ne_info$ne_index
  nonref_ne_dist <- nonref_ne_info$ne_dist
  
  post_sams <- object$post_samples
  
  DD <- cbind(1, nonref_coords)
  
  marg_family <- object$mod_spec$marg_family
  
  cop_family <- object$mod_spec$cop_family
  if (cop_family == "gaussian") {
    cop_fam <- 1
  } else if (cop_family == "gumbel") {
    cop_fam <- 2
  } else if (cop_family == "clayton") {
    cop_fam <- 3
  }    
  
  if (!is.null(nonref_covars)) {
    XX <- as.matrix(nonref_covars)
  }
  
  if (is.null(compute_control)) {
  
    lb <- -1
    ub <- max(ref_yy) + 100
    maxit <- 1000
    tol <- .Machine$double.eps^0.25
    
  } else {
    
    lb <- compute_control$lb
    ub <- compute_control$ub
    maxit <- compute_control$maxit
    tol <- compute_control$tol
    
  }

  #----------------------------------------------------
  # verbose
  #----------------------------------------------------
  if (verbose) {
     
    if (is.null(nreport)) {
      
      nreport <- floor(nrow(nonref_coords) / 10)
      
    }
    
    nsams <- length(post_sams$zeta)

    cat("--------------------------------------------------------------\n")
    cat(" << Discrete nearest-neighbor mixture process model prediction >>\n")
    cat("--------------------------------------------------------------\n")
    cat(paste0("Reference set: ", length(ref_yy), " observations\n"))
    cat(paste0("Prediction over non-reference set: ", nrow(nonref_coords), " locations\n"))
    if (is.null(compute_control)) {
      cat("Inverse cdf sampling parameters were not specified; use default values.\n")
    }
    

    cat("----------------------------------------\n")
    cat("\tNNMP model settings\n");
    cat("----------------------------------------\n")
    cat(paste0("Using ", nne, " nearest neighbors\n"))
    cat(paste0("Marginal distribution family: ", tools::toTitleCase(marg_family), "\n"))
    cat(paste0("Mixture component family: ", tools::toTitleCase(cop_family), " copula\n"))
    cat(paste0("Number of MCMC samples: ", nsams, "\n"))
    
  }

  #----------------------------------------------------
  # generate posterior predictive samples
  #----------------------------------------------------
  if (marg_family == "poisson") {
    
    marg_fam <- 3
    
    if (is.null(object$orig_data$covars)) {

      cop_param_list <- sapply(seq_along(post_sams$phi), function(x) 
        object$mod_spec$sp_func(nonref_ne_dist, post_sams$phi[x]), simplify = FALSE)
      cop_param <- simplify2array(cop_param_list)

      runtime <- system.time(
        pred_out <- predCeCopNNMP_simple(ref_yy, cop_fam, cop_param, marg_fam,
                                         cbind(post_sams$rate, post_sams$rate),
                                         nne, DD, post_sams$zeta, post_sams$ga,
                                         post_sams$kasq, post_sams$oo,
                                         nonref_ne_idx, nonref_ne_dist,
                                         probs, lb, ub, tol, maxit,
                                         nreport, verbose, sam = predict_sam)
      )

    } else {

      stop("error: for Poisson marginals, the package only supports models without covariates.")

    }
    
  }
  else if (marg_family == "negative binomial") {
    
    marg_fam <- 4
    
    if (!is.null(object$orig_data$covars)) {

      cop_param_list <- sapply(seq_along(post_sams$phi), function(x) 
        object$mod_spec$sp_func(nonref_ne_dist, post_sams$phi[x]), simplify = FALSE)
      cop_param <- simplify2array(cop_param_list)
      
      ref_XX <- object$ord_data$XX_ord
      ref_mu <- apply(post_sams$regcoef, 2, function(x) exp(ref_XX %*% x))
      nonref_mu <- apply(post_sams$regcoef, 2, function(x) exp(XX %*% x))
      nn_ref <- nrow(ref_XX)
      nn_nonref <- nrow(XX)
     
      ref_scale <- t(matrix(rep(post_sams$scale, nn_ref), ncol = nn_ref))
      nonref_scale <- t(matrix(rep(post_sams$scale, nn_nonref), ncol = nn_nonref))
      
      runtime <- system.time(
        pred_out <- predCeCopNNMP_covars(ref_yy, cop_fam, cop_param, marg_fam,
                                         ref_scale, ref_mu, nonref_scale, nonref_mu,
                                         nne, DD, post_sams$zeta, post_sams$ga,
                                         post_sams$kasq, post_sams$oo,
                                         nonref_ne_idx, nonref_ne_dist,
                                         probs, lb, ub, tol, maxit,
                                         nreport, verbose, sam = predict_sam)
      )
     
    } else {
      
      stop("error: for negative binomial marginals, the package only supports models with covariates.")
 
    }
    
  }
  else {
      
      stop("error: this family is currently not supported.")
      
  }
  
  cat("----------------------------------------\n")
  cat(paste0("Prediction running time: ", specRound(runtime[3] / 60, 2), " minutes\n\n"))
  
  pred_out
  
}

#' @title Randomized quantile residuals for model checking.
#' 
#' @description This function produces posterior samples of the randomized quantile residuals for
#' reference locations. The function currently supports only the Poisson and negative binomial NNMPs.
#' 
#' @param object an object of class dnnmp.
#'
#' @author Xiaotian Zheng \email{xiaotianzheng42@gmail.com}
#' 
#' @references 
#' Zheng, X., Kottas, A., and Sansó, B. (2023),
#' "Bayesian geostatistical modeling for discrete-valued processes".
#' Environmetrics, 34(7), e2805.
#' DOI: \href{https://doi.org/10.1002/env.2805}{https://doi.org/10.1002/env.2805}.
#' 
#' Zheng, X., Kottas, A., and Sansó, B (2023). 
#' "Nearest-neighbor mixture models for non-Gaussian spatial processes". 
#' Bayesian Anal. 18 (4) 1191 - 1222, December 2023.
#' DOI: \href{https://doi.org/10.1214/23-BA1405}{https://doi.org/10.1214/23-BA1405}.
#' 
#' @export
#' 
rqr <- function(object) {
  
  if (missing(object)) {
    stop("error: need dnnmp object for summary\n")
  }
  
  if (!inherits(object, "dnnmp")) {
    stop("error: require dnnmp object for summary\n")
  }  

  nne <- object$mod_spec$neighbor_size
  coords_ord <- object$ord_data$coords_ord
  ref_yy_ord <- object$ord_data$yy_ord
  ref_XX_ord <- object$ord_data$XX_ord

  ne_info <- neighbor(nne, coords_ord, ref_yy_ord)
  ref_ne_obs <- ne_info$ne_obs
  ref_ne_idx <- ne_info$ne_index
  ref_ne_dist <- ne_info$ne_dist
  
  marg_family <- object$mod_spec$marg_family
  
  cop_family <- object$mod_spec$cop_family
  if (!is.null(cop_family)) {
    if (cop_family == "gaussian") {
      cop_fam <- 1
    } else if (cop_family == "gumbel") {
      cop_fam <- 2
    } else if (cop_family == "clayton") {
      cop_fam <- 3
    } 
  }
  
  post_sams <- object$post_samples
  cop_param_list <- sapply(seq_along(post_sams$phi), function(x) 
    object$mod_spec$sp_func(ref_ne_dist, post_sams$phi[x]), simplify = FALSE)
  cop_param <- simplify2array(cop_param_list)
  
  if (marg_family == "poisson") {
      concump <- pCopPoNNMP_ref_simple(ref_yy_ord, cop_fam,
                                       cop_param, ref_ne_obs, ref_ne_idx,
                                       post_sams$oo, post_sams$rate, post_sams$weight)
  } 
  else if (marg_family == "negative binomial") {
      concump <- pCopNbNNMP_ref_covar(ref_yy_ord, cop_fam, cop_param, 
                                      ref_XX_ord, ref_ne_obs, ref_ne_idx,
                                      post_sams$oo, post_sams$regcoef,
                                      post_sams$scale, post_sams$weight)
  }
  else {
    stop("error: currently only the Poisson NNMP wihtout covariates and 
          the negatvie binomial with covariates are supported.")
  }
  
  qnorm(concump)
  
}

#' @title Obtain posterior estimates of the DNNMP model parameters
#'
#' @description This function calculates point and interval estimates of unknown parameters
#' using their posterior samples
#'
#' @param object an object of class dnnmp.
#' @param pt a quoted keyword that specifies point estimate.
#'           Supported keywords are \code{"mean"} and \code{"median"}.
#' @param ci a numeric vector of lower and upper probabilities with values in \eqn{[0,1]} 
#'           to generate a credible interval for the unknown parameter(s).
#' @param digit an integer indicating the number of decimal places to be used.
#' @param ... additional arguments. No additional arguments are supported currently.
#'
#' @return A vector of characters. 
#'         Each element of the vector is in the form: point estimate (lower percentile, upper percentile).
#'         
#' @export
#' 
#' @exportS3Method 
summary.dnnmp <- function(object, pt = "mean", ci = c(0.025, 0.975), digit = 2, ...) {
  
  if (missing(object)) {
    stop("error: need dnnmp object for summary\n")
  }
  
  if (!inherits(object, "dnnmp")) {
    stop("error: require dnnmp object for summary\n")
  }
  
  post_sams <- object$post_samples
  marg_family <- object$mod_spec$marg_family
  
  if (marg_family == "poisson") {
    
    if (is.null(object$orig_data$covars)) {
      
      rate_est <- summarize(post_sams$rate, pt = pt, ci = ci, k = digit)
      phi_est <- summarize(post_sams$phi, pt = pt, ci = ci, k = digit)
      zeta_est <- summarize(post_sams$zeta, pt = pt, ci = ci, k = digit)
      ga_est <- summarize(post_sams$ga, pt = pt, ci = ci, k = digit)
      kasq_est <- summarize(post_sams$kasq, pt = pt, ci = ci, k = digit)

      tbl_est <- array(NA, dim = c(7, 1))
      tbl_est[,1] <- c(rate_est, phi_est, zeta_est, ga_est, kasq_est)
      rownames(tbl_est) <- c("rate", "phi", "zeta", "gamma_0", "gamma_1", "gamma_2", "kappasq")
      colnames(tbl_est) <- "Poisson NNMP"
      
    } else {

      stop("error: for Poisson marginals, the package only supports models without covariates.")

    }
  }
  else if (marg_family == "negative binomial") {

    if (!is.null(object$orig_data$covars)) {

      bb_est <- summarize(post_sams$regcoef, pt = pt, ci = ci, k = digit)
      sz_est <- summarize(post_sams$scale, pt = pt, ci = ci, k = digit)
      phi_est <- summarize(post_sams$phi, pt = pt, ci = ci, k = digit)
      zeta_est <- summarize(post_sams$zeta, pt = pt, ci = ci, k = digit)
      ga_est <- summarize(post_sams$ga, pt = pt, ci = ci, k = digit)
      kasq_est <- summarize(post_sams$kasq, pt = pt, ci = ci, k = digit)

      tbl_est <- array(NA, dim = c(length(bb_est) + 7, 1))
      tbl_est[,1] <- c(bb_est, sz_est, phi_est, zeta_est, ga_est, kasq_est)
      rownames(tbl_est) <- c(paste0("regcoef_", seq_along(bb_est)), "scale", 
                             "phi", "zeta", "gamma_0", "gamma_1", "gamma_2", "kappasq")
      colnames(tbl_est) <- "Negative binomial NNMP"

    } else {

      stop("error: for negative binomial marginals, the package only supports models with covariates.")

    }

  }
  else {
      
      stop("error: this family is currently not supported.")
      
  }  
  
  tbl_est
  
}