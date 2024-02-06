#' @title Function for fitting a discrete nearest-neighbor mixture process model
#' 
#' @description This function fits a discrete nearest-neighbor mixture process (DNNMP) model with Poisson or 
#' negative-binomial marginal distributions via Markov chain Monte Carlo (MCMC) simulations.
#' 
#' @param response a vector of complete geostatistical data of length \eqn{n}.
#' @param covars an \eqn{n \times p} matrix of covariates. If an intercept is desired, the first 
#'               column of the matrix should be a vector of ones.
#' @param coords an \eqn{n \times 2} matrix of the observation coordinates, e.g., longitude and latitude.
#' @param neighbor_size neighborhood size (number of nearest neighbors).
#' @param marg_family a quoted keyword that specifies the family of distributions for the marginal.
#'                    Supported keywords are \code{"poisson"} and \code{"negative binomial"}.
#' @param cop_family a quoted keyword that specifies the family of copulas for the bivariate distributions
#'                   that define the first-order conditionals. 
#'                   Supported keywords are \code{"gaussian"}, \code{"gumbel"}, and \code{"clayton"}.
#' @param priors a list of priors. Each element of the list corresponds to a combination of an unknown parameter 
#'               and its prior distribution; for example, "\code{phi_invgamma}" and "\code{zeta_invgamma}".
#' @param starting a list of starting values. Each element of the list corresponds to an unknown parameter.
#' @param tuning a list of tuning parameters. Each element of the list corresponds to an unknown parameter.
#'               The value portion of each element defines the variance of the random walk Metropolis sampler 
#'               proposal distribution for the unknown parameter.
#' @param ord an index vector of length \eqn{n} to order the data. The vector defines an ordering based on which
#'            nearest neighbors are searched. If the argument is not specified, the default ordering is random ordering.
#' @param mcmc_settings a list of MCMC simulation parameters.
#'        \code{n_iter}: number of iterations; \code{n_burn}: number of burn-in samples;
#'        \code{n_thin}: thinning degree. If \code{verbose = TRUE},
#'        \code{n_report} defines the interval to report MCMC progress and Metropolis acceptance rates.
#' @param weights_samples logical; if true, samples of the mixture weights are returned. 
#' @param verbose logical; if true, model specification, MCMC progress, and Metropolis acceptance rates
#'                are printed to the screen.
#' @param neighbor_info logical; if true, a list that comprises each location's neighborhood information is returned.
#' 
#' @return
#' An object of class \code{dnnmp}. The return object is a list that comprises:
#' \tabular{ll}{
#' 
#' \code{post_samples} \tab a list of posterior samples. Each element of the list corresponds to an unknown parameter. \cr
#' \tab \cr
#' \code{orig_data} \tab a list that consists of input data 
#'         (\code{response}, \code{covars}, and \code{coords}. \cr
#' \tab \cr
#' \code{ord_data} \tab a list that consists of ordered input data
#'         (\code{ord:} index vector for ordering; 
#'          \code{yy_ord}: ordered responses;
#'          \code{XX_ord}: ordered covariates (if \code{covars = NULL}, \code{XX_ord = NULL});
#'          \code{coords_ord}: ordered observation coordinates). \cr
#' \tab \cr
#' \code{mod_spec}: \tab a list that consists of model specifications 
#'         (\code{neigbhor_size}: neighborhood size; 
#'          \code{marg_family}: family of marginal distributions;
#'          \code{cop_family}: copula family for the bivariate distributions that define the first-order conditionals;
#'          \code{sp_func}: function that introduces spatial dependence;
#'          \code{priors}: priors for unknown parameters). \cr
#' \tab \cr
#' \code{mcmc_params}: \tab a list that consists of unknown parameters' starting values and Metropolis sampler proposal distribution variances
#'         used in the MCMC simulation (
#'         \code{starting}: starting values;
#'         \code{tuning}: Metropolis sampler proposal distribution variances). \cr
#' \tab \cr
#' \code{mh_data}: \tab a list that consists of each unknown parameter's Metropolis 
#'                      sampler acceptance and rejection indicators.\cr
#' \tab \cr
#' \code{runtime}: \tab running time for MCMC simulation calculated using \code{system.time}.
#'         Running time for searching nearest neighbors is not included. 
#'}
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
dnnmp <- function(response,
                  covars = NULL,
                  coords,
                  neighbor_size = 10, 
                  marg_family, 
                  cop_family,
                  priors = NULL,
                  starting = NULL,
                  tuning = NULL,
                  ord = NULL,
                  mcmc_settings,
                  weights_samples = FALSE,
                  verbose = TRUE,
                  neighbor_info = FALSE
                  ) {

  #----------------------------------------------------
  # response
  #----------------------------------------------------
  if (missing(response)) {
    stop("error: response data must be input")
  }
  nn <- length(response)
	
  
  #----------------------------------------------------
  # nearest neighbor size
  #----------------------------------------------------
  nne <- neighbor_size
  if (nne %% 1 != 0) {
    stop("error: neighbor size must be an integer")
  }
  
  if (nne < 2) {
    stop("error: neighbor size must be greater than one")
  }

  
  #----------------------------------------------------
  # coords
  #----------------------------------------------------
  if (missing(coords)) {
    stop("error: data coordinate information must be input")
  }
  
  if (ncol(coords) != 2) {
    stop("error: the matrix must have two columns")
  }

  if (nrow(coords) != nn) {
    stop("error: the number of coordinates is not equal to the number of responses")
  }

  
  #----------------------------------------------------
  # marginal family
  #----------------------------------------------------
  if (missing(marg_family)) {
    stop("error: marginal family must be specified")
  }
  
  marg_family <- tolower(marg_family)
  marg_family_names <- c("poisson", "negative binomial")
  
  if (!marg_family %in% marg_family_names) {
    stop("error: specified family, '", cop_family, "' for the marginal distribution is not a valid option;
         available families are ", paste(marg_family_names, collapse = ", ", sep = "") ,".")
  }  
  
  
  #----------------------------------------------------
  # copula family
  #----------------------------------------------------
  if (missing(cop_family)) {
    stop("error: copula family must be specified")
  }
  cop_family <- tolower(cop_family)
  cop_family_names <- c("gaussian", "gumbel", "clayton")
  
  if (!cop_family %in% cop_family_names) {
    stop("error: specified family, '", cop_family, "' for the copula is not a valid option;
         available families are ", paste(cop_family_names, collapse = ", ", sep = "") ,".")
  }  
  
  
  #----------------------------------------------------
  # DAG
  #----------------------------------------------------
  yy <- as.numeric(response)
  coords <- as.matrix(coords)
  
  orig_data <- list(response = yy, 
                    covars = covars, 
                    coords = coords)
  
  if (is.null(ord)) {
    ord <- sample(1:nn, replace = FALSE)
  }
  
  if (is.null(covars)) {
    yy_ord <- yy[ord]
    XX_ord <- NULL
    coords_ord <- coords[ord, ]
    ord_data <- list(ord = ord,
                     yy_ord = yy_ord,
                     XX_ord = XX_ord,
                     coords_ord = coords_ord)
  } else {
    yy_ord <- yy[ord]
    XX_ord <- covars[ord, ]
    coords_ord <- coords[ord, ]
    ord_data <- list(ord = ord,
                     yy_ord = yy_ord,
                     XX_ord = XX_ord,
                     coords_ord = coords_ord) 
  }
  
  ne_info <- neighbor(nne, coords_ord, yy_ord)
  ne_info$ne_size <- nne
  
  
  #----------------------------------------------------
  # function that introduces spatial dependence
  #----------------------------------------------------
  if (cop_family == "gaussian") {
    
    sp_func <- function(dist, range_param) return(exp(-dist / range_param))
    
  } else if (cop_family == "gumbel") {
    
    sp_func <- function(dist, range_param) 1 / (1 - exp(-dist / range_param))
    
  } else if (cop_family == "clayton") {
    
    sp_func <- function(dist, range_param) {
      kk <- exp(-dist / range_param)
      2 * kk / (1 - kk)
    }
    
  } 
    
  
  #----------------------------------------------------
  # Check some priors
  #----------------------------------------------------
  if (is.null(priors$phi_invgamma)) {
    priors$phi_invgamma <- c(3, 1/3)
  }
  
  if (is.null(priors$zeta_invgamma)) {
    priors$zeta_invgamma <- c(3, 0.2)
  }
  
  if (is.null(priors$ga_gaus)) {
    priors$ga_gaus <- list(mean_vec = as.matrix(c(-1.5, 0, 0)),
                           var_mat = 2 * diag(3))
  } 
  
  if (is.null(priors$kasq_invgamma)) {
    priors$kasq_invgamma <- c(3, 1)
  }
  
  
  #----------------------------------------------------
  # Check some starting values
  #----------------------------------------------------
  if (is.null(starting$ga)) {
    starting$ga <- matrix(c(-1.5, 0, 0))
  }
  
  if (is.null(starting$kasq)) {
    starting$kasq <- 1
  }
  
  if (is.null(starting$phi)) {
    starting$phi <- 1/6
  }
  
  if (is.null(starting$zeta)) {
    starting$zeta <- 0.1
  }
  
  #----------------------------------------------
  # Check some tuning parameters 
  #----------------------------------------------
  if (is.null(tuning$phi)) {
    stop("error: random walk proposal distribution variance for log phi must be specified")
  }
  
  if (is.null(tuning$zeta)) {
    stop("error: random walk proposal distribution variance for log zeta must be specified")
  }
  
  #----------------------------------------------------
  # verbose
  #----------------------------------------------------
  if (verbose) {
    cat("--------------------------------------------------------------\n")
    cat(" << Fitting a discrete nearest-neighbor mixture process model >>\n")
    cat("--------------------------------------------------------------\n")   
    cat("----------------------------------------\n")
    cat("\t  NNMP settings\n");
    cat("----------------------------------------\n")

    cat(paste0("Number of observations: ", nn, "\n"))
    if (is.null(ord)) {
      cat(paste0("Using ", nne, " nearest neighbors with random ordering\n"))
    } else {
      cat(paste0("Using ", nne, " nearest neighbors with user-specified ordering\n"))
    }
    
    cat(paste0("Marginal distribution family: ", tools::toTitleCase(marg_family), "\n"))
    
    cat(paste0("Mixture component family: ", tools::toTitleCase(cop_family), " copula\n"))
    
    cat("----------------------------------------\n")
    cat("\t  MCMC settings\n");
    cat("----------------------------------------\n")
    cat(paste0("Number of iterations: ", mcmc_settings$n_iter, "\n"))
    cat(paste0("Number of burn-in samples: ", mcmc_settings$n_burn, "\n"))
    cat(paste0("Thinning degree: ", mcmc_settings$n_thin, "\n"))
    
  }

  
  #----------------------------------------------------
  # MCMC
  #----------------------------------------------------
  if (cop_family %in% c("gaussian", "gumbel", "clayton")) {
    
    mcmc_out <- copNNMP(yy_ord, XX_ord, coords_ord, ne_info, marg_family, cop_family, 
                        sp_func, priors, starting, tuning, mcmc_settings, weights_samples, verbose)
    
  } else {
    
    # stop(paste0("error: the current package implements only ", tools::toTitleCase(marg_family), 
    #             " NNMP with Gaussian, Gumbel, or Clayton copulas."))
    stop("error: the current package implements only models with Gaussian, Gumbel, and Clayton copulas")
    
  }
  
  post_sams <- mcmc_out$post_sams
  runtime <- mcmc_out$runtime
  mh_data <- mcmc_out$mh_acct_save
  
  cat(paste0("MCMC running time: ", specRound(runtime[3] / 60, 2), " minutes\n"))
  
  mod_spec <- list(neighbor_size = neighbor_size, 
                   marg_family = marg_family,
                   cop_family = cop_family,
                   sp_func = sp_func,
                   priors = priors)
  
  mcmc_params <- list(starting = starting, tuning = tuning)
    
  dnnmp_out <- list(
    post_samples = post_sams,
    orig_data = orig_data,
    ord_data = ord_data,
    mod_spec = mod_spec,
    mcmc_params = mcmc_params,
    mcmc_settings = mcmc_settings,
    mh_data = mh_data,
    runtime = runtime
  )
  
  
  if (neighbor_info) {
    dnnmp_out$neighbor_info = ne_info
  }
  
  
  class(dnnmp_out) <- "dnnmp"
  
  dnnmp_out
      
  
}
