#include <iostream>
#include <RcppArmadillo.h>
#include <math.h>
#include "utils.h"
#include "copulas.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

///////////////////////////////////////////////////////////////////////////////////
// Prediction functions for copula NNMPs with discrete marginals with CE
///////////////////////////////////////////////////////////////////////////////////

typedef double (*margfunc)(double, double, double);

typedef double (*ceqfunc)(double, double, double, double, double, double, int);

List predCeCop_simple(const arma::colvec& obs,
                      const int& cop_family,
                      const arma::cube& cop_param,
                      const arma::mat marg_param,
                      const int& nne,
                      const arma::mat& DD,
                      const arma::colvec zeta,
                      const arma::mat& ga,
                      const arma::colvec kasq,
                      const arma::mat oo,
                      const arma::mat grid_ne_idx,
                      const arma::mat grid_ne_dist,
                      const arma::colvec& probs,
                      const bool sam,
                      margfunc pMarg,
                      ceqfunc qMarg, 
                      double lb,
                      double ub,
                      double tol,
                      int maxit, 
                      const bool verbose,
                      const int nreport){

  int grid_size = grid_ne_dist.n_rows;
  int sample_size = zeta.n_rows;
  int probs_size = probs.n_rows;
  int label, iter_ne_label;

  double iter_ne, iter_zeta, iter_logit_mu, iter_ka,
         iter_cop_param, iter_marg_param1, iter_marg_param2,
         iter_vv, iter_uu, iter_ne_oo, y_star;

  arma::colvec label_vals = arma::linspace<arma::colvec>(0, nne - 1, nne);
  arma::colvec ipred_obs(sample_size);
  arma::colvec pred_obs_mean(grid_size);
  arma::mat pred_obs(probs_size, grid_size);
  arma::colvec ine_dist, iter_log_cutoff_rho, iter_cutoff_rho,
               iter_logit_cutoff, iter_cdf, iter_weight, iter_yy_oo;
  arma::colvec iter_cutoff(nne + 1);
  arma::rowvec iDD;
  arma::mat pred_sam(grid_size, sample_size);

  Rprintf("----------------------------------------\n");
  Rprintf("\t    Predicting\n");
  Rprintf("----------------------------------------\n");
  
  int icount = 0;
  for (int ii = 0; ii < grid_size; ++ii) {

    iDD = DD.row(ii);
    ine_dist = grid_ne_dist.row(ii).t();

    for (int iter = 0; iter < sample_size; ++iter) {
      
      // Compute weights and generate label
      iter_zeta = zeta(iter);
      iter_log_cutoff_rho = -ine_dist / iter_zeta;
      iter_cutoff_rho = arma::exp(iter_log_cutoff_rho - iter_log_cutoff_rho.max());
      iter_cutoff.fill(0);
      iter_cutoff.tail(nne) = cumsum_cpp(iter_cutoff_rho / sum(iter_cutoff_rho));
      iter_logit_cutoff = arma::log(iter_cutoff / (1.0 - iter_cutoff));
      iter_logit_mu = arma::as_scalar(iDD * ga.col(iter));
      iter_ka = sqrt(kasq(iter));
      iter_cdf = arma::normcdf(iter_logit_cutoff, iter_logit_mu, iter_ka);
      iter_weight = arma::diff(iter_cdf);
      label = arma::as_scalar(armaIntSample(label_vals, 1, 1, iter_weight));
      
      // Predic y0
      iter_ne_label = grid_ne_idx(ii, label);
      iter_ne = obs(iter_ne_label - 1);
      iter_yy_oo = oo.col(iter);
      iter_ne_oo = iter_yy_oo(iter_ne_label - 1);
      iter_cop_param = cop_param(ii, label, iter);
      iter_marg_param1 = marg_param(iter, 0);
      iter_marg_param2 = marg_param(iter, 1);
      iter_vv = pMarg(iter_ne - iter_ne_oo, iter_marg_param1, iter_marg_param2);
      iter_uu = rConCop(iter_vv, cop_family, iter_cop_param);
      y_star = qMarg(iter_uu, iter_marg_param1, iter_marg_param2, lb, ub, tol, maxit);
      ipred_obs(iter) = floor(y_star + 1.0);
      if (sam) {
        pred_sam(ii, iter) = ipred_obs(iter);
      }

    }

    pred_obs.col(ii) = arma::quantile(ipred_obs, probs);
    pred_obs_mean(ii) = mean(ipred_obs);

    icount++;
    if (verbose) {
      if (icount == nreport) {
        Rprintf("  Locations: %i/%i, %3.2f%%\n", ii+1, grid_size, 100.0*(ii+1)/grid_size);
        icount = 0;
      }
    }
    


  }

  if (sam) {
    return(List::create(Named("obs_mu") = pred_obs_mean, 
                        Named("obs_qq") = pred_obs,
                        Named("sam") = pred_sam));
  } else {
    return(List::create(Named("obs_mu") = pred_obs_mean, 
                        Named("obs_qq") = pred_obs));
  }

}


List predCeCop_covars(const arma::colvec& obs,
                      const int& cop_family,
                      const arma::cube& cop_param,
                      const arma::mat ref_marg_param1,
                      const arma::mat ref_marg_param2,
                      const arma::mat nonref_marg_param1,
                      const arma::mat nonref_marg_param2,
                      const int& nne,
                      const arma::mat& DD,
                      const arma::colvec zeta,
                      const arma::mat& ga,
                      const arma::colvec kasq,
                      const arma::mat& oo,
                      const arma::mat grid_ne_idx,
                      const arma::mat grid_ne_dist,
                      const arma::colvec& probs,
                      const bool sam,
                      margfunc pMarg,
                      ceqfunc qMarg,
                      double lb,
                      double ub,
                      double tol,
                      int maxit,
                      const bool verbose,
                      const int nreport){
  
  int grid_size = grid_ne_dist.n_rows;
  int sample_size = zeta.n_rows;
  int probs_size = probs.n_rows;
  int label, iter_ne_label;

  double iter_ne, iter_zeta, iter_logit_mu, iter_ka,
         iter_cop_param, iter_marg_param1, iter_marg_param2,
         iter_ne_marg_param1, iter_ne_marg_param2,
         iter_vv, iter_uu, iter_ne_oo, y_star;

  arma::colvec label_vals = arma::linspace<arma::colvec>(0, nne - 1, nne);
  arma::colvec ipred_obs(sample_size);
  arma::colvec pred_obs_mean(grid_size);
  arma::mat pred_obs(probs_size, grid_size);
  arma::colvec ine_dist, iter_log_cutoff_rho, iter_cutoff_rho,
               iter_logit_cutoff, iter_cdf, iter_weight, iter_yy_oo;
  arma::colvec iter_cutoff(nne + 1);
  arma::rowvec iDD;
  arma::mat pred_sam(grid_size, sample_size);

  Rprintf("----------------------------------------\n");
  Rprintf("\t  Predicting\n");
  Rprintf("----------------------------------------\n");
  
  int icount = 0;
  for (int ii = 0; ii < grid_size; ++ii) {

    iDD = DD.row(ii);
    ine_dist = grid_ne_dist.row(ii).t();

    for (int iter = 0; iter < sample_size; ++iter) {
      
      // Compute weights and generate label
      iter_zeta = zeta(iter);
      iter_log_cutoff_rho = -ine_dist / iter_zeta;
      iter_cutoff_rho = arma::exp(iter_log_cutoff_rho - iter_log_cutoff_rho.max());
      iter_cutoff.fill(0);
      iter_cutoff.tail(nne) = cumsum_cpp(iter_cutoff_rho / sum(iter_cutoff_rho));
      iter_logit_cutoff = arma::log(iter_cutoff / (1.0 - iter_cutoff));
      iter_logit_mu = arma::as_scalar(iDD * ga.col(iter));
      iter_ka = sqrt(kasq(iter));
      iter_cdf = arma::normcdf(iter_logit_cutoff, iter_logit_mu, iter_ka);
      iter_weight = arma::diff(iter_cdf);
      label = arma::as_scalar(armaIntSample(label_vals, 1, 1, iter_weight));
      
      // Predic y0
      iter_ne_label = grid_ne_idx(ii, label);
      iter_ne = obs(iter_ne_label - 1);
      iter_yy_oo = oo.col(iter);
      iter_ne_oo = iter_yy_oo(iter_ne_label - 1);
      iter_cop_param = cop_param(ii, label, iter);
      iter_marg_param1 = nonref_marg_param1(ii, iter);
      iter_marg_param2 = nonref_marg_param2(ii, iter);
      iter_ne_marg_param1 = ref_marg_param1(iter_ne_label - 1, iter);
      iter_ne_marg_param2 = ref_marg_param2(iter_ne_label - 1, iter);
      iter_vv = pMarg(iter_ne - iter_ne_oo, iter_ne_marg_param1, iter_ne_marg_param2);
      iter_uu = rConCop(iter_vv, cop_family, iter_cop_param);
      y_star = qMarg(iter_uu, iter_marg_param1, iter_marg_param2, lb, ub, tol, maxit);
      ipred_obs(iter) = floor(y_star + 1.0);
      if (sam) {
        pred_sam(ii, iter) = ipred_obs(iter);
      }

    }

    icount++;
    if (verbose) {
      if (icount == nreport) {
        Rprintf("  Locations: %i/%i, %3.2f%%\n", ii+1, grid_size, 100.0*(ii+1)/grid_size);
        icount = 0;
      }
    }    
    
    
    pred_obs.col(ii) = arma::quantile(ipred_obs, probs);
    pred_obs_mean(ii) = mean(ipred_obs);


  }

  if (sam) {
    return(List::create(Named("obs_mu") = pred_obs_mean, 
                        Named("obs_qq") = pred_obs,
                        Named("sam") = pred_sam));
  } else {
    return(List::create(Named("obs_mu") = pred_obs_mean, 
                        Named("obs_qq") = pred_obs));
  }


}


// [[Rcpp::export]]
List predCeCopNNMP_simple(const arma::colvec& obs,
                          const int& cop_family,
                          const arma::cube& cop_param,
                          const int& marg_family,
                          const arma::mat marg_param,
                          const int& nne,
                          const arma::mat& DD,
                          const arma::colvec zeta,
                          const arma::mat& ga,
                          const arma::colvec kasq,
                          const arma::mat oo,
                          const arma::mat grid_ne_idx,
                          const arma::mat grid_ne_dist,
                          const arma::colvec& probs,
                          double lb,
                          double ub,
                          double tol,
                          int maxit, 
                          int nreport,
                          const bool verbose = true,                            
                          const bool sam = false){

  List pred_out;
  
  // if (strcmp(marg_family, "poisson") == 0) 
  if (marg_family == 3)
  {
    pred_out = predCeCop_simple(obs, cop_family, cop_param, marg_param, nne, DD,
                                zeta, ga, kasq, oo, grid_ne_idx, grid_ne_dist, probs,
                                sam, pCEpois2, qCEpois, lb, ub, tol, maxit, 
                                verbose, nreport);
  } 
  else {
    stop("error: this family is currently not supported.");
  }

  return pred_out;

}


// [[Rcpp::export]]
List predCeCopNNMP_covars(const arma::colvec& obs,
                        const int& cop_family,
                        const arma::cube& cop_param,
                        const int& marg_family,
                        const arma::mat ref_marg_param1,
                        const arma::mat ref_marg_param2,
                        const arma::mat nonref_marg_param1,
                        const arma::mat nonref_marg_param2,
                        const int& nne,
                        const arma::mat& DD,
                        const arma::colvec zeta,
                        const arma::mat& ga,
                        const arma::colvec kasq,
                        const arma::mat& oo,
                        const arma::mat grid_ne_idx,
                        const arma::mat grid_ne_dist,
                        const arma::colvec& probs,
                        double lb,
                        double ub,
                        double tol,
                        int maxit,              
                        int nreport,
                        const bool verbose = true,                                 
                        const bool sam = false){

  List pred_out;
  
  // if (strcmp(marg_family, "negative binomial") == 0) 
  if (marg_family == 4)
  {
    pred_out = predCeCop_covars(obs, cop_family, cop_param, 
                              ref_marg_param1, ref_marg_param2,
                              nonref_marg_param1, nonref_marg_param2, nne, DD,
                              zeta, ga, kasq, oo, grid_ne_idx, grid_ne_dist, probs,
                              sam, pCEnbinom2, qCEnbinom, lb, ub, tol, maxit, 
                              verbose, nreport);
  } 
  else {
    stop("error: this family is currently not supported.");
  }

  return pred_out;

}

