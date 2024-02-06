#include <iostream>
#include <RcppArmadillo.h>
#include <math.h>
#include "utils.h"
#include "copulas.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//------------------------------------------------------------------------------
// Compute conditional probabilities for reference locations
//------------------------------------------------------------------------------

// Poisson NNMP
// [[Rcpp::export]]
arma::mat pCopPoNNMP_ref_simple(const arma::colvec& obs,
                                const int& cop_family,
                                const arma::cube& cop_param,
                                const arma::mat& ref_ne,
                                const arma::umat& ref_ne_idx,
                                const arma::mat& oo,
                                const arma::colvec& th,
                                const arma::cube& weight) {
  
  int data_size = obs.n_rows;
  int sample_size = th.n_rows;
  int ne_size = ref_ne.n_cols;
  arma::colvec weight_ij, oo_j, oo_ij_nes, mar_param1;
  arma::mat weight_j;
  double cop_param_j, cop_param_ikj, th_j, oo_ij, oo_ij_ne, uu_ij, vv_ij;
  
  arma::mat concump(data_size, sample_size);
  arma::colvec logcdf_ij;
  arma::uvec ine_idx;
  
  // the first obs
  for (int j = 0; j < sample_size; ++j) {
    th_j = th(j);
    concump(0, j) = pCEpois(obs(0) - oo(0,j), th_j);
  }
  
  // the second obs
  for (int j = 0; j < sample_size; ++j) {
    cop_param_j = cop_param(1, 0, j);
    th_j = th(j);
    oo_ij = oo(1, j);
    oo_ij_ne = oo(0, j);
    uu_ij = pCEpois(obs(1) - oo_ij, th_j);
    vv_ij = pCEpois(obs(0) - oo_ij_ne, th_j);
    concump(1, j) = pConCop(uu_ij, vv_ij, cop_family, cop_param_j);
  }
  
  // obs 3 - ne_size
  for (int i = 2; i < ne_size; ++i) {
    
    logcdf_ij.set_size(i);
    ine_idx = arma::trans(ref_ne_idx(i, arma::span(0, i - 1))) - 1;
    
    for (int j = 0; j < sample_size; ++j) {
      
      oo_j = oo.col(j);
      oo_ij_nes = oo_j(ine_idx);
      oo_ij = oo_j(i);
      th_j = th(j);
      
      for (int k = 0; k < i; ++k) {
        cop_param_ikj = cop_param(i, k, j);
        oo_ij_ne = oo_ij_nes(k);
        uu_ij = pCEpois(obs(i) - oo_ij, th_j);
        vv_ij = pCEpois(ref_ne(i, k) - oo_ij_ne, th_j);
        logcdf_ij(k) = pConCop(uu_ij, vv_ij, cop_family, cop_param_ikj, true);
      }
      
      weight_j = weight.slice(j);
      weight_ij = arma::trans(weight_j(i, arma::span(0, i - 1)));
      concump(i,j) = sum(arma::exp(arma::log(weight_ij) + logcdf_ij));
      
    }
    
  }
  
  // Remaining obs
  logcdf_ij.set_size(ne_size);
  for (int i = ne_size; i < data_size; ++i) {
    
    ine_idx = ref_ne_idx.row(i).t() - 1;
    
    for (int j = 0; j < sample_size; ++j) {
      
      oo_j = oo.col(j);
      oo_ij_nes = oo_j(ine_idx);
      oo_ij = oo_j(i);
      th_j = th(j);
      
      for (int k = 0; k < ne_size; ++k) {
        cop_param_ikj = cop_param(i, k, j);
        oo_ij_ne = oo_ij_nes(k);
        uu_ij = pCEpois(obs(i) - oo_ij, th_j);
        vv_ij = pCEpois(ref_ne(i, k) - oo_ij_ne, th_j);
        logcdf_ij(k) = pConCop(uu_ij, vv_ij, cop_family, cop_param_ikj, true);
      }
      
      weight_ij = weight.slice(j).row(i).t();
      concump(i,j) = sum(arma::exp(arma::log(weight_ij) + logcdf_ij));
      
    }
    
  }
  
  return concump;
  
}

// Negative binomial NNMP
// [[Rcpp::export]]
arma::mat pCopNbNNMP_ref_covar(const arma::colvec& obs,
                               const int& cop_family,
                               const arma::cube& cop_param,
                               const arma::mat& ref_XX,
                               const arma::mat& ref_ne,
                               const arma::umat& ref_ne_idx,
                               const arma::mat& oo,
                               const arma::mat& bb,
                               const arma::colvec& size,
                               const arma::cube& weight) {
  
  int data_size = obs.n_rows;
  int sample_size = size.n_rows;
  int ne_size = ref_ne.n_cols;
  arma::colvec weight_ij, oo_j, oo_ij_nes, mar_param1;
  arma::mat weight_j;
  double cop_param_j, cop_param_ikj, th_ij, th_ij_ne, size_j, oo_ij, oo_ij_ne, uu_ij, vv_ij;
  
  arma::mat concump(data_size, sample_size);
  arma::colvec logcdf_ij, bb_j, th_j, th_ij_nes;
  arma::uvec ine_idx;
  
  // the first obs
  for (int j = 0; j < sample_size; ++j) {
    bb_j = bb.col(j);
    th_ij = as_scalar(exp(ref_XX.row(0) * bb_j));
    size_j = size(j);
    concump(0, j) = pCEnbinom(obs(0) - oo(0,j), size_j, th_ij);
  }
  
  // the second obs
  for (int j = 0; j < sample_size; ++j) {
    
    cop_param_j = cop_param(1, 0, j);
    
    oo_ij = oo(1, j);
    oo_ij_ne = oo(0, j);
    
    bb_j = bb.col(j);
    th_ij = as_scalar(exp(ref_XX.row(1) * bb_j));
    th_ij_ne = as_scalar(exp(ref_XX.row(0) * bb_j));
    size_j = size(j);
    
    uu_ij = pCEnbinom(obs(1) - oo_ij, size_j, th_ij);
    vv_ij = pCEnbinom(obs(0) - oo_ij_ne, size_j, th_ij_ne);
    concump(1, j) = pConCop(uu_ij, vv_ij, cop_family, cop_param_j);
  }
  
  // obs 3 - ne_size
  for (int i = 2; i < ne_size; ++i) {
    
    logcdf_ij.set_size(i);
    ine_idx = arma::trans(ref_ne_idx(i, arma::span(0, i - 1))) - 1;
    
    for (int j = 0; j < sample_size; ++j) {
      
      oo_j = oo.col(j);
      oo_ij_nes = oo_j(ine_idx);
      oo_ij = oo_j(i);
      
      bb_j = bb.col(j);
      th_j = exp(ref_XX * bb_j);
      th_ij = th_j(i);
      th_ij_nes = th_j(ine_idx);
      size_j = size(j);
      
      for (int k = 0; k < i; ++k) {
        cop_param_ikj = cop_param(i, k, j);
        oo_ij_ne = oo_ij_nes(k);
        th_ij_ne = th_ij_nes(k);
        uu_ij = pCEnbinom(obs(i) - oo_ij, size_j, th_ij);
        vv_ij = pCEnbinom(ref_ne(i, k) - oo_ij_ne, size_j, th_ij_ne);
        logcdf_ij(k) = pConCop(uu_ij, vv_ij, cop_family, cop_param_ikj, true);
      }
      
      weight_j = weight.slice(j);
      weight_ij = arma::trans(weight_j(i, arma::span(0, i - 1)));
      concump(i,j) = sum(arma::exp(arma::log(weight_ij) + logcdf_ij));
      
    }
    
  }
  
  // Remaining obs
  logcdf_ij.set_size(ne_size);
  for (int i = ne_size; i < data_size; ++i) {
    
    ine_idx = ref_ne_idx.row(i).t() - 1;
    
    for (int j = 0; j < sample_size; ++j) {
      
      oo_j = oo.col(j);
      oo_ij_nes = oo_j(ine_idx);
      oo_ij = oo_j(i);
      
      bb_j = bb.col(j);
      th_j = exp(ref_XX * bb_j);
      th_ij = th_j(i);
      th_ij_nes = th_j(ine_idx);
      size_j = size(j);
      
      for (int k = 0; k < ne_size; ++k) {
        cop_param_ikj = cop_param(i, k, j);
        oo_ij_ne = oo_ij_nes(k);
        th_ij_ne = th_ij_nes(k);
        uu_ij = pCEnbinom(obs(i) - oo_ij, size_j, th_ij);
        vv_ij = pCEnbinom(ref_ne(i, k) - oo_ij_ne, size_j, th_ij_ne);
        logcdf_ij(k) = pConCop(uu_ij, vv_ij, cop_family, cop_param_ikj, true);
      }
      
      weight_ij = weight.slice(j).row(i).t();
      concump(i,j) = sum(arma::exp(arma::log(weight_ij) + logcdf_ij));
      
    }
    
  }
  
  return concump;
  
}
