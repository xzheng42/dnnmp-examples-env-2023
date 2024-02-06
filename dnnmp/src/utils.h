#ifndef UTILS_H
#define UTILS_H
#include <RcppArmadillo.h>

arma::colvec armaIntSample(const arma::colvec& x, const int& size, const bool& replace, const arma::colvec& prob);

arma::colvec cumsum_cpp(const arma::colvec xx);

double pCEpois(const double& xx, const double& rate);

double pCEpois2(double xx, double rate, double rate_fake);

double pCEnbinom(const double& xx, const double& size, const double& mu);

double pCEnbinom2(double xx, double size, double mu);

double qCEpois(double pp, double rate, double rate_fake,
               double lower, double upper, double tol, int maxit);

double qCEnbinom(double pp, double size, double mu,
                 double lower, double upper, double tol, int maxit);
  
#endif 