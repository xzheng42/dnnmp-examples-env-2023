#include <iostream>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>

using namespace std;
using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
// RcppArmadillo's sample function
////////////////////////////////////////////////////////////////////////////////
arma::colvec armaIntSample(const arma::colvec& x,
                           const int& size,
                           const bool& replace,
                           const arma::colvec& prob) {
  arma::colvec rsq = RcppArmadillo::sample(x, size, replace, prob);
  return(rsq);
}


////////////////////////////////////////////////////////////////////////////////
// Compute cumulative sum
////////////////////////////////////////////////////////////////////////////////
arma::colvec cumsum_cpp(const arma::colvec xx) {
  int kk = xx.n_rows;
  arma::colvec yy(kk);
  yy(0) = xx(0);
  for (int i = 1; i < kk; ++i) {
    yy(i) = sum(xx.rows(0, i));
    if (yy(i) > 1.0) {
      yy.rows(i, kk - 1).ones();
      break;
    }
  }
  return yy;
}


////////////////////////////////////////////////////////////////////////////////
// Distribution functions
////////////////////////////////////////////////////////////////////////////////
#define EPSILON DBL_EPSILON

double pCEpois(const double& xx,
               const double& rate) {
  
  int xx_int = floor(xx);
  double pp1 = R::ppois(xx_int, rate, true, false);
  double pp2 = R::dpois(xx_int + 1, rate, false);
  double pp = pp1 + (xx - xx_int) * pp2;
  
  return pp;
  
}

double pCEpois2(double xx,
                double rate,
                double rate_fake) {
  
  int xx_int = floor(xx);
  double pp1 = R::ppois(xx_int, rate, true, false);
  double pp2 = R::dpois(xx_int + 1, rate, false);
  double pp = pp1 + (xx - xx_int) * pp2;
  
  return pp;
  
}

double pCEnbinom(const double& xx,
                 const double& size,
                 const double& mu) {
  
  int xx_int = floor(xx);
  double prob = size / (size + mu);
  double pp1 = R::pnbinom(xx_int, size, prob, true, false);
  double pp2 = R::dnbinom(xx_int + 1, size, prob, false);
  double pp = pp1 + (xx - xx_int) * pp2;
  
  return pp;
  
}

double pCEnbinom2(double xx,
                  double size,
                  double mu) {
  
  int xx_int = floor(xx);
  double prob = size / (size + mu);
  double pp1 = R::pnbinom(xx_int, size, prob, true, false);
  double pp2 = R::dnbinom(xx_int + 1, size, prob, false);
  double pp = pp1 + (xx - xx_int) * pp2;
  
  return pp;
  
}

// [[Rcpp::export]]
double qCEpois(double pp,
               double rate,
               double rate_fake,
               double lower,
               double upper,
               double tol,
               int maxit) {
  
  double a, b, c, fa, fb, fc;
  int maxit_count = maxit + 1;
  
  fa = pCEpois(lower, rate) - pp;
  fb = pCEpois(upper, rate) - pp;
  
  a = lower; b = upper;
  c = a; fc = fa;
  
  if (fa == 0.0) {
    return a;
  }
  if (fb == 0.0) {
    return b;
  }
  
  while(maxit_count--) 
  {
    
    double prev_step = b - a;
    double tol_act;      // Actual tolerance		*/
    double p, q;         /* Interpolation step is calculated in the form p/q; division operations is delayed until the last moment	*/
    double new_step;     /* Step at this iteration	*/
    
    if (fabs(fc) < fabs(fb)) {
      a = b; b = c; c = a;
      fa = fb; fb = fc; fc = fa;
    }
    tol_act = 2 * EPSILON * fabs(b) + tol / 2;
    new_step = (c - b) / 2;
    
    if (fabs(new_step) <= tol_act || fb == 0.0) {
      return b;
    }
    
    if (fabs(prev_step) >= tol_act && fabs(fa) > fabs(fb)) {
      
      double t1, cb, t2;
      cb = c - b;
      
      if (a == c) {
        t1 = fb / fa;
        p = cb * t1;
        q = 1.0 - t1;
      } else {
        q = fa / fc; t1 = fb / fc; t2 = fb / fa;
        p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
        q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
      }
      
      if (p > 0.0) {
        q = -q;
      } else {
        p = -p;
      }
      
      if (p < (0.75 * cb * q - fabs(tol_act * q) /2) && p < fabs(prev_step * q / 2)) {
        new_step = p / q;
      }
      
    }
    
    if (fabs(new_step) < tol_act) {
      if (new_step > 0.0) {
        new_step = tol_act;
      } else {
        new_step = -tol_act;
      }
    }
    
    a = b; fa = fb;
    b += new_step; 
    fb = pCEpois(b, rate) - pp;
    if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
      c = a; 
      fc = fa;
    }
    
  }
  
  // fail
  warning("fail to converge");
  // stop("error: fail to converge.");
  
  return b;
  
}

// [[Rcpp::export]]
double qCEnbinom(double pp,
                 double size,
                 double mu,
                 double lower,
                 double upper,
                 double tol,
                 int maxit) {
  
  double a, b, c, fa, fb, fc;
  int maxit_count = maxit + 1;
  
  fa = pCEnbinom(lower, size, mu) - pp;
  fb = pCEnbinom(upper, size, mu) - pp;
  
  a = lower; b = upper;
  c = a; fc = fa;
  
  if (fa == 0.0) {
    return a;
  }
  if (fb == 0.0) {
    return b;
  }
  
  while(maxit_count--) 
  {
    
    double prev_step = b - a;
    double tol_act;      // Actual tolerance		*/
    double p, q;         /* Interpolation step is calculated in the form p/q; division operations is delayed until the last moment	*/
    double new_step;     /* Step at this iteration	*/
    
    if (fabs(fc) < fabs(fb)) {
      a = b; b = c; c = a;
      fa = fb; fb = fc; fc = fa;
    }
    tol_act = 2 * EPSILON * fabs(b) + tol / 2;
    new_step = (c - b) / 2;
    
    if (fabs(new_step) <= tol_act || fb == 0.0) {
      return b;
    }
    
    if (fabs(prev_step) >= tol_act && fabs(fa) > fabs(fb)) {
      
      double t1, cb, t2;
      cb = c - b;
      
      if (a == c) {
        t1 = fb / fa;
        p = cb * t1;
        q = 1.0 - t1;
      } else {
        q = fa / fc; t1 = fb / fc; t2 = fb / fa;
        p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
        q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
      }
      
      if (p > 0.0) {
        q = -q;
      } else {
        p = -p;
      }
      
      if (p < (0.75 * cb * q - fabs(tol_act * q) /2) && p < fabs(prev_step * q / 2)) {
        new_step = p / q;
      }
      
    }
    
    if (fabs(new_step) < tol_act) {
      if (new_step > 0.0) {
        new_step = tol_act;
      } else {
        new_step = -tol_act;
      }
    }
    
    a = b; fa = fb;
    b += new_step; 
    fb = pCEnbinom(b, size, mu) - pp;
    if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
      c = a; 
      fc = fa;
    }
    
  }
  
  // fail
  warning("fail to converge");
  // stop("error: fail to converge.");
  
  return b;
  
}







