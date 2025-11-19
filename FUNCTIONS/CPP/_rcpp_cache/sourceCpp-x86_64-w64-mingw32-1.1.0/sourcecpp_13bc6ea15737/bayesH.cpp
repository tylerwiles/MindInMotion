#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

double optimprice(double* x, unsigned int n);
double fmin(double* data, double ax, double bx, int n);
double logphxfunction(double H, double* x, unsigned int n);
double accrej(double* x, double logM, double add, double minu, double maxu, unsigned int n);
arma::vec acfHKp(double H, double maxlag);
arma::vec ltza(arma::vec rr, double* xr, unsigned int n);
int lev(double* r, int n, double* x, double* y, double* e, double EPS);
double levDet(int n, double* e);
double dot(int n, double* u, double* v);
double flipupdot(int n, double* u, double* v);
double sum(int n, double* u);
typedef double* VECTOR;
VECTOR cVector(long n);
void free_vector(VECTOR cVector);

//' Bayesian Inference of Hurst Exponent
//' 
//' Infer Hurst exponent of a time series based on accept-reject algorithm.
//' 
//' @param x A vector of time series
//' @param n An integer indicating the number of Hurst exponents to infer
//' 
//' @returns A probability distribution of the Hurst exponents inferred
//' 
//' @import Rcpp
//' @export
//' 
//' @details The Hurst exponent quantifies the temporal correlation among data points of a time series. This algorithm returns Hurst exponents with less variance compared to \code{dfa}. In addition, this algorithm is more robust to time series shorter than 512 data points. This algorithm estimates the Hurst exponent via Bayesian technique. Based on a predefined target distribution of Hurst exponent, an accept-reject algorithm is used to sample a posterior (probability) distribution of Hurst exponent. Common practice is to take the median of the probability distribution as the estimated Hurst exponent. For an example of using this algorithm on human movement data, refer to Likens et al. 2023.
//' 
//' @examples
//' # Generate example time series data
//' x = fgn_sim(n = 128, H = 0.9)
//' 
//' # Run BayesH
//' h.pdf = bayesH(x = x, n = 200)
//' 
//' # Take the median of the probability distribution
//' H = median(h.pdf)
//' 
//' @references 
//' Hurst, H. E. Long-Term Storage Capacity of Reservoirs. T. Am. Soc. Civ. Eng. 116, 770–799 (1951). https://doi.org/10.1061/TACEAT.0006518
//' 
//' Tyralis, H., & Koutsoyiannis, D. (2014). A Bayesian statistical model for deriving the predictive distribution of hydroclimatic variables. Climate dynamics, 42, 2867-2883. http://link.springer.com/10.1007/s00382-013-1804-y
//' 
//' Likens, A. D., Mangalam, M., Wong, A. Y., Charles, A. C., & Mills, C. (2024). Better than Detrended Fluctuation Analysis?A Bayesian Method for Estimating the Hurst Exponent in Behavioral Sciences. Preprint: https://www.researchgate.net/publication/378038184_Better_than_Detrended_Fluctuation_Analysis_A_Bayesian_method_for_estimating_the_Hurst_exponent_in_behavioral_sciences
//' 
//' Mangalam, M., Wilson, T., Sommerfeld, J., & Likens, A. D. (2024). Optimizing a Bayesian method for estimating the Hurst exponent in behavioral sciences. Preprint: https://www.researchgate.net/publication/381123281_Optimizing_a_Bayesian_method_for_estimating_the_Hurst_exponent_in_behavioral_sciences
//' 
//' 
// [[Rcpp::export]]
arma::vec bayesH(arma::vec x, unsigned int n) {
  
  double add = 0.001;
  double minu = 0.001;
  double maxu = 0.999;
  
  unsigned int nn = x.n_elem;
  
  VECTOR xx;
  xx = cVector(nn);
  for (unsigned int i = 0; i < nn; i++) {
    xx[i] = x[i];
  }
  
  double logM = optimprice(xx, nn);
  arma::vec pdf = arma::zeros(n);
  
  for (unsigned int i = 0; i < n; i ++) {
    pdf(i) = accrej(xx, logM, add, minu, maxu, nn);
  }

  return pdf;
  
}

double optimprice(double* x, unsigned int n) {
  
  double val = fmin(x, 0.00001, 0.99999, n);
  double logM = logphxfunction(val, x, n);
  
  return logM;
  
}

double fmin(double* data, double ax, double bx, int n) {
  
  const double c = (3 - sqrt(5)) * 0.5;
  double tol = pow(DBL_EPSILON,0.25);
  double u;
  double eps = DBL_EPSILON;
  double tol1 = eps + 1;
  eps = sqrt(eps);
  double a = ax;
  double b = bx;
  double v = a + c * (b - a);
  double w = v;
  double x = v;
  double d = 0;
  double e = 0;
  double fx = logphxfunction(x, data, n);
  fx = -1*fx;
  double fv = fx;
  double fw = fx;
  double tol3 = tol / 3;
  
  for(;;) {
    double xm = (a + b) * 0.5;
    double tol1 = eps * fabs(x) + tol3;
    double t2 = tol1 * 2;
    if (fabs(x - xm) <= t2 - (b - a) * 0.5) {
      break;
    }
    double p = 0;
    double q = 0;
    double r = 0;
    if (fabs(e) > tol1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2;
      if (q > 0) {
        p = -p;
      }
      else {
        q = -q;
      }
      r = e;
      e = d;
    }
    if (fabs(p) >= fabs(q * 0.5 * r) ||
        p <= q * (a - x) || p >= q * (b - x)) {
      if (x < xm) {
        e = b - x;
      }
      else {
        e = a - x;
      }
      d = c * e;
    }
    else {
      d = p / q;
      u = x + d;
      if (u - a < t2 || b - u < t2) {
        d = tol1;
        if (x >= xm) {
          d = -d;
        }
      }
    }
    if (fabs(d) >= tol1) {
      u = x + d;
    }
    else if (d > 0) {
      u = x + tol1;
    }
    else {
      u = x - tol1;
    }
    double fu = logphxfunction(u, data, n);
    fu = -1*fu;
    if (fu <= fx) {
      if (u < x) {
        b = x;
      }
      else {
        a = x;
      }
      v = w;    w = x;   x = u;
      fv = fw; fw = fx; fx = fu;
    } else {
      if (u < x) {
        a = u;
      }
      else {
        b = u;
      }
      if (fu <= fw || w == x) {
        v = w; fv = fw;
        w = u; fw = fu;
      }
      else if (fu <= fv || v == x || v == w) {
        v = u; fv = fu;
      }
    }
  }
  
  return x;
}

double accrej(double* x, double logM, double add, double minu, double maxu, unsigned int n) {
  
  double distance = maxu - minu;
  double logM1 = logM - log(distance) + add;
  
  double u = arma::randu();
  double logu = log(u) + logM1;
  
  double y = arma::randu(distr_param(minu, maxu));
  
  while (logu > (logphxfunction(y, x, n) - log(distance))) {
    
    u = arma::randu();
    
    logu = log(u) + logM1;
    
    y = arma::randu(distr_param(minu, maxu));
    
  }
  
  return y;

}

double logphxfunction(double H, double* x, unsigned int n) {
  
  unsigned int maxlag = n - 1;
  arma::vec acf = acfHKp(H, maxlag);
  arma::vec q = ltza(acf, x, n);
  double f = -0.5*q(3) - 0.5*(maxlag) * log(q(2)*q(0)-pow(q(1),2)) + (0.5*n - 1)*log(q(2));
  
  return f;
  
}

arma::vec acfHKp(double H, double maxlag) {
  
  unsigned int n = maxlag + 1;
  arma::vec acf(n);
  
  double H2 = H*2;
  
  
  acf(0) = 1;
  
  for (unsigned int i=0; i < maxlag; i++) {
    double k = i+1;
    acf(i+1) = 0.5*(pow((k+1),H2) - 2*pow(k,H2) + pow((k-1),H2));
  }    
  
  return acf;
  
}

arma::vec ltza(arma::vec rr, double* xr, unsigned int n) {
  
  
  double EPS = DBL_EPSILON;
  int fault = 1;
  int _fault1;
  int n1 = n-1;
  VECTOR r, y1, y2, e1, e2, e3;
  
  arma::vec y = zeros<arma::vec>(4);
  
  r = cVector(n);
  for (unsigned int i = 0; i < n; ++i){
    r[i] = rr[i];
  }
  
  y1 = cVector(n);
  y2 = cVector(n);
  e1 = cVector(n1);
  e2 = cVector(n1);
  e3 = cVector(n);
  for (unsigned int i = 0; i < n; ++i){
    e3[i] = 1;
  }
  
  _fault1 = lev(r,n,xr,y1,e1,EPS);
  arma::vec yy1 = zeros<arma::vec>(n);
  for (unsigned int i = 0; i < n; ++i){
    yy1[i] = y1[i];
  }
  
  if (_fault1 != 0){
    for (unsigned int i = 0; i < 2; i++){
      y(i) = 0.0;
    }
    fault = _fault1;
  }
  else
  {
    fault = 0;
    
    _fault1 = lev(r,n,e3,y2,e2,EPS);
    arma::vec yy2 = zeros<arma::vec>(n);
    for (unsigned int i = 0; i < n; ++i){
      yy2[i] = y2[i];
    }
    
    y(3) = levDet(n1,e2);
    double s1 = arma::accu(yy2);
    double s2 = arma::accu(yy1);
    double s3 = dot(n,xr,y1);
    y(0) = s3;
    y(1) = s2;
    y(2) = s1;
  }
  
  free_vector(r);
  free_vector(y1);
  free_vector(y2);
  free_vector(e1);
  free_vector(e2);
  free_vector(e3);
  
  return y;
  
}


int lev(double* r, int n, double* x, double* y, double* e, double EPS) {
  
  VECTOR v, l, b, c;
  
  int n1;
  int i, j, k, m;
  n1 = n-1;
  
  if (fabs(r[0] - 1.0) > EPS){
    return 2;
  }
  e[0] = 1.0 - r[1] * r[1];
  if (e[0] < EPS){
    return 1;
  }
  
  v = cVector(n1);
  l = cVector(n1);
  b = cVector(n);
  c = cVector(n1);
  v[0] = - r[1];
  l[0] = x[1] - r[1] * x[0];
  b[0] = - r[1];
  b[1] = 1.0;
  y[0] = (x[0] - r[1] * x[1]) / e[0];
  y[1] = l[0] / e[0];
  
  for (i = 1; i < n1; i++){
    v[i] = - dot(i + 1,r + 1,b) / e[i-1];
    e[i] = e[i-1] * (1 - v[i] * v[i]);
    l[i] = x[i+1] - flipupdot(i + 1,r + 1,y);
    
    for (k = 0; k < i + 1; k++){
      c[k] = b[i-k];
    }
    b[i + 1] = b[i];
    
    for (j = i; j > 0; j--){
      b[j] = b[j-1] + v[i] *c[j];
    }
    
    b[0] = v[i] * c[0];
    y[i+1] = (l[i] / e[i]) * b[i + 1];
    
    for (m = i; m > -1; m--){
      y[m] = y[m] + (l[i] / e[i]) * b[m];
    }
  }
  
  free_vector(v);
  free_vector(l);
  free_vector(b);
  free_vector(c);
  
  return 0;
  
}

double levDet(int n, double* e) {
  
  int i;
  double logDet;
  logDet = 0.0;
  for (i = 0; i < n; i++){
    logDet += log(e[i]);
  }
  return logDet;
  
}

double dot(int n,double* u,double* v) {
  
  double t = 0.0;
  int i;
  for (i = 0; i < n; i++)
    t += u[i] * v[i];
  return t;

}

double flipupdot(int n,double* u,double* v) {
  
  double t = 0.0;
  int i;
  for (i = 0; i < n; i++)
    t += u[n - 1 - i] * v[i];
  return t;
  
}

double sum(int n,double* u) {
  
  double t = 0.0;
  int i;
  for (i = 0; i < n; i++)
    t += u[i];
  return t;
  
}

VECTOR cVector(long n) {
  
  VECTOR vector;
  vector=(VECTOR) calloc(n, n*sizeof(double));
  
  return vector;
  
}

void free_vector(VECTOR cVector) { 
  
  free(cVector);
  
}

#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bayesH
arma::vec bayesH(arma::vec x, unsigned int n);
RcppExport SEXP sourceCpp_3_bayesH(SEXP xSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(bayesH(x, n));
    return rcpp_result_gen;
END_RCPP
}
