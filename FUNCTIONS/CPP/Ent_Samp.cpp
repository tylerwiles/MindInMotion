#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' Sample Entropy
 //' 
 //' Calculate the sample entropy of a time series.
 //' 
 //' @param x A single column time series
 //' @param m The length of the vectors to be compared for matches
 //' @param R The radius for accepting matches
 //' 
 //' @returns The output of the algorithm is a single value that reflects the entropy of the time series in bits.
 //' 
 //' @import Rcpp
 //' @export
 //'
 //' @details 
 //' Sample entropy is another of the entropy algorithms that aims to quantify the predictability of a time series. Originally developed by Richman and Moorman, sample entropy is an improvement on approximate entropy, designed to overcome the lack of consistency seen in physiological time series and to reduce the bias introduced by self-matches.
 //' 
 //' Sample entropy captures the predictability of the time series by taking the negative natural log of the conditional probability that two sequences similar for the first \eqn{m} points remain similar at the next point, within a given tolerance \eqn{r}. Unlike ApEn, SampEn does not include self-matches in this probability calculation, which helps reduce bias and dependency on the length of the time series. The full equation is as follows: 
 //' 
 //' \eqn{SE = -ln(\frac{A}{B})}
 //' 
 //' where \eqn{A} is the number of matching vectors of length \eqn{m +1} and \eqn{B} is the number of matching vectors of length m.
 //' 
 //' Like approximate entropy, sample entropy is measured in bits where lower bits (information content) indicates more predictably thus, a larger value would indicate more randomness and less predictability. While sample entropy is more robust to shorter time series and parameter selections, these things should always be kept in mind for analysis. Best practice is to use a range of parameter values as different values may reveal different results.
 //' 
 //' @examples 
 //' 
 //' x = rnorm(1000)
 //' m = 2
 //' R = 0.2
 //' 
 //' SE = Ent_Samp(x, m, R)
 //' 
 //' @references
 //' Richman, J.S., Moorman, J.R., 2000. Physiological time-series analysis using approximate entropy and sample entropy. Am. J. Physiol. Heart Circ. Physiol. 278. https://doi.org/10.1152/ajpheart.2000.278.6.H2039
 // [[Rcpp::export]]
 double Ent_Samp(arma::colvec x, int m, double R) {
   // double Ent_Samp(arma::colvec x, int m, double R) {
   
   
   double SE = 0.0;
   double Bmr, Amr;
   
   double r = R * stddev(x);
   double N = x.size();
   
   arma::mat dij = arma::mat(m + 1, N - m + 1);
   
   arma::colvec d, d1; // Initialize column vectors. Unsure if it is actually needed
   
   // Initialize row vectors
   arma::rowvec dj = arma::rowvec(N - m + 1);
   arma::rowvec dj1 = arma::rowvec(N - m);
   arma::rowvec Bm = arma::rowvec(N - m + 1);
   arma::rowvec Am = arma::rowvec(N - m + 1);
   
   for (int i = 0; i < N - m + 1; i++) {
     for (int k = 0; k < m + 1; k++) {
       if (k < m) {
         dij.row(k) = abs(x.subvec(k, N - m + k) - x(i + k)).t();
       } else {
         if (i < N - m) {
           dij.submat(m, 0, m, N - m - 1) = abs(x.subvec(m, N - 1) - x(i + k)).t();
           dij.submat(m, N - m, m, N - m) = 2*r;
         } else {
           dij.row(k) = arma::rowvec(N - m + 1, arma::fill::value(2*r));
         }
       }
     }
     dj = max(dij.rows(0, m - 1), 0);
     dj1 = max(dij.submat(0, 0, m, N - m - 1), 0);
     d = dj(find(dj <= r));
     d1 = dj1(find(dj1 <= r));
     Bm(i) = d.n_elem - 1;
     Am(i) = d1.n_elem - 1;
   }
   
   // return List::create(Named("Am") = Am, Named("Bm") = Bm);
   // return List::create(Named("dij") = dij);

   Bm /= (double)(N - m);
   Am /= (double)(N - m - 1);

   Bmr = sum(Bm) / (N - m);
   Amr = sum(Am.subvec(0,N-m-1)) / (N - m);

   SE = -log(Amr / Bmr);


   return(SE);
   
 }
 
 