#include <Rcpp.h>
using namespace Rcpp;

//' @title Newton-Raphson algorithm for cauchy(theta) using Rcpp.
//' @description Using three input to get the result of Newton-Raphson algorithm of cauchy distribution.
//' @param x the data of cauchy distribution
//' @param theta the initial value
//' @param method method=1 for one step iteration and method=2 for multi-step iteration
//' @return the estimation of param theta
//' @examples
//' \dontrun{
//' x<-rcauchy(200,location=1)
//' theta<-median(x)
//' thetahat<-NC(x,theta,1)
//' print(thetahat)
//' }
//' @export
// [[Rcpp::export]]
double NC(const NumericVector &x, double theta, int method) {
  if (method == 1) {
    return (theta + sum(2.0 * (x - theta) / (1.0 + pow(x - theta, 2))) / sum(2.0 * (1.0 - pow(x - theta, 2)) / pow(1.0 + pow(x - theta, 2), 2)));
  } else if (method == 2) {
    double t = theta + sum(2.0 * (x - theta) / (1.0 + pow(x - theta, 2))) / sum(2.0 * (1.0 - pow(x - theta, 2)) / pow(1.0 + pow(x - theta, 2), 2));
    int k = 1;
    while (std::abs(theta - t) > 1e-3) {
      theta = t;
     t = theta + sum(2.0 * (x - theta) / (1.0 + pow(x - theta, 2))) / sum(2.0 * (1.0 - pow(x - theta, 2)) / pow(1.0 + pow(x - theta, 2), 2));
     k++;
     if (k > 1000) {
       Rcpp::Rcout << "error" << std::endl;
       break;
     }
   }
   return t;
  }else{
    Rcpp::Rcout << "please input the correct method" << std::endl;
    return 0;
 }
}


//' @title Fisher scoring algorithm for cauchy(theta) using Rcpp.
//' @description Using three input to get the result of Fisher scoring algorithm of cauchy distribution.
//' @param x the data of cauchy distribution
//' @param theta the initial value
//' @param method method=1 for one step iteration and method=2 for multi-step iteration
//' @return the estimation of param theta
//' @examples
//' \dontrun{
//' x<-rcauchy(200,location=1)
//' theta<-median(x)
//' thetahat<-FSC(x,theta,1)
//' print(thetahat)
//' }
//' @export
// [[Rcpp::export]]
double FSC(const NumericVector &x, double theta, int method) {
  if (method == 1) {
    return (theta + 2.0 * mean(2.0 * (x - theta) / (1.0 + pow(x - theta, 2))));
  } else if (method == 2) {
    double t = theta + 2.0 * mean(2.0 * (x - theta) / (1.0 +pow(x - theta, 2)));
    int k = 1;
    while (std::abs(theta - t) > 1e-3) {
      theta = t;
      t = theta + 2.0 * mean(2.0 * (x - theta) / (1.0 + pow(x - theta, 2)));
      k++;
      if (k > 1000) {
        Rcpp::Rcout << "error" << std::endl;
        break;
      }
    }
    return t;
  }else{
    Rcpp::Rcout << "please input the correct method" << std::endl;
    return 0;
  }
}