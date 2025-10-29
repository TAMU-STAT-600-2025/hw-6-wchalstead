#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Soft-thresholding function, returns scalar
// [[Rcpp::export]]
double soft_c(double a, double lambda){
  // Your function code goes here
  return ((a > 0) - (a < 0)) * std::max(abs(a) - lambda, 0.0);
}

// Lasso objective function, returns scalar
// [[Rcpp::export]]
double lasso_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda){
  // Your function code goes here
  double n = Ytilde.n_rows;
  double f_obj = pow((2.0 * n), -1) * pow(arma::norm((Ytilde - Xtilde * beta)), 2) + lambda * sum(abs(beta));
  return f_obj;
}

// Lasso coordinate-descent on standardized data with one lamdba. Returns a vector beta.
// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, double lambda, const arma::colvec& beta_start, double eps = 0.001){
  // Your function code goes here
  // Initialize n and p
  double n = Xtilde.n_rows;
  double p = Xtilde.n_cols;
  
  // Initialize starting points
  arma::colvec beta = beta_start;
  double fprevious = lasso_c(Xtilde, Ytilde, beta, lambda);
  arma::colvec resid = Ytilde - Xtilde * beta;
  double error = 10000;
  
  // Main loop
  while (error > eps) {
    
    // Update beta coordinate-wise
    for (int j = 0; j < p; j++){
      double beta_old = beta(j);
      beta(j) = soft_c(as_scalar((1/n) * Xtilde.col(j).t() * (resid + Xtilde.col(j) * beta(j))), lambda);
      resid = resid - Xtilde.col(j) * (beta(j) - beta_old);
    }
    
    // Update error and previous objective function
    error = abs(fprevious - lasso_c(Xtilde, Ytilde, beta, lambda));
    fprevious = lasso_c(Xtilde, Ytilde, beta, lambda);
  }
  
  return beta;
  
}  

// Lasso coordinate-descent on standardized data with supplied lambda_seq. 
// You can assume that the supplied lambda_seq is already sorted from largest to smallest, and has no negative values.
// Returns a matrix beta (p by number of lambdas in the sequence)
// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps = 0.001){
  // Your function code goes here
  // Initialize n and p
  double n = Xtilde.n_rows;
  double p = Xtilde.n_cols;
  double nlambda = lambda_seq.n_rows;
  
  // Initialize matrix beta
  arma::mat beta(p, nlambda, arma::fill::zeros);
  
  // Main loop
  for (int j = 0; j < nlambda; j++){
    beta.col(j) = fitLASSOstandardized_c(Xtilde, Ytilde, as_scalar(lambda_seq.row(j)), beta.col(std::max(j - 1, 0)), eps);
  }
  
  return beta;
}