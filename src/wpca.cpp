// This Cpp includes  functions to impliment Weighted PCA in 
// reference Bai, J., & Liao, Y. (2013). Statistical inferences using large estimated covariances for panel data and factor models. arXiv preprint arXiv:1307.2662.
// and Inferences in panel data with interactive effects using largecovariance matrices
// It considers the heterogeneous error term in approximated factor model.
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "idrsc2.h"


#define INT_MIN (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;


/*
  * Auxiliary
*/
  //' @keywords internal
//' @noRd
//' 
//[[Rcpp::export]] 
VECTYPE calculateWeight(const MATTYPE& X, const int& nPCs){
  int n = X.n_rows;
  MATTYPE U, V;
  VECTYPE s;
  svd_econ(U, s, V, X);
  MATTYPE hX = U.cols(0, nPCs-1) *diagmat(s.subvec(0, nPCs-1)) * trans(V.cols(0,nPCs-1));
  ROWVECTYPE Lam_vec = sum((hX-X) % (hX-X))/ n;
  return Lam_vec.t();
}

//[[Rcpp::export]]   
Rcpp::List wpcaCpp(const MATTYPE&X, const int& nPCs, const bool& weighted=true){
  MATTYPE U, V;
  VECTYPE s;
  MATTYPE PCs, loadings;
  svd_econ(U, s, V, X);
  PCs = U.cols(0, nPCs-1) *diagmat(s.subvec(0, nPCs-1));
  loadings = V.cols(0, nPCs-1);
  MATTYPE dX = PCs * loadings.t() - X;
  ROWVECTYPE Lam_vec = mean(dX % dX);
  if(weighted){
    svd_econ(U, s, V, X*diagmat(1.0/ sqrt(Lam_vec)));
    // vec s2 =  s % s; // s; //
    MATTYPE loadings_unscale = diagmat(sqrt(Lam_vec)) * V.cols(0, nPCs-1);
    MATTYPE  V1;
    VECTYPE s1;
    svd_econ(loadings, s1, V1, loadings_unscale);
    PCs = U.cols(0, nPCs-1) * diagmat(s.subvec(0, nPCs-1)) * V1 * diagmat(s1);
    dX = PCs * loadings.t() - X;
    Lam_vec = mean(dX % dX);
  }
  List output = List::create(
    Rcpp::Named("PCs") = PCs,
    Rcpp::Named("loadings") = loadings,
    Rcpp::Named("Lam_vec") = Lam_vec);
  
  return output;
}
