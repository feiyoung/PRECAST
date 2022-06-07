#ifndef utilMultIntcluster_h

#define utilMultIntcluster_h

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include<ctime>
#include "idrsc2.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

sp_mat get_spNbs(ivec y, const sp_mat& Adj);
MATTYPE calYenergy2D_sp(const arma::ivec& y, const arma::sp_mat& Adj, int K, const VECTYPE alpha, const float beta);
float obj_beta(const field<ivec>& yf, const field<MATTYPE>& Rf, 
                const arma::field<sp_mat>& Adjf, int K, const VECTYPE alpha, const float beta);
float objr_beta(const ivec& y, const MATTYPE& R, 
          const sp_mat& Adj, int K, const VECTYPE alpha, const float beta);
void multi_det_SkCpp(const MATTYPE& X, const ROWVECTYPE& Lam_vec0, const MATTYPE& W0, const MATTYPE& Ck, 
                     const ROWVECTYPE Muk, 
                     float& logdSk, VECTYPE& mSk);
void multi_det_SkCpp2(const MATTYPE& X, const ROWVECTYPE& Lam_vec0,
                      const MATTYPE& W0, const MATTYPE& Ck, 
                      const ROWVECTYPE Muk, const MATTYPE& Sigmak,
                      float& logdSk, VECTYPE& mSk);
VECTYPE decomp(const MATTYPE& Cki, const MATTYPE& W0);

#endif
