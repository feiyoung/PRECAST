#ifndef idrsc2_h
#define idrsc2_h

#include <RcppArmadillo.h>
#include <math.h>
#include "PRECAST_types.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


struct Objidrsc2{
    field<ivec> yf;
    field<MATTYPE> Ezz;
    field<MATTYPE> Vf;
    field<MATTYPE> Rf;
    VECTYPE beta0;
    MATTYPE Mu0;
    CUBETYPE Sigma0;
    CUBETYPE Psi0;
    MATTYPE W0;
    MATTYPE Lam0;
    double loglik;
    VECTYPE loglik_seq;
};


Objidrsc2 idrsc2(const field<MATTYPE>& Xf, const field<sp_mat> Adjf, const field<sp_mat> Adjf_car, field<ivec> yf,
                  const MATTYPE& Mu_int, CUBETYPE Sigma0, const MATTYPE& W_int,
                  const MATTYPE& Lam_int,  const CUBETYPE& Psi_int,
                  const VECTYPE& alpha, VECTYPE beta0, const VECTYPE& beta_grid,
                  const int& maxIter_ICM, const int& maxIter, const float& epsLogLik, const bool& verbose,
                  const bool& homo, const bool& homoClust,
                  const bool& Sigma_diag, const bool& mix_prop_heter, 
                  const bool& Sp2);


#endif /* CoMM_covar_pxem_ptr_hpp */
