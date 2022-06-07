#ifndef mt_paral_job_h
#define mt_paral_job_h

#include <stdio.h>
#include <math.h>
#include <RcppArmadillo.h>
#include <thread>
#include <mutex>
#include "idrsc2.h"


using namespace std;
using namespace arma;
using namespace Rcpp;


class par_iDRSC2{
public:
    
	  field<MATTYPE> Xf;
    field<sp_mat> Adjf;
    field<sp_mat> Adjf_car;
    field<imat> yf;
    field<MATTYPE> Mu0;
    field<CUBETYPE> Sigma0;
    MATTYPE W_int;
    MATTYPE Lam0;
    CUBETYPE Psi0;
    field<VECTYPE> alpha0;
    VECTYPE beta0;
    VECTYPE beta_grid;
    int maxIter_ICM;
    int maxIter;
    float epsLogLik;
    bool verbose;
    bool homo;
    bool homoClust;
    bool Sigma_diag;
    bool mix_prop_heter;
    bool Sp2;
    int maxK, minK;
    int g;    
    int current_idx = 0;
    struct Objidrsc2 output[50];

    

	par_iDRSC2(field<MATTYPE>& Xf, field<sp_mat>& Adjf, field<sp_mat>& Adjf_car, field<imat>& yf,
                      const field<MATTYPE>& Mu0, field<CUBETYPE> Sigma0, const MATTYPE& W_int,
                      MATTYPE& Lam0, CUBETYPE& Psi0,
                      const field<VECTYPE>& alpha0, VECTYPE beta0, const VECTYPE& beta_grid,
                      const int& maxIter_ICM, const int& maxIter, const float& epsLogLik, const bool& verbose,
                      const bool& homo, const bool& homoClust, const bool& Sigma_diag, const bool& mix_prop_heter, const bool& Sp2, 
                      const int maxK, const int minK){

		this->Xf = Xf;
		this->Adjf = Adjf;
		this->Adjf_car = Adjf_car;
		this->yf = yf;
		this->Mu0 = Mu0;
		this->Sigma0 = Sigma0;
		this->W_int = W_int;
		this->Lam0 = Lam0;
        this->Psi0 = Psi0;
        this->alpha0 = alpha0;
        this->beta0 = beta0;
        this->beta_grid = beta_grid;
        this->maxIter_ICM = maxIter_ICM;
        this->maxIter = maxIter;
        this->epsLogLik = epsLogLik;
        this->verbose = verbose;
        this->homo = homo;
        this->homoClust = homoClust;
        this->Sigma_diag = Sigma_diag;
        this->mix_prop_heter = mix_prop_heter;
        this->Sp2 = Sp2;
        this->maxK = maxK;
        this->minK = minK;
        
	}

	void loop_by_K_idrsc2(int g);
	void update_by_thread_idrsc2(int thread_id);
	int  next_idrsc2();

};


#endif 

