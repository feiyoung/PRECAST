#include "mt_paral_job.h"
#include "idrsc2.h"



using namespace std;
using namespace arma;
using namespace Rcpp;


void par_iDRSC2::loop_by_K_idrsc2(int g){
    
	
  // cout << "starting parallel " << endl;
  // cout << "k is " << g << endl;    
    
  int r_max = Xf.n_rows;
    
  // cout << "checking point 1" << endl;
  MATTYPE Mu_int = Mu0(g);
  MATTYPE Lam_int = Lam0;
  CUBETYPE Psi_int = Psi0;   
  VECTYPE alpha_int = alpha0(g);
  CUBETYPE Sigma_int = Sigma0(g);
  field<ivec> yf_int(r_max);
    
    
    
  for (int r = 0; r < r_max; r++){
    ivec tmp_yf = yf(r).col(g);
    yf_int(r) = tmp_yf;
  }
  
    
  output[g] = idrsc2(Xf, Adjf, Adjf_car, yf_int,
                  Mu_int, Sigma_int,  W_int,
                  Lam_int, Psi_int,
                  alpha_int,  beta0,  beta_grid,
                  maxIter_ICM,  maxIter, epsLogLik, verbose,
                  homo, homoClust, Sigma_diag, mix_prop_heter, Sp2);  
	// reset to free memory
    
  Mu_int.reset();
  Lam_int.reset();
  Psi_int.reset();
  alpha_int.reset();
  Sigma_int.reset();
  yf_int.reset();
    


}

std::mutex _mtx22;
int par_iDRSC2::next_idrsc2(){
	std::lock_guard<std::mutex> lockGuard(_mtx22);
	if (current_idx >= maxK - minK + 1){
		return -1;
	}
	current_idx++;
	return current_idx - 1;
}

void par_iDRSC2::update_by_thread_idrsc2(int thread_id){
	while (true){
		int idx = next_idrsc2();
		if (idx == -1){
			break;
		}
		loop_by_K_idrsc2(idx);
	}
}





