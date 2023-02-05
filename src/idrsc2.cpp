// This script implement integrating any finite data samples (such as spatial transcriptomics data)
// by spatial integrative simultaneous dimension reduction and clustering with heter-variance error,
// where two Markove random fields are used to capture the spatial dependence
// updated date: 2022-03-20


#include "RcppArmadillo.h"
#include "mt_paral_job.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include<ctime>
#include"utilMultIntcluster.h"



#define INT_MIN (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;


// See "utilMultIntcluster.h" file for the following functions
// sp_mat get_spNbs(ivec y, const sp_mat& Adj) {   // ivec是索引型向量
//   // row is for pixel.
//   //output a sparse matrix, i-th row contains labels of neighbor_i. 
//   // Make const iterator
//   arma::sp_mat::const_iterator start = Adj.begin(); //构造一个sp_mat的常数迭代器,常数迭代器只可读，不可写，实现对矩阵按照列对每个非零元素进行访问。
//   //arma::sp_mat::const_iterator end   = Adj.end();
//   
//   // Calculate number of nonzero points
//   //int n = std::distance(start, end);
//   int n = Adj.n_nonzero; // 计算Adj的所有非零元的个数
//   //cout << "n=" << n << endl;
//   //cout << "n=" << Adj.n_nonzero << endl;
//   
//   sp_mat spNbs(y.n_elem, y.n_elem);    // neiborhood state matrix, matched with Adj.
//   
//   
//   arma::sp_mat::const_iterator it = start; // Note spNbs is not a symmetric matrix, the nonzero in i-th row is the class label of sample i.
//   for(int i = 0; i < n; ++i)
//   {
//     //temp(0) = it.row();
//     //temp(1) = it.col();
//     spNbs(it.row(), it.col()) = y(it.col()); // it只自加非零元个数次，得到每个i对应的邻居的状态
//     ++it; // increment
//   }
//   
//   return spNbs.t(); // return the class label of neighbor matrix, i-th column is the neighbor label of sample i
// }
// 
// 
// arma::mat calYenergy2D_sp(const arma::ivec& y, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta)	{
//   // Calculate the value of energy function of y, which is equal to negative logliklihood up to a constant
//   int n = y.n_rows;
//   arma::sp_mat spNbs_t = get_spNbs(y, Adj); // transform spNbs to iterate by column.
//   arma::mat Uy(n, K);
//   double n_sameS;
//   int i, k, nn;
//   for (k = 0; k < K; k++)
//   {
//     for (i = 0; i < n; i++)
//     {
//       arma::sp_mat col(spNbs_t.col(i)); // the class label of neighbors of i-th sample.
//       n_sameS = 0;
//       
//       nn = col.n_nonzero; // the number of neighbors of i-th sample
//       for (arma::sp_mat::iterator j = col.begin(); j != col.end(); ++j) {
//         n_sameS += ((*j) == (k+1));
//         
//       }
//       Uy(i, k) = alpha(k) + beta * (nn - n_sameS)/2;
//       
//       
//     }
//   }
//   
//   arma::mat C_mat = normalise(exp(-Uy), 1, 1); // pseudo likelihood of Y.
//   Uy = -log(C_mat); // normalized Uy, this is the energy of y.
//   return Uy;
//   
// }
// 
// //
// double obj_beta(const field<ivec>& yf, const field<mat>& Rf, 
//                 const arma::field<sp_mat>& Adjf, int K, const arma::vec alpha, const double beta)	{
//   int r, r_max = yf.n_elem;
//   double objval = 0;
//   for(r=0; r< r_max; ++r){
//     mat Uy1 = calYenergy2D_sp(yf(r), Adjf(r), K, alpha, beta); // Uy was normalized, so there is no need to normalized Uy. 
//     objval += -accu(Rf(r) % Uy1);
//   }
//   
//   return objval;
// }
//  
// //
// double objr_beta(const ivec& y, const mat& R, 
//                  const sp_mat& Adj, int K, const arma::vec alpha, const double beta)	{
//   
//   double objval = 0;
//   mat Uy1 = calYenergy2D_sp(y, Adj, K, alpha, beta); // Uy was normalized, so there is no need to normalized Uy. 
//   objval = -accu(R % Uy1);
//   
//   
//   return objval;
// }
//  
//   
// 
// void multi_det_SkCpp(const arma::mat& X, const arma::rowvec& Lam_vec0, const arma::mat& W0, const arma::mat& Ck, 
//                      const arma::rowvec Muk, 
//                      double& logdSk, arma::vec& mSk){
//   //int p = X.n_cols;
//   int n = X.n_rows;
//   // int p = X.n_cols;
//   // // mSk = zeros(n);
//   // S2k = zeros(p);
//   // dSk = 0;
//   
//   mat WC12,  tmp2;
//   vec tmp1, s, tmp3;
//   mat U, V, X_tk;
//   
//   svd(U, s, V, Ck.i());
//   
//   WC12 = W0 * (U * diagmat(sqrt(s)));
//   WC12 = sp_mat(diagmat(1/sqrt(Lam_vec0))) * WC12;  // change to sparse matrix multiplication.
//   vec d = svd(WC12);
//   //dSk = arma::as_scalar(prod(1- d % d)) / prod(Lam_vec0);
//   logdSk = accu(log(1 - d%d)) - accu(log(Lam_vec0));
//   X_tk = (X - repmat(Muk* W0.t(), n, 1)) * sp_mat(diagmat(1/sqrt(Lam_vec0))) ;  // change to sparse matrix multiplication.
//   tmp1 = sum(X_tk % X_tk, 1);
//   tmp2 = X_tk * WC12;
//   tmp3 = sum(tmp2 % tmp2, 1);
//   mSk = tmp1 - tmp3;
// }
// 
// vec decomp(const mat& Cki, const mat& W0){
//   vec s, tmp1;
//   mat U, V, WC12;
//   svd(U, s, V, Cki);
//   WC12 = W0 * (U * diagmat(sqrt(s)));
//   tmp1 = sum(WC12 % WC12, 1);
//   return tmp1;
// }  
// 


// update mean component of Gaussian mixture model: Mu0
MATTYPE update_Mu1(const arma::field<MATTYPE>& Rf, const arma::field<CUBETYPE>& Ez,
                   const CUBETYPE& Sigma, 
               const MATTYPE& Nmat){
  
  int r, k, K= Rf(0).n_cols, q= Ez(0).n_cols, r_max = Rf.n_elem;
  MATTYPE Mu(q, K, fill::zeros), A_muk;
  VECTYPE b_muk;
  
  for(k=0;k < K ; ++k){
    A_muk = zeros<MATTYPE>(q,q);
    b_muk = zeros<VECTYPE>(q,1);
    for(r=0; r< r_max; ++r){
      MATTYPE inv_sympd_Sigma = inv_sympd(Sigma.slice(k)); //yangyi
      A_muk += Nmat(r,k) * inv_sympd_Sigma; //yangyi
      b_muk += inv_sympd_Sigma * (trans(Ez(r).slice(k)) * Rf(r).col(k)); //yangyi
      
    }
    //Mu.col(k) = A_muk.i() * b_muk;
      Mu.col(k) = solve(A_muk, b_muk);
  }
  return Mu.t();
}
  

// update covariance component of GMM: Sigma0
CUBETYPE update_Sigma1(const field<MATTYPE>& Rf, const field<CUBETYPE>& Ez, const field<CUBETYPE>& Ci,
                         const MATTYPE& Mu, const MATTYPE& Nmat, 
                         const bool& homoClust=false, const bool& Sigma_diag=false){
  int r,k, K= Mu.n_rows, q=Mu.n_cols, r_max=Rf.n_elem;
  CUBETYPE Sigma0(q,q , K);
  MATTYPE Smat(q, q, fill::zeros);
  float n_sum;
  if(homoClust){// Sigma is not related to  k.
    Smat = zeros<MATTYPE>(q,q);
    n_sum = 0;
    for(k = 0; k<K; ++k){  
      for(r=0; r< r_max; ++r){
          Smat += (trans(Ez(r).slice(k) - repmat(Mu.row(k), Rf(r).n_rows, 1)) % trans(repmat(Rf(r).col(k), 1, q))) * 
            (Ez(r).slice(k) - repmat(Mu.row(k), Rf(r).n_rows, 1)) + Nmat(r,k)*  (Ci(r).slice(k)); //yangyi n-by-n matrix to n-by-q matrix
          n_sum  += Nmat(r,k);
      }
    }
    for(k = 0; k<K; ++k){ 
        if(Sigma_diag){
          Sigma0.slice(k) = diagmat(Smat) / n_sum;
        }else{
          Sigma0.slice(k) = Smat / n_sum;
        }
    }
    }else{// Sigma is  related to  k.
      for(k = 0; k<K; ++k){  
        Smat = zeros<MATTYPE>(q,q);
        n_sum = 0;
        for(r=0; r< r_max; ++r){
          Smat += (trans(Ez(r).slice(k) - repmat(Mu.row(k), Rf(r).n_rows, 1)) % trans(repmat(Rf(r).col(k), 1, q))) * 
            (Ez(r).slice(k) - repmat(Mu.row(k), Rf(r).n_rows, 1)) + Nmat(r,k)*  (Ci(r).slice(k)); //yangyi n-by-n matrix to n-by-q matrix
          n_sum  += Nmat(r,k);
        }
        if(Sigma_diag){
          Sigma0.slice(k) = diagmat(Smat) / n_sum;
        }else{
          Sigma0.slice(k) = Smat / n_sum;
        }
      }
    }   

  return(Sigma0);
}
  
  

  
  
// update the loading matrix in probabilistic PCA(PPCA) model: W
MATTYPE update_W1(const field<MATTYPE>& Xf, const field<MATTYPE>& Rf, const field<CUBETYPE>& Ez,  
              const MATTYPE& Lam, const field<CUBETYPE>& Ci){
  
  int  r, j, k, K= Rf(0).n_cols, q= Ez(0).n_cols, p = Xf(0).n_cols, r_max=Rf.n_elem;
  MATTYPE  A_w(p,q, fill::zeros);
  for(k=0; k<K; ++k){
    for(r=0; r<r_max; ++r){
      A_w += trans(Xf(r).each_row() % (1.0/Lam.row(r))) *  (Ez(r).slice(k).each_col()% Rf(r).col(k));
      // A_w += (sp_mat(diagmat(1.0/Lam.row(r))) * Xf(r).t()) * sp_mat(diagmat(Rf(r).col(k))) * Ez(r).slice(k);
      
    }
  }
  
  CUBETYPE B_arr(q,q,r_max, fill::zeros);
  for(r=0; r<r_max; ++r){
    for(k=0; k<K; ++k){
      B_arr.slice(r) += (Ez(r).slice(k).t()%trans(repmat(Rf(r).col(k), 1, q)))*  Ez(r).slice(k) + sum(Rf(r).col(k))*Ci(r).slice(k);
      
    }
  }
  
  MATTYPE W(p, q, fill::zeros), tmpMat;
  for(j=0; j<p; ++j){
    tmpMat = zeros<MATTYPE>(q,q);
    for(r=0; r<r_max; ++r){
      tmpMat += 1.0/Lam(r,j)*B_arr.slice(r);
      
    }
    W.row(j) += A_w.row(j)* tmpMat.i();
  }
  
  return W;
}
  
  
// update the covariance matrix of residual in PPCA model: Lambda0
MATTYPE update_Lam1(const field<MATTYPE>& Rf, const field<MATTYPE>& Xf, const MATTYPE& W, const field<CUBETYPE>& Ez, 
                const field<CUBETYPE>& Ci,
                const MATTYPE& Nmat, const bool& homo = false){
  int r, k, K= Rf(0).n_cols, p = Xf(0).n_cols, r_max=Rf.n_elem;
  VECTYPE Lsum(p,fill::zeros);
  MATTYPE Lam(r_max, p);
  for(r=0; r< r_max; ++r){
    Lsum = zeros<VECTYPE>(p,1);
    MATTYPE tmpXk;
    for(k=0; k<K; ++k){
      tmpXk = (Xf(r) - Ez(r).slice(k) * W.t() );
      Lsum += trans(sum(tmpXk % tmpXk % repmat(Rf(r).col(k), 1, p)));
      Lsum += Nmat(r,k) * decomp(Ci(r).slice(k), W); // fast SVD method
    }
    if(homo){
      Lam.row(r) = mean(Lsum)* ones<MATTYPE>(1, p) / (Xf(r).n_rows*1.0);
    }else{
      Lam.row(r) = Lsum.t()/(Xf(r).n_rows*1.0);
    }
    
  }
  // replace little values to increase the stability of the algorithm.
  uvec id1 =  find(Lam <1e-7);
  Lam(id1) = 1e-7*ones<VECTYPE>(id1.n_elem, 1);
  return Lam; 
}

// evaluate the mean of spatial component of each spot in latent features.
MATTYPE get_Vmean(const MATTYPE& V, const arma::sp_mat& Adj){
  int i, n = V.n_rows, q= V.n_cols;
  VECTYPE m(n);
  MATTYPE Uv(n, q, fill::zeros);
  for (i = 0; i < n; i++)
  {
    arma::vec col(Adj.col(i)); // the class label of neighbors of i-th sample.
    uvec q1 = find(col > 0);
    // cout<<q1.n_rows<<endl;
    if( q1.n_rows>0){
      Uv.row(i) = mean(V.rows(q1));
    }
    
  }

  return Uv;
}
  

// Conduct ICM step in ICM-EM algorithm
void runICM_sp1(const arma::field<MATTYPE>& Xf, arma::field<MATTYPE>& Vf, arma::field<ivec>& yf, const MATTYPE& W0, 
                       const MATTYPE& Lam0, const MATTYPE& Mu0,
                       const CUBETYPE& Sigma0, CUBETYPE& Psi0,const arma::field<sp_mat>& Adjf, const arma::field<sp_mat>& Adjf_car,
                       const VECTYPE& alpha, const VECTYPE& beta_grid,
                       VECTYPE& beta, int maxIter_ICM, const bool& mix_prop_heter, const bool& Sp2, 
                       field<MATTYPE>& Rf, field<MATTYPE>& Muv, field<CUBETYPE>& Ezv, field<CUBETYPE>& Varzv, 
                       field<CUBETYPE>& Ez, field<CUBETYPE>& Varz, field<CUBETYPE>& Ev,  field<CUBETYPE>& Varv, float& loglik){
  // C++ can only set the last parameters to have default values.
  // Target: predict class labels Y and spatial components of latent features V,
  // posterior probability of Y: R, posterior expectation and covariance of V, Z and (V+Z):
  // Ev, Ez, Evz, Varv, Varz, Varvz, and update covariance of V (Psi_r) and 
  // smoothing parameter beta_r by using grid search.
  
  // basic info.
  int r, r_max = Xf.n_elem, K = Mu0.n_rows, q= Mu0.n_cols;
  int i, iter, k, n;
  
  // two cached objects used for parameters update.
  field<MATTYPE> Ux(r_max);
  field<CUBETYPE> WtSW(r_max), XSW(r_max); // q*q*K, n*q*K
  
  float  logdSk;
  // evaluate energy of x, Ux
  MATTYPE Crk, bCrk, WtLrW, WtbSrkW;
  
  for(r = 0; r < r_max; ++r){
    // evaluate energy of x, Ux
    
    //
    MATTYPE Lam0_rW0 = trans(repmat(1.0/ Lam0.row(r), q, 1))%W0; //yangyi
    WtLrW = W0.t() * Lam0_rW0; // O(p^2 q)
    MATTYPE XrLrW = Xf(r)* Lam0_rW0; // O(npq + p^2 q)
    MATTYPE  XrbSrkW_uv;
    n = Xf(r).n_rows; // tmp object
    Ux(r) = zeros<MATTYPE>(n, K);
    if(Sp2){
      WtSW(r) = zeros<CUBETYPE>(q,q,K);
      XSW(r) = zeros<CUBETYPE>(n, q,K);
    }
    // Rprintf("evaluate posterior! \n");
    VECTYPE mSk;
    for (k = 0; k < K; k++)	{
      
      if(Sp2){
        // save some objects for predict V.
        Crk = WtLrW +  inv_sympd(Sigma0.slice(k));
        WtSW(r).slice(k) = WtLrW -WtLrW*Crk.i()*WtLrW; // save some objects for update v
        XSW(r).slice(k) = XrLrW * (diagmat(ones<VECTYPE>(q,1))- Crk.i()*WtLrW); // O(nq^2)
        
      }
      
      
      // save some objects for posterior of z.
      bCrk = WtLrW +  inv_sympd(Sigma0.slice(k) + Psi0.slice(r));
      WtbSrkW = WtLrW -WtLrW*bCrk.i()*WtLrW; // O(q^3)
      XrbSrkW_uv = XrLrW * (diagmat(ones<VECTYPE>(q,1))- bCrk.i()*WtLrW) -
        (Muv(r)+repmat(Mu0.row(k), n, 1))* WtbSrkW; // O(pq^2+nq^2)
      
      
      Ez(r).slice(k) = XrbSrkW_uv*Sigma0.slice(k) + repmat(Mu0.row(k), n, 1); // complexity: O(nq^2)
      Varz(r).slice(k) = Sigma0.slice(k) - Sigma0.slice(k)*WtbSrkW*Sigma0.slice(k); // O(q^3)
      
      // save some objects for posterior of v
      Ev(r).slice(k) = XrbSrkW_uv *Psi0.slice(r) + Muv(r); // 
      Varv(r).slice(k) = Psi0.slice(r) - Psi0.slice(r)* WtbSrkW*Psi0.slice(r);
      
      // Rprintf("evaluate posterior! \n");
      // save some objects for posterior of z+v
      Varzv(r).slice(k) = Varz(r).slice(k)+ Varv(r).slice(k) - Psi0.slice(r)*WtbSrkW*Sigma0.slice(k) -
        Sigma0.slice(k)*WtbSrkW*Psi0.slice(r);
      Ezv(r) = Ez(r) + Ev(r); // use the addition of two cubes.
      
      // evaluate energy of x, Ux for updating y and caculating responsibility R
      multi_det_SkCpp(Xf(r) - Muv(r)*W0.t(), Lam0.row(r),W0, bCrk, Mu0.row(k), // Use SVD to speed up.
                      logdSk, mSk);
      
      Ux(r).col(k) = -0.5*logdSk  + 0.5 * mSk; // calculate energy by column.
      
    }
  }
  
  // Estimate Y by ICM
  VECTYPE Energy(maxIter_ICM);
  //--------------------------------------------------------------------------------	
  // ICM algrithm to predict Y and V
  //--------------------------------------------------------------------------------
  
  field<MATTYPE> U(r_max);
  Rprintf("predict Y and V! \n");
  for(r= 0; r< r_max; r++){
    
    Energy(0) = INFINITY;
    MATTYPE Uy1, U1;
    VECTYPE U1min;
    uvec y1_u;
    // ivec y = conv_to< ivec >::from(yf(r));
    n = Xf(r).n_rows; // tmp object
    
    // ivec y1;
    for (iter = 1; iter < maxIter_ICM;  ++iter ) {
      
      Uy1 = calYenergy2D_sp(yf(r), Adjf(r), K, alpha, beta(r));
      
      U1 = Uy1 + Ux(r); // log likelihood of (x, y).
      U1min = min(U1, 1);
      y1_u = index_min(U1, 1);
      yf(r) = conv_to< ivec >::from(y1_u) + 1;
      // Rprintf("predict Y and V! \n");
      if(Sp2){
        Muv(r) = get_Vmean(Vf(r),  Adjf_car(r));
        // cout<<"iter="<< iter <<endl;
        // for(i=0; i<n; i++){  // O(nq^2*maxIter_ICM)
        // 
        //   Vf(r).row(i)= (XSW(r).slice(y1_u(i)).row(i) - (Mu0.row(y1_u(i)))* WtSW(r).slice(y1_u(i)) +  Muv(r).row(i)* (Psi0.slice(r)).i())*
        //     inv_sympd(WtSW(r).slice(y1_u(i)) + (Psi0.slice(r)).i());
        // 
        // }
        
        for(k = 0; k<K; ++k){
          uvec index_k = find(y1_u == k);
          int nk = index_k.n_elem;
          // Rprintf("k= %d,  nk = %d ! \n", k, nk);
          if(nk > 0){// if the number of spots whose cluster is k is greater than 0
            Vf(r).rows(index_k) = (XSW(r).slice(k).rows(index_k)- repmat(Mu0.row(k), nk,1) * WtSW(r).slice(k) + Muv(r).rows(index_k)*Psi0.slice(r).i()) *
              inv_sympd(WtSW(r).slice(k) + Psi0.slice(r).i());
          }

        }
        
        
        // Since energy_V is computationally high, we do not caculate it.
        // Energy(iter) = energy_V(X, V, W0, Lam_vec0, Muv, Mu0, Sigma0,Psi0,y, Cki) + sum(Umin); // 
        
      }
      
      //Rprintf("predict Y and V! \n");
      Energy(iter) = sum(U1min);
      if (Energy(iter) - Energy(iter - 1) > 1e-5) {
        // cout << "diff Energy = " << Energy(iter) - Energy(iter - 1)  << endl;
        Rprintf("diff Energy = %4f \n", Energy(iter) - Energy(iter - 1));
        break;
      }
      
      if(!Sp2){
        if (Energy(iter-1) - Energy(iter) < 1e-5)
        {
          Rprintf("ICM Converged at Iteration = %d \n", iter);
          break;
        }
      }
    }
    U(r) = U1;
    
  }
  
  
  // calculate R and pseudo observed loglikelihood
  loglik=0;
  for(r=0; r< r_max; ++r){
    VECTYPE maxA1 = max(-U(r), 1);
    //cout<<"good3"<<endl;
    MATTYPE negU = (-U(r) - repmat(maxA1, 1, K));
    VECTYPE loglik_more_vec = sum(exp(negU),1);
    loglik += sum(log(loglik_more_vec) + maxA1);  
    Rf(r) = exp(negU) / repmat(loglik_more_vec, 1, K);
  }
  
  
  // Rprintf("update beta! \n");
  // update beta: grid search.
  int ng_beta = beta_grid.n_elem;
  VECTYPE objBetaVec(ng_beta);
  if(mix_prop_heter){
    for(r=0; r < r_max; ++r){
      objBetaVec = zeros<VECTYPE>(ng_beta,1);
      for(k = 0; k < ng_beta; ++k){ // each sample has a smoothing parameter
        objBetaVec(k) =  objr_beta(yf(r), Rf(r), Adjf(r), K, alpha, beta_grid(k)); 
      }
      beta(r) = beta_grid(index_max(objBetaVec));
    }
  }else{
    for(k=0; k < ng_beta; ++k){ // all samples have a same smoothing parameter.
      objBetaVec(k) = obj_beta(yf, Rf, Adjf, K, alpha, beta_grid(k));
    }
    beta = ones<VECTYPE>(r_max, 1) * beta_grid(index_max(objBetaVec));
  }
  
  // Rprintf("update Psi! \n");
  // update Psi0
  if(Sp2){
    VECTYPE N;
    for(r=0; r < r_max; ++r){
      N = arma::sum(Rf(r).t(), 1);
      // update Psi0
      MATTYPE Psi_new(q,q, fill::zeros);
      for(k=0; k< K; ++k){
        Psi_new += (trans(Muv(r)- Ev(r).slice(k)) % trans(repmat(Rf(r).col(k), 1, q))) * (Muv(r)- Ev(r).slice(k)) + N(k) * Varv(r).slice(k); //yangyi
      }
      Psi0.slice(r) = Psi_new / Rf(r).n_rows;
      
    }
    
  }
    
    
  
}
    

  
Objidrsc2 idrsc2(const field<MATTYPE>& Xf, const field<sp_mat> Adjf, const field<sp_mat> Adjf_car, field<ivec> yf,
                               const MATTYPE& Mu_int, CUBETYPE Sigma0, const MATTYPE& W_int,
                               const MATTYPE& Lam_int,  const CUBETYPE& Psi_int,
                               const VECTYPE& alpha, VECTYPE beta0, const VECTYPE& beta_grid,
                               const int& maxIter_ICM, const int& maxIter, const float& epsLogLik, const bool& verbose,
                               const bool& homo, const bool& homoClust,
                               const bool& Sigma_diag, const bool& mix_prop_heter, 
                               const bool& Sp2){
      
      
  int r_max = Xf.n_rows;
  int K = Mu_int.n_rows, q= Mu_int.n_cols;
  int r;
      
      
  // Initialize the other iterative parameters
  MATTYPE Mu0(Mu_int), W0(W_int), Lam0(Lam_int);
  CUBETYPE Psi0(Psi_int);
 
  
  // If p is sufficient large, loglik can not be computed.
  // But this can be solved by some programming tricks.
  VECTYPE loglik(maxIter);
  loglik(0) = INT_MIN;
  MATTYPE Nmat(r_max, K);
  
  // Initiailize some objects that will be used in algorithm
  // Posterior of y
  field<MATTYPE> Rf(r_max),  Muv(r_max), Vf(r_max);
  // Posterior of z, v and z+v
  field<CUBETYPE> Ezv(r_max), Varzv(r_max), Ez(r_max), Varz(r_max), Ev(r_max),  Varv(r_max);
  int n;
  for(r = 0; r < r_max; ++r){// initialize the shape of cube in these fields.
    n = Xf(r).n_rows; // tmp object
    Muv(r) = zeros<MATTYPE>(n, q);
    Vf(r) = Muv(r);
    
    Ezv(r) = zeros<CUBETYPE>(n, q, K);
    Ez(r) = Ezv(r);
    Ev(r) = Ezv(r);
    Varzv(r) = zeros<CUBETYPE>(q,q,K);
    Varz(r) = Varzv(r);
    Varv(r) = Varzv(r);
  }
  // pseudo obserbed loglikelihood.
  float loglikVal;
  int k, iter;
  
  Rprintf("variable initialize finish! \n");
  
  
  // begin ICM-EM algorithm
  for(iter = 1; iter < maxIter; iter++){
    //clock_t start1, finish1; // compute the running time of codes.
    //start1 = clock();
    // cache some objects
    // predict y and V, update beta and Psi, and cache some object
    // Rprintf("Satrt ICM step! \n");
    runICM_sp1(Xf, Vf, yf, W0, Lam0, Mu0,Sigma0, Psi0,
                         Adjf, Adjf_car, alpha, beta_grid,beta0, maxIter_ICM, mix_prop_heter, Sp2,
                         Rf, Muv, Ezv, Varzv, Ez, Varz, Ev, Varv, loglikVal);
    loglik(iter) = loglikVal; // obtain the pseudo observed log-likelihood.
    Rprintf("Finish ICM step! \n");
    
    // compute N
    for(r=0; r< r_max; ++r){
      Nmat.row(r) = sum(Rf(r));
    }
    
    // double Q1 = Q_fun1(Xf, Rf, Ez,Ci_ara, W0, Mu0, Sigma0,
    //                                      Lam0, tau0);
    
    
    // update Mu0
    Mu0 = update_Mu1(Rf,  Ez,  Sigma0,  Nmat);
    // cout<<"good!"<< Mu0.row(0)<<endl;
    // double Q2 =  Q_fun1(Xf, Rf, Ez,Ci_ara, W0, Mu0, Sigma0, Lam0, tau0);
    // cout<<"dQ_Mu="<< Q2-Q1<<endl;
    
    // cout<<"tau0="<< tau0.row(1)<<endl;
    // double Q21 = Q_fun1(Xf, Rf, Ez,Ci_ara, W0, Mu0, Sigma0, Lam0, tau0);
    //cout<<"dQ_tau20="<< Q21-Q2<<endl;
    // update Sigma0
    //clock_t start, finish;
    //start = clock();
    Sigma0 = update_Sigma1(Rf, Ez, Varz, Mu0, Nmat, homoClust, Sigma_diag); // 
    // cout<<"Sigma0="<< Sigma0.slice(0).row(3)<<endl;
    //finish = clock();
    //cout << finish - start << "/" << CLOCKS_PER_SEC << " (s) " << endl;
    
    // double Q3 =  Q_fun1(Xf, Rf, Ez,Ci_ara, W0, Mu0, Sigma0, Lam0, tau0);
    // cout<<"dQ_Sigma0="<< Q3-Q21<<endl;
    // update W
    // start = clock();
    W0 = update_W1(Xf, Rf,  Ezv,  Lam0,  Varzv);
    // cout<<"W05="<< W0.row(5)<<endl;
    // double Q4 =  Q_fun1(Xf, Rf, Ez,Ci_ara, W0, Mu0, Sigma0, Lam0, tau0);
    // cout<<"dQ_W="<< Q4-Q3<<endl;
    //finish = clock();
    //cout << finish - start << "/" << CLOCKS_PER_SEC << " (s) " << endl;
    
    // update  Lambda
    // start = clock();
    
    Lam0 = update_Lam1(Rf, Xf, W0, Ezv, Varzv, Nmat, homo);
    // cout<<"Lam1="<< Lam0.row(0).subvec(0,3)<<endl;
    // double Q5 =  Q_fun1(Xf, Rf, Ez,Ci_ara, W0, Mu0, Sigma0, Lam0, tau0);
    // cout<<"dQ_Lambda="<< Q5-Q4<<endl;
    //finish = clock();
    //cout << finish - start << "/" << CLOCKS_PER_SEC << " (s) " << endl;
    
    // calculate loglikelihood
    if(loglik(iter)  - loglik(iter-1)   < -1e-7){
      // perror("The likelihood failed to increase!");
      break;
    }
    
    //finish1 = clock();
    //cout << finish1 - start1 << "/" << CLOCKS_PER_SEC << " (s) " << endl;
    // output algorithm info.
    if(verbose){
      Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", 
              iter +1, loglik(iter), (loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1)));
    }
    if(abs((loglik(iter)  - loglik(iter-1))/ loglik(iter-1)) < epsLogLik) break;
    // if(abs(Q  - tmp_Q) < epsLogLik) break;
    
  }
  
  field<MATTYPE> Ezz(r_max);
  for(r=0; r< r_max; ++r){
    Ezz(r) = zeros<MATTYPE>(Rf(r).n_rows, q);
    for(k=0; k<K; k++){
      
      Ezz(r) +=   Ez(r).slice(k) % repmat(Rf(r).col(k), 1, q);
    }
  }
    
    Objidrsc2 output;
    
    output.yf = yf;
    output.Ezz = Ezz;
    output.Vf = Vf;
    output.Rf = Rf;
    output.beta0 = beta0;
    output.Mu0 = Mu0;
    output.Sigma0 = Sigma0;
    output.Psi0 = Psi0;
    output.W0 = W0;
    output.Lam0 = Lam0;
    output.loglik= loglik(iter-1);
    output.loglik_seq = loglik;
    
    return output;
    
}   
    
    

//' @title
//' idrsc2Cpp
//' @description
//' the main function of ICM-EM algorithm for integrative dimension reduction and clustering.
//'
//' @param Xlist is a r_max-length Rcpp list, including log-normalized gene expression matrix of each sample.
//' @param Adjlist is a r_max-length Rcpp list, including adjoint matrix (obtained by spatial coordinates) of each sample.
//' @param y_intlist is a r_max-length Rcpp list, including initial cluster labels (vector) of each sample.
//' @param Mu_int is a K-by-q matrix, where each row is an intial mean component of gaussian mixture model shared by all samples.
//' @param Sigma_int is a q-by-q-by-K array, where each slice is an initial covariance component of GMM that is sample-specified.
//' @param W_int is a p-by-q matrix that is the intial loading matrix in the probabilistic PCA model.
//' @param Lam_int is a r_max-by-p matrix, where each row is an initial variance vector of error in PPCA model for each sample.
//' @param Psi_int is a q-by-q-by-r_max array, where each slice is an initialized covariance of CAR model for each sample.
//' @param alpha is a K-length vector that equals to zero in Potts model.
//' @param beta_int is a positive real (0~5) that is the initial value of smoothing parameter in Potts model
//' @param beta_grid is a vector that is a grid for searching smoothing parameter beta in the ICM step.
//' @param maxIter_ICM is a positive integer that specifies the maximum iterations in ICM step.
//' @param maxIter is a positive integer that specifies the maximum iterations in this ICM-EM algorithm.
//' @param epsLogLik is a postive real that specifies the relative tolerance value of pseudo observed log-likelihood.
//' @param verbose is a bool indicating whether output the information of algorithm
//' @param homo is a bool indicating whether the homogeneous error covariance in PPCA model is used.
//' @param homoClust is a bool indicating whehter mixture covariance only related to r, or related to r and k.
//' @param Sigma_diag is a bool indicating whether the mixture covariance is a diagonal matrix.
//' @param mix_prop_heter is a bool indicating whether the smoothing parameter beta_r is different among samples or shared by all samples
//' @param Sp2 is a bool indicating whether the CAR model is used.
//' @return A list includes cluster (cluster labels list), hZ (the posterior estimation of Z), hV (predictions of V)
//'
// [[Rcpp::export]]
Rcpp:: List idrsc2Cpp(const Rcpp::List& Xlist, const Rcpp::List& Adjlist, const Rcpp::List& Adjlist_car, const MATTYPE hZ, const arma::imat& ymat,
                               const Rcpp::List& Mu_intList, const Rcpp::List& Sigma_intList, const MATTYPE& W_int,
                               const Rcpp::List& alpha_intList, const float& beta_int, const VECTYPE& beta_grid,
                               const int& maxIter_ICM,const int& maxIter, const float& epsLogLik, const bool& verbose,
                               const bool& homo = false, const bool& homoClust= false,
                               const bool& Sigma_diag=false, const bool& mix_prop_heter = false, const bool& Sp2=true,
                               const int maxK = 16, const int minK = 4, const int& coreNum = 1){
  
  // homo denotes error covariance is a scalar * identity matrix.
  // homoClust denotes mixture covariance related to r and k, or only related to r.
  // Sp2 control whether add CAR model V
  // basic info
  //     cout << "checking point 1" << endl;
  int r,r_max = Xlist.length(); // get the number of data source
  
  // transfer list to field.
  field<MATTYPE> Xf(r_max);
  field<sp_mat> Adjf(r_max);
  field<sp_mat> Adjf_car(r_max);
  VECTYPE beta0(r_max);
  field<imat> yf(r_max);
  for(r=0; r < r_max; ++r){ 
    MATTYPE Xtmp = Xlist[r]; // enforce to become a matrix.
    Xf(r) = Xtmp;
    sp_mat Adjtmp = Adjlist[r];
    Adjf(r) = Adjtmp;
    sp_mat Adjtmp2 = Adjlist_car[r];
    Adjf_car(r) = Adjtmp2;
    beta0(r) = beta_int;
  }    
 
  
    // transfer list to 
    int lengthK = maxK - minK + 1;
    field<MATTYPE> Mu0(lengthK);
    field<CUBETYPE> Sigma0(lengthK);
    field<VECTYPE> alpha0(lengthK);
    
    for (int i = 0; i < lengthK; i++){
        
        MATTYPE tmp_Mu0 = Mu_intList[i];
        Mu0(i) = tmp_Mu0;
        
        CUBETYPE tmp_Sigma0 = Sigma_intList[i];
        Sigma0(i) = tmp_Sigma0;
        
        VECTYPE tmp_alpha0 = alpha_intList[i];
        alpha0(i) = tmp_alpha0;
    }
    
    
    // cout << "checking point 2" << endl;
    // transfer list to 
    int p = Xf(0).n_cols;
    int q = W_int.n_cols;
    MATTYPE Lam0 = zeros<MATTYPE>(r_max, p);
    CUBETYPE Psi0 = zeros<CUBETYPE>(q, q, r_max);
    
    
    uvec nf_cumsum = zeros<uvec>(r_max);
    for (int r = 0; r < r_max; r++){
        nf_cumsum(r) = Xf(r).n_rows;
    }
    nf_cumsum = cumsum(nf_cumsum);
    // cout << "checking point 4" << endl;
    //cout << nf_cumsum << endl;
    
    
    yf(0) = ymat.rows(0,nf_cumsum(0)-1);
    if(!homo){
      Lam0.row(0) = var(Xf(0) - hZ.rows(0,nf_cumsum(0)-1) * trans(W_int), 0, 0);
    }else{
      Lam0.row(0) = ones<ROWVECTYPE>(1, p) * mean(var(Xf(0) - hZ.rows(0,nf_cumsum(0)-1) * trans(W_int), 0, 0));
    }
    
    if(Sp2){
      Psi0.slice(0) = cov(hZ.rows(0,nf_cumsum(0)-1));
    }
    
    // cout << "checking point 5" << endl;
    if(r_max>1){
      for(int r = 1; r < r_max; r++){
        yf(r) = ymat.rows(nf_cumsum(r-1), nf_cumsum(r)-1);
        uvec idx = linspace<uvec>(nf_cumsum[r-1], nf_cumsum[r]-1, nf_cumsum[r] - nf_cumsum[r-1]);
        if(Sp2){
          Psi0.slice(r) = cov(hZ.rows(idx));
        }
        if(!homo){
          Lam0.row(r) = var(Xf(r) - hZ.rows(nf_cumsum(r-1), nf_cumsum(r)-1) * trans(W_int), 0, 0);
        }else{
          Lam0.row(r) = ones<ROWVECTYPE>(1, p) * mean(var(Xf(r) - hZ.rows(nf_cumsum(r-1), nf_cumsum(r)-1) * trans(W_int), 0, 0));
        }
      }
    }
    
    
    if(lengthK==1){
      field<ivec> yf_int(r_max);
      for (int r = 0; r < r_max; r++){
        ivec tmp_yf = yf(r).col(0);
        yf_int(r) = tmp_yf;
      }
      MATTYPE Mu_int = Mu0(0);
      MATTYPE Lam_int = Lam0;
      CUBETYPE Psi_int = Psi0;   
      VECTYPE alpha_int = alpha0(0);
      CUBETYPE Sigma_int = Sigma0(0);
      Objidrsc2 output = idrsc2(Xf, Adjf, Adjf_car, yf_int,
                         Mu_int, Sigma_int,  W_int,
                         Lam_int, Psi_int,
                         alpha_int,  beta0,  beta_grid,
                         maxIter_ICM,  maxIter, epsLogLik, verbose,
                         homo, homoClust, Sigma_diag, mix_prop_heter, Sp2); 
      
      List Objidrsc2Rcpp(lengthK);
      
        // output return value
        Objidrsc2Rcpp[0] = List::create(
          Rcpp::Named("cluster") = output.yf,
          Rcpp::Named("hZ") = output.Ezz,
          Rcpp::Named("hV") = output.Vf,
          Rcpp::Named("Rf") = output.Rf,
          Rcpp::Named("beta") = output.beta0,
          Rcpp::Named("Mu") = output.Mu0,
          Rcpp::Named("Sigma") = output.Sigma0,
          Rcpp::Named("Psi") = output.Psi0,
          Rcpp::Named("W") = output.W0,
          Rcpp::Named("Lam") = output.Lam0,
          Rcpp::Named("loglik") = output.loglik,
          Rcpp::Named("loglik_seq") = output.loglik_seq);
      
      
      return(Objidrsc2Rcpp);
      
      
    }else{
      //set parallel structure object
      par_iDRSC2 parObj(Xf, Adjf, Adjf_car, yf,  Mu0, Sigma0, W_int,  Lam0, Psi0,
                        alpha0, beta0, beta_grid, maxIter_ICM, maxIter, epsLogLik, verbose,
                        homo, homoClust, Sigma_diag, mix_prop_heter, Sp2, maxK, minK);
      
      const int n_thread = coreNum;
      std::vector<std::thread> threads(n_thread);
      
      for (int i_thread = 0; i_thread < n_thread; i_thread++){
        threads[i_thread] = std::thread(&par_iDRSC2::update_by_thread_idrsc2, &parObj, i_thread);
      }
      for (int i = 0; i < n_thread; i++){
        threads[i].join();
      }
      
      
      List Objidrsc2Rcpp(maxK-minK+1);
      
      
      for (int k = 0; k < maxK - minK + 1; k++){
        // output return value
        Objidrsc2Rcpp[k] = List::create(
          Rcpp::Named("cluster") = parObj.output[k].yf,
          Rcpp::Named("hZ") = parObj.output[k].Ezz,
          Rcpp::Named("hV") = parObj.output[k].Vf,
          Rcpp::Named("Rf") = parObj.output[k].Rf,
          Rcpp::Named("beta") = parObj.output[k].beta0,
          Rcpp::Named("Mu") = parObj.output[k].Mu0,
          Rcpp::Named("Sigma") = parObj.output[k].Sigma0,
          Rcpp::Named("Psi") = parObj.output[k].Psi0,
          Rcpp::Named("W") = parObj.output[k].W0,
          Rcpp::Named("Lam") = parObj.output[k].Lam0,
          Rcpp::Named("loglik") = parObj.output[k].loglik,
          Rcpp::Named("loglik_seq") = parObj.output[k].loglik_seq);
      }
      
      
      return(Objidrsc2Rcpp);
    }
    
    
    
}








