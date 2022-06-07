// This scripts includes the all common functions used in multInteCluster.cpp and SISDRclusterv2.cpp
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include<ctime>
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
//' Evaluate class label of neighbor matrix.

sp_mat get_spNbs(ivec y, const sp_mat& Adj) {   // ivec是索引型向量
  // row is for pixel.
  //output a sparse matrix, i-th row contains labels of neighbor_i. 
  // Make const iterator
  arma::sp_mat::const_iterator start = Adj.begin(); //构造一个sp_mat的常数迭代器,常数迭代器只可读，不可写，实现对矩阵按照列对每个非零元素进行访问。
  // The nonzero elements of i-th column in Adj is the neighbors of spot i.
  
  // Calculate number of nonzero points
  //int n = std::distance(start, end);
  int n = Adj.n_nonzero; // 计算Adj的所有非零元的个数
  //cout << "n=" << n << endl;
  //cout << "n=" << Adj.n_nonzero << endl;
  
  sp_mat spNbs(y.n_elem, y.n_elem);    // neiborhood state matrix, matched with Adj.
  
  
  arma::sp_mat::const_iterator it = start; // Note spNbs is not a symmetric matrix, the nonzero in i-th row is the class label of sample i.
  for(int i = 0; i < n; ++i)
  {
    //temp(0) = it.row();
    //temp(1) = it.col();
    spNbs(it.row(), it.col()) = y(it.row()); // it只自加非零元个数次，得到每个i对应的邻居的状态
    ++it; // increment
  }
  
  return spNbs; // return the class label of neighbor matrix, i-th column is the neighbor label of sample i
}

// Calculate the value of energy function of y, which is equal to negative logliklihood up to a constant
MATTYPE calYenergy2D_sp(const arma::ivec& y, const arma::sp_mat& Adj, int K, const VECTYPE alpha, const float beta)	{
  
  int n = y.n_rows;
  arma::sp_mat spNbs = get_spNbs(y, Adj); // transform spNbs to iterate by column.
  MATTYPE Uy(n, K);
  double n_sameS;
  int i, k, nn;
  for (k = 0; k < K; k++)
  {
    for (i = 0; i < n; i++)
    {
      arma::sp_mat col(spNbs.col(i)); // the class label of neighbors of i-th sample.
      n_sameS = 0;
      
      nn = col.n_nonzero; // the number of neighbors of i-th sample
      for (arma::sp_mat::iterator j = col.begin(); j != col.end(); ++j) {
        n_sameS += ((*j) == (k+1));
        
      }
      Uy(i, k) = alpha(k) + beta * (nn - (float)n_sameS)/2;
      
      
    }
  }
  
  MATTYPE C_mat = normalise(exp(-Uy), 1, 1); // pseudo likelihood of Y.
  Uy = -log(C_mat); // normalized Uy, this is the energy of y.
  return Uy;
  
}

// When smoothing parameter is shared by all sample, this function is used to evaluate the objective function of beta.
float obj_beta(const field<ivec>& yf, const field<MATTYPE>& Rf, 
                  const arma::field<sp_mat>& Adjf, int K, const VECTYPE alpha, const float beta)	{
    int r, r_max = yf.n_elem;
    double objval = 0;
    for(r=0; r< r_max; ++r){
      MATTYPE Uy1 = calYenergy2D_sp(yf(r), Adjf(r), K, alpha, beta); // Uy was normalized, so there is no need to normalized Uy. 
      objval += -accu(Rf(r) % Uy1);
    }
    
    return objval;
}
// 
// When smoothing parameter is sample-specified, this function is used to evaluate the objective function of beta_r for each sample.
float objr_beta(const ivec& y, const MATTYPE& R, 
                const sp_mat& Adj, int K, const VECTYPE alpha, const float beta)	{
  
  float objval = 0;
  MATTYPE Uy1 = calYenergy2D_sp(y, Adj, K, alpha, beta); // Uy was normalized, so there is no need to normalized Uy. 
  objval = -accu(R % Uy1);
  
  
  return objval;
}
  

// Evaluate the log-determinant of  $\bar S_{rk}$   and (x_{ri}-W(\mu_k +\tau_r + \mu_{v_{ri}}))^T * \bar S_{rk} *(x_{ri}-W(\mu_k +\tau_r + \mu_{v_{ri}}))
// in an efficient method used in spatial data
void multi_det_SkCpp(const MATTYPE& X, const ROWVECTYPE& Lam_vec0, const MATTYPE& W0, const MATTYPE& Ck, 
                       const ROWVECTYPE Muk, 
                       float& logdSk, VECTYPE& mSk){
    //int p = X.n_cols;
    int n = X.n_rows;
    // int p = X.n_cols;
    // // mSk = zeros(n);
    // S2k = zeros(p);
    // dSk = 0;
    
    MATTYPE WC12,  tmp2;
    VECTYPE tmp1, s, tmp3;
    MATTYPE U, V, X_tk;
    
    svd(U, s, V, Ck.i());
    
    WC12 = W0 * (U * diagmat(sqrt(s)));
    WC12 = diagmat(1/sqrt(Lam_vec0)) * WC12;  // change to sparse matrix multiplication.
    VECTYPE d = svd(WC12);
    //dSk = arma::as_scalar(prod(1- d % d)) / prod(Lam_vec0);
    logdSk = accu(log(1 - d%d)) - accu(log(Lam_vec0));
    X_tk = (X - repmat(Muk* W0.t(), n, 1)) * diagmat(1/sqrt(Lam_vec0)) ;  // change to sparse matrix multiplication.
    tmp1 = sum(X_tk % X_tk, 1);
    tmp2 = X_tk * WC12;
    tmp3 = sum(tmp2 % tmp2, 1);
    mSk = tmp1 - tmp3;
}

// used for non-spatial data
void multi_det_SkCpp2(const MATTYPE& X, const ROWVECTYPE& Lam_vec0,
                        const MATTYPE& W0, const MATTYPE& Ck, 
                        const ROWVECTYPE Muk, const MATTYPE& Sigmak,
                        float& logdSk, VECTYPE& mSk){
    //int p = X.n_cols;
    int n = X.n_rows;
    
    
    MATTYPE WC12,  tmp2;
    VECTYPE tmp1, s, tmp3;
    MATTYPE U, V, X_tk;
    
    // method B: use SVD to compute |Srk.i()|
    svd(U, s, V, Sigmak);
    WC12 = W0 * (U * diagmat(sqrt(s)));
    WC12 = diagmat(1.0/sqrt(Lam_vec0)) * WC12;  // change to sparse matrix multiplication.
    VECTYPE d = svd(WC12);
    logdSk = -accu(log(1 +  d%d)) - accu(log(Lam_vec0));
    // method A: directly compuate log|Srki|
    // mat Srki = W0*Sigmak*W0.t()+ sp_mat(diagmat(Lam_vec0));
    // s = eig_sym(Srki);
    // //svd(U,s, V, Srki);
    // logdSk = -accu(log(s));
    // WC12 = U* diagmat(sqrt(1.0/ s));
    // X_tk = (X - repmat(Muk* W0.t(), n, 1));  // change to sparse matrix multiplication.
    // tmp2 = X_tk * WC12;
    // mSk  = sum(tmp2 % tmp2, 1);
    svd(U, s, V, Ck.i());
    WC12 = W0 * (U * diagmat(sqrt(s)));
    WC12 = diagmat(1.0/sqrt(Lam_vec0)) * WC12;  
    X_tk = (X - repmat(Muk* W0.t(), n, 1)) * diagmat(1/sqrt(Lam_vec0)) ;  // change to sparse matrix multiplication.
    tmp1 = sum(X_tk % X_tk, 1);
    tmp2 = X_tk * WC12;
    tmp3 = sum(tmp2 % tmp2, 1);
    mSk = tmp1 - tmp3;
    
  }
  
  
// This function is to compute diag(W0* Cki^(-1) *W0^T)  \in R^p  efficiently  
VECTYPE decomp(const MATTYPE& Cki, const MATTYPE& W0){
    VECTYPE s, tmp1;
    MATTYPE U, V, WC12;
    svd(U, s, V, Cki);
    WC12 = W0 * (U * diagmat(sqrt(s)));
    tmp1 = sum(WC12 % WC12, 1);
    return tmp1;
}  

