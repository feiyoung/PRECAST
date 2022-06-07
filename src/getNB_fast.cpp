#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "idrsc2.h"

using namespace Rcpp;
using namespace arma;
using namespace std;



// [[Rcpp::export]]
arma::sp_umat get_fixedNumber_neighbors(const MATTYPE  x, int  number=6)	{
  // Get the fixed number of neighbors for each spot.
  int i,j, jj, N = x.n_rows;
  int number4 = 4* number+1;
  arma::sp_umat D(N, N);
  //float dis;
  uvec idx, idx2, idx_all;
  for (j = 0; j < N; ++j)
  {   
    //Rprintf("good R!\n");
    idx_all = sort_index(abs(x(j,0) - x.col(0)));
    idx = idx_all.subvec(0, number4);
    //Rprintf("good R %d !\n",j);
    idx2 = find(idx != j);
    int p = idx2.n_elem;
    vec dis(p);
    
    for (i = 0; i < p; ++i)
    {
      dis(i) = norm(x.row(idx(idx2(i))) - x.row(j), 2);
    }
    uvec idx3 = sort_index(dis);
    for(i=0; i<number; ++i){
      jj = idx(idx2(idx3(i)));
      // D(jj,j) = 1;
      D(jj,j) = 1;
    }
  }
  return D;
}


//' @title
//' getneighborhood_fast
//' @description
//' an efficient function to find the neighborhood based on the matrix of position and a pre-defined cutoff
//'
//' @param x is a n-by-2 matrix of position.
//' @param radius is a threashold of Euclidean distance to decide whether a spot is an neighborhood of another spot. For example, if the Euclidean distance between spot A and B is less than cutoff, then A is taken as the neighbourhood of B. 
//' @return A sparse matrix containing the neighbourhood
//'
//' @export
// [[Rcpp::export]]
arma::sp_umat getneighborhood_fast(const MATTYPE x, float radius)	{
  int N = x.n_rows;
  arma::sp_umat D(N, N);
  float dis;
  uvec idx, idx2;
  for (int j = 0; j < N-1; ++j)
  {    
    idx = find(abs(x(j,0) - x.col(0))<radius); 
    idx2 = find(idx>j);
    int p = idx2.n_elem;
    for (int i = 0; i < p; ++i)
    {
      dis = norm(x.row(idx(idx2(i))) - x.row(j), 2);
      if (dis < radius){
        D(idx(idx2(i)),j) = 1;
        D(j,idx(idx2(i))) = 1;
      }
    }
  }
  return D;
}


//' Calculate column-wise or row-wise mean
//' @param sp_data A sparse matrix
//' @param rowMeans A boolean value, whether to calculate row-wise mean
//' @return A n x 1 or p x 1 matrix 
//' @export
// [[Rcpp::export]]
arma::vec sp_means_Rcpp(arma::sp_mat sp_data, bool rowMeans = false) {
  
  arma::sp_mat norm_col_sums;
  mat tmp_mat;
  mat tmp_mat2;
  
  if (rowMeans) {
    norm_col_sums = arma::mean(sp_data, 1);
    tmp_mat = arma::conv_to< mat >::from(norm_col_sums.col(0));
    tmp_mat2 =  arma::conv_to< mat >::from(tmp_mat);
  }
  else {
    
    norm_col_sums = arma::mean(sp_data, 0);
    
    tmp_mat = arma::conv_to< mat >::from(norm_col_sums.row(0).t());
    tmp_mat2 =  arma::conv_to< mat >::from(tmp_mat);
  }
  
  return tmp_mat2;
}




//' Calculate column-wise or row-wise sum
//' @param sp_data A sparse matrix
//' @param rowSums A boolean value, whether to calculate row-wise sum
//' @return A n x 1 or p x 1 matrix 
//' @export
// [[Rcpp::export]]
arma::vec sp_sums_Rcpp(arma::sp_mat sp_data, bool rowSums = false) {
  mat tmp_mat;
  mat tmp_mat2;
  
  arma::sp_mat norm_col_sums;
  
  if (rowSums) {
    norm_col_sums = arma::sum(sp_data, 1);
    tmp_mat = arma::conv_to< mat >::from(norm_col_sums.col(0));
    tmp_mat2 =  arma::conv_to< mat >::from(tmp_mat);
  }
  else {
    norm_col_sums = arma::sum(sp_data, 0);
    tmp_mat = arma::conv_to< mat >::from(norm_col_sums.row(0).t());
    tmp_mat2 =  arma::conv_to< mat >::from(tmp_mat);
      
  }
  
  return tmp_mat2;
}