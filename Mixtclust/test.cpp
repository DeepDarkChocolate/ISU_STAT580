#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppDist.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

arma::vec mahadist(const arma::mat& x, const arma::vec& mu,
                   const arma::mat& S){
  arma::uword n = x.n_rows, m = x.n_cols;
  arma::mat S_inv = S.i();
  arma::vec result(n);
  arma::vec X(m);
  for ( arma::uword i = 0; i < n; ++i ) {
    X = x.row(i).t() - mu;
    arma::mat tmp = X.t() * S_inv * X;
    result(i) = tmp(0,0);
  }
  return result;

}

// [[Rcpp::export]]
arma::mat up_Z(arma::cube y, arma::mat mus, arma::cube sigmas, arma::vec nus, arma::vec pis){
  int J = y.n_rows, n = y.n_slices, K = mus.n_rows;
  arma::mat ans(n, K);
  for(int k = 0; k < K; k++){
	arma::vec tmp(n, arma::fill::ones);
	for(int i = 0; i < n; i++){
	  arma::vec tmpmat(J);
	  tmpmat = dmvt(y.slice(i), mus.row(k).t(), sigmas.slice(k), nus(k), false);
	  //tmpmat = dmvtInt(y.slice(i), mus.row(k).t(), sigmas.slice(k), nus(k), false);
	  //tmpmat = dMVT(y.slice(i), mus.row(k).t(), sigmas.slice(k), nus(k), false, false);
	  //if(i == 0 & k == 0){
		//Rcout << "y.slice(i)" << y.slice(i);
		//Rcout << "mus.row(k).t()" << mus.row(k).t();
		//Rcout << "sigmas.slice(k)" <<sigmas.slice(k);
		//Rcout << "nus(k)" << nus(k);
		//Rcout << "tmpmat" << tmpmat;
	  //}
	  for(int j = 0; j < J; j++){
		if(R_IsNaN(tmpmat(j))) tmpmat(j) = 1;
	  }
	  tmp(i) = prod(tmpmat);
	}
	ans.col(k) = tmp * pis(k);
  }
  for (int i=0; i<n; i++) {
	double rowsum = 0.0;
	for (int k=0; k<K; k++) {
	  rowsum += ans(i,k);
	}
	ans.row(i) = ans.row(i) / rowsum;
  }
  return ans;
}

// [[Rcpp::export]]
arma::cube up_U(arma::cube y, arma::mat mus, arma::cube sigmas, arma::vec nus, arma::vec pis){
  int J = y.n_rows, p = y.n_cols, n = y.n_slices, K= mus.n_rows;
  arma::cube ans(n, J, K);
  for(int i = 0; i < n; i++){
	arma::mat y_tmp = y.slice(i);
	for(int k = 0; k < K; k++){
	  arma::vec maha = mahadist(y_tmp, mus.row(k).t(), sigmas.slice(k));
	  //arma::vec maha = mahalanobis(y_tmp, mus.row(k).t(), sigmas.slice(k), false);
	  //arma::vec maha = mahaInt(y_tmp, mus.row(k).t(), sigmas.slice(k), false);
	  for(int j = 0; j < J; j++){
		ans(i, j, k) = (nus(k) + p) / (nus(k) + maha(j));
	  }
	}
  }
  return ans;
}

// [[Rcpp::export]]
arma::vec up_pi(arma::mat z) {
  int n = z.n_rows, K = z.n_cols;
  arma::vec ans(K); 
  ans.zeros();
  for (int k=0; k<K; k++) {
	for (int i=0; i<n; i++) {
	  ans(k) += z(i,k);
	}
	ans(k) = ans(k) / n;
  }
  return ans;
}

// [[Rcpp::export]]
arma::mat up_mu(arma::cube y, arma::mat z, arma::cube u){
  int J = y.n_rows, p = y.n_cols, n = y.n_slices, K= z.n_cols;
  arma::mat num(K, p), den(K, p);
  num.zeros(); den.zeros();
  for (int k = 0; k < K; k++){
	for (int i = 0; i < n; i++){
	  for(int j = 0; j < J; j++){
		double utmp;
		arma::rowvec ytmp(p);
		if(R_IsNA(u(i, j, k))){
		  utmp = 0;
		  ytmp.zeros();
		}
		else{
		  utmp = u(i, j, k);
		  ytmp = y.slice(i).row(j);
		}
		num.row(k) += z(i, k) * utmp * ytmp;
		den.row(k) += z(i, k) * utmp;
	  }
	}
  }
  arma::mat ans = num / den;
  return ans;
}

// [[Rcpp::export]]
arma::vec up_nk(arma::cube y, arma::mat z){
  int J = y.n_rows, n = y.n_slices, K = z.n_cols;
  arma::vec ans(K);
  int J_i; double sum;
  for(int k = 0; k < K; k++){
	sum = 0;
	for(int i = 0; i < n; i++){
	  J_i = 0;
	  for(int j = 0; j < J; j++){
		if(!y.slice(i).row(j).has_nan()) J_i += 1;
	  }
	  sum += z(i, k) * J_i;
	}
	ans(k) = sum;
  }
  return ans;
}
			
// [[Rcpp::export]]
arma::cube up_Sk(arma::cube y, arma::cube u, arma::mat z, arma::vec nk, arma::mat mus){
  int J = y.n_rows, p = y.n_cols, n = y.n_slices, K = mus.n_rows;
  arma::cube ans(p, p, K);
  arma::mat sum(p,p);
  for(int k = 0; k < K; k++){
	arma::mat sum(p,p, arma::fill::zeros);
	for(int i = 0; i < n; i++){
	  arma::mat sumi(p, p, arma::fill::zeros);
	  for(int j = 0; j <J; j++){
		arma::rowvec tmp(p, arma::fill::zeros);
		tmp = y.slice(i).row(j) - mus.row(k);
		if(!R_IsNA(u(i,j,k))) sumi = sumi + u(i,j,k) * tmp.t() * tmp;
	  }
	  sum = sum + z(i, k) * sumi;
	}
	ans.slice(k) = sum / nk(k);
  }
  return ans;
}

// [[Rcpp::export]]
arma::vec up_lambdak(arma::cube Sk, arma::mat C){
  int K = Sk.n_slices, p = Sk.n_cols;
  arma::vec ans(K);
  for(int k = 0; k < K; k++){
    arma::mat tmp = solve(C, Sk.slice(k));
    ans(k) = trace(tmp) / p;
  }
  return ans;
}

// [[Rcpp::export]]
arma::mat up_C(arma::cube Sk, arma::vec nk, arma::vec lambdak){
  int K = Sk.n_slices, p = Sk.n_cols;
  arma::mat ans(p, p);
  for(int k = 0; k < K; k++){
    ans = ans + nk(k) / lambdak(k) * Sk.slice(k);
  }
  double tmp = pow(det(ans), (double)1 / (double)p);
  ans = ans / tmp;
  return ans;
}

// [[Rcpp::export]]
arma::cube up_Sigma(arma::mat C, arma::vec lambdak){
  int K = lambdak.n_elem, p = C.n_cols;
  arma::cube ans(p, p, K);
  for(int k = 0; k < K; k++){
    ans.slice(k) = lambdak(k) * C;
  }
  return ans;
}

// [[Rcpp::export]]
arma::vec approx_nu(arma::mat z, arma::cube u, arma::vec nus, int p, bool constraint) {
  int K = z.n_cols, n = z.n_rows, J = u.n_cols;
  arma::vec ans(K);
  if(constraint == false){
    for (int k = 0; k < K; k++){
      double zsum = 0, asum = 0;
      for (int i = 0; i < n; i++){
        for(int j = 0; j < J; j++){
          if(!R_IsNA(u(i, j, k))){
            asum += z(i, k) * ( log(u(i, j, k)) - u(i, j, k));
            zsum += z(i, k);
          }
        }
      }
      double cc = 1.0 + R::digamma((nus(k) + p)/2.0) - log((nus(k) + p)/2.0) + (asum/zsum);
      cc = -cc;
      ans(k) = (-exp(cc)+2.0*(exp(cc))*(exp(R::digamma(nus(k)/2.0)) - (nus(k)/2.0-0.5)))/(1.0-exp(cc));
      
    }
  }
  else{
    double zsum = 0, asum = 0;
    for (int k = 0; k < K; k++){
      for (int i = 0; i < n; i++){
        for(int j = 0; j < J; j++){
          if(!R_IsNA(u(i, j, k))){
            asum += z(i, k) * ( log(u(i, j, k)) - u(i, j, k));
            zsum += z(i, k);
          }
        }
      }
    }
    double cc = 1.0 + R::digamma((nus(1) + p)/2.0) - log((nus(1) + p)/2.0) + (asum/zsum);
    cc = -cc;
    double tmp = (-exp(cc)+2.0*(exp(cc))*(exp(R::digamma(nus(1)/2.0)) - (nus(1)/2.0-0.5)))/(1.0-exp(cc));
    for (int k = 0; k < K; k++){ 
      ans(k) = tmp;
    }
  }
  return ans;
}

// [[Rcpp::export]]
double loglik(arma::cube y, arma::mat z, arma::mat mus, arma::cube sigmas,
   	arma::vec pis, arma::vec nus){
  int J = y.n_rows, n = y.n_slices, K = z.n_cols;
  double ans;
  for(int i = 0; i < n; i++){
	for(int k = 0; k < K; k++){
	  ans += z(i, k) * log(pis(k));
	  for(int j = 0; j < J; j++){
		arma::vec tmpmat(J);
		tmpmat = dmvt(y.slice(i), mus.row(k).t(), sigmas.slice(k), nus(k), true);
		//tmpmat = dmvtInt(y.slice(i), mus.row(k).t(), sigmas.slice(k), nus(k), true);
		if(!R_IsNaN(tmpmat(j))) ans += tmpmat(j) * z(i, k);
	  }
	}
  }
  return ans;
}
