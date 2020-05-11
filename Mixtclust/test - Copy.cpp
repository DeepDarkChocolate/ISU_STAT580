#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppDist.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
/* https://github.com/mfasiolo/mvnfast */
/*
arma::vec mahaInt(arma::mat X,  
                  arma::vec mu,  
                  arma::mat sigma,
                  bool isChol = false)
{
  using namespace arma;
  
  // Some sanity checks 
  if(mu.n_elem != sigma.n_cols) Rcpp::stop("The mean vector has a different dimensions from the covariance matrix.");
  if(X.n_cols != sigma.n_cols)  Rcpp::stop("The number of columns of X is different from the dimension of the covariance matrix.");
  
  // Calculate transposed cholesky dec. unless sigma is alread a cholesky dec.
  mat cholDec;
  if( isChol == false ) {
    cholDec = trimatl(chol(symmatu(sigma)).t());
  }
  else{
    cholDec = trimatl(sigma.t()); 
    if(any(cholDec.diag() <= 0.0))  Rcpp::stop("The supplied cholesky decomposition has values <= 0.0 on the main diagonal.");
  }
  
  vec D = cholDec.diag();
  
  vec out(X.n_rows);
  
  // Declaring some private variables
  uint32_t d = X.n_cols;
  uint32_t n = X.n_rows;
  
  vec tmp(d);  
  
  double acc;
  uint32_t icol, irow, ii;  
  
  // For each of the "n" random vectors, forwardsolve the corresponding linear system.
  // Forwardsolve because I'm using the lower triangle Cholesky.

  for(icol = 0; icol < n; icol++)
  {
    
    for(irow = 0; irow < d; irow++)
    {
      acc = 0.0;
      
      for(ii = 0; ii < irow; ii++) acc += tmp.at(ii) * cholDec.at(irow, ii);
      
      tmp.at(irow) = ( X.at(icol, irow) - mu.at(irow) - acc ) / D.at(irow);
    }
    
    out.at(icol) = sum(square(tmp)); 
  }

return out;
}

arma::vec dmvtInt( arma::mat X, arma::vec mu, arma::mat Sigma, double df, bool log)
{
  using namespace arma;
  
  unsigned int d = X.n_cols;
  arma::mat cholDec = chol(symmatu(Sigma));
  
  vec out = mahaInt(X, mu, cholDec, true);
  
  if( df <= 0.0 ){ // Multivariate normal density OR ...
    
    out = - 0.5 * out - ( (d / 2.0) * std::log(2.0 * M_PI) + sum(arma::log(cholDec.diag())) );
    
  } else { // ... multivariate Student-t density
  
  uint32_t ii;  
  uint32_t n = X.n_rows;  
  double logDet = sum(arma::log(cholDec.diag())); 
  double c = lgamma((d + df)/2.0) - (lgamma(df/2.0) + logDet + d/2.0 * std::log(M_PI * df));

  for(ii = 0; ii < n; ii++)
  {
    out.at(ii) = c - 0.5 * (df + d) * log1p(out.at(ii)/df);
  }

  }
  
  if (log == false) out = exp(out);
  
  return( out );
}
*/

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

/*
arma::vec mahalanobis(arma::mat x, arma::vec mu, arma::mat sigma, bool ischol = false) {
    if (mu.n_elem != sigma.n_cols) {
	      Rcpp::stop("The supplied mean vector and covariance matrix have incompatible dimensions.");
		    }
	  if (x.n_cols != sigma.n_cols)  {
		    Rcpp::stop("The supplied data matrix and covariance matrix have incompatible dimensions.");
			  }
	    int n = x.n_rows, p = x.n_cols;
		  arma::mat A(p,p);
		    if (!ischol) {
			      arma::mat Atmp(p,p);
				      bool success = arma::chol(Atmp, sigma);
					      if (!success) {
							      Atmp = arma::chol(sigma);
								      }
						      A = arma::trimatl(Atmp.t());
							    } else {
								      A = arma::trimatl(sigma.t());
									    }
			  arma::vec D = A.diag();
			    arma::vec ans(n), tmp(p);
				  double s;
				    for (int i = 0; i < n; i++) {
					      for (int j = 0; j < p; j++) {
							      s = 0.0;
								        for (int k = 0; k < j; k++) {
										          s += tmp(k) * A(j, k);
												        }
										      tmp(j) = ( x(i, j) - mu(j) - s ) / D(j);
											      }
						      ans.at(i) = sum(square(tmp)); 
							    }
					  return ans;
}
*/

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

arma::vec getEigenValues(arma::mat M) {
return arma::eig_sym(M);
}
