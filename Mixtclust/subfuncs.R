EM_iter <- function(oldpars, y, sigma.constr = "VEE", df.constr = FALSE) {
  pis <- oldpars$pi
  nus <- oldpars$nu
  mus <- oldpars$mu
  Sigmas <- oldpars$Sigmas
  lambdas <- oldpars$lambdas
  C <- oldpars$C
  # 1st cycle: (pi, mu, Sigma)
  z <- up_Z(y, mus, Sigmas, nus, pis)
  u <- up_U(y, mus, Sigmas, nus, pis)
  pis <- up_pi(z)
  mus <- up_mu(y, z, u)
  # 2ns cycle: nu, Sigma
  z <- up_Z(y, mus, Sigmas, nus, pis)
  u <- up_U(y, mus, Sigmas, nus, pis)
  p <- ncol(mus)
  nus <- approx_nu(z, u, nus, p, df.constr)
  nk <- up_nk(y, z)
  Sk <- up_Sk(y, u, z, nk, mus)
  lambdas <- drop(up_lambdak(Sk, C))
  C <- up_C(Sk, nk, lambdas)
  Sigmas <- up_Sigma(C, lambdas)
  newpars <- list(pis, nus, mus, Sigmas, lambdas, C)
  names(newpars) <- c("pi", "nu", "mu", "Sigmas", "labmdas", "C")
  return(newpars)
}

getlambda = function(sigmas, K, p){
  lambdas <- rep(0, K)
  for(k in 1:K){
    lambdas[k] <- det(sigmas[,,k])^(1/p)
  }
  lambdas
}

getC = function(sigmas, K, p){
  lsigma <- lapply(1:K, function(x) sigmas[,,x] / det(sigmas[,,x])^(1/p))
  C <- matrix(0L, nr = p, nc = p)
  for(k in 1:K){
    C = C + lsigma[[k]]
  }
  C = C / K
}

skewtmix_cpp <- function(X,
                         initial.values, df.constr = FALSE, sigma.constr = "VEE",
                         model = "EEV", # constraints on df, lambda, Sigma.
                         max.iter = 1000,
                         tol = 1e-3) {
  K <- length(initial.values$pi)
  J <- dim(X)[1]; p <- dim(X)[2]; n <- dim(X)[3]
  old <- initial.values
  old$lambdas = getlambda(initial.values$Sigma, K, p)
  old$C = getC(initial.values$Sigma, K, p)
  parnames <- c("pi", "nu", "mu", "Sigmas", "labmdas", "C")
  names(old) <- parnames
  del <- 1e6 # Initialize change in loglik holder.
  iter <- 0
  LLs <- rep(NA, max.iter) # Store loglikelihood at each iteration.
  converged <- FALSE
  while (del > tol) {
    iter <- iter + 1
    if (iter > max.iter) {
      break
    }
    new <- EM_iter(old, X, sigma.constr, df.constr)
    names(new) <- parnames
    Zs <- up_Z(X, new$mu, new$Sigmas, new$nu, new$pi)
    oldLik <- loglik(X, Zs, old$mu, old$Sigmas, old$pi, old$nu)
    newLik <- loglik(X, Zs, new$mu, new$Sigmas, new$pi, new$nu)
    del <- newLik - oldLik # Change in loglikelihoods.
    LLs[iter] <- newLik
    old <- new
  }
  if(del < tol) converged <- TRUE
  #Zs <- up_Z(X, new$mu, new$Sigma, new$nu, new$pi)
  classification <- apply(Zs, 1, which.max)
  if(sigma.constr == "VEE"){
    sigmanpar = K - 1 + p*(p+1)/2
  }
  npar <- (K-1) + K*p + sigmanpar + ifelse(df.constr, 1, K)
  BIC <- - 2*LLs[iter] + npar*log(n)
  res <- list(new, iter, Zs, classification, LLs[1:iter], BIC, converged)
  names(res) <- c("estimates", "iterations", "Zs", "classification", "loglik", "bic", "converged")
  return(res)
}


