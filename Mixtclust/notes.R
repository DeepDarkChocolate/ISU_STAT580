Z = up_Z(y, mus, sigmas, nus, pis)
Z2 <- matrix(0L, n, K)
for(i in 1:n){
  for(k in 1:K){
    tmp <- sapply(1:J, function(x) dmvt(y[x,,i], mus[k,], sigmas[,,k], nus[k], log = FALSE))
    Z2[i, k] <- pis[k] * prod(tmp)
  }
}
all.equal(Z, Z2)

U = up_U(y, mus, sigmas, nus, pis)
U2 <- array(0L, dim = c(n, J, K))
for(i in 1:n){
  for(k in 1:K){
    for(j in 1:J){
      U2[i,j,k] <- (nus[k] + p) / (nus[k] + mahalanobis(y[j,,i], mus[k,], sigmas[,,k]))
    }
  }
}

all.equal(U, U2)

mu <- up_mu(y, Z, U)
mu2 <- matrix(0L, nr = K, nc = p)
for(k in 1:K){
  num <- rep(0L, p)
  denum <- 0
  for(i in 1:n){
    for(j in 1:J){
      num<- num+Z[i, k] * U[i,j,k] * y[j,,i]
      denum<- denum+Z[i, k] * U[i,j,k]
    }
  }
  mu2[k,] <- num / denum
}
all.equal(mu, mu2)

nk <- up_nk(y, Z)
Ji <- apply(y, 3, function(x){sum(!is.na(x[,1]))})
nk2 <- t(Z) %*% Ji
all.equal(nk, nk2)

Sk <- up_Sk(y, U, Z, nk, mus)

sum(diag(solve(Cinit, Sk[,,1])))

lambdak <- drop(up_lambdak(Sk, Cinit))
lambdak2 <- rep(0, K)
for(k in 1:K){
  lambdak2[k] <- sum(diag(solve(Cinit, Sk[,,k])))/p # det(sigmas[,,k])^(1/p)
}
lambdak2
all.equal(lambdak, lambdak2)

C <- up_C(Sk, nk, lambdak)
sigmas[,,1] / det(sigmas[,,1])^(1/p)

up_Sigma(C, lambdak)

i = 1; k = 1
sapply(1:J, function(x) dmvt(y[x,,i], mus[k,], sigmas[,,k], nus[k], log = FALSE))
y[,,i]
mvnfast::dmvt(y[x,,i], mus[k,], sigmas[,,k], nus[k])
sapply(1:J, function(x) mvnfast::dmvt(y[x,,i], mus[k,], sigmas[,,k], nus[k]))
MixtClust:::dMVT(Y, mus[k,], Sigmas[,,k], nus[k])

Z2 / apply(Z2, 1, sum)
Z
matrix(c(1,2,3,4), nr = 2) / c(-1,2)


sapply(1:J, function(x) dmvt(y[x,,i], mus[k,], sigmas[,,k], nus[k], log = FALSE))

apply(rmvt(1000, sigma = sigmas[,,k], df = nus[k], delta = mus[k,]), 2, mean)
