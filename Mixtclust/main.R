library(Rcpp)
library(mvtnorm)
library(RcppDist)

#sourceCpp("test.cpp")
source("subfuncs.R")

library(devtools)
find_rtools()

n = 400
p = 3
J = 10
K = 4

set.seed(1)

mus = matrix(rnorm(K*p, sd = 4), nr = K, nc = p)
sigmas = array(diag(1:p), dim = c(p,p,K))
for(k in 1:K){
  sigmas[,,k] <- k*sigmas[,,k] + k*rep(1, p) %*% t(rep(1,p))
}
nus = rep(10, K)
pis = runif(K)+0.5
pis = pis / sum(pis)
ztmp <- rmultinom(n, 1, pis)
y = array(0L, dim = c(J,p,n))
for(i in 1:n){
  k = drop((1:K) %*% ztmp[,i])
  y[,,i] <- rmvt(J, sigma = sigmas[,,k], df = nus[k], delta = mus[k,])
}

#y[3,,1] <- NA
for(j in 2:J){
  for(i in 1:n){
    if(rbinom(1,1,0.5)) y[j,,i] <- NA
  }
}

y_naomit <- apply(y, 3, function(x){data.frame(na.omit(x))})

y_numJ <- sapply(y_naomit, nrow)

#Z = up_Z(y, mus, sigmas, nus, pis)
#U = up_U(y, mus, sigmas, nus, pis)

#up_pi(Z)
#up_mu(y, Z, U)
#up_nk(y, Z)
#loglik(y,Z, mus, sigmas, pis, nus)

#library(ClusterR)
#GMM1 <- GMM(t(sapply(y_naomit, colMeans)), K, dist_mode = "maha_dist",
#            km_iter = 100, em_iter = 100)
#mus2 <- GMM1$centroids
#nus2 <- rpois(K, 10) + 1
#pis2 <- GMM1$weights
#sigmas2 <- array(0L, dim = c(p,p,K))
#for(k in 1:K){
#  sigmas2[,,k] <- diag(GMM1$covariance_matrices[k,])
#}

K = 4
kmeans1 <- kmeans(t(sapply(y_naomit, colMeans)), centers = K, iter.max = 100, nstart = 10)
mus2 <- kmeans1$centers
mus2 <- unname(mus2)
nus2 <- rpois(K, 10) + 1
pis2 <- kmeans1$size / sum(kmeans1$size)
zvec2 <- kmeans1$cluster
sigmas2 <- unclass(by(t(sapply(y_naomit, colMeans)), kmeans1$cluster, var))
attributes(sigmas2) <- NULL
sigmas2 <- array(as.numeric(unlist(sigmas2)), dim=c(p,p,K))

perm1 <- permutations(K,K)
permlen <- 1:nrow(perm1)
norm1 <- sapply(permlen, function(x){norm(mus2[perm1[x,],] - mus, "F")})
#norm1 <- apply(perm1, 1, function(x){norm(mus2[x,] - mus, "F")})
perm2 <- perm1[which.min(norm1),]

mus2 <- mus2[perm2,]
pis2 <- pis2[perm2]
sigmas2 <- sigmas2[,,perm2]

#pis2 <- rmultinom(1, 10000, pis) / 10000
#nus2 <- 5 * nus
#mus2 <- mus * c(0.9, 1.5, 0.5, 1.1, 0.9)
#sigmas2 <- 0.25 * sigmas

#EM_iter(oldpars, y)

oldpars<- list(pi = pis2, nu = nus2, mu = mus2, Sigma = sigmas2)

res <- skewtmix_cpp(y, oldpars, df.constr = TRUE)

res$estimates$mu
mus

drop(res$estimates$pi)
pis

res$estimates$Sigmas
sigmas

sum(zvec2 == perm2[res$classification]) / length(res$classification) 
sum(drop((1:K) %*% ztmp) == res$classification) / length(res$classification)

res$estimates$nu
nus

res$bic
#k = 2 49829.14
#k = 3 46861.53
#k = 4 43842.76
#k = 5 44076.86
#k = 6 44202.41

library(scatterplot3d)
#scatterplot3d(apply(y[,,], 2, rbind)[,1:3], color = rep(apply(ztmp, 2, function(x){return((1:K)%*%x)}), each = J))
#scatterplot3d(apply(y[1,,], 1, rbind)[,1:3], color = res$classification)


par(mfrow = c(2,2))
scatterplot3d(t(sapply(y_naomit, colMeans)), color = drop((1:K) %*% ztmp),
              main = "True classification; mean y_1, ... y_n")

scatterplot3d(do.call(rbind, y_naomit)[,1:3], 
              color = rep(drop((1:K) %*% ztmp), y_numJ),
              main = "True classification; y_11, ..., y_nn_J")

scatterplot3d(t(sapply(y_naomit, colMeans)), color = res$classification,
              main = "model based classification; mean y_1, ... y_n")

scatterplot3d(do.call(rbind, y_naomit)[,1:3], 
              color = rep(res$classification, y_numJ),
              main = "model based classification; y_11, ..., y_n{J_n}")


BICres <- c()
Kvec <- 3:7
for(K in Kvec){
  kmeans1 <- kmeans(t(sapply(y_naomit, colMeans)), centers = K, iter.max = 100, nstart = 10)
  mus2 <- kmeans1$centers
  mus2 <- unname(mus2)
  nus2 <- rpois(K, 10) + 1
  pis2 <- kmeans1$size / sum(kmeans1$size)
  zvec2 <- kmeans1$cluster
  sigmas2 <- unclass(by(t(sapply(y_naomit, colMeans)), kmeans1$cluster, var))
  attributes(sigmas2) <- NULL
  sigmas2 <- array(as.numeric(unlist(sigmas2)), dim=c(p,p,K))
  
  oldpars<- list(pi = pis2, nu = nus2, mu = mus2, Sigma = sigmas2)
  
  res <- skewtmix_cpp(y, oldpars, df.constr = TRUE)
  BICres <- c(BICres, res$bic)
}

BICres
plot(Kvec, BICres, type = "l", main = "BIC", xlab= "number of clusters: K")



