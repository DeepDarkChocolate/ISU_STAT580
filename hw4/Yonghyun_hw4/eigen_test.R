mat <- diag(2:10)
mat[1,-1] <- 1
mat[-1,1] <- 1
mat

eigen1 <- eigen(mat)
lambda <- eigen1$values[9:1]
P <- eigen1$vectors[,9:1]
Lambda <- diag(lambda)
P %*% Lambda %*% t(P)
t(P) %*% P
all.equal(P %*% Lambda %*% t(P), mat)
all.equal(t(P) %*% P, diag(9))
