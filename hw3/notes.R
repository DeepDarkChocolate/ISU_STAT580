# 1

mat = matrix(c(4,3,2,1,3,4,3,2,2,3,4,3,1,2,3,4), nr = 4)
eigen(mat)
e <- eigen(mat)
sum(e$vectors[,1]^2)
install.packages("matlib")
library(matlib)



mat %*% rep(0.5, 4) / sqrt(122)

res = powerMethod(mat)
vec1 = res$vector
lambda1 = res$value
mat2 = mat - lambda1 * vec1 %*% t(vec1)


res2 = powerMethod(mat2, eps = 1e-20)
vec2 = res2$vector
lambda2 = res2$value
mat3 = mat2 - lambda2 * vec2 %*% t(vec2)
powerMethod(mat3, eps = 1e-20)


n = 20
mat = matrix(nrow = n, ncol = n)
for(i in 1:n){
  for(j in 1:n){
    mat[i, j] <- 1/(i+j-1)
  }
}

eigen(mat)

n = 10
x = rep(1:n, each = n)
y = rep(1:n, n)
plot(x*(1+ 1/(abs(x-y) + 1)), y*(1+1/(abs(x-y) + 1)))

pts = matrix(c(1,1,0,0,1,0,1,0),nc = 2)
for(i in 1:10) pts = rbind(pts, matrix(c(0, 1/2^i, 1/2^i, 1 -1/2^i, 1 -1/2^i, 1,
                                         1/2^i, 0, 1/2^i, 1 -1/2^i, 1, 1 -1/2^i), nc = 2))
plot(pts)


mat = matrix(nr = 10, nc = 5);
for(i in 1:10){
  for(j in 1:5){
    mat[i, j] = 1/(abs(i - j)+1)
  }
}
cor(mat)
cor(mat[1,], mat[2,])
cormat = cor(t(mat))
cormat[,5]

cormat[,4]
cormat[,6]

cormat[,1]
cormat[,2]
cormat[,3]
cormat[,4]
cormat[,5]
cormat[,6]
cormat[,7]
cormat[,8]
cormat[,9]
cormat[,10]

mat <- diag(10)
mat[2,1] <- 1
mat[6,1] <- 1
mat[1,2] <- 1
mat[3,2] <- 1
mat[2,3] <- 1
mat[4,2] <- 1
mat[4,3] <- 1
mat[7,3] <- 1
mat[2,4] <- 1
mat[3,4] <- 1
mat[7,4] <- 1
mat[1,6] <- 1
mat[3,7] <- 1
mat[4,7] <- 1


mat
g <- apply(mat, 1, sum)
mat <- diag(g^{-1/2}) %*% mat %*% diag(g^{-1/2})
eigen(mat)
mat %*% c(0, 0, 0, 0, 1/2, 0, 0, 1/2, 1/2, 1/2)

mat <- diag(1:10)
mat[1,] <- 1
mat
eigen(mat)


mat = matrix(c(3,2,1,0,2,3,2,0,1,2,3,0,0,0,0,3), nr =4)
eigen(mat)

mat <- matrix(c(3, 0, 1, 0, 2, 0, 1,0,1), nr = 3)
eigen(mat)

mat %*% rep(1/sqrt(3), 3)

t(evec2) %*% evec1

evec1 <- eigen(mat)$vectors[,1]
evec2 <- matrix(rep(1/sqrt(3), 3), nc = 1)
for(i in 1:100){
  evec2 <- evec2 / norm(evec2, type = "F")
  evec2 <- evec2 - drop(t(evec2) %*% evec1) * evec1
  evec2 <- mat %*% evec2
}
evec2 <- evec2 / norm(evec2, type = "F")

evec2

sum(evec1)
rep(1/sqrt(3), 3)  - 1/sqrt(3)*sum(evec1)*evec1
tmp =rep(1/sqrt(3), 3)  - 1/sqrt(3)*sum(evec1)*evec1

(tmp <- mat %*% tmp)

v <- c(-0.123064,  0.991089,  -0.050975)
t(v) %*% mat %*% v


mat <- diag(1:1000)
mat[1,] <- 1
eigen(mat)$values

P <- eigen(mat)$vectors
evec1 <- eigen(mat)$vectors[,1]
lambda1 <- eigen(mat)$values[1]
evec2 <- eigen(mat)$vectors[,2]
lambda2 <- eigen(mat)$values[2]

eigen(mat - lambda1 * evec1 %*% t(evec1))

v <- matrix(rep(1/sqrt(10), 10), nc = 1)
for(i in 1:1000){
  v <- mat %*% v - lambda1 * sum(evec1 * v) * evec1
  v <- v / norm(v, type = "F")
}
evec2 <- v
lambda2 <- drop(t(v) %*% mat %*% v)

for(i in 1:10000){
  v <- mat %*% v - lambda1 * sum(evec1 * v) * evec1 - lambda2 * sum(evec2 * v) * evec2
  v <- v / norm(v, type = "F")
}

v

sum(evec1 * evec2)
sum(evec2 * evec2)

round(v, 3)
lamb <- drop(t(v) %*% mat %*% v)
round(mat %*% v, 3)
round(drop(t(v) %*% mat %*% v) * v, 3)

mat %*% P[,2]
eigen(mat)$values[2] %*% P[,2]

sum(eigen(mat)$vectors[,1])
# 4

diurnal = read.csv("diurnaldata.csv", header = TRUE)
mat = read.csv("evecs-3-20-0.99.txt", header = FALSE)
mat = as.matrix(mat)
par(mar=c(0, 0, 0, 0))
image(t(mat), useRaster=TRUE, axes=FALSE)



which.max(mat[,2])
order(mat[,2], decreasing = TRUE)[1:5]
sort(mat[,2], decreasing = TRUE)[1:5]
diurnal[7741, ]
diurnal[4589, ]
cor(t(drop(diurnal[7741, 2:12])),
    t(drop(diurnal[4589, 2:12])))
cor(t(diurnal[order(mat[,2], decreasing = TRUE)[1:5], 2:12]))

which
sort(c(4,1,2), decreasing = TRUE)

diurnal[order(mat[,2], decreasing = TRUE)[1:5], 1]
