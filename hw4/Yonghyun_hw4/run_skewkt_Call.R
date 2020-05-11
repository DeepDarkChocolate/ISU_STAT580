# source this file in R to test the skewkt functions
# gcc -shared -fPIC -o Rskewkt.so skewkt.c
# Rscript run_skewkt.R

library(mvtnorm)
source("skewkt_Call.R")

# randomly generate pos. def. sigma using multivariate normal theory
Sigma <- function(p) var(matrix(rnorm(p*(p + 1)), ncol = p))

# generate data
n <- 1000
X <- rmvt(n = n, sigma = Sigma(10), df = 1)

sk.c <- skewkt(X, scale = T)
cat("Skew estimated in C:", sk.c, "\n")
sk.r <- Rskewkt(X, scale = T)
cat("Skew estimated in R:", sk.r, "\n")
