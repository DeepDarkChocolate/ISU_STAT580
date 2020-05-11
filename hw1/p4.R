rtnorm <- function(n, d,..., rate = FALSE){
  res <- c() # random deviates
  if(rate == FALSE){ # Only generate random deviates
    for(i in 1:n){
      while(TRUE){
        y = sqrt(rexp(n = 1, rate = 1/2) + d^2)
        u = runif(1)
        if(u < d / y){
          res <- c(res, y)
          break
        }
      }
    }
    return(res)
  }else if(rate == TRUE){ # Also compute the acceptance rate
    cnt = 0 # number of all trials
    for(i in 1:n){
      while(TRUE){
        cnt = cnt + 1
        y = sqrt(rexp(n = 1, rate = 1/2) + d^2)
        u = runif(1)
        if(u < d / y){
          res <- c(res, y)
          break
        }
      }
    }
    return(list(res = res, rate = n / cnt))
  }
}

acceptrate <- function(d){
  return(1/(1/pnorm(-d)/sqrt(2*pi)*exp(-d^2/2)/d))
}

dvec <- c(1,3,5,7)
sample1 <- sapply(dvec, function(x){rtnorm(10,x)})
colnames(sample1) <- dvec
sample1

sapply(dvec, function(x){mean(rtnorm(1000,x))})
sapply(dvec, function(x){dnorm(x)/pnorm(-x)})
# mean(rtnorm(1000, 5))
# dnorm(5)/pnorm(-5)


hist(rtnorm(10000, 5), breaks = 100)

windows()
par(mfrow = c(2,2))
for(i in dvec){
  hist(rtnorm(10000, i), xlim = c(0, 10), freq = FALSE, ylim = c(0, 5.5), xlab = NULL,
       main = paste("d = ", i), breaks = 80)
}

table1 <- sapply(dvec, function(x){
  simul.rr = 1 - rtnorm(10000, x, rate = TRUE)$rate
  theo.rr = 1 - acceptrate(x)
  sd = sqrt(simul.rr * (1 - simul.rr) / 10000)
  upper = simul.rr + sd*qnorm(0.975)
  lower = simul.rr - sd*qnorm(0.975)
  return(c(theo.rr = theo.rr, simul.rr = simul.rr,
           upperCI = upper, lowerCI = lower
           ))})
colnames(table1) <- dvec
t(table1)

windows()
par(mfrow = c(1,1))
curve(acceptrate, xlim = c(0, 5), col = 4, ylab = "p", xlab = "d", main = "Acceptance rate")
seq1 <- exp(seq(from = 0.01, to = log(6), length.out = 15)) - 1
lines(x = seq1, y = sapply(seq1, function(x){rtnorm(10000, x, rate = TRUE)$rate}), col = 2)
legend("bottomright", c("Theoretical value", "Simulated value"), col = c(4, 2), lty = 1)
