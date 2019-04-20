set.seed(7937)
library(quantreg)
beta.true = rep(1,3)
n= 5000


# Y is the response, X is the covariates, B is the bootstrapping times, 
# quantile is the desired quantile
RW_bootstrapping <- function(Y, X, B = 300, quantile = 0.5){
  n = length(Y)
  p = dim(X)[2] - 1
  Y_star = Y; X_star = X
  beta.result = c()
  i = 1
  while(i <= B){
    fit = rq(Y_star ~ X_star[,-1], tau = quantile)
    beta.hat = fit$coefficients
    beta.result = rbind(beta.result,beta.hat)
    weight = rexp(n,1)
    boot.index = sample.int(n, n, replace = TRUE,prob = weight)
    Y_star = Y_star[boot.index]
    X_star = X_star[boot.index,]
    # error = Y_star - predict(fit, newdata = as.data.frame(X_star))
    # weight = rexp(n,1)
    # weight = sample(c(-1,1),size = n,replace=TRUE,prob = c(0.5,0.5))
    # error_star = weight * error
    # Y_star = X_star %*% beta.hat + error_star
    # i = i+1
    # Y_transform = X_star %*% beta.hat + error_star

     i = i+1
  }
  se.beta = apply(beta.result,2,sd)
  CI_length = 2*qnorm(0.95)*se.beta
  coverage = numeric(3)
  for(i in 1:3){
    min = mean(beta.result[,i]) - qnorm(0.95)*se.beta[i]
    max = mean(beta.result[,i]) + qnorm(0.95)*se.beta[i]
    coverage[i] = ifelse(((min <= 1) & (1<= max)),1,0)
  }
  return(list(beta.result = beta.result, se.beta = se.beta,
              CI_length = CI_length, coverage = coverage, 
              CI = c(mean(beta.result[,1]) - qnorm(0.95)*se.beta[1], 
                     mean(beta.result[,1]) + qnorm(0.95)*se.beta[1],
                     mean(beta.result[,2]) - qnorm(0.95)*se.beta[2],
                     mean(beta.result[,2]) + qnorm(0.95)*se.beta[2],
                     mean(beta.result[,3]) - qnorm(0.95)*se.beta[3],
                     mean(beta.result[,3]) + qnorm(0.95)*se.beta[3])))
}

CI_result = c()
for(L in 1:1000){
  error <- rt(n, 3)
  x1 <- rlnorm(n, 0, 1)
  x2 <- c(rep(1,n*0.8), rep(0, n*0.2))
  X = cbind(rep(1,n),x1, x2)
  coeff <- function(a, b){
    return(((1 + (a-8)^2 + b)/10 + 2) * 3^(-0.5))
  }
  
  Y <- numeric(n)
  
  
  for(i in 1:n){
    Y[i] = sum(beta.true * X[i,]) + coeff(X[i,1], X[i,2]) * error[i]
  }
  
  CI_result = rbind(CI_result, RW_bootstrapping(Y, X)$CI_length)
}
