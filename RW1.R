#using rq
set.seed(793)
L = 10000
n = 100
B = 1000
C1_total = 0
C2_total = 0
C3_total = 0
length = c()
for(i in 1:L){
  print(i)
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
  beta.result = c()
  for(b in 1:B){
    fit = rq(Y ~ X[,-1], weights = rexp(n,1))
    beta.hat = fit$coefficients
    beta.result = rbind(beta.result,beta.hat)
  }
  sort1 = sort(beta.result[,1])
  sort2 = sort(beta.result[,2])
  sort3 = sort(beta.result[,3])
  index.min = floor(B*0.05)
  index.max = floor(B*0.95)
  CI1 = c(sort1[index.min],sort1[index.max])
  CI2 = c(sort2[index.min],sort2[index.max])
  CI3 = c(sort3[index.min],sort3[index.max])
  length1 = sort1[index.max]-sort1[index.min]
  length2 = sort2[index.max]-sort2[index.min]
  length3 = sort3[index.max]-sort3[index.min]
  coverage1 = ifelse(((sort1[index.max] >= 1) &(sort1[index.min] <= 1)),1,0)
  coverage2 = ifelse(((sort2[index.max] >= 1) &(sort2[index.min] <= 1)),1,0)
  coverage3 = ifelse(((sort3[index.max] >= 1) &(sort3[index.min] <= 1)),1,0)
  C1_total = C1_total + coverage1
  C2_total = C2_total + coverage2
  C3_total = C3_total + coverage3
  length = rbind(length,c(length1,length2,length3))
}


