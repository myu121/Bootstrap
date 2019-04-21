library(quantreg)

set.seed(793)
numit = 10000
n=50

#Rk
l1.final = matrix(rep(0,3*numit),nrow = numit, ncol = 3)
covq.final = matrix(rep(0,3*numit),nrow = numit, ncol = 3)

for(num in 1:numit){
  
  covq = c(0,0,0)
  
  x0 = rep(1,n)
  x1 = rlnorm(n, meanlog= 0 , sdlog = 1)
  x2= c(rep(1,n*0.8), rep(0, n*0.2))
  e = rt(n, df=3)
  b = c(1,1,1)
  y = b[1] + b[2]*x1 + b[3]*x2 + 3^(-0.5)*(2+(1+(x1-8)^2+x2)/10)*e
  
  m = rq(y ~ x1+x2, tau = 0.5)
  up1 = summary(m)$coefficients[,3]
  lo1 = summary(m)$coefficients[,2]
  l = up1 - lo1
  
  for(i in 1:3){
    if (up1[i] >= 1 && lo1[i] <= 1){
      covq[i] = 1
    }
  }
  
  l1.final[num,] = l
  covq.final[num,] = covq
}

colSums(covq.final)/numit
colSums(l1.final)/numit
sd(l1.final[,1])/sqrt(numit)
sd(l1.final[,3])/sqrt(numit)


