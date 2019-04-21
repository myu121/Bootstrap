#Benchmark SE Values
library(quantreg)
set.seed(793)
numit = 50
n=50
beta.final = matrix(rep(0,3*numit),nrow = 5000, ncol = 3)

for(it in 1:5000){
  x0 = rep(1,n)
  x1 = rlnorm(n, meanlog= 0 , sdlog = 1)
  x2= c(rep(1,n*0.8), rep(0, n*0.2))
  e = rt(n, df=3)
  b = c(1,1,1)
  y = b[1] + b[2]*x1 + b[3]*x2 + 3^(-0.5)*(2+(1+(x1-8)^2+x2)/10)*e
  
  me = rq(y ~ x1+x2, tau = 0.5)
  beta.final[it,] = summary(me)$coefficients[,1]}

b0sd = sd(beta.final[,1])
b1sd = sd(beta.final[,2])
b2sd = sd(beta.final[,3])