set.seed(793)
numit = 10000
n=50


#Linear Regression

l2.final = matrix(rep(0,3*numit),nrow = numit, ncol = 3)
covr.final = matrix(rep(0,3*numit),nrow = numit, ncol = 3)
ser.final = matrix(rep(0,3*numit),nrow = numit, ncol = 3)



for(num in 1:numit){
  
  covr = c(0,0,0)
  
  x0 = rep(1,n)
  x1 = rlnorm(n, meanlog= 0 , sdlog = 1)
  x2= c(rep(1,n*0.8), rep(0, n*0.2))
  e = rt(n, df=3)
  b = c(1,1,1)
  y = b[1] + b[2]*x1 + b[3]*x2 + 3^(-0.5)*(2+(1+(x1-8)^2+x2)/10)*e
  
  m2 = lm(y~ x1+x2)
  se2 = summary(m2)$coefficients[,2]
  creg = confint(m2,level = 0.9)
  l2 = c(creg[1,2]-creg[1,1],creg[2,2]-creg[2,1],creg[3,2]-creg[3,1])
  
  for(j in 1:3){
    if (creg[j,2] >=1 && creg[j,1] <=1){
      covr[j] = 1
    }
  }
  

  l2.final[num,] = l2
  covr.final[num,] = covr
  ser.final[num,] = se2
}


colSums(covr.final)/numit
colSums(l2.final)/numit
colMeans(ser.final)
sd(l2.final[,1])/sqrt(numit)
sd(l2.final[,2])/sqrt(numit)
sd(l2.final[,3])/sqrt(numit)



