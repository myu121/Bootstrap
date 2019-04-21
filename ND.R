library(quantreg)

set.seed(793)
numit = 10000
n=50

#Nd
lqs.final = matrix(rep(0,3*numit),nrow = numit, ncol = 3)
covqs.final = matrix(rep(0,3*numit),nrow = numit, ncol = 3)
seqs.final = matrix(rep(0,3*numit),nrow = numit, ncol = 3)


for(num in 1:numit){

  covqs = c(0,0,0)
  
  x0 = rep(1,n)
  x1 = rlnorm(n, meanlog= 0 , sdlog = 1)
  x2= c(rep(1,n*0.8), rep(0, n*0.2))
  e = rt(n, df=3)
  b = c(1,1,1)
  y = b[1] + b[2]*x1 + b[3]*x2 + 3^(-0.5)*(2+(1+(x1-8)^2+x2)/10)*e
  

  mqs = rqss(y ~ x1+x2, tau = 0.5)
  seqs = summary(mqs)$coef[,2]
  paraqs = summary(mqs)$coef[,1]
  upqs = paraqs + seqs*qt(0.95, df=48)
  loqs = paraqs - seqs*qt(0.95, df=48)
  
  for(k in 1:3){
    if (upqs[k] >= 1 && loqs[k] <= 1){
      covqs[k] = 1
    }
  }
  
  lqs = upqs - loqs

  lqs.final[num,] = lqs
  covqs.final[num,] = covqs
  seqs.final[num,] = seqs

}


colSums(covqs.final)/numit
colSums(lqs.final)/numit
colMeans(seqs.final)
sd(lqs.final[,1])/sqrt(numit)
sd(lqs.final[,2])/sqrt(numit)
sd(lqs.final[,3])/sqrt(numit)

