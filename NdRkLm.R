library(quantreg)

#Estimate Quantile Regressions se
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


#Various types of estimators
l1.final = matrix(rep(0,3*numit),nrow = numit, ncol = 3)
lqs.final = matrix(rep(0,3*numit),nrow = numit, ncol = 3)
l2.final = matrix(rep(0,3*numit),nrow = numit, ncol = 3)
covq.final = matrix(rep(0,3*numit),nrow = numit, ncol = 3)
covqs.final = matrix(rep(0,3*numit),nrow = numit, ncol = 3)
covr.final = matrix(rep(0,3*numit),nrow = numit, ncol = 3)
ser.final = matrix(rep(0,3*numit),nrow = numit, ncol = 3)
seqs.final = matrix(rep(0,3*numit),nrow = numit, ncol = 3)


for(num in 1:numit){
  
  covq = c(0,0,0)
  covqs = c(0,0,0)
  covr = c(0,0,0)
  
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
  
  m2 = lm(y~ x1+x2)
  se2 = summary(m2)$coefficients[,2]
  creg = confint(m2,level = 0.9)
  l2 = c(creg[1,2]-creg[1,1],creg[2,2]-creg[2,1],creg[3,2]-creg[3,1])
  
  for(j in 1:3){
    if (creg[j,2] >=1 && creg[j,1] <=1){
      covr[j] = 1
    }
  }
  
  l1.final[num,] = l
  covq.final[num,] = covq
  lqs.final[num,] = lqs
  covqs.final[num,] = covqs
  seqs.final[num,] = seqs
  l2.final[num,] = l2
  covr.final[num,] = covr
  ser.final[num,] = se2
}

colSums(covq.final)/numit
colSums(covqs.final)/numit
colSums(covr.final)/numit
colSums(l1.final)/numit
colSums(lqs.final)/numit
colSums(l2.final)/numit
colMeans(ser.final)
colMeans(seqs.final)



