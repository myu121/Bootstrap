library(boot)
library(quantreg)
library(sn)
test_fun = function(x){
  return(qsn(0.9,alpha=x))
}
uniroot(test_fun,c(-20,0),tol=1e-8)

slant = -3.077684


set.seed(793)
MC_num = 100
n = 100
boot_num = 1000
beta <- c(1,1,1)
qreg_fun <- function(df,ind){
  #df$y <- jitter(y,1e-8)
  fit <- summary(rq(y~x1+x2,tau=0.9,data=df[ind,]))
  boot_coef = as.numeric(fit$coefficients[,1])
  return(boot_coef)
}
MC.per.b0 <- matrix(0,ncol=2,nrow=MC_num)
MC.per.b1 <- matrix(0,ncol=2,nrow=MC_num)
MC.per.b2 <- matrix(0,ncol=2,nrow=MC_num)
se.beta <- matrix(0,ncol=3,nrow=100)
for(i in 1:MC_num){
  x1 <- rlnorm(n)
  x2 <- c(rep(1,n*0.8),rep(0,n*0.2))
  y <- numeric(n)
  for(j in 1:n){
    y[j]<- beta[1]+beta[2]*x1[j]+beta[3]*x2[j]+1/sqrt(3)*(2+(1+(x1[j]-8)^2+x2[j])/10)*rsn(1,alpha=slant)
  }
  boot_df <- data.frame(y=y,x1=x1,x2=x2)
  options(warn=-1)
  boot.out <- boot(data=boot_df,qreg_fun,R=999)
  options(warn=0)
  se.beta[i,1] <- sd(boot.out$t[,1])
  se.beta[i,2] <- sd(boot.out$t[,2])
  se.beta[i,3] <- sd(boot.out$t[,3])
  MC.per.b0[i,] <- as.numeric(quantile(boot.out$t[,1],c(0.05,0.95)))
  MC.per.b1[i,] <- as.numeric(quantile(boot.out$t[,2],c(0.05,0.95)))
  MC.per.b2[i,] <- as.numeric(quantile(boot.out$t[,3],c(0.05,0.95)))
  cat("Iteration finished:",i,"\n")
}
colnames(se.beta) <- c("b0","b1","b2")
head(se.beta)
coverage.out = c(length(which((MC.per.b0[,1]<beta[1])&(MC.per.b0[,2]>beta[1]))),
                 length(which((MC.per.b1[,1]<beta[2])&(MC.per.b1[,2]>beta[2]))),length(which((MC.per.b2[,1]<beta[3])&(MC.per.b2[,2]>beta[3]))))/MC_num


set.seed(793)
numit = 100
n=100
beta.final = matrix(rep(0,3*numit),nrow = 100000, ncol = 3)

for(it in 1:100000){
  x0 = rep(1,n)
  x1 = rlnorm(n, meanlog= 0 , sdlog = 1)
  x2= c(rep(1,n*0.2), rep(0, n*0.8))
  e = rsn(n, alpha=slant)
  b = c(1,1,1)
  y = b[1] + b[2]*x1 + b[3]*x2 + 3^(-0.5)*(2+(1+(x1-8)^2+x2)/10)*e
  
  me = rq(y ~ x1+x2, tau = 0.9)
  beta.final[it,] = summary(me)$coefficients[,1]
  if(it%%1e4==0) cat("Iteration finished:",it,"\n")
}

b0sd = sd(beta.final[,1])
b1sd = sd(beta.final[,2])
b2sd = sd(beta.final[,3])
b0.real.9 = 0.66
b1.real.9 = 0.125
b2.real.9 = 0.66
