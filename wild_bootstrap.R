library(quantreg)

set.seed(793)


#Probability Integral Transform for weight distribution
finv=function(u){
  if(u<=0.5){
    w=-sqrt(25/16-2*u)
  }
  else{
    w=sqrt(2*u-7/16)
  }
  return(w)
}


Beta_hat={}
Length={}
L_Beta={}
U_Beta={}

#Replicate
for (L in 1:1000){
  #Date generate
  n=50
  x0=rep(1,n)
  x1=rlnorm(n,meanlog=0,sdlog=1)
  x2=c(rep(1, n*0.8),rep(0, n*0.2))
  X=cbind(x0,x1,x2)
  epsilon=rt(n, df=3)
  beta=c(1,1,1)
  y=beta[1]+beta[2]*x1+beta[3]*x2+3^(-0.5)*(2+(1+(x1-8)^2+x2)/10)*epsilon
  fit=rq(y~x1+x2, tau=0.5)
  beta_hat=fit$coefficients
  Beta_hat=rbind(Beta_hat,beta_hat)
  ehat=y-X%*%beta_hat
  
  #Wild Bootstrap
  B=1000
  Beta_star={}
  
  for(l in 1:B){
    u=runif(n)
    w=sapply(u,finv)
    estar=w*abs(ehat)
    ystar=X%*%beta_hat+estar
    fit_star=rq(ystar~x1+x2, tau=0.5)
    beta_star=fit_star$coefficients
    Beta_star=rbind(Beta_star, beta_star)
  }
  
  
  sd=apply(Beta_star,2,sd)
  #CI Length
  length=2*qnorm(0.975)*sd
  Length=rbind(Length,length)
  
  #Construct Confidence Interval
  L_beta=beta_hat-0.5*length
  U_beta=beta_hat+0.5*length
  L_Beta=rbind(L_Beta,L_beta)
  U_Beta=rbind(U_Beta,U_beta)
}


#Average Length
apply(Length,2,mean)
#Coverage Probability
length(which((L_Beta[,1]<1)&(U_Beta[,1]>1)))/1000
length(which((L_Beta[,2]<1)&(U_Beta[,2]>1)))/1000
length(which((L_Beta[,3]<1)&(U_Beta[,3]>1)))/1000


#Paired Bootstrap
library(boot)
library(quantreg)
set.seed(793)
#Monte Carlo sample size
MC_num = 100
#Sample size
n = 100
#Bootstrapping sample size
boot_num = 999
#True beta value
beta <- c(1,1,1)
#Confidence level
perc <- 0.95
qreg_fun <- function(df,ind){
  fit <- summary(rq(y~x1+x2,tau=0.5,data=df[ind,]),se="nid")
  boot_coef = as.numeric(fit$coefficients[,1])
  boot_se = as.numeric(fit$coefficients[,2])
  return(c(boot_coef,boot_se))
}
MC.ci.b0 <- matrix(0,ncol=2,nrow=MC_num)
MC.ci.b1 <- matrix(0,ncol=2,nrow=MC_num)
MC.ci.b2 <- matrix(0,ncol=2,nrow=MC_num)
tic("Paired Bootstrap")
for(i in 1:MC_num){
  x1 <- rlnorm(n)
  x2 <- c(rep(1,n*0.8),rep(0,n*0.2))
  y <- numeric(n)
  for(j in 1:n){
    y[j]<- beta[1]+beta[2]*x1[j]+beta[3]*x2[j]+1/sqrt(3)*(2+(1+(x1[j]-8)^2+x2[j])/10)*rt(1,3)
  }
  boot_df <- data.frame(y=y,x1=x1,x2=x2)
  options(warn=-1)
  boot.out <- boot(data=boot_df,qreg_fun,R=boot_num)
  options(warn=0)
  #Obtain Bootstrap-t CI (default 95%)
  MC.ci.b0[i,] <- boot.ci(boot.out,type="stud",index=c(1,4),conf=perc)$student[4:5]
  MC.ci.b1[i,] <- boot.ci(boot.out,type="stud",index=c(2,5),conf=perc)$student[4:5]
  MC.ci.b2[i,] <- boot.ci(boot.out,type="stud",index=c(3,6),conf=perc)$student[4:5]
}
toc()

#Coverage
length(which((MC.ci.b0[,1]<beta[1])&(MC.ci.b0[,2]>beta[1])))/MC_num
length(which((MC.ci.b1[,1]<beta[2])&(MC.ci.b1[,2]>beta[2])))/MC_num
length(which((MC.ci.b2[,1]<beta[3])&(MC.ci.b2[,2]>beta[3])))/MC_num
#Length
mean(MC.ci.b0[,2]-MC.ci.b0[,1])
mean(MC.ci.b1[,2]-MC.ci.b1[,1])
mean(MC.ci.b2[,2]-MC.ci.b2[,1])
