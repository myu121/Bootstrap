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

M=10000         #MC-samples
B=1000        #Bootstrap replicates
n=100         #sample size
tau=0.5       #quantile level
alpha=0.1     #1 minus theoretical coverage probability
beta=c(1,1,1) #true parameter

Beta_hat={}
SD={}
Per.Length={}
Per.L={}
Per.U={}
Bt.Length={}
Bt.L={}
Bt.U={}

#Replicate
for (L in 1:M){
  #Date generate
  x0=rep(1,n)
  x1=rlnorm(n,meanlog=0,sdlog=1)
  x2=c(rep(1, n*0.8),rep(0, n*0.2))
  X=cbind(x0,x1,x2)
  epsilon=rt(n, df=3)
  y=beta[1]+beta[2]*x1+beta[3]*x2+3^(-0.5)*(2+(1+(x1-8)^2+x2)/10)*epsilon
  fit=rq(y~x1+x2, tau=0.5)
  beta_hat=fit$coefficients
  Beta_hat=rbind(Beta_hat,beta_hat)
  ehat=y-X%*%beta_hat
  
  #correction for ehat
  fhat_0=akj(ehat,0)$dens
  H=diag(X%*%solve(t(X)%*%X)%*%t(X))
  Phi=tau-(ehat<0)
  ehat=ehat+1/fhat_0*H*Phi
  
  #Wild Bootstrap
  Beta_star={}
  Z_B={}
  
  for(l in 1:B){
    u=runif(n)
    w=sapply(u,finv)
    estar=w*abs(ehat)
    ystar=X%*%beta_hat+estar
    fit_star=rq(ystar~x1+x2, tau=tau)
    beta_star=fit_star$coefficients
    Beta_star=rbind(Beta_star, beta_star)
    se_b=summary(fit_star, se="nid")$coefficients[,2]
    z_b=(beta_star-beta_hat)/se_b
    Z_B=rbind(Z_B,z_b)
  }
  
  #compute standard error
  sd=apply(Beta_star,2,sd)
  #compute bootstrap-t CI and length
  q_u=apply(Z_B,2,quantile,probs=alpha/2)
  q_l=apply(Z_B,2,quantile,probs=1-alpha/2)
  bt.L=beta_hat-q_l*sd
  bt.U=beta_hat-q_u*sd
  bt.length=(q_l-q_u)*sd
  #compute percentile CI and length
  per.U=apply(Beta_star,2,quantile,probs=1-alpha/2)
  per.L=apply(Beta_star,2,quantile,probs=alpha/2)
  per.length=per.U-per.L
  #bind results
  SD=rbind(SD,sd)
  Per.Length=rbind(Per.Length,per.length)
  Per.L=rbind(Per.L,per.L)
  Per.U=rbind(Per.U,per.U)
  Bt.Length=rbind(Bt.Length,bt.length)
  Bt.L=rbind(Bt.L,bt.L)
  Bt.U=rbind(Bt.U,bt.U)
}

#Average Length and Coverage Probability for percentile CI
apply(Per.Length,2,mean)
apply(Per.Length,2,sd)/sqrt(M)
length(which((Per.L[,1]<beta[1])&(Per.U[,1]>beta[1])))/M
length(which((Per.L[,2]<beta[2])&(Per.U[,2]>beta[2])))/M
length(which((Per.L[,3]<beta[3])&(Per.U[,3]>beta[3])))/M

#Average Length and Coverage Probability for bootstrap-t CI
apply(Bt.Length,2,mean)
apply(Bt.Length,2,sd)/sqrt(M)
length(which((Bt.L[,1]<beta[1])&(Bt.U[,1]>beta[1])))/M
length(which((Bt.L[,2]<beta[2])&(Bt.U[,2]>beta[2])))/M
length(which((Bt.L[,3]<beta[3])&(Bt.U[,3]>beta[3])))/M

#Save estimated standard error
write.csv(SD, file="SE_for_wild_bootstrap.csv")