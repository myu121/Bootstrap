library(boot)
library(quantreg)
#seed for hard fitting 35
set.seed(793)
MC_num = 10000
n = 100
boot_num = 1000
beta <- c(1,1,1)
qreg_fun <- function(df,ind){
  while(sum(df[ind,]$x2==1) == n){
    ind = sample(ind,replace=TRUE)
  }
  fit <- summary(rq(y~x1+x2,tau=0.5,data=df[ind,]),se="nid")
  boot_coef = as.numeric(fit$coefficients[,1])
  boot_se = as.numeric(fit$coefficients[,2])
  return(c(boot_coef,boot_se))
}
MC.ci.b0 <- matrix(0,ncol=2,nrow=MC_num)
MC.ci.b1 <- matrix(0,ncol=2,nrow=MC_num)
MC.ci.b2 <- matrix(0,ncol=2,nrow=MC_num)
MC.per.b0 <- matrix(0,ncol=2,nrow=MC_num)
MC.per.b1 <- matrix(0,ncol=2,nrow=MC_num)
MC.per.b2 <- matrix(0,ncol=2,nrow=MC_num)
for(i in 1:MC_num){
  x1 <- rlnorm(n)
  x2 <- c(rep(1,n*0.8),rep(0,n*0.2))
  y <- numeric(n)
  for(j in 1:n){
    y[j]<- beta[1]+beta[2]*x1[j]+beta[3]*x2[j]+1/sqrt(3)*(2+(1+(x1[j]-8)^2+x2[j])/10)*rt(1,3)
  }
  boot_df <- data.frame(y=y,x1=x1,x2=x2)
  options(warn=-1)
  boot.out <- boot(data=boot_df,qreg_fun,R=999)
  options(warn=0)
  MC.per.b0[i,] <- as.numeric(quantile(boot.out$t[,1],c(0.05,0.95)))
  MC.per.b1[i,] <- as.numeric(quantile(boot.out$t[,2],c(0.05,0.95)))
  MC.per.b2[i,] <- as.numeric(quantile(boot.out$t[,3],c(0.05,0.95)))
  MC.ci.b0[i,] <- boot.ci(boot.out,type="stud",index=c(1,4),conf = 0.9)$student[4:5]
  MC.ci.b1[i,] <- boot.ci(boot.out,type="stud",index=c(2,5),conf = 0.9)$student[4:5]
  MC.ci.b2[i,] <- boot.ci(boot.out,type="stud",index=c(3,6),conf = 0.9)$student[4:5]
}

#Coverage
coverage.out = c(length(which((MC.ci.b0[,1]<beta[1])&(MC.ci.b0[,2]>beta[1]))),length(which((MC.ci.b1[,1]<beta[2])&(MC.ci.b1[,2]>beta[2]))),
                 length(which((MC.ci.b2[,1]<beta[3])&(MC.ci.b2[,2]>beta[3]))),length(which((MC.per.b0[,1]<beta[1])&(MC.per.b0[,2]>beta[1]))),
                 length(which((MC.per.b1[,1]<beta[2])&(MC.per.b1[,2]>beta[2]))),length(which((MC.per.b2[,1]<beta[3])&(MC.per.b2[,2]>beta[3]))))/MC_num
#Mean length
ave.length = c(mean(MC.per.b0[,2]- MC.per.b0[,1]),mean(MC.per.b1[,2]- MC.per.b1[,1]),mean(MC.per.b2[,2]- MC.per.b2[,1]), 
               mean(MC.ci.b0[,2]- MC.ci.b0[,1]),mean(MC.ci.b1[,2]- MC.ci.b1[,1]),mean(MC.ci.b2[,2]- MC.ci.b2[,1]))
#SD length
sd.length = c(sd(MC.per.b0[,2]- MC.per.b0[,1]),sd(MC.per.b1[,2]- MC.per.b1[,1]),sd(MC.per.b2[,2]- MC.per.b2[,1]), 
              sd(MC.ci.b0[,2]- MC.ci.b0[,1]),sd(MC.ci.b1[,2]- MC.ci.b1[,1]),sd(MC.ci.b2[,2]- MC.ci.b2[,1]))
summary.out = rbind(coverage.out,ave.length,sd.length)
colnames(summary.out)=c("NP.b0","NP.b1","NP.b2","BT.b0","BT.b1","BT.b2")
rownames(summary.out)=c("Coverage","Mean","SD")
summary.out
