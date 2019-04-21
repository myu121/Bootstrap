library(foreach)
library(parallel)
library(doParallel)
library(tictoc)
set.seed(793)
cl <-makeCluster(detectCores())
registerDoParallel(cl)
MC_num = 10000
n = 100
boot_num = 1000
beta <- c(1,1,1)
qreg_fun <- function(df,ind){
  fit <- summary(rq(y~x1+x2,tau=0.5,data=df[ind,]),se="nid")
  boot_coef = as.numeric(fit$coefficients[,1])
  boot_se = as.numeric(fit$coefficients[,2])
  return(c(boot_coef,boot_se))
}
tic("Paired Bootstrap")
output<-foreach(i = 1:MC_num,.combine = append, .packages = c("parallel","doParallel","boot","quantreg")) %dopar%{
  x1 <- rlnorm(n)
  x2 <- c(rep(1,n*0.8),rep(0,n*0.2))
  y <- numeric(n)
  for(j in 1:n){
    y[j]<- beta[1]+beta[2]*x1[j]+beta[3]*x2[j]+1/sqrt(3)*(2+(1+(x1[j]-8)^2+x2[j])/10)*rt(1,3)
  }
  boot_df <- data.frame(y=y,x1=x1,x2=x2)
  options(warn=-1)
  boot.out <- boot(data=boot_df,qreg_fun,R=999,parallel = "multicore",ncpus=detectCores())
  options(warn=0)
  par_results <- list(
    Per.b0 =  as.numeric(quantile(boot.out$t[,1],c(0.05,0.95))),
    Per.b1 =  as.numeric(quantile(boot.out$t[,2],c(0.05,0.95))),
    Per.b2 =  as.numeric(quantile(boot.out$t[,3],c(0.05,0.95))),
    Boot.t.b0 = boot.ci(boot.out,type="stud",index=c(1,4),conf = 0.9)$student[4:5],
    Boot.t.b1 = boot.ci(boot.out,type="stud",index=c(2,5),conf = 0.9)$student[4:5],
    Boot.t.b2 = boot.ci(boot.out,type="stud",index=c(3,6),conf = 0.9)$student[4:5]
  )
  par_results <- list(par_results)
  names(par_results) <- paste0("step_", i)
  return(par_results)
}
toc()
stopCluster(cl)

coverage.prob = rep(0,6)
names(coverage.prob)=c("NP.b0","NP.b1","NP.b2","BT.b0","BT.b1","BT.b2")
for(i in 1:MC_num){
  if((output[[i]]$Per.b0[1]<beta[1])&&(output[[i]]$Per.b0[2]>beta[1])) coverage.prob[1]=coverage.prob[1]+1
  if((output[[i]]$Per.b1[1]<beta[2])&&(output[[i]]$Boot.t.b1[2]>beta[2])) coverage.prob[2]=coverage.prob[2]+1
  if((output[[i]]$Per.b2[1]<beta[3])&&(output[[i]]$Boot.t.b2[2]>beta[3])) coverage.prob[3]=coverage.prob[3]+1
  if((output[[i]]$Boot.t.b0[1]<beta[1])&&(output[[i]]$Boot.t.b0[2]>beta[1])) coverage.prob[4]=coverage.prob[4]+1
  if((output[[i]]$Boot.t.b1[1]<beta[2])&&(output[[i]]$Boot.t.b1[2]>beta[2])) coverage.prob[5]=coverage.prob[5]+1
  if((output[[i]]$Boot.t.b2[1]<beta[3])&&(output[[i]]$Boot.t.b2[2]>beta[3])) coverage.prob[6]=coverage.prob[6]+1
}
coverage.prob/MC_num
ci.length = matrix(0,ncol=2,nrow=6)
colnames(ci.length) <- c("5%","95%")
for(i in 1:MC_num){
  ci.length[1,1] = ci.length[1,1]+output[[i]]$Per.b0[1]
  ci.length[1,2] = ci.length[1,2]+output[[i]]$Per.b0[2]
  ci.length[2,1] = ci.length[2,1]+output[[i]]$Per.b1[1]
  ci.length[2,2] = ci.length[2,2]+output[[i]]$Per.b1[2]
  ci.length[3,1] = ci.length[3,1]+output[[i]]$Per.b2[1]
  ci.length[3,2] = ci.length[3,2]+output[[i]]$Per.b2[2]
  ci.length[4,1] = ci.length[4,1]+output[[i]]$Boot.t.b0[1]
  ci.length[4,2] = ci.length[4,2]+output[[i]]$Boot.t.b0[2]
  ci.length[5,1] = ci.length[5,1]+output[[i]]$Boot.t.b1[1]
  ci.length[5,2] = ci.length[5,2]+output[[i]]$Boot.t.b1[2]
  ci.length[6,1] = ci.length[6,1]+output[[i]]$Boot.t.b2[1]
  ci.length[6,2] = ci.length[6,2]+output[[i]]$Boot.t.b2[2]
}
ci.length = ci.length/MC_num
ave.length = ci.length[,2]-ci.length[,1]
names(ave.length) = names(coverage.prob)
