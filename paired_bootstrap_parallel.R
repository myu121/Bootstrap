library(foreach)
library(parallel)
library(doParallel)
library(tictoc)
library(doRNG)
out.array = array(0,c(3,6,5))
for(m in 1:5){
  cl <-makeCluster(detectCores())
  registerDoParallel(cl)
  MC_num = 10000
  n = m*100
  boot_num = 1000
  beta <- c(1,1,1)
  qreg_fun <- function(df,ind){
    while((sum(df[ind,]$x2==1) == n)||(sum(df[ind,]$x2==1) == 0)){
      ind = sample(ind,replace=TRUE)
    }
    fit <- summary(rq(y~x1+x2,tau=0.5,data=df[ind,]),se="nid")
    boot_coef = as.numeric(fit$coefficients[,1])
    boot_se = as.numeric(fit$coefficients[,2])
    return(c(boot_coef,boot_se))
  }
  cat("Iteration started:",m,"\n")
  tic("Paired Bootstrap")
  output<-foreach(i = 1:MC_num,.combine = append, .options.RNG=793,.packages = c("parallel","doParallel","boot","quantreg")) %dorng%{
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
    if((output[[i]]$Per.b1[1]<beta[2])&&(output[[i]]$Per.b1[2]>beta[2])) coverage.prob[2]=coverage.prob[2]+1
    if((output[[i]]$Per.b2[1]<beta[3])&&(output[[i]]$Per.b2[2]>beta[3])) coverage.prob[3]=coverage.prob[3]+1
    if((output[[i]]$Boot.t.b0[1]<beta[1])&&(output[[i]]$Boot.t.b0[2]>beta[1])) coverage.prob[4]=coverage.prob[4]+1
    if((output[[i]]$Boot.t.b1[1]<beta[2])&&(output[[i]]$Boot.t.b1[2]>beta[2])) coverage.prob[5]=coverage.prob[5]+1
    if((output[[i]]$Boot.t.b2[1]<beta[3])&&(output[[i]]$Boot.t.b2[2]>beta[3])) coverage.prob[6]=coverage.prob[6]+1
  }
  Coverage = coverage.prob/MC_num
  cat("The coverage probablities for the 90% confidence intervals are:","\n")
  Coverage
  ci.out = matrix(0,ncol=12,nrow=MC_num)
  for(i in 1:MC_num){
    ci.out[i,1] = output[[i]]$Per.b0[1]
    ci.out[i,2] = output[[i]]$Per.b0[2]
    ci.out[i,3] = output[[i]]$Per.b1[1]
    ci.out[i,4] = output[[i]]$Per.b1[2]
    ci.out[i,5] = output[[i]]$Per.b2[1]
    ci.out[i,6] = output[[i]]$Per.b2[2]
    ci.out[i,7] = output[[i]]$Boot.t.b0[1]
    ci.out[i,8] = output[[i]]$Boot.t.b0[2]
    ci.out[i,9] = output[[i]]$Boot.t.b1[1]
    ci.out[i,10] = output[[i]]$Boot.t.b1[2]
    ci.out[i,11] = output[[i]]$Boot.t.b2[1]
    ci.out[i,12] = output[[i]]$Boot.t.b2[2]
  }
  length.out = matrix(0,ncol=6,nrow=MC_num)
  for(i in 1:6){
    length.out[,i] <- ci.out[,2*i]-ci.out[,2*i-1] 
  }
  ave.length.out = colMeans(length.out)
  sd.length.out = apply(length.out,2,sd)
  summary.out = rbind(ave.length.out,sd.length.out)
  colnames(summary.out)=c("NP.b0","NP.b1","NP.b2","BT.b0","BT.b1","BT.b2")
  rownames(summary.out)=c("Average","SD")
  cat("The means and standard deviations of the 90% confidence intervals are:","\n")
  summary.out
  sim.out <- rbind(Coverage,summary.out)
  out.array[,,m]<-sim.out
  cat("Iterations finished:",m,"\n")
}
  write.csv(out.array,file="Paired_Bootstrap_MC2.csv")
