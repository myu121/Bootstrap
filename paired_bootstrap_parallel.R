library(foreach)
library(parallel)
library(doParallel)
library(tictoc)
library(doRNG)
sample.size = c(100,200,400,800,1600,3200,5000)
out.array = array(0,c(3,3,2))
for(m in 1:7){
  cl <-makeCluster(detectCores())
  registerDoParallel(cl)
  MC_num = 2500
  n = sample.size[m]
  boot_num = 1000
  beta <- c(1,1,1)
  qreg_fun <- function(df,ind){
    df$y <- jitter(df$y,1e-8)
    while((sum(df[ind,]$x2==1) == n)||(sum(df[ind,]$x2==1) == 0)){
      ind = sample(ind,replace=TRUE)
    }
    fit <- summary(rq(y~x1+x2,tau=0.5,data=df[ind,]))
    boot_coef = as.numeric(fit$coefficients[,1])
    return(boot_coef)
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
      Per.b2 =  as.numeric(quantile(boot.out$t[,3],c(0.05,0.95)))
    )
    par_results <- list(par_results)
    names(par_results) <- paste0("step_", i)
    return(par_results)
  }
  toc()
  stopCluster(cl)
  coverage.prob = rep(0,3)
  names(coverage.prob)=c("NP.b0","NP.b1","NP.b2")
  for(i in 1:MC_num){
    if((output[[i]]$Per.b0[1]<beta[1])&&(output[[i]]$Per.b0[2]>beta[1])) coverage.prob[1]=coverage.prob[1]+1
    if((output[[i]]$Per.b1[1]<beta[2])&&(output[[i]]$Per.b1[2]>beta[2])) coverage.prob[2]=coverage.prob[2]+1
    if((output[[i]]$Per.b2[1]<beta[3])&&(output[[i]]$Per.b2[2]>beta[3])) coverage.prob[3]=coverage.prob[3]+1
  }
  Coverage = coverage.prob/MC_num
  cat("The coverage probablities for the 90% confidence intervals are:","\n")
  Coverage
  ci.out = matrix(0,ncol=6,nrow=MC_num)
  for(i in 1:MC_num){
    ci.out[i,1] = output[[i]]$Per.b0[1]
    ci.out[i,2] = output[[i]]$Per.b0[2]
    ci.out[i,3] = output[[i]]$Per.b1[1]
    ci.out[i,4] = output[[i]]$Per.b1[2]
    ci.out[i,5] = output[[i]]$Per.b2[1]
    ci.out[i,6] = output[[i]]$Per.b2[2]
  }
  length.out = matrix(0,ncol=3,nrow=MC_num)
  for(i in 1:3){
    length.out[,i] <- ci.out[,2*i]-ci.out[,2*i-1] 
  }
  ave.length.out = colMeans(length.out)
  sd.length.out = apply(length.out,2,sd)
  summary.out = rbind(ave.length.out,sd.length.out)
  colnames(summary.out)=c("NP.b0","NP.b1","NP.b2")
  rownames(summary.out)=c("Average","SD")
  cat("The means and standard deviations of the 90% confidence intervals are:","\n")
  summary.out
  sim.out <- rbind(Coverage,summary.out)
  out.array[,,m]<-sim.out
  cat("Iterations finished:",m,"\n")
}
write.csv(out.array,file="Paired_Bootstrap_MC2.csv")
