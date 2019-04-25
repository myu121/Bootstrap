library(ggplot2)
library(Cairo)
library(gridExtra)
library(vioplot)
b0.sd.true = 1.21
b1.sd.true = 0.28
b2.sd.true = 1.18
#b0.sd.true = 0.66     (tau=0.9)
#b1.sd.true = 0.125
#b2.sd.true = 0.66
wb1.se = read.csv("SE1.csv")[,-1]
wb2.se = read.csv("SE2.csv")[,-1]
rw.se = read.csv("rw_se.csv")
pb.se = read.csv("Paired_Bootstrap_SE.csv")[,5:7]
lm.se = read.csv("100SdEstimates.csv")[,2:4]
nd.se = read.csv("100SdEstimates.csv")[,5:7]
b0.se <- data.frame(WB1 = wb1.se[,1],WB2 = wb2.se[,1],PB = pb.se[,1],RW = rw.se[,1])
b1.se <- data.frame(WB1 = wb1.se[,2],WB2 = wb2.se[,2],PB = pb.se[,2],RW = rw.se[,2])
b2.se <- data.frame(WB1 = wb1.se[,3],WB2 = wb2.se[,3],PB = pb.se[,3],RW = rw.se[,3])
b0.ratio <- b0.se/b0.sd.true
b1.ratio <- b1.se/b1.sd.true
b2.ratio <- b2.se/b2.sd.true
b0.ratio <- data.frame(b0.ratio, LM = lm.se[,1], ND = nd.se[,1])
b1.ratio <- data.frame(b1.ratio, LM = lm.se[,2], ND = nd.se[,2])
b2.ratio <- data.frame(b2.ratio, LM = lm.se[,3], ND = nd.se[,3])
cbp1 <- c("#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")

b0.name <- c(rep("WB+G",100),rep("WB+PM",100),rep("PB",100),rep("RW",100),rep("LM",100),rep("ND",100))
b0.name <- factor(b0.name,levels = c("WB+G","WB+PM","PB","RW","LM","ND"))

b0.gg.df <- data.frame(Method=b0.name,se=c(b0.ratio[,1],b0.ratio[,2],b0.ratio[,3],b0.ratio[,4],b0.ratio[,5],b0.ratio[,6]))
p1<-ggplot(data = b0.gg.df,aes(x=Method,y=se))+scale_alpha(guide = 'none')+ylim(0,3.5)+
  geom_violin(trim=TRUE,alpha=0.5,aes(fill=Method))+geom_boxplot(width=0.1)+scale_color_manual(values = cbp1) + geom_hline(yintercept = 1,linetype=2)+
  theme_bw()+labs(x="Method",y="SE Ratio",title=TeX('$\\hat{\\beta}_0(0.5)$'))+ 
  guides(fill="none")+theme(axis.text=element_text(size=15),plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=15,face="bold"))


b1.gg.df <- data.frame(Method=b0.name,se=c(b1.ratio[,1],b1.ratio[,2],b1.ratio[,3],b1.ratio[,4],b1.ratio[,5],b1.ratio[,6]))
p2<-ggplot(data = b1.gg.df,aes(x=Method,y=se))+scale_alpha(guide = 'none')+ylim(0,3.5)+
  geom_violin(trim=TRUE,alpha=0.5,aes(fill=Method))+geom_boxplot(width=0.1)+scale_color_manual(values = cbp1) + geom_hline(yintercept = 1,linetype=2)+
  theme_bw()+labs(x="Method",y="SE Ratio",title=TeX('$\\hat{\\beta}_1(0.5)$'))+ 
  guides(fill="none")+theme(axis.text=element_text(size=15),plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=15,face="bold"))


b2.gg.df <- data.frame(Method=b0.name,se=c(b2.ratio[,1],b2.ratio[,2],b2.ratio[,3],b2.ratio[,4],b2.ratio[,5],b2.ratio[,6]))
p3<-ggplot(data = b2.gg.df,aes(x=Method,y=se))+scale_alpha(guide = 'none')+ylim(0,3.5)+
  geom_violin(trim=TRUE,alpha=0.5,aes(fill=Method))+geom_boxplot(width=0.1)+scale_color_manual(values = cbp1) + geom_hline(yintercept = 1,linetype=2)+
  theme_bw()+labs(x="Method",y="SE Ratio",title=TeX('$\\hat{\\beta}_2(0.5)$'))+ 
  guides(fill="none")+theme(axis.text=element_text(size=15),plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=15,face="bold"))



Cairo(file="se_comparison.png", 
      type="png",
      width=640,height=1280,pointsize = 14*1.3,dpi=72*1.3)
grid.arrange(p1, p2,p3, ncol=1)
dev.off()


Cairo(file="se_comparison_t.png", 
      type="png",
      width=1720,height=480,pointsize = 14*1.3,dpi=72*1.3)
grid.arrange(p1, p2,p3, ncol=3)
dev.off()


ci.coverage <- read.csv("CI_cov_com.csv")
ci.df <- data.frame(n=as.numeric(ci.coverage[,1]),b0=ci.coverage$b0,b1=ci.coverage$b1,b2=ci.coverage$b2,Method=ci.coverage$Method)
p4<-ggplot(data=ci.df,aes(x=n,y=b0))+geom_line(aes(color=Method,linetype=Method,alpha=0.9),lwd=1.2)+
  scale_linetype_manual(values=c(rep("solid",5),"dashed","dashed"))+
  geom_hline(yintercept=0.9,linetype="twodash",lwd=1.2,alpha=0.5)+scale_alpha(guide = 'none')+labs(x="Sample Size",y="Coverage")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"))
p4
p5<-ggplot(data=ci.df,aes(x=n,y=b1))+geom_line(aes(color=Method,linetype=Method,alpha=0.9),lwd=1.2)+
  scale_linetype_manual(values=c(rep("solid",5),"dashed","dashed"))+
  geom_hline(yintercept=0.9,linetype="twodash",lwd=1.2,alpha=0.5)+scale_alpha(guide = 'none')+labs(x="Sample Size",y="Coverage")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"))
p5
p6<-ggplot(data=ci.df,aes(x=n,y=b2))+geom_line(aes(color=Method,linetype=Method,alpha=0.9),lwd=1.2)+
  scale_linetype_manual(values=c(rep("solid",5),"dashed","dashed"))+
  geom_hline(yintercept=0.9,linetype="twodash",lwd=1.2,alpha=0.5)+scale_alpha(guide = 'none')+labs(x="Sample Size",y="Coverage")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"))
p6


Cairo(file="coverage_trend_b2.png", 
      type="png",
      width=840,height=360,pointsize = 14*1.3,dpi=72*1.3)
p6
dev.off()
