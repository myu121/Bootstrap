b0.sd.true = 1.21
b1.sd.true = 0.28
b2.sd.true = 1.18
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
library(ggplot2)
b0.name <- c(rep("WB+G",100),rep("WB+PM",100),rep("PB",100),rep("RW",100),rep("LM",100),rep("ND",100))
b0.name <- factor(b0.name,levels = c("WB+G","WB+PM","PB","RW","LM","ND"))
                  
b0.gg.df <- data.frame(name=b0.name,se=c(b0.ratio[,1],b0.ratio[,2],b0.ratio[,3],b0.ratio[,4],b0.ratio[,5],b0.ratio[,6]))
ggplot(data = b0.gg.df,aes(x=name,y=se))+scale_alpha(guide = 'none')+ylim(0,3)+
  geom_violin(trim=TRUE,alpha=0.5,aes(fill=name),lwd=1)+geom_boxplot(width=0.1,lwd=1)+scale_color_manual(values = cbp1) + geom_hline(yintercept = 1,linetype=2)+
  theme_bw()+xlab("Method")+ylab("Standard Error")+ggtitle("Variance Estimation of QR Estimators")


b1.gg.df <- data.frame(name=b0.name,se=c(b1.ratio[,1],b1.ratio[,2],b1.ratio[,3],b1.ratio[,4],b1.ratio[,5],b1.ratio[,6]))
ggplot(data = b1.gg.df,aes(x=name,y=se))+scale_alpha(guide = 'none')+ylim(0,4)+
  geom_violin(trim=TRUE,alpha=0.5,aes(fill=name),lwd=1)+geom_boxplot(width=0.1,lwd=1)+scale_color_manual(values = cbp1) + geom_hline(yintercept = 1,linetype=2)+
  theme_bw()+xlab("Method")+ylab("Standard Error")+ggtitle("Variance Estimation of QR Estimators")


b2.gg.df <- data.frame(name=b0.name,se=c(b2.ratio[,1],b2.ratio[,2],b2.ratio[,3],b2.ratio[,4],b2.ratio[,5],b2.ratio[,6]))
ggplot(data = b2.gg.df,aes(x=name,y=se))+scale_alpha(guide = 'none')+ylim(0,3)+
  geom_violin(trim=TRUE,alpha=0.5,aes(fill=name),lwd=1)+geom_boxplot(width=0.1,lwd=1)+scale_color_manual(values = cbp1) + geom_hline(yintercept = 1,linetype=2)+
  theme_bw()+xlab("Method")+ylab("Standard Error")+ggtitle("Variance Estimation of QR Estimators")


