## ---------------------------------------------------------------- ##
## By Yuzheng Zhang (yzhang@fredhutch.org) ##
## Data: CAZymes species RNA/DNA ratios ##
## Fit Table3 & Figure3B, Mixed-effect model "species  ~ diet + age + sex + bodyfat + kcal + diet.seq + (1|ID) ##
## Use files: "Cazymes_organized_filtered_logratio.csv" and "BactoCarb_Clinical_info_3timepoint.csv" ##
## Output file: "CAZymes_RNADNAratio_MixedEffectModel_dietLGLvsHGL.csv" ##
## ---------------------------------------------------------------- ##
rm(list=ls())
library(lme4) ## for lme
library(lmerTest) ## for lme p value
graphics.off()
setwd("Your_working_directory")
## CAZymes species RNA/DNA ratios ##
data.logratio <- read.csv("Cazymes_organized_filtered_logratio.csv")
rownames(data.logratio) <- data.logratio[,1]
data.logratio <- data.logratio[,-1]
colnames(data.logratio) <- substr(colnames(data.logratio),2,nchar(colnames(data.logratio)))

info <- read.csv("BactoCARB_Clinical_Info_3timepoint.csv")
info <- info[!info$diet2%in%"B",]
info$Index <- paste(info$PptID,info$diet2,sep="_")
table(info$Index%in%colnames(data.logratio))

info <- info[info$Index%in%colnames(data.logratio),]
m <- match(colnames(data.logratio),info$Index)
info <- info[m,]
table(info$Index==colnames(data.logratio))

## Genewise standardize to HGL mean 0 sd 1 ##
norm <- t(apply(data.logratio,1,function(x) (x-mean(x[info$diet2=="HGL"],na.rm=T))/sd(x[info$diet2=="HGL"],na.rm=T)))
mean(norm[1,info$diet2=="HGL"])
data.logratio <- norm

info$kcalcat <- factor(info$kcalcat)
pvalue <- coef <- NA
i <- 1
for(i in 1:nrow(data.logratio)){
    m <- as.numeric(data.logratio[i,])
    mylm <- lmer(m ~ diet2 + Diet_sequence + female + age + X.CL_DXA_FAT_CAT_0low_1high + kcalcat+ (1|PptID), data=info)
    summary(mylm)$coef
    pvalue[i] <- summary(mylm)$coef["diet2LGL",5]
    coef[i] <- summary(mylm)$coef["diet2LGL",1]
}
qvalue <- p.adjust(pvalue,method="BH")
result <- data.frame(Name=rownames(data.logratio),coef=coef,pvalue=pvalue,qvalue=qvalue)
result <- result[order(result$pvalue),]
write.csv(result,row.names=F,file="../Result/CAZymes_RNADNAratio_MixedEffectModel_dietLGLvsHGL.csv")

## Figure 3B ##
par(mfrow=c(2,2),font=2)
cols <- ifelse(result$pvalue<0.05&result$coef>0,"darkorange","black")
cols <- ifelse(result$pvalue<0.05&result$coef<0,"blue",as.character(cols))
cexs <- ifelse(result$pvalue<0.05,1.6,1)
plot(as.numeric(as.character(result$coef)),-1*log10(as.numeric(as.character(result$pvalue))),cex=cexs,pch=20,col=cols,xlab="coef",ylab="-1*log10(p)",main="LGL vs. HGL")
abline(h=-1*log10(0.05),col="blue")
legend("top",bty='n',fill=c("darkorange","blue"),legend=c(paste(sum(as.numeric(as.character(result$coef))>0&as.numeric(as.character(result$pvalue))<0.05,na.rm=T)," LGL",sep=""),paste(sum(as.numeric(as.character(result$coef))<0&as.numeric(as.character(result$pvalue))<0.05,na.rm=T)," HGL",sep="")))
hist(as.numeric(as.character(result$pvalue)),breaks=20,col='lightblue',xlab="p",main="Species  ~ diet + age + sex + \nbodyfat + kcal +  diet.seq + (1|ID)")
legend("topright",bty='n',legend=paste(sum(as.numeric(as.character(result$pvalue))<0.05,na.rm=T)," sp p<0.05 over ",nrow(result)," (",round(sum(as.numeric(as.character(result$pvalue))<0.05,na.rm=T)/nrow(result)*100,1),"%)",sep=""))
abline(v=0.05,col='blue',lwd=2)

png(file="../Result/Figure3_B.png",width=600,height=600)
par(mfrow=c(1,1),cex=2,font=2,font.lab=2,font.axis=2)
cols <- ifelse(result$pvalue<0.05&result$coef>0,"darkorange","black")
cols <- ifelse(result$pvalue<0.05&result$coef<0,"blue",as.character(cols))
cexs <- ifelse(result$pvalue<0.05,1.6,1)
plot(as.numeric(as.character(result$coef)),-1*log10(as.numeric(as.character(result$pvalue))),cex=cexs,pch=20,xlim=c(-1.2,1.7),col=cols,xlab="coef",ylab="-1*log10(p)",main="")
abline(h=-1*log10(0.05),col="blue")
legend("top",bty='n',fill=c("darkorange","blue"),legend=c(paste(sum(as.numeric(as.character(result$coef))>0&as.numeric(as.character(result$pvalue))<0.05,na.rm=T)," LGL",sep=""),paste(sum(as.numeric(as.character(result$coef))<0&as.numeric(as.character(result$pvalue))<0.05,na.rm=T)," HGL",sep="")))
dev.off()
