## ---------------------------------------------------------------- ##
## By Yuzheng Zhang (yzhang@fredhutch.org) ##
## Fit Table2 & Figure2A, Mixed-effect model "pathway  ~ diet + age + sex + bodyfat + kcal + diet.seq + (1|ID) ##
## Use files: "Pathway_organized_filtered_logratio.csv" and "BactoCarb_Clinical_info_3timepoint.csv" ##
## Output file: "MetaOmicsPathway_RNADNAratio_MixedEffectModel_dietLGLvsHGL.csv" ##
## ---------------------------------------------------------------- ##
rm(list=ls())
library(lme4) ## for lme
library(lmerTest) ## for lme p value
library(Rcpp)
graphics.off()
setwd("Your_working_directory")

data.logratio <- read.csv("Pathway_organized_filtered_logratio.csv")
rownames(data.logratio) <- data.logratio[,1]
data.logratio <- data.logratio[,-1]
colnames(data.logratio) <- substr(colnames(data.logratio),2,nchar(colnames(data.logratio)))

info <- read.csv("BactoCARB_Clinical_Info_3timepoint.csv")
info <- info[!info$diet2%in%"B",]
info$Index <- paste(info$PptID,info$diet2,sep=".")
table(info$Index%in%colnames(data.logratio))

info <- info[info$Index%in%colnames(data.logratio),]
dim(info)
m <- match(colnames(data.logratio),info$Index)
info <- info[m,]
table(info$Index==colnames(data.logratio))

## standardize to HGL mean 0 sd 1 ##
norm <- t(apply(data.logratio,1,function(x) (x-mean(x[info$diet2=="HGL"],na.rm=T))/sd(x[info$diet2=="HGL"],na.rm=T)))
mean(norm[1,info$diet2=="HGL"],na.rm=T)
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
write.csv(result,row.names=F,file="../Result/MetaOmicsPathway_RNADNAratio_MixedEffectModel_dietLGLvsHGL.csv")

par(mfrow=c(2,2),font=2)
cols <- ifelse(result$pvalue<0.05&result$coef>0,"darkorange","black")
cols <- ifelse(result$pvalue<0.05&result$coef<0,"blue",as.character(cols))
cexs <- ifelse(result$pvalue<0.05,1.6,1)
plot(as.numeric(as.character(result$coef)),-1*log10(as.numeric(as.character(result$pvalue))),cex=cexs,pch=20,col=cols,xlab="coef",ylab="-1*log10(p)",main="LGL vs. HGL")
abline(h=-1*log10(0.05),col="blue")
legend("top",bty='n',fill=c("darkorange","blue"),legend=c(paste(sum(as.numeric(as.character(result$coef))>0&as.numeric(as.character(result$pvalue))<0.05,na.rm=T)," LGL",sep=""),paste(sum(as.numeric(as.character(result$coef))<0&as.numeric(as.character(result$pvalue))<0.05,na.rm=T)," HGL",sep="")))
hist(as.numeric(as.character(result$pvalue)),breaks=20,col='lightblue',xlab="p",main="Pathway  ~ diet + age + sex + \nbodyfat + kcal +  diet.seq + (1|ID)")
legend("topright",bty='n',legend=paste(sum(as.numeric(as.character(result$pvalue))<0.05,na.rm=T)," pw p<0.05 over ",nrow(result)," (",round(sum(as.numeric(as.character(result$pvalue))<0.05,na.rm=T)/nrow(result)*100,1),"%)",sep=""))
abline(v=0.05,col='blue',lwd=2)


png(file="../Result/Figure2_A.png",width=600,height=600)
par(mfrow=c(1,1),cex=2,font=2,font.lab=2,font.axis=2)
cols <- ifelse(result$pvalue<0.05&result$coef>0,"darkorange","black")
cols <- ifelse(result$pvalue<0.05&result$coef<0,"blue",as.character(cols))
cexs <- ifelse(result$pvalue<0.05,1.6,1)
plot(as.numeric(as.character(result$coef)),-1*log10(as.numeric(as.character(result$pvalue))),cex=cexs,pch=20,col=cols,xlab="coef",ylab="-1*log10(p)",main="")
abline(h=-1*log10(0.05),col="blue")
legend("top",bty='n',fill=c("darkorange","blue"),legend=c(paste(sum(as.numeric(as.character(result$coef))>0&as.numeric(as.character(result$pvalue))<0.05,na.rm=T)," LGL",sep=""),paste(sum(as.numeric(as.character(result$coef))<0&as.numeric(as.character(result$pvalue))<0.05,na.rm=T)," HGL",sep="")))
dev.off()


## ---------------------------------------------------------------- ##
## By Yuzheng Zhang (yzhang@fredhutch.org) ##
## Fit Table4, HOMA  ~ pathway + (age + sex + bodyfat + kcal + diet.seq + (1|ID)) ##
## Use files: "Pathway_organized_filtered_logratio.csv" and "BactoCarb_Clinical_info_3timepoint.csv" ##
## Output file: "MetaOmics_RNADNAratio_MixedEffectModel_HOMA.csv" ##
## Output variables:
## coef.mLGL:	coef of pathway for LGL (beta1+beta3);
## coef.mHGL:	coef of pathway for HGL (beta1);
## coef.mdiet:	coef for interactive term "diet*path";
## p.mdiet:	p value for interactive term "diet*path";
## q.mdiet:	BH adjustment of "p.mdiet";
## coef.m:	coef for term "pathway";
## p.m:	p value for term "pathway";
## coef.diet:	coef of "diet";
## p.diet:	p value of "diet";
## ---------------------------------------------------------------- ##
rm(list=ls())
library(lme4) ## for lme
library(lmerTest) ## for lme p value
graphics.off()

data.logratio <- read.csv("Pathway_organized_filtered_logratio.csv")
rownames(data.logratio) <- data.logratio[,1]
data.logratio <- data.logratio[,-1]
colnames(data.logratio) <- substr(colnames(data.logratio),2,nchar(colnames(data.logratio)))

info <- read.csv("BactoCARB_Clinical_Info_3timepoint.csv")
info <- info[!info$diet2%in%"B",]
info$Index <- paste(info$PptID,info$diet2,sep=".")
table(info$Index%in%colnames(data.logratio))

info <- info[info$Index%in%colnames(data.logratio),]
dim(info)
m <- match(colnames(data.logratio),info$Index)
info <- info[m,]
table(info$Index==colnames(data.logratio))

info$kcalcat <- factor(info$kcalcat)
info$logHOMA <- log2(info$HOMA)
## interact term ##
p.diet <- coef.diet <- coef.mHGL <- coef.mLGL <- p.m <- p.mdiet <- coef.m <- coef.mdiet <- NA
i <- 1
for(i in 1:nrow(data.logratio)){
    m <- as.numeric(data.logratio[i,])
    mylm <- lmer(logHOMA ~ m+diet2+diet2*m + Diet_sequence + female + age + X.CL_DXA_FAT_CAT_0low_1high + kcalcat+(1|PptID) , data=info)
    summary(mylm)$coef
    p.m[i] <- summary(mylm)$coef["m",5]
    coef.m[i] <- summary(mylm)$coef["m",1]
    p.diet[i] <- summary(mylm)$coef["diet2LGL",5]
    coef.diet[i] <- summary(mylm)$coef["diet2LGL",1]
    coef.mdiet[i] <- summary(mylm)$coef["m:diet2LGL",1]
    p.mdiet[i] <- summary(mylm)$coef["m:diet2LGL",5]
    coef.mLGL[i] <- coef.m[i]+coef.mdiet[i]
    coef.mHGL[i] <- coef.m[i]
}
q.mdiet <- p.adjust(p.mdiet,"BH")
q.m <- p.adjust(p.m,"BH")

result <- data.frame(Name=rownames(data.logratio),coef.mLGL,coef.mHGL,coef.mdiet,p.mdiet,q.mdiet,coef.m,p.m,q.m,coef.diet,p.diet)
result <- result[order(result$p.m),]
write.csv(result,row.names=F,file="../Result/MetaOmicsPathway_RNADNAratio_MixedEffectModel_HOMA.csv")
