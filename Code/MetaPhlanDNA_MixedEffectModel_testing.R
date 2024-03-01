## ---------------------------------------------------------------- ##
## By Yuzheng Zhang (yzhang@fredhutch.org) ##
## Fit Figure1C, Mixed-effect model "species  ~ diet + age + sex + bodyfat + kcal + diet.seq + (1|ID) on Phlan species DNA data ##
## Use files: "MetaPhlanDNA_organized_filtered_log2scale.csv" and "BactoCarb_Clinical_info_3timepoint.csv" ##
## Output file: "PhlanSpecies_DNA_MixedEffectModel_dietLGLvsHGL.csv" ##
## ---------------------------------------------------------------- ##
rm(list=ls())
library(lme4) ## for lme
library(lmerTest) ## for lme p value
library(Rcpp)
graphics.off()
setwd("Your_working_directory")

## Phlan species RNA/DNA ratios ##
data.log <- read.csv("MetaPhlanDNA_organized_filtered_log2scale.csv")
rownames(data.log) <- data.log[,1]
colnames(data.log) <- substr(colnames(data.log),2,nchar(colnames(data.log)))
data.log <- data.log[,-1]
species.info <- data.frame(meta.name=paste("Meta",1:ncol(data.log),sep=""),name=colnames(data.log))

info <- read.csv("BactoCarb_Clinical_info_3timepoint.csv")
info <- info[!info$diet2%in%"B",]
info$Index <- paste(info$PptID,info$diet2,sep=".")
table(info$Index%in%colnames(data.log))

info <- info[info$Index%in%colnames(data.log),]
dim(info)
m <- match(colnames(data.log),info$Index)
info <- info[m,]
table(info$Index==colnames(data.log))

info$kcalcat <- factor(info$kcalcat)
pvalue <- coef <- NA
i <- 1
for(i in 1:nrow(data.log)){
    m <- as.numeric(data.log[i,])
    mylm <- lmer(m ~ diet2 + Diet_sequence + female + age + X.CL_DXA_FAT_CAT_0low_1high + kcalcat+ (1|PptID), data=info)
    summary(mylm)$coef
    pvalue[i] <- summary(mylm)$coef["diet2LGL",5]
    coef[i] <- summary(mylm)$coef["diet2LGL",1]
}
qvalue <- p.adjust(pvalue,method="BH")

table(species.info$meta.name==rownames(data.log))
result <- data.frame(Name=rownames(data.log),coef=coef,pvalue=pvalue,qvalue=qvalue)
result <- result[order(result$pvalue),]

write.csv(result,row.names=F,file="../Result/PhlanSpecies_DNA_MixedEffectModel_dietLGLvsHGL.csv")

par(mfrow=c(2,2),font=2)
cols <- ifelse(result$pvalue<0.05&result$coef>0,"darkorange","black")
cols <- ifelse(result$pvalue<0.05&result$coef<0,"blue",as.character(cols))
cexs <- ifelse(result$pvalue<0.05,1.6,1)
plot(as.numeric(as.character(result$coef)),-1*log10(as.numeric(as.character(result$pvalue))),cex=cexs,pch=20,xlim=c(-1,1.5),col=cols,xlab="coef",ylab="-1*log10(p)",main="LGL vs. HGL")
abline(h=-1*log10(0.05),col="blue")
legend("top",bty='n',fill=c("darkorange","blue"),legend=c(paste(sum(as.numeric(as.character(result$coef))>0&as.numeric(as.character(result$pvalue))<0.05,na.rm=T)," LGL",sep=""),paste(sum(as.numeric(as.character(result$coef))<0&as.numeric(as.character(result$pvalue))<0.05,na.rm=T)," HGL",sep="")))
hist(as.numeric(as.character(result$pvalue)),breaks=20,col='lightblue',xlab="p",main="Species  ~ diet + age + sex + \nbodyfat + kcal +  diet.seq + (1|ID)")
legend("topright",bty='n',legend=paste(sum(as.numeric(as.character(result$pvalue))<0.05,na.rm=T)," sp p<0.05 over ",nrow(result)," (",round(sum(as.numeric(as.character(result$pvalue))<0.05,na.rm=T)/nrow(result)*100,1),"%)",sep=""))
abline(v=0.05,col='blue',lwd=2)

png(file="../Result/Figure1_C.png",width=600,height=600)
par(mfrow=c(1,1),cex=2,font=2,font.lab=2,font.axis=2)
cols <- ifelse(result$pvalue<0.05&result$coef>0,"darkorange","black")
cols <- ifelse(result$pvalue<0.05&result$coef<0,"blue",as.character(cols))
cexs <- ifelse(result$pvalue<0.05,1.6,1)
plot(as.numeric(as.character(result$coef)),-1*log10(as.numeric(as.character(result$pvalue))),cex=cexs,pch=20,xlim=c(-1,1.5),col=cols,xlab="coef",ylab="-1*log10(p)",main="")
abline(h=-1*log10(0.05),col="blue")
legend("top",bty='n',fill=c("darkorange","blue"),legend=c(paste(sum(as.numeric(as.character(result$coef))>0&as.numeric(as.character(result$pvalue))<0.05,na.rm=T)," LGL",sep=""),paste(sum(as.numeric(as.character(result$coef))<0&as.numeric(as.character(result$pvalue))<0.05,na.rm=T)," HGL",sep="")))
dev.off()
