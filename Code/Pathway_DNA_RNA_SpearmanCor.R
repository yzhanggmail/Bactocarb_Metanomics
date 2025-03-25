## ---------------------------------------------------------------- ##
## By Yuzheng Zhang (yzhang@fredhutch.org) ##
## For Table S6, Table S7 and Supplemental Figure2 ##
## Output file: "Metagenomics_Metatranscriptomics_cor.csv", "Correlation_RNAvsDNA_LGL.png" and "Correlation_RNAvsDNA_HGL.png" ##
## ---------------------------------------------------------------- ##
rm(list=ls())
graphics.off()
library(plyr)  ## for rbind.fill ##
##setwd("Your_working_directory")
setwd("S:/yzhang/projects/Tim/Bactocarb2018/code_3timepoints/data_tobesent/Data")
setwd("S:/yzhang/projects/Tim/Bactocarb2018/data")
load("./Data_organized/Microbiome_Metagenomics_Metatranscpt_pathway_organized.Rdata")

dna.data <- data.metagn.tp
dna.pathinfo <- info.metagn.tp
rna.data <- data.metatrns.tp
rna.pathinfo <- info.metatrns.tp

## pathway info ##
path.info <- read.csv("MetaOmics_filtered_pathway_finallist_052021.csv")
## choose only core microbiome ##
path.info <- path.info[which(path.info$filter.final==0),]

table(path.info$partitioned_within_metatranscriptomics)
m <- match(path.info$pathway,colnames(dna.data))
dna.data <- dna.data[,m]
rna.data <- rna.data[,m]
table(path.info$pathway==colnames(dna.data))
table(path.info$pathway==colnames(rna.data))

dna.hgl <- dna.data[grep("HGL",rownames(dna.data)),]
dna.lgl <- dna.data[grep("LGL",rownames(dna.data)),]
rna.hgl <- rna.data[grep("HGL",rownames(rna.data)),]
rna.lgl <- rna.data[grep("LGL",rownames(rna.data)),]

select <- which(rownames(dna.hgl)%in%rownames(rna.hgl))
dna.hgl <- dna.hgl[select,]
select <- which(rownames(dna.lgl)%in%rownames(rna.lgl))
dna.lgl <- dna.lgl[select,]

## Weighted Spearman Rank Correlation ##
cor.hgl <- sapply(1:ncol(dna.hgl),function(x) weightedCorr(dna.hgl[!is.na(dna.hgl[,x])&!is.na(rna.hgl[,x]),x],rna.hgl[!is.na(dna.hgl[,x])&!is.na(rna.hgl[,x]),x],method="Spearman"))
cor.lgl <- sapply(1:ncol(dna.lgl),function(x) weightedCorr(dna.lgl[!is.na(dna.lgl[,x])&!is.na(rna.lgl[,x]),x],rna.lgl[!is.na(dna.lgl[,x])&!is.na(rna.lgl[,x]),x],method="Spearman"))

table(cor.hgl=="NaN")
select <- which(cor.hgl!="NaN")
dna.hgl <- dna.hgl[,select]
cor.hgl <- cor.hgl[select]
select <- which(cor.lgl!="NaN")
dna.lgl <- dna.lgl[,select]
cor.lgl <- cor.lgl[select]
path.info <- path.info[select,]
table(path.info$pathway==colnames(dna.hgl))
table(path.info$pathway==colnames(dna.lgl))
output <- data.frame(PHWY=path.info$PHWY.Name,cor.hgl,cor.lgl,pathway.info=path.info$pathway)
output <- output[order(output$PHWY),]
write.csv(output,row.names=F,file="Metagenomics_Metatranscriptomics_cor.csv")

vcol<-c("tomato", "royalblue", "mediumpurple", "lightpink4", "tan",
          "maroon", "turquoise", "firebrick", "gold", "purple",
          "cadetblue2", "wheat4", "orangered", "blue4","darkorchid",
          "skyblue", "brown2", "blue1", "seagreen", "chartreuse",
          "violetred", "darkseagreen4", "hotpink","cyan","lightyellow","lightpink")

cols <- factor(output$PHWY)
levels(cols) <- vcol
names(cols) <- output$PHWY
cbind(levels(factor(output$PHWY)),levels(cols))

png("../graphs/Correlation_RNAvsDNA_LGL.png",height=800,width=1500)
par(omd=c(0,1,0.2,1),cex=1.6,font=2,font.lab=2,font.axis=2)
barplot(output$cor.lgl,xaxt="n",ylim=c(-0.3,1),col=as.character(cols),ylab="Weighted Spearman Correlation",main="LGL")
dev.off()

png("../graphs/Correlation_RNAvsDNA_HGL.png",height=800,width=1500)
par(omd=c(0,1,0.2,1),cex=1.6,font=2,font.lab=2,font.axis=2)
barplot(output$cor.hgl,las=3,ylim=c(-0.3,1),names=output$PHWY,col=as.character(cols),ylab="Weighted Spearman Correlation",main="HGL")
dev.off()
