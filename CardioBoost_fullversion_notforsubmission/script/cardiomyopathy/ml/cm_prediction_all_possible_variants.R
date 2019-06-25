## generate the prediction for all possible variants in cardiomyopathy-associated genes

rm(list=ls())
load("../../../../data/cardiomyopathy/cm_all_rare_mutation.RData")
cm_all_rare$pathogenic<-as.factor(cm_all_rare$pathogenic)

load("../../../../data/cardiomyopathy/preprocess.RData")
source("../../../src/preprocess_test.R")
source("../../../src/predict.R")
load("../../../../data/cardiomyopathy/ml/train_ada.RData")
load("../../../../data/cardiomyopathy/train_scaled.RData")
#preprocess the file

cm_prediction<-predict.cardioboost(train_ada,cm_all_rare)

#select the needed columns
cm_prediction<-subset(cm_prediction,select=c(CHROM,POS,REF,ALT,gene,HGVSc,HGVSp,AF_Adj,gnomAD_AF,MCAP,REVEL,pathogenicity,key))

train_scaled$key<-paste(train_scaled$CHROM,train_scaled$POS,train_scaled$REF,train_scaled$ALT,sep="_")
cm_prediction$isTrain<-is.element(cm_prediction$key,train_scaled$key)

cm_prediction$cDot<-sapply(as.character(cm_prediction$HGVSc),function(x){strsplit(x,split=":")[[1]][2]})
cm_prediction$pDot<-sapply(as.character(cm_prediction$HGVSp),function(x){strsplit(x,split=":")[[1]][2]})

#rearrange the columns
cm_prediction<-cm_prediction[,c("CHROM","POS","REF","ALT","gene","HGVSc","HGVSp","cDot","pDot","AF_Adj","gnomAD_AF","pathogenicity","isTrain","key")]
save(cm_prediction,file="../../../../data/cardiomyopathy/prediction/cm_prediction.RData")

