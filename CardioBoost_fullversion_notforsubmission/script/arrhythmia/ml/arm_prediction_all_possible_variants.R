## generate the prediction for all possible variants in cardiomyopathy-associated genes
rm(list=ls())
load("../../../data/arrhythmia/arm_all_rare_mutation.RData")
arm_all_rare$pathogenic<-as.factor(arm_all_rare$pathogenic)

load("../../../data/arrhythmia/preprocess.RData")
source("../../src/preprocess_test.R")
source("../../src/predict.R")
load("../../../data/arrhythmia/ml/train_ada.RData")
load("../../../data/arrhythmia/train_scaled.RData")
#preprocess the file

arm_prediction<-predict.cardioboost(train_ada,arm_all_rare)

#select the needed columns
arm_prediction<-subset(arm_prediction,select=c(CHROM,POS,REF,ALT,gene,HGVSc,HGVSp,AF_Adj,gnomAD_AF,MCAP,REVEL,pathogenicity,key))

train_scaled$key<-paste(train_scaled$CHROM,train_scaled$POS,train_scaled$REF,train_scaled$ALT,sep="_")
arm_prediction$isTrain<-is.element(arm_prediction$key,train_scaled$key)

arm_prediction$cDot<-sapply(as.character(arm_prediction$HGVSc),function(x){strsplit(x,split=":")[[1]][2]})
arm_prediction$pDot<-sapply(as.character(arm_prediction$HGVSp),function(x){strsplit(x,split=":")[[1]][2]})

#rearrange the columns
arm_prediction<-arm_prediction[,c("CHROM","POS","REF","ALT","gene","HGVSc","HGVSp","cDot","pDot","AF_Adj","gnomAD_AF","MCAP","REVEL","pathogenicity","isTrain","key")]
save(arm_prediction,file="../../../data/arrhythmia/prediction/arm_prediction.RData")



