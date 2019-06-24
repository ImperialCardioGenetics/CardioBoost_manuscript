## generate the prediction for all possible variants in myh7-associated genes

rm(list=ls())
load("../../../data/myh7_dcm/cm_all_rare_mutation.RData")
cm_all_rare$pathogenic<-as.factor(cm_all_rare$pathogenic)
cm_all_rare<-subset(cm_all_rare,gene=="MYH7")

load("../../../data/myh7_dcm/preprocess.RData")
source("../../src/myh7_preprocess_test.R")
source("../../src/predict.R")
load("../../../data/myh7_dcm/ml/train_ada.RData")
load("../../../data/myh7_dcm/train_scaled.RData")
#preprocess the file

##TODO: bug here
cm_prediction<-predict.cardioboost(train_ada,cm_all_rare)

cm_prediction<-data_imputed

#select the needed columns
cm_prediction<-subset(cm_prediction,select=c(CHROM,POS,REF,ALT,gene,HGVSc,HGVSp,AF_Adj,gnomAD_AF,MCAP,REVEL,pathogenicity,key))

train_scaled$key<-paste(train_scaled$CHROM,train_scaled$POS,train_scaled$REF,train_scaled$ALT,sep="_")
cm_prediction$isTrain<-is.element(cm_prediction$key,train_scaled$key)

cm_prediction$cDot<-sapply(as.character(cm_prediction$HGVSc),function(x){strsplit(x,split=":")[[1]][2]})
cm_prediction$pDot<-sapply(as.character(cm_prediction$HGVSp),function(x){strsplit(x,split=":")[[1]][2]})

#rearrange the columns
cm_prediction<-cm_prediction[,c("CHROM","POS","REF","ALT","gene","HGVSc","HGVSp","cDot","pDot","AF_Adj","gnomAD_AF","MCAP","REVEL","pathogenicity","isTrain","key")]
save(cm_prediction,file="../../../data/myh7_dcm/prediction/cm_prediction.RData")

