setwd("./")
rm(list = ls())
require(mlr)
set.seed(123) 
load("../../data/processed/train_scaled.RData")
load("../../data/processed/tunecontrol.RData")
#inter outer - resampling desc is loaded from fixed file

source("./prauc.R")
self_measure<-makeMeasure(id="imbalance_measure_1",minimize = FALSE,properties = c("classif"),fun=prauc)


cart<-makeLearner("classif.rpart",predict.type="prob")
meas = list(mmce, brier,auc,fdr,tpr,tnr,ppv,mcc)
tuned_learners_full<-list(cart)

bmr_cart_full <- benchmark(tuned_learners_full,task,outer,measures=list(mmce, brier,auc,fdr,tpr,tnr,ppv,mcc,self_measure))

save(bmr_cart_full,file="../../data/ml/bmr_cart.RData")
