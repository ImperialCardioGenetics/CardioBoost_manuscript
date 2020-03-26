setwd("./")
rm(list = ls())

require(mlr)
set.seed(123)
load("../../../data/arrhythmia/train_scaled.RData")
load("../../../data/arrhythmia/tunecontrol.RData")
#inter outer - resampling desc is loaded from fixed file

source("../../src/prauc.R")
self_measure<-makeMeasure(id="imbalance_measure_1",minimize = FALSE,properties = c("classif"),fun=prauc)

#SVM tunning using the pratical guide to support vector classification
tune_control_svm <- makeTuneControlRandom(maxit=200)

svm<-makeLearner("classif.ksvm",predict.type="prob",kernel="rbfdot")
svm_grid <- makeParamSet(makeDiscreteParam("C", values=2^seq(-3,3,1)),
                         makeDiscreteParam("sigma", values=2^seq(-8,-1,0.5)))

svm_tuned <- makeTuneWrapper(svm, resampling = inter,
                             par.set = svm_grid, control = tune_control_svm,measures=list(self_measure))

meas = list(mmce,brier,auc,fdr,tpr,tnr,ppv,mcc,self_measure)
tuned_learners_full<-list(svm_tuned)

bmr_svm_full <- benchmark(tuned_learners_full,task,outer,meas)
save(bmr_svm_full,file="../../../data/arrhythmia/ml/bmr_svm_RBF.RData")