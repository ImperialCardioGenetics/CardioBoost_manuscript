setwd("./")
rm(list = ls())
require(mlr)
set.seed(123) 
load("../../data/processed/train_scaled.RData")
load("../../data/processed/tunecontrol.RData")
#inter outer - resampling desc is loaded from fixed file

source("./prauc.R")
self_measure<-makeMeasure(id="imbalance_measure_1",minimize = FALSE,properties = c("classif"),fun=prauc)


#SVM tunning using the pratical guide to support vector classification: https://www.csie.ntu.edu.tw/~cjlin/papers/guide/guide.pdf
tune_control_svm <- makeTuneControlGrid()

svm<-makeLearner("classif.ksvm",predict.type="prob",kernel="rbfdot")
svm_grid <- makeParamSet(makeDiscreteParam("C", values=2^seq(-10,10,1)),
                         makeDiscreteParam("sigma", values=2^seq(-10,10,1)))

svm_tuned <- makeTuneWrapper(svm, resampling = inter,
                            par.set = svm_grid, control = tune_control_svm,measures=list(self_measure))

meas = list(mmce,brier,auc,fdr,tpr,tnr,ppv,mcc,self_measure)
tuned_learners_full<-list(svm_tuned)

bmr_svm_full <- benchmark(tuned_learners_full,task,outer,meas)
save(bmr_svm_full,file="../../data/ml/bmr_svm_RBF.RData")
