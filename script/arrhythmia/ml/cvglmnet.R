setwd("./")
rm(list = ls())

require(mlr)
set.seed(123)
load("../../../data/arrhythmia/train_scaled.RData")
load("../../../data/arrhythmia/tunecontrol.RData")
#inter outer - resampling desc is loaded from fixed file

source("../../src/prauc.R")
self_measure<-makeMeasure(id="imbalance_measure_1",minimize = FALSE,properties = c("classif"),fun=prauc)

tune_control_glmnet <- makeTuneControlRandom(maxit = 100)

cvglmnet<-makeLearner("classif.cvglmnet",predict.type="prob")
cvglmnet_grid<-makeParamSet(makeNumericParam("alpha",lower=0,upper=1))
cvglmnet_tuned <- makeTuneWrapper(cvglmnet, resampling = inter,
                                  par.set = cvglmnet_grid, control = tune_control_glmnet,measures=list(self_measure))

meas = list(mmce, brier,auc,fdr,tpr,tnr,ppv,mcc,fpr,self_measure)
tuned_learners_full<-list(cvglmnet_tuned)

bmr_cvglmnet_full <- benchmark(tuned_learners_full,task,outer,meas)
save(bmr_cvglmnet_full,file="../../../data/arrhythmia/ml/bmr_cvglmnet.RData")