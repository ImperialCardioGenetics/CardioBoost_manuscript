setwd("./")
rm(list = ls())

require(mlr)
set.seed(123)
load("../../../data/arrhythmia/train_scaled.RData")
load("../../../data/arrhythmia/tunecontrol.RData")
#inter outer - resampling desc is loaded from fixed file

source("../../src/prauc.R")
self_measure<-makeMeasure(id="imbalance_measure_1",minimize = FALSE,properties = c("classif"),fun=prauc)

tune_control_xgboost <- makeTuneControlRandom(maxit = 200)
xgboost<-makeLearner("classif.xgboost",predict.type="prob",max_depth=6,nthread=10,nrounds=10,objective="binary:logistic")
xgboost_grid<-makeParamSet(makeNumericParam("eta",lower=0.2,upper=0.75))
xgboost_tuned <- makeTuneWrapper(xgboost, resampling = inter,
                                 par.set = xgboost_grid, control = tune_control_xgboost,measures=list(self_measure))
#meas = list(mmce, brier,auc,fdr,tpr,tnr,ppv,mcc)
meas<-list(brier,auc)
bmr_xgboost_full <- mlr::benchmark(xgboost_tuned,task,outer,meas)
save(bmr_xgboost_full,file="../../../data/arrhythmia/ml/bmr_xgboost.RData")