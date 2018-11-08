setwd("./")
rm(list = ls())
require(mlr)
set.seed(123)
load("../../data/processed/train_scaled.RData")
load("../../data/processed/tunecontrol.RData")
#inter outer - resampling desc is loaded from fixed file

source("./prauc.R")
self_measure<-makeMeasure(id="imbalance_measure_1",minimize = FALSE,properties = c("classif"),fun=prauc)


tune_control_knn <- makeTuneControlGrid(resolution=20L)
knn <- makeLearner("classif.kknn",predict.type = "prob")
knn_grid<-makeParamSet(makeIntegerParam("k",lower=1,upper=60))
knn_tuned <- makeTuneWrapper(knn, resampling = inter,
                                  par.set = knn_grid, control = tune_control_knn,measures=list(self_measure))

meas = list(mmce, brier,auc,fdr,tpr,tnr,ppv,mcc,self_measure)
tuned_learners_full<-list(knn_tuned)

bmr_knn_full <- benchmark(tuned_learners_full,task,outer,meas)
save(bmr_knn_full,file="../../data/ml/bmr_knn.RData")
