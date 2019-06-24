setwd("./")
rm(list = ls())
require(mlr)
set.seed(123)
load("../../../data/cardiomyopathy/train_scaled.RData")
load("../../../data/cardiomyopathy/tunecontrol.RData")
#inter outer - resampling desc is loaded from fixed file

source("../../src/prauc.R")
self_measure<-makeMeasure(id="imbalance_measure_1",minimize = FALSE,properties = c("classif"),fun=prauc)

tune_control_knn <- makeTuneControlRandom(maxit=100)
knn <- makeLearner("classif.kknn",predict.type = "prob")
knn_grid<-makeParamSet(makeIntegerParam("k",lower=1,upper=100))
knn_tuned <- makeTuneWrapper(knn, resampling = inter,
                             par.set = knn_grid, control = tune_control_knn,measures=list(self_measure))

meas = list(mmce, brier,auc,fdr,tpr,tnr,ppv,mcc,self_measure)
tuned_learners_full<-list(knn_tuned)

bmr_knn_full <- benchmark(tuned_learners_full,task,outer,meas)

save(bmr_knn_full,file="../../../data/cardiomyopathy/ml/bmr_knn.RData")
