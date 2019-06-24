setwd("./")
rm(list = ls())

require(mlr)
set.seed(123)
load("../../../data/arrhythmia/train_scaled.RData")
load("../../../data/arrhythmia/tunecontrol.RData")
#inter outer - resampling desc is loaded from fixed file

source("../../src/prauc.R")
self_measure<-makeMeasure(id="imbalance_measure_1",minimize = FALSE,properties = c("classif"),fun=prauc)

tune_control_rf <- makeTuneControlRandom(maxit = 200)
rf<-makeLearner("classif.ranger",predict.type="prob")
rf_grid <- makeParamSet(makeIntegerParam("mtry", lower=5, upper=50),
                        makeLogicalParam("replace"),
                        makeDiscreteParam("num.trees", values=250:2500))
rf_tuned <- makeTuneWrapper(rf, resampling = inter,
                            par.set = rf_grid, control = tune_control_rf,measures=list(self_measure))

meas = list(mmce, brier,auc,fdr,tpr,tnr,ppv,mcc,self_measure)
tuned_learners_full<-list(rf_tuned)

bmr_rf_full <- benchmark(tuned_learners_full,task,outer,meas)
save(bmr_rf_full,file="../../../data/arrhythmia/ml/bmr_rf.RData")