setwd("./")


rm(list = ls())
set.seed(123)
load("../../../data/arrhythmia/train_scaled.RData")
load("../../../data/arrhythmia/tunecontrol.RData")
#inter outer - resampling desc is loaded from fixed file

library(mlr)
source("../../src/prauc.R")
self_measure<-makeMeasure(id="imbalance_measure_1",minimize = FALSE,properties = c("classif"),fun=prauc)

options(java.parameters = "-Xmx10g")
library(rJava)
library(bartMachine)

tune_control_bart <- makeTuneControlGrid(resolution = 50L)
bartMachine <- makeLearner("classif.bartMachine",predict.type = "prob",mem_cache_for_speed=FALSE)
bartMachine_grid <- makeParamSet(makeDiscreteParam("num_trees",values=25*seq(2,100,by=2)))
bartMachine_tuned <- makeTuneWrapper(bartMachine, resampling = inter,
                                     par.set = bartMachine_grid, control = tune_control_bart,measures=list(self_measure))
meas = list(mmce, brier,auc,fdr,tpr,tnr,ppv,mcc,self_measure)
tuned_learners_full<-list(bartMachine_tuned)

bmr_bartMachine_full <- benchmark(tuned_learners_full,task,outer,meas)

save(bmr_bartMachine_full,file="./bmr_bartMachine.RData")

save(bmr_bartMachine_full,file="../../../data/arrhythmia/ml/bmr_bartMachine.RData")
