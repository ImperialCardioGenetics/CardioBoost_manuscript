setwd("./")
rm(list = ls())
options(java.parameters = "-Xmx20g")
require(mlr) 
require(rJava)
require(bartMachine)
set.seed(123)
load("../../../data/cardiomyopathy/train_scaled.RData")
load("../../../data/cardiomyopathy/tunecontrol.RData")
#inter outer - resampling desc is loaded from fixed file

source("../../src/prauc.R")
self_measure<-makeMeasure(id="imbalance_measure_1",minimize = FALSE,properties = c("classif"),fun=prauc)

set_bart_machine_num_cores(4)

#options(java.parameters = "-Xmx10g")


tune_control_bart <- makeTuneControlRandom(maxit = 50)
bartMachine <- makeLearner("classif.bartMachine",predict.type = "prob",mem_cache_for_speed=FALSE)
bartMachine_grid <- makeParamSet(makeDiscreteParam("num_trees",values=25*seq(2,100,by=2)))
bartMachine_tuned <- makeTuneWrapper(bartMachine, resampling = inter,
                                     par.set = bartMachine_grid, control = tune_control_bart,measures=list(self_measure))
meas = list(mmce, brier,auc,fdr,tpr,tnr,ppv,mcc,self_measure)
tuned_learners_full<-list(bartMachine_tuned)

bmr_bartMachine_full <- benchmark(tuned_learners_full,task,outer,meas)


save(bmr_bartMachine_full,file="../../../data/cardiomyopathy/ml/bmr_bartMachine.RData")
