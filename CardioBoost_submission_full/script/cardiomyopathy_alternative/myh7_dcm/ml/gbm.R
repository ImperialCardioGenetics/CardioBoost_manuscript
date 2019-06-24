setwd("./")
rm(list = ls())
require(mlr)
set.seed(123) 
load("../../../data/cardiomyopathy/train_scaled.RData")
load("../../../data/cardiomyopathy/tunecontrol.RData")
#inter outer - resampling desc is loaded from fixed file

source("../../src/prauc.R")
self_measure<-makeMeasure(id="imbalance_measure_1",minimize = FALSE,properties = c("classif"),fun=prauc)

tune_control_gbm <- makeTuneControlRandom(maxit = 200)

gbm<-makeLearner("classif.gbm",predict.type = "prob")
gbm_grid <- makeParamSet(makeIntegerParam("n.trees",200,2000),
                         makeNumericParam("shrinkage",lower=0.01,upper=0.07),
                         makeDiscreteParam("distribution","bernoulli"))
gbm_tuned <- makeTuneWrapper(gbm, resampling = inter,
                             par.set = gbm_grid, control = tune_control_gbm,measures=list(self_measure))

meas = list(mmce,brier,auc,fdr,tpr,tnr,ppv,mcc,self_measure)
tuned_learners_full<-list(gbm_tuned)

bmr_gbm_full<- benchmark(tuned_learners_full,task,outer,meas)
save(bmr_gbm_full,file="../../../data/cardiomyopathy/ml/bmr_gbm.RData")
