setwd("./")
rm(list = ls())
require(mlr) 
set.seed(123)
load("../../../../data/cardiomyopathy_alternative/myh7_dcm/train_scaled.RData")
load("../../../../data/cardiomyopathy_alternative/myh7_dcm/tunecontrol.RData")
#inter outer - resampling desc is loaded from fixed file

source("../../../src/prauc.R")
self_measure<-makeMeasure(id="imbalance_measure_1",minimize = FALSE,properties = c("classif"),fun=prauc)

tune_control_ada <- makeTuneControlRandom(maxit = 200)
ada<-makeLearner("classif.ada",predict.type = "prob")
ada_grid <- makeParamSet(makeDiscreteParam("loss",values=c("exponential","logistic")),
                         makeNumericParam("nu",lower=0.01,upper=1))

ada_tuned <- makeTuneWrapper(ada, resampling = inter,
                             par.set = ada_grid, control = tune_control_ada,measures=list(self_measure))

meas = list(mmce, brier,auc,fdr,tpr,tnr,ppv,mcc,self_measure)
tuned_learners_full<-list(ada_tuned)

bmr_ada_full<- benchmark(tuned_learners_full,task,outer,meas)
save(bmr_ada_full,file="../../../../data/cardiomyopathy_alternative/myh7_dcm/ml/bmr_ada.RData")
