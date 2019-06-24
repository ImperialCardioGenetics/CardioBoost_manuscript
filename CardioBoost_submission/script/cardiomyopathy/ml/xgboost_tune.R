setwd("./")
rm(list = ls())
require(mlr) 
set.seed(123)
load("../../../data/cardiomyopathy/train_scaled.RData")
load("../../../data/cardiomyopathy/tunecontrol.RData")
#inter outer - resampling desc is loaded from fixed file

source("../../src/prauc.R")
self_measure<-makeMeasure(id="imbalance_measure_1",minimize = FALSE,properties = c("classif"),fun=prauc)

tune_control_xgboost <- makeTuneControlRandom(maxit=200)
xgboost<-makeLearner("classif.xgboost",predict.type="prob",max_depth=6,nthread=10,nrounds=10,objective="binary:logistic")
xgboost_grid<-makeParamSet(makeNumericParam("eta",lower=0.1,upper=0.8))
res<-tuneParams(xgboost,task=task,resampling = makeResampleDesc("CV",iters=5,stratify = TRUE),par.set =xgboost_grid,control = tune_control_xgboost,measures=list(self_measure))

save(res,file="../../../data/cardiomyopathy/ml/xgboost_tune.RData")