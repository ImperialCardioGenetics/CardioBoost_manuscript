setwd("./")
rm(list = ls())
require(mlr) 
set.seed(2)
load("../../../data/cardiomyopathy/train_scaled.RData")
load("../../../data/cardiomyopathy/tunecontrol.RData")
#inter outer - resampling desc is loaded from fixed file

source("../../src/prauc.R")
self_measure<-makeMeasure(id="imbalance_measure_1",minimize = FALSE,properties = c("classif"),fun=prauc)


tune_control_ada <- makeTuneControlRandom(maxit = 500)
ada<-makeLearner("classif.ada",predict.type = "prob")
ada_grid <- makeParamSet(makeDiscreteParam("loss",values=c("exponential")),
                         makeNumericParam("nu",lower=0.1,upper=0.9))

res<-tuneParams(ada,task=task,resampling = makeResampleDesc("CV",iters=5,stratify = TRUE),par.set =ada_grid,control = tune_control_ada,measures=list(self_measure))

save(res,file="../../../data/cardiomyopathy/ml/ada_tune.RData")