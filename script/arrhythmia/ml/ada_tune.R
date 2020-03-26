setwd("./")
rm(list = ls())

require(mlr)
set.seed(123)
load("../../../data/arrhythmia/train_scaled.RData")
load("../../../data/arrhythmia/tunecontrol.RData")
#inter outer - resampling desc is loaded from fixed file

source("../../src/prauc.R")
self_measure<-makeMeasure(id="imbalance_measure_1",minimize = FALSE,properties = c("classif"),fun=prauc)

# #Adaboos: ada
tune_control_ada <- makeTuneControlRandom(maxit = 200)
ada<-makeLearner("classif.ada",predict.type = "prob")
ada_grid <- makeParamSet(makeDiscreteParam("loss",values=c("exponential")),
                         makeNumericParam("nu",lower=0.1,upper=0.9))

res<-tuneParams(ada,task=task,resampling = makeResampleDesc("CV",iters=5,stratify = TRUE),par.set =ada_grid,control = tune_control_ada,measures=list(self_measure))

save(res,file="../../../data/arrhythmia/ml/ada_tune.RData")
