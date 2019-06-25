setwd("./")
rm(list = ls())
require(mlr) 
set.seed(123)
load("../../../../data/cardiomyopathy_alternative/mybpc3_hcm/train_scaled.RData")
#inter outer - resampling desc is loaded from fixed file

source("../../src/prauc.R")
self_measure<-makeMeasure(id="imbalance_measure_1",minimize = FALSE,properties = c("classif"),fun=prauc)

## Define task
task <- makeClassifTask(id = "train", data = train_scaled[,c(4,9:72)],
                        target = "pathogenic")


tune_control_ada <- makeTuneControlRandom(maxit = 200)
ada<-makeLearner("classif.ada",predict.type = "prob")
ada_grid <- makeParamSet(makeDiscreteParam("loss",values=c("exponential","logistic")),
                         makeNumericParam("nu",lower=0.01,upper=1))

res<-tuneParams(ada,task=task,resampling = makeResampleDesc("CV",iters=5,stratify = TRUE),par.set =ada_grid,control = tune_control_ada,measures=list(self_measure))

save(res,file="../../../../data/cardiomyopathy_alternative/mybpc3_hcm/ml/ada_tune.RData")
