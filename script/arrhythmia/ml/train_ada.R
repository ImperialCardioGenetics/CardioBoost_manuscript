rm(list=ls())
gc(reset=TRUE)
set.seed(42) #From random.org

setwd("./")
load("../../../data/arrhythmia/train_scaled.RData")
library(mlr)
load("../../../data/arrhythmia/ml/ada_tune.RData")

train_task <- makeClassifTask(id = "train", data = train_scaled[,c(5,11:78,80,81)],
                        target = "pathogenic")
train_task <- removeConstantFeatures(train_task)
ada<-makeLearner("classif.ada",predict.type = "prob",loss=res$x$loss,nu=res$x$nu)

train_ada<-train(ada,train_task)
save(train_ada,file="../../../data/arrhythmia/ml/train_ada.RData")
