rm(list=ls())
gc(reset=TRUE)
set.seed(42) #From random.org

setwd("./")
load("../../../data/arrhythmia/train_scaled.RData")
library(mlr)
outer<-makeResampleDesc("CV",iters=10)

train_task <- makeClassifTask(id = "train", data = train_scaled[,c(6,11:80,82,83)],
                        target = "pathogenic")
train_task <- removeConstantFeatures(train_task)
ada<-makeLearner("classif.ada",predict.type = "prob",loss="exponential",nu=0.6461263)

train_ada<-train(ada,train_task)
save(train_ada,file="../../../data/arrhythmia/ml/train_ada.RData")
