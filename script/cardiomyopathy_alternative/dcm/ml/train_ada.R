##Train AdaBoost on the whole training data

rm(list=ls())
gc(reset=TRUE)
set.seed(42) #From random.org

load("../../../../data/cardiomyopathy_alternative/dcm/train_scaled.RData")
load("../../../../data/cardiomyopathy_alternative/dcm/ml/ada_tune.RData")

library(mlr)

train_task <- makeClassifTask(id = "train", data = train_scaled[,c(6,11:80,82,83)],
                              target = "pathogenic")
train_task <- removeConstantFeatures(train_task)
ada<-makeLearner("classif.ada",predict.type = "prob",loss=res$x$loss,nu=res$x$nu)

train_ada<-train(ada,train_task)
save(train_ada,file="../../../../data/cardiomyopathy_alternative/dcm/ml/train_ada.RData")
