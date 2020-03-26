##Train AdaBoost on the whole training data

rm(list=ls())
gc(reset=TRUE)


load("../../../data/cardiomyopathy/train_scaled.RData")

library(mlr)
load("../../../data/cardiomyopathy/ml/ada_tune.RData")

train_task <- makeClassifTask(id = "train", data = train_scaled[,c(6,11:80,82,83)],
                              target = "pathogenic")
train_task <- removeConstantFeatures(train_task)
ada<-makeLearner("classif.ada",predict.type = "prob",loss="exponential",nu=res$x$nu)

train_ada<-train(ada,train_task)
save(train_ada,file="../../../data/cardiomyopathy/ml/train_ada.RData")
