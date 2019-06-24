##Train AdaBoost on the whole training data

rm(list=ls())
gc(reset=TRUE)
set.seed(42) #From random.org


load("../../../data/mybpc3/train_scaled.RData")
load("../../../data/mybpc3/ml/ada_tune.RData")

library(mlr)
outer<-makeResampleDesc("CV",iters=10)

train_task <- makeClassifTask(id = "train", data = train_scaled[,c(4,9:72)],
                              target = "pathogenic")
train_task <- removeConstantFeatures(train_task)
ada<-makeLearner("classif.ada",predict.type = "prob",loss=res$x$loss,nu=res$x$nu)

train_ada<-train(ada,train_task)
save(train_ada,file="../../../data/mybpc3/ml/train_ada.RData")
