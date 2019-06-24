setwd("./")

load("../../../data/arrhythmia/train_scaled.RData")
load("../../../data/arrhythmia/test_data_leakage.RData")

library(dplyr)
library(mlr)

#to ensure the results are reproducible for each child process
set.seed(3) 

train_scaled<-sample_n(train_scaled,nrow(train_scaled))
train_filtered_scaled<-sample_n(train_filtered_scaled,nrow(train_filtered_scaled))
train_used_scaled<-sample_n(train_used_scaled,nrow(train_used_scaled))

## Define task
task_1 <- makeClassifTask(id = "train", data = train_scaled[,c(6,11:80,82,83)],
                          target = "pathogenic")
task_2 <- makeClassifTask(id = "train_filtered", data = train_filtered_scaled[,c(6,11:79,81)],
                          target = "pathogenic")
task_3 <- makeClassifTask(id="train_used_scaled",data=train_used_scaled[,c(6,11:79,81,82)],target = "pathogenic")

#preprocessing
task_1 <- removeConstantFeatures(task_1)
task_2 <- removeConstantFeatures(task_2)
task_3 <- removeConstantFeatures(task_3)

## Make tune wrapper
task<-list(task_1,task_2,task_3)

library(PRROC)
f=function(task,model,pred,feats,extra.args){
  prob_positive<-(1-getPredictionProbabilities(pred))[pred$data$truth==1]
  prob_negative<-(1-getPredictionProbabilities(pred))[pred$data$truth==0]
  pr <- pr.curve(prob_positive,prob_negative,  curve = F)
  pr$auc.integral
}


self_measure<-makeMeasure(id="prauc",minimize = FALSE,properties = c("classif"),fun=f)

tune_control <- makeTuneControlRandom(maxit = 20)

inter<-makeResampleDesc("CV",iters=5,stratify = TRUE)
outer<-makeResampleDesc("CV",iters=10,stratify = TRUE)

tune_control_ada <- makeTuneControlGrid(resolution = 100L)
ada<-makeLearner("classif.ada",predict.type="prob")
ada_grid <- makeParamSet(makeDiscreteParam("loss",values=c("exponential","logistic")),
                         makeNumericParam("nu",lower=0.2,upper=1))

ada_tuned <- makeTuneWrapper(ada, resampling = inter,
                             par.set = ada_grid, control = tune_control_ada,measures=list(auc))

meas<-list(mmce,tpr,ppv,tnr,auc,self_measure)
bmr_ada_test_filtered <- mlr::benchmark(ada_tuned,task,outer,meas)
save(bmr_ada_test_filtered,file="../../../data/arrhythmia/ml/bmr_ada_test_filtered.RData")