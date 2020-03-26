##Split the samples into inner CV and outer CV loops so that they are fixed in different model training. 
setwd("./")
load("../../../data/arrhythmia/train_scaled.RData")

library(dplyr)
library(mlr)

#to ensure the results are reproducible for each child process
set.seed(7) 

## Define task
task <- makeClassifTask(id = "train", data = train_scaled[,c(5,11:78,80,81)],
                        target = "pathogenic")

task <- removeConstantFeatures(task)

tune_control <- makeTuneControlRandom(maxit = 20)

inter<-makeResampleDesc("CV",iters=5,stratify = TRUE)
outer<-makeResampleDesc("CV",iters=10,stratify = TRUE)

save(task,inter,outer,tune_control,file="../../../data/arrhythmia/tunecontrol.RData")

