##Set up the number of inner and outer CV loops


load("../../data/processed/train_scaled.RData")

library(dplyr)
library(mlr)

#to ensure the results are reproducible for each child process
set.seed(12) 



## Define task
task <- makeClassifTask(id = "train", data = train_scaled[,c(6,11:80,82,83)],
                          target = "pathogenic")

task <- removeConstantFeatures(task)

tune_control <- makeTuneControlRandom(maxit = 20)

inter<-makeResampleDesc("CV",iters=5,stratify = TRUE)
outer<-makeResampleDesc("CV",iters=10,stratify = TRUE)

save(task,inter,outer,tune_control,file="../../data/processed/tunecontrol.RData")

