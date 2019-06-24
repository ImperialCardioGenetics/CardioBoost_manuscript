#calculate PRAUC given a mlr object
library(PRROC)
prauc=function(task,model,pred,feats,extra.args){
  prob_positive<-(1-getPredictionProbabilities(pred))[pred$data$truth==1]
  prob_negative<-(1-getPredictionProbabilities(pred))[pred$data$truth==0]
  pr <- pr.curve(prob_positive,prob_negative,  curve = F)
  pr$auc.integral
}