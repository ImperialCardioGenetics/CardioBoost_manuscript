library(mlr)
mcc<-function(task,model,pred,feats,extra.args){
  measureMCC(pred$data$truth,ifelse((1-getPredictionProbabilities(pred)>=0.9),1,0),0,1)
  
}
high_confidence_mcc=function(task,model,pred,feats,extra.args){
  TP = sum((pred$data$truth==1) & (1-getPredictionProbabilities(pred)>=0.90))
  TN = sum((pred$data$truth==0) & (getPredictionProbabilities(pred)>=0.90))
  
  FP = sum((pred$data$truth==0) & (1-getPredictionProbabilities(pred)>=0.90))
  FN = sum((pred$data$truth==1) & (getPredictionProbabilities(pred)>=0.90))

  mcc = (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
  mcc
  
}

f1_measure<-function(task,model,pred,feats,extra.args){
  f1=measureF1(pred$data$truth, ifelse((1-getPredictionProbabilities(pred)>=0.9),1,0), 1)
  f1
}