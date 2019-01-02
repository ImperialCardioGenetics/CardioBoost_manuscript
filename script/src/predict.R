##Given a test set of variants, preprocess its annotation, get the features of the variants and make pathogenic prediction on them
library(mlr)



imputed.cardioboost<-function(data){
  data<-subset(data,is.element(gene,disease_gene))
  data<-preprocess(data)
  data_rare<-extract_rare(data)
  data_imputed<-imputation(impute_mean,data_rare)
  data_imputed
}

scaled.cardioboost<-function(data_imputed){
  data_scaled<-normalization(data_imputed,scale_mean,scale_std)
  data_scaled<-select_feature(data_scaled)
  data_scaled  
}


predict.cardioboost<-function(train_model,data){
  data_imputed<-imputed.cardioboost(data)
  data_scaled<-scaled.cardioboost(data_imputed)
  test_task<-makeClassifTask(id = "test", data = data_scaled[,c(train_model$features,"pathogenic")],
                             target = "pathogenic")
  test_predict=predict(train_model,test_task)
  data_imputed$pathogenicity<-(test_predict$data$prob.1)
  data_imputed
  }
    