## preprocess script: preprocess sample to generate predictions

#preprocess
preprocess<-function(data_table){
  data_table<-data_table[!duplicated(data_table),]
  if(is.element("integrated_confidence_value",colnames(data_table))==TRUE){
    data_table<-subset(data_table,select=-c(integrated_confidence_value))}
  if(is.element("Eigen_coding_or_noncoding",colnames(data_table))==TRUE){
      data_table<-subset(data_table,select=-c(Eigen_coding_or_noncoding))}

  data_table[data_table=="."] <- NA
  
  data_table[,11:42]<-data.matrix(data_table[,11:42])
  data_table[,11:42]<-data.frame(data_table[,11:42])
  return(data_table)
  
}


#keeping the rare variants
extract_rare<-function(data_table){
  data_table<-subset(data_table,gnomAD_AF<=0.001 | is.na(gnomAD_AF))
  return(data_table)
}



#imputation
imputation<-function(impute_mean,data_table){
  data_imputed<-data_table
  names<-impute_features
  data_imputed[,14:39]<-mean_imputation(data_imputed[,14:39],impute_mean,names)
  
  ##For Column PARASCOREzSCORE, mis_badness and MPC (missing values means not applicable, replace with mean and add one more indicator variable to indicate the missingness
  data_imputed$PARASCORE_exist<-ifelse(is.na(data_imputed$PARASCOREzSCORE),0,1)
  data_imputed$PARASCOREzSCORE[which(is.na(data_imputed$PARASCOREzSCORE)==TRUE)]<-0
  
  ##missense badness [0,1]
  data_imputed$mis_badness_exist<-ifelse(is.na(data_imputed$mis_badness),0,1)
  data_imputed$mis_badness[is.na(data_imputed$mis_badness)]<--1
  ##MPC >0
  data_imputed$MPC[is.na(data_imputed$MPC)]<--1
  return(data_imputed)
}

z_normalization<-function(data,mean_feature,std_feature){
  for(i in 1:ncol(data)){
    data[,i]=(data[,i]-mean_feature[i])/std_feature[i]
  }
  return(data)
}

#feature normalization using the mean and std from the training data
normalization<-function(data_table,scale_mean,scale_std){
  data_table[,11:39]<-z_normalization(data_table[,11:39],scale_mean,scale_std)
  return(data_table)
}


#selecting wanted features
select_feature<-function(data_table){
  data_table<-data_table[,feature_column]
  return(data_table)
}

