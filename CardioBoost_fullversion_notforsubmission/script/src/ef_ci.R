check_input<-function(case_freq,control_freq,predict=NA){
  #Pathogenic variants: case_freq could be 0 since rare variants on the given gene in a patient cohort could be not disease-causing.
  #Benign variants: control_freq could not be 0 since rare variants on the control population on the given gene should not be all disease-causing. 
  
  
    if(is.na(case_freq)==TRUE){case_freq=0}
    if(is.na(control_freq)==TRUE){control_freq=0}
    
    #if case_freq == 0: OR=0
    #if control_freq == 0: OR is undefined
  
    if(predict!="benign" || control_freq!=0){
    if((case_freq==0 || is.na(case_freq)==TRUE)||(control_freq==0 || is.na(control_freq)==TRUE)){
    case_freq=case_freq+0.5
    control_freq=control_freq+0.5
  }
    }
  if(case_freq==0 && control_freq ==0){
    return(c(0.0,0.0))
  }
  return(c(case_freq,control_freq))
}

or<-function(case_freq,control_freq,predict=NA){
  if(predict=="benign" && control_freq==0){
    or<-NA
  }else{
  or<-(case_freq/control_freq)/((1-case_freq)/(1-control_freq))}
  return(or)
}



############################
##Situations that the calculation of EF is problematic, should be noted and have a check on the data. 
##1) For predicted pathogenic variants,
##    a classifier doesn't output any pathogenic predictions in a patient cohort
##    here case_freq = 0

##2) For predicted benign variants,
##    a classifier doesn't output any benign predictions in a control population
##    here control_freq = 0


ef<-function(case_freq,control_freq,predict=NA){
  or<-or(case_freq,control_freq)
  ef_tilde<-1-(1/or)
  if(predict=="pathogenic"){
  ef_tilde<-ifelse(case_freq==0.5,NA,ef_tilde)}#control_freq ==0.5 means no variants in control population is predicted to be benign
  if(predict=="benign"){
  ef_tilde<-ifelse(control_freq==0.5,NA,ef_tilde)
  }  
  ef_tilde<-ifelse(ef_tilde<0,0,ef_tilde)
}


se_log_or<-function(case_freq,control_freq,patient_number,gnomAD_number,predict=NA){
  if(predict=="benign" && control_freq==0){se_or=NA}
  else{
  value<-1/((case_freq*patient_number))+(1/(control_freq*gnomAD_number))+(1/((1-case_freq)*patient_number))+(1/((1-control_freq)*gnomAD_number))
  se_or=sqrt(value)
  }
  se_or
}

or_ci_lower<-function(case_freq,control_freq,patient_number,gnomAD_number,predict=NA){
  or<-or(case_freq,control_freq,predict=predict)
  se_log_or<-se_log_or(case_freq,control_freq,patient_number,gnomAD_number,predict=predict)
  or_ci_lower<-or-1.96*se_log_or
  or_ci_lower<-ifelse(or_ci_lower<0,0,or_ci_lower)
}

or_ci_upper<-function(case_freq,control_freq,patient_number,gnomAD_number,predict=NA){
  or<-or(case_freq,control_freq,predict=predict)
  se_log_or<-se_log_or(case_freq,control_freq,patient_number,gnomAD_number,predict=predict)
  or_ci_upper<-or+1.96*se_log_or
  
}

#EF range [0,1]
ef_ci_lower<-function(case_freq,control_freq,patient_number,gnomAD_number,predict){
  ef<-ef(case_freq,control_freq,predict)
  log_or<-log(or(case_freq,control_freq))
  se_log_or<-se_log_or(case_freq,control_freq,patient_number,gnomAD_number)
  lower<-1-1/(exp(log_or-1.96*se_log_or))
  lower<-ifelse(lower<0,0,lower)
  lower<-ifelse(is.na(ef),NA,lower)
}


ef_ci_upper<-function(case_freq,control_freq,patient_number,gnomAD_number,predict){
  ef<-ef(case_freq,control_freq,predict)
  log_or<-log(or(case_freq,control_freq))
  se_log_or<-se_log_or(case_freq,control_freq,patient_number,gnomAD_number)
  upper<-1-1/(exp(log_or+1.96*se_log_or))
  upper<-ifelse(upper>1,1,upper)
  upper<-ifelse(is.na(ef),NA,upper)
}
