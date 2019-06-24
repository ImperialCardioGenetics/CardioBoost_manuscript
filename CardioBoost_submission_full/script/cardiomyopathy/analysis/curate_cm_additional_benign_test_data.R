load("../../../data/cardiomyopathy/prediction/cm_prediction_full.RData")
load("../../../data/cardiomyopathy/cm_holdout_test.RData")
cm_gnomad_rare_noexac<-subset(cm_prediction,gnomAD_AF>=0)
cm_gnomad_rare_noexac<-subset(cm_gnomad_rare_noexac,is.na(AF_Adj)==TRUE)
cm_gnomad_rare_noexac<-subset(cm_gnomad_rare_noexac,isTrain==FALSE)
cm_gnomad_rare_noexac<-subset(cm_gnomad_rare_noexac,!is.element(key,cm_patho_test$key))
cm_gnomad_rare_noexac<-subset(cm_gnomad_rare_noexac,!is.element(key,test$key))

save(cm_gnomad_rare_noexac,file="../../../data/cardiomyopathy/cm_additional_test_benign.RData")
