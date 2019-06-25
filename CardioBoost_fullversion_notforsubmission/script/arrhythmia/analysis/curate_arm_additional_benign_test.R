load("../../../data/arrhythmia/prediction/arm_prediction.RData")
load("../../../data/arrhythmia/arm_holdout_test.RData")
load("../../../data/arrhythmia/arm_additional_patho_test.RData")
arm_gnomad_rare_noexac<-subset(arm_prediction,gnomAD_AF>=0)
arm_gnomad_rare_noexac<-subset(arm_gnomad_rare_noexac,is.na(AF_Adj)==TRUE)
arm_gnomad_rare_noexac<-subset(arm_gnomad_rare_noexac,isTrain==FALSE)

arm_gnomad_rare_noexac<-subset(arm_gnomad_rare_noexac,!is.element(key,arm_patho_test$key))
arm_gnomad_rare_noexac<-subset(arm_gnomad_rare_noexac,!is.element(key,test$key))
arm_gnomad_rare_noexac$pathogenic=0
save(arm_gnomad_rare_noexac,file="../../../data/arrhythmia/arm_additional_test_benign.RData")
