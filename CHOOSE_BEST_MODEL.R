library(reshape)
inputDir<-crop_dir
mx_out <- paste(inputDir,"/maxent_modeling",sep="")
summ_dir<-paste0(mx_out,"/","summary-files"); if (!file.exists(summ_dir)) {dir.create(summ_dir)}

asd <- read.csv(paste(crop_dir,"/maxent_modeling/summary-files/ASD15.csv",sep=""))

########################################################################################
calc_cauc <- function(x){
  AUC  <- newpriFile$auc[newpriFile$Taxon==x]
  nAUC <- nullModel$AUC[nullModel==x]
  cAUC <- AUC+.5-max(c(.5,nAUC,na.rm=T))
  return(cAUC)
}
#calc_cauc <- Vectorize(FUN=calc_cauc, vectorize.args=x)

# newpriFile<-read.csv(paste0(summ_dir,"/","modelsMets.csv"),header=T)
# cauc <- calc_cauc(as.character(newpriFile$TAXON))


########################################################################################

ENM_PA_MEDIAN <-read.csv(paste0(summ_dir,"/","ensemble_metrics_PA_MEDIAN.csv"),header = T)
ENM_PA_MEDIAN<-ENM_PA_MEDIAN[,1:11]
ENM_PA_SUM    <-read.csv(paste0(summ_dir,"/","ensemble_metrics_PA_SUM.csv"),header = T)
ENM_PA_SUM<-ENM_PA_SUM[,1:11]

ENM_PA_MEAN   <-read.csv(paste0(summ_dir,"/","ensemble_metrics_PA_MEAN.csv"),header = T)
ENM_PA_MEAN<-ENM_PA_MEAN[,1:11]

########################################################################################
ENM_PA_MEDIAN$asd15<-NA
ENM_PA_SUM$asd15<-NA
ENM_PA_MEAN$asd15<-NA


########################################################################################
ENM_PA_MEDIAN$cAUC<-NA
ENM_PA_SUM$cAUC<-NA
ENM_PA_MEAN$cAUC<-NA

########################################################################################
ENM_PA_MEDIAN$nAUC<-NA
ENM_PA_SUM$nAUC<-NA
ENM_PA_MEAN$nAUC<-NA

########################################################################################
ENM_PA_MEDIAN$ValidModel<-NA
ENM_PA_SUM$ValidModel<-NA
ENM_PA_MEAN$ValidModel<-NA
########################################################################################
 nullModel <- read.csv(paste(crop_dir, "/maxent_modeling/summary-files/nullModel_auc.csv", sep=""))
 nullModel <- nullModel[order(nullModel$SPID),]
 rownames(nullModel) <- 1:nrow(nullModel)
 
 newpriFile<-ENM_PA_MEAN
  ENM_PA_MEAN$cAUC <- calc_cauc(as.character(newpriFile$Taxon))
 ENM_PA_MEAN$nAUC <- merge(ENM_PA_MEAN,nullModel,by.x="Taxon",by.y="SPID")$AUC
 
 
 newpriFile<-ENM_PA_MEDIAN
 ENM_PA_MEDIAN$cAUC <- calc_cauc(as.character(newpriFile$Taxon))
 ENM_PA_MEDIAN$nAUC <- merge(ENM_PA_MEDIAN,nullModel,by.x="Taxon",by.y="SPID")$AUC
 
 newpriFile<-ENM_PA_SUM
 ENM_PA_SUM$cAUC <- calc_cauc(as.character(newpriFile$Taxon))
 ENM_PA_SUM$nAUC <- merge(ENM_PA_SUM,nullModel,by.x="Taxon",by.y="SPID")$AUC
# newpriFile$cAUC <- cauc
# newpriFile$nAUC<-nullModel$AUC
########################################################################################
for (i in 1:nrow(ENM_PA_MEAN)) {
spp<-as.character(ENM_PA_MEAN$Taxon[[i]])
    cat("Processing taxon",paste(spp),"\n")
  
  #getting the quality metrics
#  atauc <- acc$TestAUC[which(acc$SPID==spp)]
 # stauc <- acc$TestAUCSD[which(acc$SPID==spp)]

  ENM_PA_MEDIAN[i,12]<-asd$rateThresholded[which(asd$taxon==paste(spp))]
  ENM_PA_SUM[i,12]<-asd$rateThresholded[which(asd$taxon==paste(spp))]
  ENM_PA_MEAN[i,12]<-asd$rateThresholded[which(asd$taxon==paste(spp))]
  
  };rm(i)
########################################################################################

for (i in 1:nrow(ENM_PA_MEAN)) {
  
  if (is.na(ENM_PA_MEDIAN[i,2])) {ENM_PA_MEDIAN[i,2] <- 0}
  if (is.na(ENM_PA_SUM[i,2])) {ENM_PA_SUM[i,2] <- 0}
  if (is.na(ENM_PA_MEAN[i,2])) {ENM_PA_MEAN[i,2] <- 0}
  
  if (is.na(ENM_PA_MEDIAN[i,3])) {ENM_PA_MEDIAN[i,3] <- 1}
  if (is.na(ENM_PA_SUM[i,3])) {ENM_PA_SUM[i,3] <- 1}
  if (is.na(ENM_PA_MEAN[i,3])) {ENM_PA_MEAN[i,3] <- 1}
 
  if (is.na(ENM_PA_MEDIAN[i,10])) {ENM_PA_MEDIAN[i,10] <- 100}
  if (is.na(ENM_PA_SUM[i,10])) {ENM_PA_SUM[i,10] <- 100}
  if (is.na(ENM_PA_MEAN[i,10])) {ENM_PA_MEAN[i,10] <- 100}

};rm(i)


########################################################################################
for (i in 1:nrow(ENM_PA_MEAN)) {

#if (ENM_PA_MEAN[i,2]>=0.7 & ENM_PA_MEAN[i,3]<=0.15 & ENM_PA_MEAN[i,10]<=10 & ENM_PA_MEAN[i,9]>=0.7) {#using ASD15
#  if (ENM_PA_MEAN[i,2]>=0.7 & ENM_PA_MEAN[i,3]<=0.15 & ENM_PA_MEAN[i,5]>=0.6 & ENM_PA_MEAN[i,7]>=0.7 & ENM_PA_MEAN[i,9]>=0.6) {#using PCC
 # if (ENM_PA_MEAN[i,2]>=0.7 & ENM_PA_MEAN[i,3]<=0.15  & ENM_PA_MEAN[i,7]>=0.7 & ENM_PA_MEAN[i,9]>=0.4) { # using sensitivity
 # if (ENM_PA_MEAN[i,2]>=0.7 & ENM_PA_MEAN[i,3]<=0.15   & ENM_PA_MEAN[i,9]>=0.4) { # using TSS
    if (ENM_PA_MEAN[i,"auc"]>=0.7 &
        ENM_PA_MEAN[i,"sd_auc"]<=0.15  & 
      #  ENM_PA_MEAN[i,"TSS"]>=0.4 &
        ENM_PA_MEAN[i,"cAUC"]>=0.4 &
        ENM_PA_MEAN[i,"asd15"]<=10 
                ) { # using TSS and sensitivity
      
  ENM_PA_MEAN[i,"ValidModel"] <- 1
} else {
  ENM_PA_MEAN[i,"ValidModel"]  <- 0
}
###
#  if (ENM_PA_SUM[i,2]>=0.7 & ENM_PA_SUM[i,3]<=0.15 & ENM_PA_SUM[i,10]<=10 & ENM_PA_MEAN[i,9]>=0.7) { #using ASD15
#  if (ENM_PA_SUM[i,2]>=0.7 & ENM_PA_SUM[i,3]<=0.15  &  ENM_PA_SUM[i,5]>=0.6 & ENM_PA_SUM[i,7]>=0.7 & ENM_PA_SUM[i,9]>=0.6) { #using PCC
 # if (ENM_PA_SUM[i,2]>=0.7 & ENM_PA_SUM[i,3]<=0.15 & ENM_PA_SUM[i,7]>=0.7 & ENM_PA_SUM[i,9]>=0.4) { # using sensitivity
 # if (ENM_PA_SUM[i,2]>=0.7 & ENM_PA_SUM[i,3]<=0.15 & ENM_PA_SUM[i,9]>=0.4) { # using TSS
  if (
    ENM_PA_SUM[i,"auc"]>=0.7 &
    ENM_PA_SUM[i,"sd_auc"]<=0.15  & 
    #ENM_PA_SUM[i,"TSS"]>=0.4 &
    ENM_PA_SUM[i,"cAUC"]>=0.4 &
    ENM_PA_SUM[i,"asd15"]<=10
    
  ){
    ENM_PA_SUM[i,"ValidModel"] <- 1
  } else {
    ENM_PA_SUM[i,"ValidModel"]  <- 0
  }  
  ###
#  if (ENM_PA_MEDIAN[i,2]>=0.7 & ENM_PA_MEDIAN[i,3]<=0.15 & ENM_PA_MEDIAN[i,10]<=10 & ENM_PA_MEAN[i,9]>=0.7) {#using ASD15
 # if (ENM_PA_MEDIAN[i,2]>=0.7 & ENM_PA_MEDIAN[i,3]<=0.15 & ENM_PA_MEDIAN[i,5]>=0.6 & ENM_PA_MEDIAN[i,7]>=0.7 & ENM_PA_MEDIAN[i,9]>=0.6) {  #using PCC
 # if (ENM_PA_MEDIAN[i,2]>=0.7 & ENM_PA_MEDIAN[i,3]<=0.15 & ENM_PA_MEDIAN[i,7]>=0.7 & ENM_PA_MEDIAN[i,9]>=0.4) { # using sensitivity
 # if (ENM_PA_MEDIAN[i,2]>=0.7 & ENM_PA_MEDIAN[i,3]<=0.15  & ENM_PA_MEDIAN[i,9]>=0.4) { # using TSS
  if (
    ENM_PA_MEDIAN[i,"auc"]>=0.7 &
    ENM_PA_MEDIAN[i,"sd_auc"]<=0.15  & 
   # ENM_PA_MEDIAN[i,"TSS"]>=0.4 &
    ENM_PA_MEDIAN[i,"cAUC"]>=0.4 &
    ENM_PA_MEDIAN[i,"asd15"]<=10
    
  ) { # using TSS and sensitivity
    
    ENM_PA_MEDIAN[i,"ValidModel"] <- 1
  } else {
    ENM_PA_MEDIAN[i,"ValidModel"]  <- 0
    }  
};rm(i)
########################################################################################
 
 mxe_out <- paste(inputDir,"/maxent_modeling/models",sep="")
 
 spList2<-list.dirs(mxe_out,full.names = FALSE, recursive = FALSE)
 
 points_info<-lapply(1:length(spList2), function(i){
   spID<-spList2[[i]]
   outFolder <- paste(mxe_out,"/",spID,sep="")
   crossval_dir<- paste(outFolder, "/", "crossval",sep="")  
   x<-read.csv(paste0(crossval_dir,"/","maxentResults.csv"),header=T)
   x<-cbind(spID,as.numeric(x[nrow(x),2]+x[nrow(x),7]))
 return(x)
   })
 points_info<-as.data.frame(do.call(rbind,points_info))
 points_info$V2<-as.numeric(as.character(points_info$V2))
 colnames(points_info)<-c("Taxon","Occurrences")
 
 points_info$POINTS_CRITERIA<-NA
 points_info$POINTS_CRITERIA[which(points_info$Occurrences<20)]<-1
 points_info$POINTS_CRITERIA[which(points_info$Occurrences>=20)]<-0

 #######################################################################################
VALID_MODELS<-as.data.frame(matrix(nrow = nrow(ENM_PA_MEAN),ncol=2))

VALID_MODELS[,1]<-ENM_PA_MEAN$Taxon
VALID_MODELS[,2]<-ENM_PA_MEAN$ValidModel
# VALID_MODELS[,3]<-ENM_PA_MEDIAN$ValidModel
# VALID_MODELS[,4]<-ENM_PA_SUM$ValidModel
# VALID_MODELS[,5]<-ENM_PA_MEAN$TSS
# VALID_MODELS[,6]<-ENM_PA_MEDIAN$TSS
# VALID_MODELS[,7]<-ENM_PA_SUM$TSS
colnames(VALID_MODELS)<-c("Taxon","BEST_MODEL")

VALID_MODELS<- merge(VALID_MODELS,points_info,by="Taxon")
VALID_MODELS<- merge(VALID_MODELS,ENM_PA_MEAN,by="Taxon")
#VALID_MODELS<-VALID_MODELS[,-c(5:12,14:16)]
#colnames(VALID_MODELS)<-c("Taxon","ENM_PA_MEAN","ENM_PA_MEDIAN","ENM_PA_SUM","ENM_PA_MEAN_TSS","ENM_PA_MEDIAN_TSS","ENM_PA_SUM_TSS","BEST_MODEL")



for (i in 1:nrow(VALID_MODELS)) {
  ###############
  PA_MEAN<-NA;if(ENM_PA_MEAN[i,"ValidModel"]==0){PA_MEAN<-NA}else{PA_MEAN<-ENM_PA_MEAN[i,"ValidModel"]}
  PA_MEDIAN<-NA;if(ENM_PA_MEDIAN[i,"ValidModel"]==0){PA_MEDIAN<-NA}else{PA_MEDIAN<-ENM_PA_MEDIAN[i,"ValidModel"]}
  PA_SUM<-NA;if(ENM_PA_SUM[i,"ValidModel"]==0){PA_SUM<-NA}else{PA_SUM<-ENM_PA_SUM[i,"ValidModel"]}
  ###############
  PA_MEAN_TSS<-ENM_PA_MEAN[i,"TSS"];if(is.na(PA_MEAN)){PA_MEAN_TSS<--100}
  PA_MEDIAN_TSS<-ENM_PA_MEDIAN[i,"TSS"];if(is.na(PA_MEDIAN)){PA_MEDIAN_TSS<--100}
  PA_SUM_TSS<-ENM_PA_SUM[i,"TSS"];if(is.na(PA_SUM)){PA_SUM_TSS<--100}
  ###############
  ###############
  PA_TSS<-as.data.frame(cbind(PA_MEAN_TSS,PA_MEDIAN_TSS,PA_SUM_TSS));colnames(PA_TSS)<-c("ENM_PA_MEAN","ENM_PA_MEDIAN","ENM_PA_SUM")
  
  if(length(PA_TSS[which(PA_TSS==max(PA_TSS,na.rm = T))])>1){
  PA_TSS$BEST_MODELS1<-colnames(PA_TSS)[max.col(matrix(PA_TSS[,1:3],nrow=1),"first")] 
  PA_TSS$BEST_MODELS2<-colnames(PA_TSS)[max.col(matrix(PA_TSS[,1:3],nrow=1),"last")]  
  }else{
  PA_TSS$BEST_MODELS1<-colnames(PA_TSS)[max.col(matrix(PA_TSS[,1:3],nrow=1),"first")]   
  PA_TSS$BEST_MODELS2<-NA  
  }
  
  if(PA_TSS[,1]==-100 & PA_TSS[,2]==-100 & PA_TSS[,3]==-100){
    PA_TSS[,1]<-NA
    PA_TSS[,2]<-NA
    PA_TSS[,3]<-NA
    PA_TSS[,4]<-NA
    PA_TSS[,5]<-NA
  }
  
  ###############
  ###############
  PA_MEAN_AUC<-ENM_PA_MEAN[i,"auc"];if(is.na(PA_MEAN)){PA_MEAN_AUC<--100}
  PA_MEDIAN_AUC<-ENM_PA_MEDIAN[i,"auc"];if(is.na(PA_MEDIAN)){PA_MEDIAN_AUC<--100}
  PA_SUM_AUC<-ENM_PA_SUM[i,"auc"];if(is.na(PA_SUM)){PA_SUM_AUC<--100}
  ###############
  ###############
  PA_AUC<-as.data.frame(cbind(PA_MEAN_AUC,PA_MEDIAN_AUC,PA_SUM_AUC));colnames(PA_AUC)<-c("ENM_PA_MEAN","ENM_PA_MEDIAN","ENM_PA_SUM")
  if(length(PA_AUC[which(PA_AUC==max(PA_AUC,na.rm = T))])>1){
    PA_AUC$BEST_MODELS1<-colnames(PA_AUC)[max.col(matrix(PA_AUC[,1:3],nrow=1),"first")] 
    PA_AUC$BEST_MODELS2<-colnames(PA_AUC)[max.col(matrix(PA_AUC[,1:3],nrow=1),"last")]  
  }else{
    PA_AUC$BEST_MODELS1<-colnames(PA_AUC)[max.col(matrix(PA_AUC[,1:3],nrow=1),"first")]   
    PA_AUC$BEST_MODELS2<-NA  
  }

  if(PA_AUC[,1]==-100 & PA_AUC[,2]==-100 & PA_AUC[,3]==-100){
    PA_AUC[,1]<-NA
    PA_AUC[,2]<-NA
    PA_AUC[,3]<-NA
    PA_AUC[,4]<-NA
    PA_AUC[,5]<-NA
    
  }
  
  ###############
  ###############
  PA_MEAN_SENS<-ENM_PA_MEAN[i,"Sensitivity"];if(is.na(PA_MEAN)){PA_MEAN_SENS<--100}
  PA_MEDIAN_SENS<-ENM_PA_MEDIAN[i,"Sensitivity"];if(is.na(PA_MEDIAN)){PA_MEDIAN_SENS<--100}
  PA_SUM_SENS<-ENM_PA_SUM[i,"Sensitivity"];if(is.na(PA_SUM)){PA_SUM_SENS<--100}
  ###############
  ###############
  PA_SENS<-as.data.frame(cbind(PA_MEAN_SENS,PA_MEDIAN_SENS,PA_SUM_SENS));colnames(PA_SENS)<-c("ENM_PA_MEAN","ENM_PA_MEDIAN","ENM_PA_SUM")
  if(length(PA_SENS[which(PA_SENS==max(PA_SENS,na.rm = T))])>1){
    PA_SENS$BEST_MODELS1<-colnames(PA_SENS)[max.col(matrix(PA_SENS[,1:3],nrow=1),"first")] 
    PA_SENS$BEST_MODELS2<-colnames(PA_SENS)[max.col(matrix(PA_SENS[,1:3],nrow=1),"last")]  
  }else{
    PA_SENS$BEST_MODELS1<-colnames(PA_SENS)[max.col(matrix(PA_SENS[,1:3],nrow=1),"first")]   
    PA_SENS$BEST_MODELS2<-NA  
  }
  
  if(PA_SENS[,1]==-100 & PA_SENS[,2]==-100 & PA_SENS[,3]==-100){
    PA_SENS[,1]<-NA
    PA_SENS[,2]<-NA
    PA_SENS[,3]<-NA
    PA_SENS[,4]<-NA
    PA_SENS[,5]<-NA
    
  }
  ###############
  ###############
  PA_MEAN_SPEC<-ENM_PA_MEAN[i,"Specificity"];if(is.na(PA_MEAN)){PA_MEAN_SPEC<--100}
  PA_MEDIAN_SPEC<-ENM_PA_MEDIAN[i,"Specificity"];if(is.na(PA_MEDIAN)){PA_MEDIAN_SPEC<--100}
  PA_SUM_SPEC<-ENM_PA_SUM[i,"Specificity"];if(is.na(PA_SUM)){PA_SUM_SPEC<--100}
  ###############
  ###############
  PA_SPEC<-as.data.frame(cbind(PA_MEAN_SPEC,PA_MEDIAN_SPEC,PA_SUM_SPEC));colnames(PA_SPEC)<-c("ENM_PA_MEAN","ENM_PA_MEDIAN","ENM_PA_SUM")
  if(length(PA_SPEC[which(PA_SPEC==max(PA_SPEC,na.rm = T))])>1){
    PA_SPEC$BEST_MODELS1<-colnames(PA_SPEC)[max.col(matrix(PA_SPEC[,1:3],nrow=1),"first")] 
    PA_SPEC$BEST_MODELS2<-colnames(PA_SPEC)[max.col(matrix(PA_SPEC[,1:3],nrow=1),"last")]  
  }else{
    PA_SPEC$BEST_MODELS1<-colnames(PA_SPEC)[max.col(matrix(PA_SPEC[,1:3],nrow=1),"first")]   
    PA_SPEC$BEST_MODELS2<-NA  
  }
  
  if(PA_SPEC[,1] ==-100 & PA_SPEC[,2] ==-100 & PA_SPEC[,3]==-100){
    PA_SPEC[,1]<-NA
    PA_SPEC[,2]<-NA
    PA_SPEC[,3]<-NA
    PA_SPEC[,4]<-NA
    PA_SPEC[,5]<-NA
    
  }
  
  ###############
  ###############
  FINAL_BM<-as.data.frame(cbind(PA_TSS$BEST_MODELS1,
                                PA_TSS$BEST_MODELS2,
                                PA_AUC$BEST_MODELS1,
                                PA_AUC$BEST_MODELS2,
                                PA_SENS$BEST_MODELS1,
                                PA_SENS$BEST_MODELS2,
                                PA_SPEC$BEST_MODELS1,
                                PA_SPEC$BEST_MODELS2
                                ))
  colnames(FINAL_BM)<-c("TSS_1","TSS_2","AUC_1","AUC_2","SENS_1","SENS_2","SPEC_1","SPEC_2")
  FINAL_BM<-as.data.frame(t(FINAL_BM))
  colnames(FINAL_BM)<-"MODEL"
  s_BM<-as.data.frame(t(tapply(FINAL_BM$MODEL,FINAL_BM$MODEL,length)))
  colnames(s_BM)[which.max(s_BM)]
  if(ncol(s_BM)>0){
  VALID_MODELS[i,2]<- colnames(s_BM)[which.max(s_BM)]
  }else{VALID_MODELS[i,2]<-NA
  }
   };rm(i,
        PA_MEAN,
        PA_MEDIAN,
        PA_SUM,
        PA_MEAN_TSS,
        PA_MEDIAN_TSS,
        PA_SUM_TSS,
        PA_TSS,
        PA_MEAN_AUC,
        PA_MEDIAN_AUC,
        PA_SUM_AUC,
        PA_AUC,
        PA_MEAN_SENS,
        PA_MEDIAN_SENS,
        PA_SUM_SENS,
        PA_MEAN_SPEC,
        PA_MEDIAN_SPEC,
        PA_SUM_SPEC,
        PA_SPEC,
        FINAL_BM,
        asd
   ) 



for(i in 1:nrow(VALID_MODELS)){
  cat("processing ",as.character(VALID_MODELS$Taxon[[i]])," | ",i,"\n")
  if(VALID_MODELS$POINTS_CRITERIA[[i]]==1 & VALID_MODELS$asd15[[i]]>0.10 | is.na(VALID_MODELS$asd15[[i]])){
    VALID_MODELS$FILENAME[[i]]<- paste0(VALID_MODELS$Taxon[[i]],"_worldclim2_5_EMN_PA_SUM_ensemble.asc.gz")
    VALID_MODELS$BEST_MODEL[[i]]<-"ENM_PA_SUM"
  }else{
    if(is.na(VALID_MODELS$BEST_MODEL[[i]])){
      VALID_MODELS$BEST_MODEL[[i]]<-"ENM_PA_SUM"
      VALID_MODELS$FILENAME[[i]]<- paste0(VALID_MODELS$Taxon[[i]],"_worldclim2_5_EMN_PA_SUM_ensemble.asc.gz")
      
    }else if(VALID_MODELS$BEST_MODEL[[i]]=="ENM_PA_MEAN"){
      VALID_MODELS$FILENAME[[i]]<- paste0(VALID_MODELS$Taxon[[i]],"_worldclim2_5_EMN_PA_mean.asc.gz")
    }else if(VALID_MODELS$BEST_MODEL[[i]]=="ENM_PA_MEDIAN"){
      VALID_MODELS$FILENAME[[i]]<- paste0(VALID_MODELS$Taxon[[i]],"_worldclim2_5_EMN_PR_THR_ensemble.asc.gz")
    }else if(VALID_MODELS$BEST_MODEL[[i]]=="ENM_PA_SUM"){
      VALID_MODELS$FILENAME[[i]]<- paste0(VALID_MODELS$Taxon[[i]],"_worldclim2_5_EMN_PA_SUM_ensemble.asc.gz")
    }
  }
}#;rm(i)  

########################################################################################

ModelsMets<-as.data.frame(matrix(nrow = nrow(VALID_MODELS),ncol=18))
colnames(ModelsMets)<-c(
  "TAXON",
  "APPROACH",
  "UNIQUE_POINTS",
  "POINT_CRITERIA",
  "THRESHOLD",
  "OBS_PREV",
  "PRED_PREV",
  "PCC",
  "ATAUC",
  "STAUC",
  "ASD15",
  "TSS",
  "SENSITIVITY",
  "SPECIFICITY",
  "cAUC",
  "nAUC",
  "ValidModel_APPROACH",
  "ValidModel"
)

ModelsMets[,1]<-VALID_MODELS$Taxon
ModelsMets[,2]<-VALID_MODELS$BEST_MODEL

for(i in 1:nrow(VALID_MODELS)){
  if(VALID_MODELS[i,2]=="ENM_PA_MEAN"){
    ModelsMets[i,3]<-VALID_MODELS[which(VALID_MODELS$Taxon==ModelsMets[i,1]),"Occurrences"]
    ModelsMets[i,4]<-VALID_MODELS[which(VALID_MODELS$Taxon==ModelsMets[i,1]),"POINTS_CRITERIA"]
    ModelsMets[i,5]<-ENM_PA_MEAN[which(ENM_PA_MEAN$Taxon==ModelsMets[i,1]),"Threshold"]
    ModelsMets[i,6]<-ENM_PA_MEAN[which(ENM_PA_MEAN$Taxon==ModelsMets[i,1]),"Obs_Prev"]
    ModelsMets[i,7]<-ENM_PA_MEAN[which(ENM_PA_MEAN$Taxon==ModelsMets[i,1]),"Pred_Prev"]
    ModelsMets[i,8]<-ENM_PA_MEAN[which(ENM_PA_MEAN$Taxon==ModelsMets[i,1]),"pcc"]
    ModelsMets[i,9]<-ENM_PA_MEAN[which(ENM_PA_MEAN$Taxon==ModelsMets[i,1]),"auc"]
    ModelsMets[i,10]<-ENM_PA_MEAN[which(ENM_PA_MEAN$Taxon==ModelsMets[i,1]),"sd_auc"]
    ModelsMets[i,11]<-ENM_PA_MEAN[which(ENM_PA_MEAN$Taxon==ModelsMets[i,1]),"asd15"]
    ModelsMets[i,12]<-ENM_PA_MEAN[which(ENM_PA_MEAN$Taxon==ModelsMets[i,1]),"TSS"]
    ModelsMets[i,13]<-ENM_PA_MEAN[which(ENM_PA_MEAN$Taxon==ModelsMets[i,1]),"Sensitivity"]
    ModelsMets[i,14]<-ENM_PA_MEAN[which(ENM_PA_MEAN$Taxon==ModelsMets[i,1]),"Specificity"]
    ModelsMets[i,15]<-ENM_PA_MEAN[which(ENM_PA_MEAN$Taxon==ModelsMets[i,1]),"cAUC"]
    ModelsMets[i,16]<-ENM_PA_MEAN[which(ENM_PA_MEAN$Taxon==ModelsMets[i,1]),"nAUC"]
    ModelsMets[i,17]<-ENM_PA_MEAN[which(ENM_PA_MEAN$Taxon==ModelsMets[i,1]),"ValidModel"]
    
  }else if(VALID_MODELS[i,2]=="ENM_PA_MEDIAN"){
    ModelsMets[i,3]<-VALID_MODELS[which(VALID_MODELS$Taxon==ModelsMets[i,1]),"Occurrences"]
    ModelsMets[i,4]<-VALID_MODELS[which(VALID_MODELS$Taxon==ModelsMets[i,1]),"POINTS_CRITERIA"]
    ModelsMets[i,5]<-ENM_PA_MEDIAN[which(ENM_PA_MEDIAN$Taxon==ModelsMets[i,1]),"Threshold"]
    ModelsMets[i,6]<-ENM_PA_MEDIAN[which(ENM_PA_MEDIAN$Taxon==ModelsMets[i,1]),"Obs_Prev"]
    ModelsMets[i,7]<-ENM_PA_MEDIAN[which(ENM_PA_MEDIAN$Taxon==ModelsMets[i,1]),"Pred_Prev"]
    ModelsMets[i,8]<-ENM_PA_MEDIAN[which(ENM_PA_MEDIAN$Taxon==ModelsMets[i,1]),"pcc"]
    ModelsMets[i,9]<-ENM_PA_MEDIAN[which(ENM_PA_MEDIAN$Taxon==ModelsMets[i,1]),"auc"]
    ModelsMets[i,10]<-ENM_PA_MEDIAN[which(ENM_PA_MEDIAN$Taxon==ModelsMets[i,1]),"sd_auc"]
    ModelsMets[i,11]<-ENM_PA_MEDIAN[which(ENM_PA_MEDIAN$Taxon==ModelsMets[i,1]),"asd15"]
    ModelsMets[i,12]<-ENM_PA_MEDIAN[which(ENM_PA_MEDIAN$Taxon==ModelsMets[i,1]),"TSS"]
    ModelsMets[i,13]<-ENM_PA_MEDIAN[which(ENM_PA_MEDIAN$Taxon==ModelsMets[i,1]),"Sensitivity"]
    ModelsMets[i,14]<-ENM_PA_MEDIAN[which(ENM_PA_MEDIAN$Taxon==ModelsMets[i,1]),"Specificity"]
    ModelsMets[i,15]<-ENM_PA_MEDIAN[which(ENM_PA_MEDIAN$Taxon==ModelsMets[i,1]),"cAUC"]
    ModelsMets[i,16]<-ENM_PA_MEDIAN[which(ENM_PA_MEDIAN$Taxon==ModelsMets[i,1]),"nAUC"]
    ModelsMets[i,17]<-ENM_PA_MEDIAN[which(ENM_PA_MEDIAN$Taxon==ModelsMets[i,1]),"ValidModel"]
    
    
  }else if(VALID_MODELS[i,2]=="ENM_PA_SUM"){
    ModelsMets[i,3]<-VALID_MODELS[which(VALID_MODELS$Taxon==ModelsMets[i,1]),"Occurrences"]
    ModelsMets[i,4]<-VALID_MODELS[which(VALID_MODELS$Taxon==ModelsMets[i,1]),"POINTS_CRITERIA"]
    ModelsMets[i,5]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"Threshold"]
    ModelsMets[i,6]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"Obs_Prev"]
    ModelsMets[i,7]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"Pred_Prev"]
    ModelsMets[i,8]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"pcc"]
    ModelsMets[i,9]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"auc"]
    ModelsMets[i,10]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"sd_auc"]
    ModelsMets[i,11]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"asd15"]
    ModelsMets[i,12]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"TSS"]
    ModelsMets[i,13]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"Sensitivity"]
    ModelsMets[i,14]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"Specificity"]
    ModelsMets[i,15]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"cAUC"]
    ModelsMets[i,16]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"nAUC"]
    ModelsMets[i,17]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"ValidModel"]
    
  }else if(is.na(VALID_MODELS[i,2])){
    ModelsMets[i,3]<-VALID_MODELS[which(VALID_MODELS$Taxon==ModelsMets[i,1]),"Occurrences"]
    ModelsMets[i,4]<-VALID_MODELS[which(VALID_MODELS$Taxon==ModelsMets[i,1]),"POINTS_CRITERIA"]
    ModelsMets[i,5]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"Threshold"]
    ModelsMets[i,6]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"Obs_Prev"]
    ModelsMets[i,7]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"Pred_Prev"]
    ModelsMets[i,8]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"pcc"]
    ModelsMets[i,9]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"auc"]
    ModelsMets[i,10]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"sd_auc"]
    ModelsMets[i,11]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"asd15"]
    ModelsMets[i,12]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"TSS"]
    ModelsMets[i,13]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"Sensitivity"]
    ModelsMets[i,14]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"Specificity"]
    ModelsMets[i,15]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"cAUC"]
    ModelsMets[i,16]<-ENM_PA_SUM[which(ENM_PA_SUM$Taxon==ModelsMets[i,1]),"nAUC"]
    ModelsMets[i,17]<-NA
    
  }
};rm(i)
ModelsMets$ValidModel<-1
  
  
  
  
########################################################################################

write.csv(ENM_PA_SUM, paste0(summ_dir,"/","ensemble_metrics_PA_SUM.csv"),row.names=F,quote = F)
write.csv(ENM_PA_MEAN, paste0(summ_dir,"/","ensemble_metrics_PA_MEAN.csv"),row.names=F,quote = F)
write.csv(ENM_PA_MEDIAN, paste0(summ_dir,"/","ensemble_metrics_PA_MEDIAN.csv"),row.names=F,quote = F)
write.csv(VALID_MODELS, paste0(summ_dir,"/","VALID_MODELS.csv"),row.names=F,quote = F)
write.csv(ModelsMets, paste0(summ_dir,"/","modelsMets.csv"),row.names=F,quote = F)

############################################

for(i in 1:nrow(VALID_MODELS)){
cat("copying ",paste0(mx_out,"/","models","/",as.character(VALID_MODELS[i,1]),"/","projections","/",VALID_MODELS[i,"FILENAME"]),"\n")
  
file.copy(from=paste0(mx_out,"/","models","/",as.character(VALID_MODELS[i,1]),"/","projections","/",VALID_MODELS[i,"FILENAME"]),
          to=paste0(mx_out,"/","models","/",as.character(VALID_MODELS[i,1]),"/","projections","/",as.character(VALID_MODELS[i,1]),"_worldclim2_5_EMN_PA.asc.gz"),
          overwrite =T)
cat("NEW FILE: ",paste0(mx_out,"/","models","/",as.character(VALID_MODELS[i,1]),"/","projections","/",VALID_MODELS[i,"FILENAME"]),"\n")

};rm(i)
