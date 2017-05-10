################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
###LOADING FUNCTIONS###
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

###LOADING PACKAGES AND SOURCES###

require(PresenceAbsence);require(raster);require(SDMTools);require(dismo);require(parallel);require(geosphere);require(ape)
# src.dir <- paste("D:/SDMs/A/gap_analysis/gap_",crop,"/_scripts",sep="") # !!! change accordingly !!!
source(paste(src.dir,"/mclapply2.R",sep=""))


inputDir=crop_dir
county=F



get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}


################################################################################################################################################################
####COUNTIES ADM_AREA AS NAREA########################################################################################################
################################################################################################################################################################

COUNTIES_TILE_FUNCTION<-function(county,spList){
  
if(county==F){
  
  cat("NO COUNTIES REQUIRED...","OMITTING COUNTY APPROACH...","\n")
  counties<-list()
  
}else{
 cat("loading GADM shapefile","\n")
  adm_shp<-shapefile(paste0(geo_data_dir,"/","GAUL_2014/G2014_2013_2.shp"))
 # inputDir<-crop_dir
  mx_out <- paste(inputDir,"/maxent_modeling",sep="")
  summ_dir<-paste0(mx_out,"/","summary-files"); if (!file.exists(summ_dir)) {dir.create(summ_dir)}
  
counties<-lapply(1:length(spList), function(i){
  
  spID<-spList[[i]]
  cat("Processing ",spID,"\n")
  
  occ_points<-read.csv(paste0(inputDir,"/","occurrence_files","/",spID,".csv"),header = T)
  occ_points2<- unique(cbind(occ_points$lon,occ_points$lat))
  
  mxe_out <- paste(inputDir,"/maxent_modeling/models",sep="")
  outFolder <- paste(mxe_out,"/",spID,sep="")
  proj_dir<-paste(outFolder, "/", "projections",sep="")  
  crossval_dir<- paste(outFolder, "/", "crossval",sep="")  
  
  ENM_PA<-raster(paste(outFolder, "/", "projections","/",spID,"_worldclim2_5_EMN_PR_MEAN.asc",sep=""))
  N_AR<-ENM_PA
  N_AR[which(N_AR[]<1)]<-1
  
  ngrids<-  (cellFromXY(N_AR,cbind(occ_points2[,1],occ_points2[,2])))
  ngrids2<-matrix(nrow=length(ngrids),ncol=2)
  ngrids2[,1]<-ngrids
  for(i in 1:length(ngrids)){
    ngrids2[i,2]<- N_AR[][[ngrids[i]]]
    
  }
  
  ngrids<-as.data.frame(cbind(occ_points2,ngrids2))
  ngrids<-ngrids[complete.cases(ngrids),]
  colnames(ngrids)<-c("lon","lat","cell","NA_P")
  coordinates(ngrids)<-~lon+lat
  crs(ngrids) <- crs(adm_shp)
  scoun<-over(ngrids,adm_shp)
  shp_NA2 <- subset(adm_shp, ADM2_CODE %in% c(scoun$ADM2_CODE))
  
  return(shp_NA2)
})


names(counties)<-spList
cat("returning shapefiles","\n")

return(counties)
#rm(adm_shp);
gc()
  }
}

cat("SUBSETTING COUNTIES","\n")  
counties<-COUNTIES_TILE_FUNCTION(county=county,spList = spList);gc()  

#counties<-COUNTIES_TILE_FUNCTION(county=T,spList = c("Helianthus_annuus","Helianthus_niveus_subsp._tephrodes"));gc()

###############################################################################
#######COMPRESSING ASCII FILES##################################################  
################################################################################


compress_function<-function(spID,inputDir,OSys){
 # inputDir<-crop_dir
  mxe_out <- paste(inputDir,"/maxent_modeling/models",sep="")
  outFolder <- paste(mxe_out,"/",spID,sep="")
  
  ftoZIP <- list.files(paste(outFolder, "/projections/", sep=""), pattern=".asc$")
  cat("Compressing... \n")
  for (fz in ftoZIP) {
    fName <- paste(outFolder, "/projections/", fz, sep="")
    if (OSys == "linux") {
      system(paste("gzip", fName))
    } else {
      setwd("C:/Program Files/7-Zip")
      system(paste('7z.exe a -tgzip "',gsub(paste(fName,".gz",sep=""),pattern="/",replacement="\\\\"),'" "',gsub(fName,pattern="/",replacement="\\\\"),'"',sep=''),wait=T)
      file.remove(fName)
    }
  }
}
#################
# compress_function<-function(spID,inputDir){
#   cat("Compressing ASCII files", "\n")
#    inputDir<-crop_dir
#   mxe_out <- paste(inputDir,"/maxent_modeling/models",sep="")
#   outFolder <- paste(mxe_out,"/",spID,sep="")
#   proj_dir<-paste(outFolder, "/", "projections",sep="")
#   #proj_files<-list.files(proj_dir,pattern = ".asc$",full.names = T)
#   proj_files_T<-list.files(proj_dir,pattern = ".asc$",full.names = T)

# lapply(1:length(proj_files_T), function(k){
# 
#     write.asc.gz(x=read.asc(proj_files_T[[k]]),file=paste0(proj_files_T[[k]]))
#     file.remove(proj_files_T[[k]])
#   })
# }



#RUNNING COMPRESSING FUNCTION

# mclapply2(1:length(spList),fun=function(i){
#   compress_function(spID=spList[[i]],inputDir=crop_dir,OSys="linux")
#   
# })


################################################################################################################################################################
################################################################################################################################################################
##############################################################################################################################################################

test_ensemble_function<-function(spID,inputDir,county){
  require(PresenceAbsence);require(raster);require(SDMTools);require(dismo);require(parallel)
  #####CALLING EMSEMBLE FUNCTION#####
  
  quantile_ensemble <- function(x, quantile)
  {
    tot_ensemble <- sum(x)
    tot_ensemble[which(tot_ensemble[] < nlayers(x)*quantile)] <- NA
    tot_ensemble[which(!is.na(tot_ensemble[]))] <- 1
    return(tot_ensemble)
  }
  ####CALLING OCCURRENCES DATA####
  
  #####CALLING NATIVE AREA#####
  mxe_out <- paste(inputDir,"/maxent_modeling/models",sep="")
  outFolder <- paste(mxe_out,"/",spID,sep="")
  proj_dir<-paste(outFolder, "/", "projections",sep="")  
  crossval_dir<- paste(outFolder, "/", "crossval",sep="")  
  
 if(file.exists(paste0(proj_dir,"/",spID,"_worldclim2_5_EMN_PR_ensemble.asc"))& file.exists(paste0(proj_dir,"/",spID,"_worldclim2_5_EMN_PA_SUM_ensemble.asc"))){
   cat("omitting ", as.character(spID),"\n")
    
 }else{
   cat("processing ensemble PR median and PA Sum for ", as.character(spID),"\n")
    
    ENM_PA<-raster(paste(outFolder, "/", "projections","/",spID,"_worldclim2_5_EMN_PR_MEAN.asc",sep=""))
    N_AR<-ENM_PA
    N_AR[which(N_AR[]<1)]<-1
    #####CALLING MODELS AND THRESHOLDS#####
    maxent_results<-as.data.frame(read.csv(paste0(crossval_dir,"/","UpperLeftROC.csv"),header=T))
    #maxent_th<-maxent_results[1:10,"Maximum.training.sensitivity.plus.specificity.logistic.threshold"]
    maxent_th<-maxent_results[1:nrow(maxent_results),"UpperLeftROC"]
    
    ####################################################################################
    
    
    #####THRESHOLDING MODELS AND FITTING TO NATIVE AREA#####
    #models_f<-paste0(proj_dir,"/",spID,"_worldclim2_5_f",1:10,".asc")
    models_f<-list.files(proj_dir,pattern = "_worldclim2_5_f",full.names = T)
    
    models_f<-lapply(models_f,raster)
    models_f<-stack(models_f)
    
    #models<-paste0(proj_dir,"/",spID,"_worldclim2_5_f",1:10,".asc")
    models<-list.files(proj_dir,pattern = "_worldclim2_5_f",full.names = T)
    
    models<-lapply(1:length(models),function(i){
      mod<-raster(models[[i]])
      mod[which(mod[] < maxent_th[[i]])] <- 0
      mod[which(mod[] != 0)] <- 1
      mod <- mod * N_AR
      return(mod)
    })
    
    models<-stack(models)
    
    result_PR <- calc(models_f, fun=median)
    result_PR<-result_PR*N_AR
    if(county==F){
      cat("using county approach for ",as.character(spID),"\n")
      writeRaster(result_PR, filename=paste0(proj_dir,"/",spID,"_worldclim2_5_EMN_PR_ensemble.asc"), overwrite=T)
    }else{
      result_PR<-mask(result_PR,counties[[spID]]) #COUNTY APPROACH
      writeRaster(result_PR, filename=paste0(proj_dir,"/",spID,"_worldclim2_5_EMN_PR_ensemble.asc"), overwrite=T)
    }
    
    
    #####WRITING MODELS AND PERFORMING QUANTILE ENSEMBLE APPROACH#####
    
    result <- quantile_ensemble(models, quantile=0.5)
    
    #
    if(county==F){
      cat("using county approach for ",as.character(spID),"\n")
      writeRaster(result, filename=paste0(proj_dir,"/",spID,"_worldclim2_5_EMN_PA_SUM_ensemble.asc"), overwrite=T)
      
    }else{
      result<-mask(result,counties[[spID]]) #COUNTY APPROACH
      writeRaster(result, filename=paste0(proj_dir,"/",spID,"_worldclim2_5_EMN_PA_SUM_ensemble.asc"), overwrite=T)
      
    }

    sum_met<-sum(models)/nlayers(models)
    
    if(county==F){
      cat("using county approach for ",as.character(spID),"\n")
      writeRaster(sum_met, filename=paste0(proj_dir,"/",spID,"_worldclim2_5_EMN_PA_MET_ensemble.asc"), overwrite=T)
      
    }else{
        sum_met<-mask(sum_met,counties[[spID]]) #COUNTY APPROACH
        writeRaster(sum_met, filename=paste0(proj_dir,"/",spID,"_worldclim2_5_EMN_PA_MET_ensemble.asc"), overwrite=T)
        
    }
  }
}


cat("CALCULATING PA_SUM AND MEDIANS","\n")  

# mclapply2(1:length(spList),fun=function(i){
#   test_ensemble_function(spID=spList[[i]],inputDir=inputDir,county=F)
#   
# });gc()

lapply(1:length(spList),function(i){
  test_ensemble_function(spID=spList[[i]],inputDir=inputDir,county=F)

});gc()


################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

################################################################################
#########ENSEMBLE FUNCTION PR  MEDIAN###########################################  
################################################################################

cat("THRESHOLDING MEDIANS","\n")  

#MEDIAN_ENSEMBLE_APPROACH<-function(inputDir,spList,county){
  
  predictions <- lapply(1:length(spList), function(i){
    spID<-spList[[i]]
   #spID<-spList[[23]]
   # inputDir<-crop_dir
    mxe_out <- paste(inputDir,"/maxent_modeling/models",sep="")
    outFolder <- paste(mxe_out,"/",spID,sep="")
    proj_dir<-paste(outFolder, "/", "projections",sep="")  
    model <- raster(paste0(proj_dir,"/",spID,"_worldclim2_5_EMN_PR_ensemble.asc"))
    
    
    occVal <- read.csv(paste(inputDir, "/occurrence_files","/",spID,".csv", sep=""))
    occVal <- occVal[,c("lon","lat")]
    occVal$observed <- 1  
    #
    ENM_PA<-raster(paste(outFolder, "/", "projections","/",spID,"_worldclim2_5_EMN_PR_MEAN.asc",sep=""))
    N_AR<-ENM_PA
    N_AR[which(N_AR[]<1)]<-1
    if(county==T){
      cat("using county approach for ",as.character(spID),"\n")
      N_AR<-mask(N_AR,counties[[spID]]) #COUNTY APPROACH
    }else{
      N_AR<-N_AR
    }
    
    model<-model*N_AR
    # Background unicamente en área nativa
    bckVal<-read.csv(paste0(inputDir,"/","maxent_modeling/background","/",spID,".csv"),header=T)
    bckVal <- bckVal[,c("lon","lat")]
    bckVal <- as.data.frame(cbind(bckVal,extract(N_AR,bckVal)))
    names(bckVal) <- c("lon","lat","ver")
    bckVal <- bckVal[which(bckVal$ver==1),]
    bckVal$ver <- NULL
    bckVal$observed <- 0
    
    allVal <- as.data.frame(rbind(occVal, bckVal))
    allVal <- as.data.frame(cbind(allVal, extract(model,allVal[,c("lon","lat")])))
    names(allVal)[4] <- "predicted"
    
    allVal$spID <- spList[[i]]
    
    return(allVal)
    
  })
  
  ####CALCULATING METRICS FOR PR_ENSEMBLE#########################################
  
  predictions3<-as.data.frame(matrix(nrow = length(spList),ncol=11))
  
  for(i in 1:length(spList)){
    z <- predictions[[i]]
    z <- rbind(z[which(z$observed==1 & !is.na(z$predicted)),], z[which(z$observed==0),])
    z<-z[complete.cases(z),]
    z <- as.data.frame(z)
    spIDM<-spList[[i]]
   # spIDM<-spList[[23]]
    z <- z[,c("spID","observed","predicted")]
    
    predictions3[i,1]<-spIDM
    predictions3[i,2]<-PresenceAbsence::auc(z)[1]
    predictions3[i,3]<-PresenceAbsence::auc(z,st.dev = T)[2]
    predictions3[i,4]<-PresenceAbsence::optimal.thresholds(z,opt.methods = 9)[2]
    predictions3[i,5]<-PresenceAbsence::pcc(cmx(z))[1]
    predictions3[i,6]<-PresenceAbsence::Kappa(cmx(z))[1]
    predictions3[i,7]<-PresenceAbsence::sensitivity(cmx(z))[1]
    predictions3[i,8]<-PresenceAbsence::specificity(cmx(z))[1]
    predictions3[i,9]<-(PresenceAbsence::sensitivity(cmx(z))[1]+PresenceAbsence::specificity(cmx(z))[1])-1
    predictions3[i,10]<-PresenceAbsence::predicted.prevalence(z,threshold =PresenceAbsence::optimal.thresholds(z,opt.methods = 9)[2])[2]
    predictions3[i,11]<-PresenceAbsence::predicted.prevalence(z,threshold =PresenceAbsence::optimal.thresholds(z,opt.methods = 9)[2])[3]
  };rm(i)
  
  colnames(predictions3)<-c("Taxon","auc","sd_auc","Threshold","pcc","Kappa","Sensitivity","Specificity","TSS","Obs_Prev","Pred_Prev")
  
  results <- predictions3
  cat("Writting ensemble metrics csv file for PR_ENSEMBLE PR", "\n")
  
  mx_out <- paste(inputDir,"/maxent_modeling",sep="")
  summ_dir<-paste0(mx_out,"/","summary-files"); if (!file.exists(summ_dir)) {dir.create(summ_dir)}
  
  write.csv(results, paste0(summ_dir,"/","ensemble_metrics_PA_MEDIAN.csv"),row.names=F,quote = F)
  
  #######THRESHOLDING USING PR_ENSEMBLE###########################################################################################################################  
  
 lapply(1:length(spList),function(i){
  #i=1
    spID<-spList[[i]]
    cat("Thresholding PR Median raster for", as.character(spID),"\n")
    #inputDir<-crop_dir
    mxe_out <- paste(inputDir,"/maxent_modeling/models",sep="")
    outFolder <- paste(mxe_out,"/",spID,sep="")
    proj_dir<-paste(outFolder, "/", "projections",sep="")  
    PR_MODEL <- raster(paste0(proj_dir,"/",spID,"_worldclim2_5_EMN_PR_ensemble.asc"))
    thr<-predictions3[which(predictions3$Taxon==spID),]
    thr<-thr[1,4]
    PR_MODEL[which(PR_MODEL[] < thr)] <- 0
    PR_MODEL[which(PR_MODEL[] != 0)] <- 1
    
    if(county==T){
      PR_MODEL<-mask(PR_MODEL,counties[[spID]]) #COUNTY APPROACH
    }else{
      PR_MODEL<-PR_MODEL
    }
    
    writeRaster(PR_MODEL, filename=paste0(proj_dir,"/",spID,"_worldclim2_5_EMN_PR_THR_ensemble.asc"), overwrite=T)
  })
#}


#######################################################################
############METRICS USING PA_SUM###########################################  
################################################################################

  cat("STATISTICS PA_SUM","\n")  
  
  
#PA_SUM_ENSEMBLE_APPROACH<-function(inputDir,spList,county){
  
  predictions_PA_SUM <- lapply(1:length(spList), function(i){
    spID<-spList[[i]]
    cat("PA SUM for ",as.character(spID),"\n")
    #inputDir<-crop_dir
    mxe_out <- paste(inputDir,"/maxent_modeling/models",sep="")
    outFolder <- paste(mxe_out,"/",spID,sep="")
    proj_dir<-paste(outFolder, "/", "projections",sep="")  
    model <- raster(paste0(proj_dir,"/",spID,"_worldclim2_5_EMN_PA_MET_ensemble.asc"))
    
    
    occVal <- read.csv(paste(inputDir, "/occurrence_files","/",spID,".csv", sep=""))
    occVal <- occVal[,c("lon","lat")]
    occVal$observed <- 1  
    #
    ENM_PA<-raster(paste(outFolder, "/", "projections","/",spID,"_worldclim2_5_EMN_PR_MEAN.asc",sep=""))
    N_AR<-ENM_PA
    N_AR[which(N_AR[]<1)]<-1
    
    if(county==T){
      cat("using county approach for ",as.character(spID),"\n")
      N_AR<-mask(N_AR,counties[[spID]]) #COUNTY APPROACH
    }else{
      N_AR<-N_AR
    }
    # Background unicamente en área nativa
    bckVal<-read.csv(paste0(inputDir,"/","maxent_modeling/background","/",spID,".csv"),header=T)
    bckVal <- bckVal[,c("lon","lat")]
    bckVal <- as.data.frame(cbind(bckVal,extract(N_AR,bckVal)))
    names(bckVal) <- c("lon","lat","ver")
    bckVal <- bckVal[which(bckVal$ver==1),]
    bckVal$ver <- NULL
    bckVal$observed <- 0
    
    allVal <- as.data.frame(rbind(occVal, bckVal))
    allVal <- as.data.frame(cbind(allVal, extract(model,allVal[,c("lon","lat")])))
    names(allVal)[4] <- "predicted"
    
    allVal$spID <- spList[[i]]
    
    return(allVal)
    
  })
  
  ####METRICS USING PREDICTIONS PA
  
  
  predictions4<-as.data.frame(matrix(nrow = length(spList),ncol=11))
  
  for(i in 1:length(spList)){
    z <- predictions_PA_SUM[[i]]
    z <- rbind(z[which(z$observed==1 & !is.na(z$predicted)),], z[which(z$observed==0),])
    z<-z[complete.cases(z),]
    z <- as.data.frame(z)
    spIDM<-spList[[i]]
    z <- z[,c("spID","observed","predicted")]
    predictions4[i,1]<-spIDM
    predictions4[i,2]<-PresenceAbsence::auc(z)[1]
    predictions4[i,3]<-PresenceAbsence::auc(z,st.dev = T)[2]
    predictions4[i,4]<-NA
    predictions4[i,5]<-PresenceAbsence::pcc(cmx(z))[1]
    predictions4[i,6]<-PresenceAbsence::Kappa(cmx(z))[1]
    predictions4[i,7]<-PresenceAbsence::sensitivity(cmx(z))[1]
    predictions4[i,8]<-PresenceAbsence::specificity(cmx(z))[1]
    predictions4[i,9]<-(PresenceAbsence::sensitivity(cmx(z))[1]+PresenceAbsence::specificity(cmx(z))[1])-1
    predictions4[i,10]<-PresenceAbsence::predicted.prevalence(z,threshold =PresenceAbsence::optimal.thresholds(z,opt.methods = 9)[2])[2]
    predictions4[i,11]<-PresenceAbsence::predicted.prevalence(z,threshold =PresenceAbsence::optimal.thresholds(z,opt.methods = 9)[2])[3]
    
  };rm(i)
  
  colnames(predictions4)<-c("Taxon","auc","sd_auc","Threshold","pcc","Kappa","Sensitivity","Specificity","TSS","Obs_Prev","Pred_Prev")
  
  results2 <- predictions4
  cat("Writting ensemble metrics csv file for PA_ENSEMBLE SUM", "\n")
  
  mx_out <- paste(inputDir,"/maxent_modeling",sep="")
  summ_dir<-paste0(mx_out,"/","summary-files"); if (!file.exists(summ_dir)) {dir.create(summ_dir)}
  
  write.csv(results2, paste0(summ_dir,"/","ensemble_metrics_PA_SUM.csv"),row.names=F,quote = F)
  
#}


################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
##METRICS PA ENSEMBLE MEAN###############
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################


  cat("THRESHOLDING MEANS","\n")  
  
#MEAN_ENSEMBLE_APPROACH<-function(inputDir,spList,county){
  
  
  predictions_PA_MEAN <- lapply(1:length(spList), function(i){
    spID<-spList[[i]]
    #inputDir<-crop_dir
    mxe_out <- paste(inputDir,"/maxent_modeling/models",sep="")
    outFolder <- paste(mxe_out,"/",spID,sep="")
    proj_dir<-paste(outFolder, "/", "projections",sep="")  
    model <- raster(paste0(proj_dir,"/",spID,"_worldclim2_5_EMN_PR_MEAN.asc"))
    
    
    occVal <- read.csv(paste(inputDir, "/occurrence_files","/",spID,".csv", sep=""))
    occVal <- occVal[,c("lon","lat")]
    occVal$observed <- 1  
    #
    ENM_PA<-raster(paste(outFolder, "/", "projections","/",spID,"_worldclim2_5_EMN_PR_MEAN.asc",sep=""))
    N_AR<-ENM_PA
    N_AR[which(N_AR[]<1)]<-1
    
    if(county==T){
      cat("using county approach for ",as.character(spID),"\n")
      N_AR<-mask(N_AR,counties[[spID]]) #COUNTY APPROACH
    }else{
      N_AR<-N_AR
    }
    
    model<-model*N_AR
    # Background unicamente en área nativa
    bckVal<-read.csv(paste0(inputDir,"/","maxent_modeling/background","/",spID,".csv"),header=T)
    bckVal <- bckVal[,c("lon","lat")]
    bckVal <- as.data.frame(cbind(bckVal,extract(N_AR,bckVal)))
    names(bckVal) <- c("lon","lat","ver")
    bckVal <- bckVal[which(bckVal$ver==1),]
    bckVal$ver <- NULL
    bckVal$observed <- 0
    
    allVal <- as.data.frame(rbind(occVal, bckVal))
    allVal <- as.data.frame(cbind(allVal, extract(model,allVal[,c("lon","lat")])))
    names(allVal)[4] <- "predicted"
    
    allVal$spID <- spList[[i]]
    
    return(allVal)
    
  })
  
  
  ####METRICS USING PREDICTIONS PA_MEAN
  
  
  predictions5<-as.data.frame(matrix(nrow = length(spList),ncol=11))
  
  for(i in 1:length(spList)){
    z <- predictions_PA_MEAN[[i]]
    z <- rbind(z[which(z$observed==1 & !is.na(z$predicted)),], z[which(z$observed==0),])
    z<-z[complete.cases(z),]
    z <- as.data.frame(z)
    spIDM<-spList[[i]]
    z <- z[,c("spID","observed","predicted")]
    predictions5[i,1]<-spIDM
    predictions5[i,2]<-PresenceAbsence::auc(z)[1]
    predictions5[i,3]<-PresenceAbsence::auc(z,st.dev = T)[2]
    predictions5[i,4]<-PresenceAbsence::optimal.thresholds(z,opt.methods = 9)[2]
    predictions5[i,5]<-PresenceAbsence::pcc(cmx(z))[1]
    predictions5[i,6]<-PresenceAbsence::Kappa(cmx(z))[1]
    predictions5[i,7]<-PresenceAbsence::sensitivity(cmx(z))[1]
    predictions5[i,8]<-PresenceAbsence::specificity(cmx(z))[1]
    predictions5[i,9]<-(PresenceAbsence::sensitivity(cmx(z))[1]+PresenceAbsence::specificity(cmx(z))[1])-1
    predictions5[i,10]<-PresenceAbsence::predicted.prevalence(z,threshold =PresenceAbsence::optimal.thresholds(z,opt.methods = 9)[2])[2]
    predictions5[i,11]<-PresenceAbsence::predicted.prevalence(z,threshold =PresenceAbsence::optimal.thresholds(z,opt.methods = 9)[2])[3]
    
    
  };rm(i)
  
  colnames(predictions5)<-c("Taxon","auc","sd_auc","Threshold","pcc","Kappa","Sensitivity","Specificity","TSS","Obs_Prev","Pred_Prev")
  
  results3 <- predictions5
  cat("Writting ensemble metrics csv file for PA_ENSEMBLE MEAN", "\n")
  
  mx_out <- paste(inputDir,"/maxent_modeling",sep="")
  summ_dir<-paste0(mx_out,"/","summary-files"); if (!file.exists(summ_dir)) {dir.create(summ_dir)}
  
  write.csv(results3, paste0(summ_dir,"/","ensemble_metrics_PA_MEAN.csv"),row.names=F,quote = F)
  
  
  #######THRESHOLDING USING PA MEAN_ENSEMBLE###########################################################################################################################  
  
  lapply(1:length(spList),function(i){
    spID<-spList[[i]]
    cat("Thresholding PA Mean raster for", as.character(spID),"\n")
    #inputDir<-crop_dir
    mxe_out <- paste(inputDir,"/maxent_modeling/models",sep="")
    outFolder <- paste(mxe_out,"/",spID,sep="")
    proj_dir<-paste(outFolder, "/", "projections",sep="")  
    PR_MODEL <- raster(paste0(proj_dir,"/",spID,"_worldclim2_5_EMN_PR_MEAN.asc"))
    thr<-predictions5[which(predictions5$Taxon==spID),]
    thr<-thr[1,4]
    PR_MODEL[which(PR_MODEL[] < thr)] <- 0
    PR_MODEL[which(PR_MODEL[] != 0)] <- 1
    
    if(county==T){
      PR_MODEL<-mask(PR_MODEL,counties[[spID]]) #COUNTY APPROACH
    }else{
      PR_MODEL<-PR_MODEL
    }
    writeRaster(PR_MODEL, filename=paste0(proj_dir,"/",spID,"_worldclim2_5_EMN_PA_mean.asc"), overwrite=T)
  })
  
#}

  cat("COMPRESSING FILES IN GZ FORMAT","\n")  
  
  
  mclapply2(1:length(spList),fun=function(i){
    compress_function(spID=spList[[i]],inputDir=inputDir,OSys=as.character(get_os()))
    
  });gc()

  
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################################################################################################

#spID<-"Helianthus_niveus_subsp._tephrodes"
#spID<-"Helianthus_annuus"


#ENSEMBLE_FUNCTION_RASTER<-function(inputDir,county,spList){
  
 # 

  
###PA_SUM AND MEDIAN CALCULATIONS###
  


# 
# lapply(1:length(spList),function(i){
#   test_ensemble_function(spID=spList[[i]],inputDir=crop_dir,county=F)
#   
# })#;gc()
# # ###MEDIAN CALCULATIONS AND THRESHOLD###
  
#cat("MEDIAN CALCULATIONS AND THRESHOLDING","\n")  


# x<-MEDIAN_ENSEMBLE_APPROACH(inputDir,spList,county);gc()
#   
# ###PA_SUM CALCULATIONS###
# 
# cat("PA_SUM CALCULATIONS","\n")  
# 
# x<-PA_SUM_ENSEMBLE_APPROACH(inputDir,spList,county);gc()
# 
# ###MEAN CALCULATIONS AND THRESHOLD###
# 
# cat("MEAN CALCULATIONS AND THRESHOLDING","\n")  
# 
# 
# x<-MEAN_ENSEMBLE_APPROACH(inputDir,spList,county);gc()
#   
#RUNNING COMPRESSING FUNCTION


# 
# cat("     ", "\n")
# cat("DONE!", "\n")
# cat("     ", "\n")
# 
# }


#spList <- list.files(paste(crop_dir, "/occurrence_files", sep=""),pattern=".csv")
#spList<-sub(".csv","",spList)
#spList<-c("Helianthus_annuus","Helianthus_niveus_subsp._tephrodes")

  

