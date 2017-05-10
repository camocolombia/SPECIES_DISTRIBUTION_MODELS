require(rgdal)
require(raster)

source(paste(src.dir,"/000.getMetrics.R",sep=""))
source(paste(src.dir,"/000.zipRead.R",sep=""))
source(paste(src.dir,"/000.zipWrite.R",sep=""))
source(paste(src.dir,"/000.createChullBuffer.R",sep=""))

nullModel_calculate <- function(spID, inputDir) {
  
  mxe_out <- paste(inputDir,"/maxent_modeling/models",sep="")
  if (!file.exists(mxe_out)) {dir.create(mxe_out)}
  spDir <- paste(mxe_out,"/",spID,sep="")
  
  inProjClimDir <- paste(inputDir, "/biomod_modeling/current-clim", sep="")
  maxentApp <- paste(inputDir, "/maxent_modeling/lib/maxent.jar", sep="")
  mskDir <- paste(inputDir, "/masks", sep="")
  backoutdir <- paste(inputDir, "/maxent_modeling/background", sep="")
  NADir <- paste(inputDir, "/biomod_modeling/native-areas/asciigrids", sep="")
  
  cat("Taxon ", spID, "\n")
  
  #1. Load species data
  
  occFile <- paste(inputDir, "/occurrence_files/", spID, ".csv", sep="")
  
  if (file.exists(occFile)) {
    
    #1.1 Load the data
    
    library(dismo)
    #1.1a Load the data and fitting to a Native area
    
    
    NAREA_COUNTY_DIR<-paste0(inputDir, "/NAREAS/ASCII", sep="")
    
    if(file.exists(paste0(NAREA_COUNTY_DIR,"/",spID,"/","narea.asc.gz"))){
      N_AR<-zipRead(path=paste0(NAREA_COUNTY_DIR,"/",spID),fname="narea.asc.gz")
    }else{
      
      N_AR<-zipRead(path=paste0(crop_dir,"/","biomod_modeling/native-areas/asciigrids","/",spID),fname="narea.asc.gz")
    }
    
    NAREA<-N_AR
    #NAREA<-zipRead(paste0(NADir,"/",as.character(spID)),"narea.asc.gz")
      #raster(paste0(NADir,"/",as.character(spID),"/","narea.asc"))
    NAREA[which(NAREA[]<1)]<-NA
    
    inData <- read.csv(occFile)
    
    inData<-cbind(inData,raster::extract(NAREA,inData[,2:3]))
    inData<-inData[!is.na(inData[,ncol(inData)]),]
    
    nOcc <- nrow(inData)
    
    env_data <- list.files(inProjClimDir, pattern=".asc$", full.names=T)
    env_data <- lapply(env_data, raster)
    
    if(nOcc==0){
      
      cat("NO OCCS for",as.character(spID)," SKIPPING...","\n")
      AUC<-NA
    }else if(nOcc > 2){
      cat(nOcc," records for ",as.character(spID),"\n")
      
      # Step 2: training and testing sets
      set.seed(1235)
      training <- inData[which(kfold(inData, k=3)!=3),c("lon", "lat")]
      rownames(training) <- 1:nrow(training)
      set.seed(1235)
      testing  <- inData[which(kfold(inData, k=3)==3),c("lon", "lat")]
      rownames(testing) <- 1:nrow(testing)
      
      # Step 3: Geographic distance model
      Null_model <- geoDist(training, lonlat=T)
      
      # Step 4: evaluate model
      pabs <- randomPoints(mask=env_data[[1]], n=10000)
      pabs <- as.data.frame(pabs)
      names(pabs) <- c("lon", "lat")
      eval <- dismo::evaluate(p=testing, a=pabs, model=Null_model)
      AUC <- eval@auc
      
    }else if(nOcc ==2){
      cat("WARNING...",nOcc," records for ",as.character(spID),"\n")
    set.seed(1235)
    training <- inData[which(kfold(inData, k=2)!=2),c("lon", "lat")]
    rownames(training) <- 1:nrow(training)
    set.seed(1235)
    testing  <- inData[which(kfold(inData, k=2)==2),c("lon", "lat")]
    rownames(testing) <- 1:nrow(testing)
    
    # Step 3: Geographic distance model
    Null_model <- geoDist(training, lonlat=T)
    
    # Step 4: evaluate model
    pabs <- randomPoints(mask=env_data[[1]], n=10000)
    pabs <- as.data.frame(pabs)
    names(pabs) <- c("lon", "lat")
    eval <- dismo::evaluate(p=testing, a=pabs, model=Null_model)
    AUC <- eval@auc
    }else if(nOcc ==1){
      cat("WARNING...",nOcc," records for ",as.character(spID),"\n")
      set.seed(1235)
      training <- inData[which(kfold(inData, k=1)==1),c("lon", "lat")]
      rownames(training) <- 1:nrow(training)
      set.seed(1235)
      testing  <- inData[which(kfold(inData, k=1)==1),c("lon", "lat")]
      rownames(testing) <- 1:nrow(testing)
      
      # Step 3: Geographic distance model
      Null_model <- geoDist(training, lonlat=T)
      
      # Step 4: evaluate model
      pabs <- randomPoints(mask=env_data[[1]], n=10000)
      pabs <- as.data.frame(pabs)
      names(pabs) <- c("lon", "lat")
      eval <- dismo::evaluate(p=testing, a=pabs, model=Null_model)
      AUC <- eval@auc
      
      
      
      #AUC<-0.5
      }
  }  
  
  return(AUC)
  
}


spList <- list.dirs(paste(crop_dir, "/maxent_modeling/models", sep=""),recursive = F)
spList <- gsub(pattern=paste(crop_dir, "/maxent_modeling/models/", sep=""), replacement="", x=spList)

#spList <- list.files(paste(crop_dir, "/occurrence_files", sep=""),pattern=".csv")
#spList <- gsub(pattern=".csv", replacement="", x=spList)

if(!file.exists(paste(crop_dir, "/maxent_modeling/summary-files/nullModel_auc_NAREA.csv", sep="")))
{
  results <- list(0)
  for(i in 1:length(spList))
 # for(i in 73:length(spList))
  {
    auc_vec <- nullModel_calculate(spID=spList[i], inputDir=crop_dir)
    results[[i]] <- data.frame(SPID=paste(spList[i]), AUC=auc_vec)
  }
  
  results <- Reduce(function(...) rbind(..., deparse.level=1), results)
  results <- results[order(results$SPID),]
  rownames(results) <- 1:nrow(results)
  
  write.csv(results, paste(crop_dir, "/maxent_modeling/summary-files/nullModel_auc_NAREA.csv", sep=""), row.names=F)
} else {
cat("NULL MODEL PROCESSED!","\n")  
}



  # nullModel <- read.csv(paste(crop_dir, "/maxent_modeling/summary-files/nullModel_auc.csv", sep=""))
  # nullModel <- nullModel[order(nullModel$SPID),]
  # rownames(nullModel) <- 1:nrow(nullModel)
  # 
  # calc_cauc <- function(x){
  #   AUC  <- newpriFile$ATAUC[newpriFile$TAXON==x]
  #   nAUC <- nullModel$AUC[nullModel==x]
  #   cAUC <- AUC+.5-max(c(.5,nAUC,na.rm=T))
  #   return(cAUC)
  # }
  # calc_cauc <- Vectorize(FUN=calc_cauc, vectorize.args=x)
  # inputDir<-crop_dir
  # mx_out <- paste(inputDir,"/maxent_modeling",sep="")
  # summ_dir<-paste0(mx_out,"/","summary-files"); if (!file.exists(summ_dir)) {dir.create(summ_dir)}
  # newpriFile<-read.csv(paste0(summ_dir,"/","modelsMets.csv"),header=T)
  # cauc <- calc_cauc(as.character(newpriFile$TAXON))
  # newpriFile$cAUC <- cauc
  # newpriFile$nAUC<-nullModel$AUC
  # 
 # newpriFile <- newpriFile[,c("TAXON","TOTAL_RP","TOTAL","HS","HS_RP","GS","GS_RP","SRS","ATAUC","cAUC","STAUC","ASD15","IS_VALID","GRS","ERS","FPS","FPCAT","MAP_AVAILABLE")]
  
#}



