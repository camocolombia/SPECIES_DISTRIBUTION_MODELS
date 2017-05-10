
#crop<-"Vegies_root_Kantar"
#crop<-"Ornamental_medicinal_McCoy"
crop<-"Industrial_indust_Dierig"


#basic stuff - where is the code
#src.dir <- paste("/curie_data/ncastaneda/gap-analysis/gap_",crop,"/_scripts",sep="") # !!! change accordingly !!!
src.dir <- paste("/mnt/workspace_cluster_6/ccsosa/GAP_ANALYSIS_US_2016/gap_analysis/gap_",crop,"/_scripts",sep="") # !!! change accordingly !!!

#src.dir <- paste("//dapadfs/workspace_cluster_6/CWR/CWR_PROJECT_CC_BD/ccsosa/GAP_ANALYSIS_US_2016/gap_analysis/gap_",crop,"/_scripts",sep="") # !!! change accordingly !!!

#gap.dir <-"/curie_data/ncastaneda/gap-analysis" # !!! change accordingly !!!

#gap.dir <-"//dapadfs/workspace_cluster_6/CWR/CWR_PROJECT_CC_BD/ccsosa/GAP_ANALYSIS_US_2016/gap_analysis" # !!! change accordingly !!!
gap.dir <-"/mnt/workspace_cluster_6/ccsosa/GAP_ANALYSIS_US_2016/gap_analysis" # !!! change accordingly !!!

#crop details
crop_dir <- paste(gap.dir,"/gap_",crop,sep="")

if (!file.exists(crop_dir)) {dir.create(crop_dir)}
setwd(crop_dir)

source(paste(src.dir,"/mclapply2.R",sep=""))
require(parallel)

spList <- list.dirs(paste(crop_dir, "/maxent_modeling/models", sep=""),recursive = F)
spList<-sub(paste(crop_dir, "/maxent_modeling/models/", sep=""),"",spList)

inputDir<-crop_dir
uncompress_function<-function(spID,inputDir){
  
  library(R.utils)
  
  # inputDir<-crop_dir
  mxe_out <- paste(inputDir,"/maxent_modeling/models",sep="")
  outFolder <- paste(mxe_out,"/",spID,sep="")
  
  ftoZIP <- list.files(paste(outFolder, "/projections/", sep=""), pattern=".gz$",full.names = T)
  if(length(ftoZIP)==0){
    cat("skipping ",spID,"\n")
  }else{
  cat("Uncompressing ",spID,"\n")
  for (i in 1:length(ftoZIP)){
   fz<-ftoZIP[[i]] 
    x<-gunzip(fz)
      file.remove(fz)
    }
  }
}  


lapply(1:length(spList),function(i){
  x<-uncompress_function(spID=spList[[i]],inputDir=inputDir)
})

# mclapply2(1:length(spList),fun=function(i){
#   uncompress_function(spID=spList[[i]],inputDir=inputDir)
#   
# });gc()


##########

uncompress_function_NAREA<-function(spID,inputDir){
  
  library(R.utils)
  
  # inputDir<-crop_dir
  mxe_out <- paste(inputDir,"/maxent_modeling/models",sep="")
  outFolder <- paste(mxe_out,"/",spID,sep="")
  
  ftoZIP <- list.files(paste(outFolder, "/projections/NAREA_APPROACH", sep=""), pattern=".gz$",full.names = T)
  if(length(ftoZIP)==0){
    cat("skipping ",spID,"\n")
  }else{
    cat("Uncompressing ",spID,"\n")
    for (i in 1:length(ftoZIP)){
      fz<-ftoZIP[[i]] 
      x<-gunzip(fz)
      file.remove(fz)
    }
  }
}  

lapply(1:length(spList),function(i){
  x<-uncompress_function_NAREA(spID=spList[[i]],inputDir=inputDir)
})