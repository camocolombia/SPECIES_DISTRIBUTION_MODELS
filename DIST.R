require(stringdist);require(psych);library(geosphere)

#####################################################
DATASETS_PATH<- "D:/Dropbox/Dropbox/NPGS georeferencing project/DATASETS/DATASETS_GEOREF"



IDs<-as.data.frame(read.csv(paste0("D:/Dropbox/Dropbox/NPGS georeferencing project/DATASETS/","acc_id.csv"),sep="|",header=T,na.strings = ""));gc()


IDs$lat_RAW<-NA
IDs$lon_RAW<-NA
###h
IDs$lat_RDGI<-NA
IDs$lon_RDGI<-NA
###
IDs$lat_IR<-NA
IDs$lon_IR<-NA
###
IDs$DIST_RAW_RDGI<-NA #DIST RAW - IRRI IMPROVED GEOGRAPHY
IDs$THR_RAW_RDGI<-NA
###
IDs$DIST_RAW_IR<-NA #DIST RAW -  HIJMANS
IDs$THR_RAW_IR<-NA
###
IDs$DIST_RDGI_IR<-NA #DIST IRRI IMPROVED GEOGRAPHY -  HIJMANS
IDs$THR_RDGI_IR<-NA
###
gc()
#####################################################
# IR<-read.csv(paste0(DATASETS_PATH,"/","TO_USDA_NPGS_GRIN_IRRI_FF.csv"),sep="|",header=T)     ## HIJMANS RESULTS
# RDGI<-read.csv(paste0(DATASETS_PATH,"/","raw_USDA_NPGS_GRIN_IRRI_FF.csv"),sep="|",header=T)  ## IRRI IMPROVED GEOGRAPHY
# RAW<-read.csv(paste0(DATASETS_PATH,"/","GRIN_GLOBAL_ORI_DATA_FF.csv"),sep="|",header=T)      ## RAW DATA FROM GRIN GLOBAL
# gc()
#####################################################

is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

#####################################################

COORDS<-as.data.frame(read.csv(paste0("D:/Dropbox/Dropbox/NPGS georeferencing project/DATASETS/","acc_id_coords.csv"),sep="|",header=T,
                               dec = ".",na.strings = ""));gc()




COORDS$raw_lat<-as.numeric(as.character(COORDS$raw_lat))
COORDS$raw_lon<-as.numeric(as.character(COORDS$raw_lon))
COORDS$GI_lat<-as.numeric(as.character(COORDS$GI_lat))
COORDS$GI_lon<-as.numeric(as.character(COORDS$GI_lon))
COORDS$IRRI_lat<-as.numeric(as.character(COORDS$IRRI_lat))
COORDS$IRRI_lon<-as.numeric(as.character(COORDS$IRRI_lon))




RAW<-as.data.frame(cbind(COORDS$ACID,COORDS$raw_lon,COORDS$raw_lat));colnames(RAW)<-c("ACID","LON","LAT");gc()## RAW DATA FROM GRIN GLOBAL
RDGI<-as.data.frame(cbind(COORDS$ACID,COORDS$GI_lon,COORDS$GI_lat));colnames(RDGI)<-c("ACID","LON","LAT");gc()## IRRI IMPROVED GEOGRAPHY
IR<-as.data.frame(cbind(COORDS$ACID,COORDS$IRRI_lon,COORDS$IRRI_lat));colnames(IR)<-c("ACID","LON","LAT");gc() ## HIJMANS RESULTS



IDs$DIST_RAW_RDGI<-NA #DIST RAW - IRRI IMPROVED GEOGRAPHY
IDs$THR_RAW_RDGI<-NA
###
IDs$DIST_RAW_IR<-NA #DIST RAW -  HIJMANS
IDs$THR_RAW_IR<-NA
###
IDs$DIST_RDGI_IR<-NA #DIST IRRI IMPROVED GEOGRAPHY -  HIJMANS
IDs$THR_RDGI_IR<-NA




for(i in 1:nrow(IDs)){
  cat(as.character(IDs$ACID[[i]])," | ",i,"\n")
  
  IDs$lon_RAW[[i]]<-RAW$LON[[i]]
  IDs$lat_RAW[[i]]<-RAW$LAT[[i]]
  ###
  IDs$lon_RDGI[[i]]<-RDGI$LON[[i]]
  IDs$lat_RDGI[[i]]<-RDGI$LAT[[i]]
  ###
  IDs$lon_IR[[i]]<-IR$LON[[i]]
  IDs$lat_IR[[i]]<-IR$LAT[[i]]
  
  

IDs$DIST_RAW_RDGI[[i]]<-distm(c(RAW$LON[[i]],RAW$LAT[[i]]),c(RDGI$LON[[i]],RDGI$LAT[[i]]),distHaversine)/1000

if(is.na(IDs$DIST_RAW_RDGI[[i]])){
  IDs$THR_RAW_RDGI[[i]]<-NA
}else if(IDs$DIST_RAW_RDGI[[i]]>10){
  IDs$THR_RAW_RDGI[[i]]<-0
}else if(IDs$DIST_RAW_RDGI[[i]]<=10){
  IDs$THR_RAW_RDGI[[i]]<-1
}

IDs$DIST_RAW_IR[[i]]<-distm(c(RAW$LON[[i]],RAW$LAT[[i]]),c(IR$LON[[i]],IR$LAT[[i]]),distHaversine)/1000

if(is.na(IDs$DIST_RAW_IR[[i]])){
  IDs$THR_RAW_IR[[i]]<-NA
}else if(IDs$DIST_RAW_IR[[i]]>10){
  IDs$THR_RAW_IR[[i]]<-0
}else if(IDs$DIST_RAW_IR[[i]]<=10){
  IDs$THR_RAW_IR[[i]]<-1
}

IDs$DIST_RDGI_IR[[i]]<-distm(c(RDGI$LON[[i]],RDGI$LAT[[i]]),c(IR$LON[[i]],IR$LAT[[i]]),distHaversine)/1000
  
if(is.na(IDs$DIST_RDGI_IR[[i]])){
  IDs$THR_RDGI_IR[[i]]<-NA
    }else if(IDs$DIST_RDGI_IR[[i]]>10){
  IDs$THR_RDGI_IR[[i]]<-0
    }else if(IDs$DIST_RDGI_IR[[i]]<=10){
  IDs$THR_RDGI_IR[[i]]<-1
  }

}



write.table(IDs,"D:/Dropbox/Dropbox/NPGS georeferencing project/DATASETS/acc_id_coords_dist.csv",sep="|",quote =F,row.names = F,na= "");gc()

#IDs<-read.csv("D:/Dropbox/Dropbox/NPGS georeferencing project/DATASETS/acc_id_coords_dist.csv",sep="|",na.strings = "")
