require(shapefiles);require(psych);library(geosphere);require(raster)

#####################################################
DATASETS_PATH<- "D:/Dropbox/Dropbox/NPGS georeferencing project/DATASETS/DATASETS_GEOREF"

IDs<-read.csv("D:/Dropbox/Dropbox/NPGS georeferencing project/DATASETS/acc_id_coords_dist.csv",sep="|",na.strings = "");gc()

adm<-shapefile("D:/ADMIN/GAUL_2014/gisdata/G2014_2013_2.shp");gc()
gc()
#####################################################

RAW<-as.data.frame(cbind(IDs$ACID,IDs$lon_RAW,IDs$lat_RAW));colnames(RAW)<-c("ACID","LON","LAT");gc()## RAW DATA FROM GRIN GLOBAL
RDGI<-as.data.frame(cbind(IDs$ACID,IDs$lon_RDGI,IDs$lat_RDGI));colnames(RDGI)<-c("ACID","LON","LAT");gc()## IRRI IMPROVED GEOGRAPHY
IR<-as.data.frame(cbind(IDs$ACID,IDs$lon_IR,IDs$lat_IR));colnames(IR)<-c("ACID","LON","LAT");gc() ## HIJMANS RESULTS

RAW<-RAW[complete.cases(RAW),];gc()
RDGI<-RAW[complete.cases(RDGI),];gc()
IR<-RAW[complete.cases(IR),];gc()

RAW<-RAW[which(!is.na(RAW$LON)),];gc()
RDGI<-RDGI[which(!is.na(RDGI$LON)),];gc()
IR<-IR[which(!is.na(IR$LON)),];gc()



coordinates(RAW) <- ~LON+LAT;gc()
coordinates(RDGI)<- ~LON+LAT;gc()
coordinates(IR)  <- ~LON+LAT;gc()


crs(RAW) <- crs(adm);crs(RDGI) <- crs(adm);crs(IR) <- crs(adm)

ovr_RAW <- over(RAW, adm);gc()
ovr_RDGI <- over(RDGI, adm);gc()
ovr_IR <- over(IR, adm);gc()

ovr_RAW<-cbind(RAW$ACID,ovr_RAW);gc()
ovr_RDGI<-cbind(RDGI$ACID,ovr_RDGI);gc()
ovr_IR<-cbind(IR$ACID,ovr_IR);gc()

colnames(ovr_RAW)[[1]]<-"ACID";colnames(ovr_RDGI)[[1]]<-"ACID";colnames(ovr_IR)[[1]]<-"ACID"


o_R<-merge(IDs,ovr_RAW,by.x="ACID", by.y="ACID", incomparables = NA,all.x=T)
o_RD<-merge(IDs,ovr_RDGI,by.x="ACID",by.y="ACID",  incomparables = NA,all.x=T)
o_I<-merge(IDs,ovr_IR,by.x="ACID",by.y="ACID",  incomparables = NA,all.x=T)

IDs$RAW_ADM0_CODE<-NA;IDs$RAW_ADM0_CODE<-o_R$ADM0_CODE#;rm(o_R)
IDs$RDGI_ADM0_CODE<-NA;IDs$RDGI_ADM0_CODE<-o_RD$ADM0_CODE#;rm(o_RD)
IDs$IR_ADM0_CODE<-NA;IDs$IR_ADM0_CODE<-o_I$ADM0_CODE#;rm(o_I)


IDs$RAW_ADM0_NAME<-NA;IDs$RAW_ADM0_NAME<-o_R$ADM0_NAME#;rm(o_R)
IDs$RDGI_ADM0_NAME<-NA;IDs$RDGI_ADM0_NAME<-o_RD$ADM0_NAME#;rm(o_RD)
IDs$IR_ADM0_NAME<-NA;IDs$IR_ADM0_NAME<-o_I$ADM0_NAME#;rm(o_I)




write.table(IDs,"D:/Dropbox/Dropbox/NPGS georeferencing project/DATASETS/acc_id_coords_dist_coun_COUN.csv",sep="|",na="",quote = F,row.names = F);gc()


