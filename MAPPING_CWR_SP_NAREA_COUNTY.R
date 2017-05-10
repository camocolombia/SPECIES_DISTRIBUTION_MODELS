
MAPPING_CWR_FUNCTION_NAREA<-function(spID){
  ###############################
  #                             #
  #     LOADING LIBRARIES       #             
  #                             #
  ###############################
  
  
  require(rgeos);
  require(rgdal);
  require(ffbase);
  require(sp);
  require(ggplot2);
  require(ggmap);
  library(RgoogleMaps)
  require(utils)
  require(colorRamps)
  require(dismo)
  require(rasterVis)
  require(ggthemes)
  require(maptools)
  require(plyr)
  ###############################
  
  source(paste(src.dir,"/000.zipRead.R",sep=""))
  
  ###############################
  
  out_dir<-paste0(crop_dir,"/","figures","/","SDMs");if (!file.exists(out_dir)) {dir.create(out_dir)}
  
  ###############################
  
  spID<-sub(".csv","",spID)
  
  tax_rich_table<-read.csv(paste0(crop_dir,"/","summary-files","/","taxaForRichness_NAREA.csv"),header = T)
  tax_rich_table<-tax_rich_table[which(!is.na(tax_rich_table$IS_VALID)),]
  
  
Tax_status<-tax_rich_table[which(tax_rich_table$TAXON==spID),]

if(nrow(Tax_status)>0){

  ###############################
  #                             #
  #     LOADING SDM DATA        #             
  #                             #
  ###############################

if(Tax_status$IS_VALID==1){

  spFolder <- paste(crop_dir, "/maxent_modeling/models/", spID, sep="")
  
  projFolder <- paste(spFolder, "/projections", sep="")
  projFolder2 <- paste0(projFolder,"/","NAREA_APPROACH")
  
  ENM_PA <- paste(spID, "_worldclim2_5_EMN_PA.asc.gz", sep="")
  
  dumm <- zipRead(projFolder2, ENM_PA)
  dumm[which(dumm[]<0)]<-NA
  dumm[which(dumm[]<1)]<-NA
}else{
  
  spFolder <- paste(crop_dir, "/samples_calculations/", spID, sep="")
  projFolder <- spFolder
  projFolder2 <- spFolder
  
  
  ENM_PA <- paste(spID, "_samples-buffer-na.asc.gz", sep="")
  dumm <- zipRead(projFolder2, ENM_PA)
  dumm[which(dumm[]<0)]<-NA
  dumm[which(dumm[]<1)]<-NA
  
}

  ###############################
  #                             #
  # LOADING NAREA DATA EXTENT   #             
  #                             #
  ###############################
  
  
if(file.exists(paste0(crop_dir,"/","NAREAS","/","SHP","/",spID,"/","narea.shp"))){
  
  narea_folder<-paste0(crop_dir,"/","NAREAS","/","SHP")
  
}else{
  narea_folder<-paste0(crop_dir,"/","biomod_modeling/native-areas/polyshps")
  
}

  
if(file.exists(paste0(narea_folder,"/",spID,"/","narea.shp"))){
cat("Loading NAREA Shapefile for ",as.character(spID),"\n")
sp_NA<-extent(shapefile(paste0(narea_folder,"/",spID,"/","narea.shp")))

}else{
  cat("Using Convex Hull as NAREA for ",as.character(spID),"\n")
  
sp_NA<-ENM_PA;  sp_NA <- zipRead(projFolder2, ENM_PA)
sp_NA[which(sp_NA[]<0)]<-NA
sp_NA[which(sp_NA[]<1)]<-1
sp_NA_NaM <- is.na(as.matrix(sp_NA))
colNotNA <- which(colSums(sp_NA_NaM) != nrow(sp_NA))
rowNotNA <- which(rowSums(sp_NA_NaM) != ncol(sp_NA))

r3Extent <- extent(sp_NA, rowNotNA[1], rowNotNA[length(rowNotNA)],
                   colNotNA[1], colNotNA[length(colNotNA)])

sp_NA <- extent(crop(sp_NA, r3Extent))

}

  ###############################
  #                             #
  #  LOADING LON LAT IN A FF    #             
  #                             #
  ###############################

occ_Narea_Dir<- paste(crop_dir, "/occurrence_files_narea",sep="")
occ_NA_file<-read.csv(paste0(occ_Narea_Dir,"/",spID,".csv_NAREA.csv"),header=T)
  
 
  ######################################
  
  #                                    #
  # RANGING COORDINATES TO GET MAP     #             
  #                                    #
  ######################################
  

xmin<-round(sp_NA[1]) ;if(xmin< -180){xmin=-180}
xmax<-round(sp_NA[2]) 
ymin<-round(sp_NA[3]) 
ymax<-round(sp_NA[4]) ;if(ymax>90){ymax=90}



  ###############################
  #                             #
  #     TESTING RANKING         #             
  #                             #
  ###############################
  
taxon3<-spID

  
  cat("Testing ranking such as:  subsp. ,var., f. for ", as.character(taxon3),"\n")
  
  taxon3<-sub("_"," ",taxon3)
  taxon3<-sub("_subsp._"," subsp. ",taxon3)
  taxon3<-sub("_var._"," var. ",taxon3)
  taxon3<-sub("_f._","f. ",taxon3)
  taxon3<-sub("_nothosubsp._"," nothosubsp. ",taxon3)
  split_text <- unlist(strsplit(taxon3, split=' '))
  cond <- c('subsp.', 'var.', 'f.')
  
  selected_text <- split_text[split_text %in% cond]
  
  
  cat("Savig file for ", as.character(taxon3),"\n")
  
  
  if(length(selected_text)==0){
    
    NAME<-substitute(expr=italic(taxon3),env=list(taxon3=taxon3))
  }else{
    a<-split_text[[1]]
    b<-split_text[[2]]
    c<-split_text[[4]]
    
    NAME<-substitute(expr=paste(italic(a)," ",italic(b)," ",selected_text," ",italic(c)),
                     env=list(a=a,b=b,selected_text=selected_text,c=c))
  }

  ###############################
  #                             #
  #     CONVERTING SDMs         #             
  #                             #
  ###############################

  map.p <- rasterToPoints(dumm)
  df <- data.frame(map.p)
  colnames(df) <- c("Longitude", "Latitude", "MAP")
  
  cat("Loading  geographical information","\n")
  
  usa <- map_data("world")
  usa <- subset(usa, region %in% c("Canada", "Mexico", "USA"))

  ###############################
  #                             #
  #     MAKING MAPS             #             
  #                             #
  ###############################

  cat("Saving map for ",as.character(spID),"\n")
  
 map2<-ggplot() + 
    geom_polygon(data=usa,aes(x=long,y=lat,group=group),fill="gray70",na.rm = T, show.legend = F) +
      geom_raster(data=df,aes(y=Latitude, x=Longitude,fill=MAP),show.legend=F,interpolate=T,na.rm=T) +
   geom_polygon(data=usa,aes(x=long,y=lat,group=group),colour="gainsboro",fill=NA,na.rm = T, show.legend = F, size = 0.5) +
   coord_fixed(1.3)+
       geom_path(color="white") +
    coord_equal() +
     geom_point(aes(x=lon, y=lat,col="Taxon"),size=0.8,stroke = 1,data=occ_NA_file, show.legend = TRUE,shape=17)+
     coord_equal() + 
          ggtitle(NAME)+
   labs(x ="Longitude", y="Latitude")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10),
          legend.position="bottom",
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey40"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.title=element_blank())+
    scale_colour_manual(name="Taxon occurrences",values="black")+
    guides(fill=FALSE)
   #coord_map("albers",  at0 = 45.5, lat1 = 29.5)
 
    
 map2<-map2+coord_fixed(xlim = c(xmin,xmax),  ylim = c(ymin,ymax), ratio = 1.3)
 
ggsave(paste0(out_dir,"/",as.character(spID),"_",Sys.Date(),"_NAREA.pdf"),units="in",width=4,height=4.8,scale=2,dpi=600)
 
}else{
      cat("skipping","\n")
      }
    }
  
##############################################################

spList <- list.files(paste(crop_dir, "/occurrence_files", sep=""))

lapply(1:length(spList),function(i){
  spID<-spList[[i]]
  x<-MAPPING_CWR_FUNCTION_NAREA(spID)
  })
