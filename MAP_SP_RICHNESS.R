MAPPING_CWR_FUNCTION_SP_RICH<-function(sp_rich_file,crop){
  ###############################
  #                             #
  #     LOADING LIBRARIES       #             
  #                             #
  ###############################
  
  require(raster);
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
  crop<-crop
  src.dir <- paste("X:/GAP_ANALYSIS_US_2016/gap_analysis/gap_",crop,"/_scripts",sep="") # !!! change accordingly !!!
  
  #gap.dir <-"/curie_data/ncastaneda/gap-analysis" # !!! change accordingly !!!
  
  gap.dir <-"X:/GAP_ANALYSIS_US_2016/gap_analysis" # !!! change accordingly !!!
  #gap.dir <-"/mnt/workspace_cluster_6/ccsosa/GAP_ANALYSIS_US_2016/gap_analysis" # !!! change accordingly !!!
  
  #crop details
  crop_dir <- paste(gap.dir,"/gap_",crop,sep="")
  
  if (!file.exists(crop_dir)) {dir.create(crop_dir)}
  setwd(crop_dir)
  
  ###############################
  source(paste(src.dir,"/000.zipRead.R",sep=""))
  
  ###############################

  
  sp_rich_file<-sp_rich_file
  out_dir<-crop_dir
  
  ###############################
 
    ###############################
    #                             #
    #     LOADING SDM DATA        #             
    #                             #
    ###############################
    
       projFolder2 <- crop_dir
#       c
      ENM_PA <- sp_rich_file
      
      dumm <- zipRead(projFolder2, ENM_PA)
      dumm[which(dumm[]<0)]<-NA
      dumm[which(dumm[]<1)]<-NA
    
    
    ###############################
    #                             #
    # LOADING NAREA DATA EXTENT   #             
    #                             #
    ###############################
    
   
    
      sp_NA <- extent(dumm)
      
    
    
  
    
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
    
    taxon3<-sp_rich_file
    taxon3<-sub("_sp_richness.asc.gz","",taxon3)
    
    cat("Testing spaces for ", as.character(taxon3),"\n")
    
    taxon3<-sub("_"," ",taxon3)
#     taxon3<-sub("_subsp._"," subsp. ",taxon3)
#     taxon3<-sub("_var._"," var. ",taxon3)
#     taxon3<-sub("_f._","f. ",taxon3)
#     taxon3<-sub("_nothosubsp._"," nothosubsp. ",taxon3)
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
    
    cat("Saving map for ",as.character(taxon3),"\n")
    
    map2<-ggplot() + 
      geom_polygon(data=usa,aes(x=long,y=lat,group=group),fill="gray70",na.rm = T, show.legend = F) +
      geom_raster(data=df,aes(y=Latitude, x=Longitude,fill=MAP),show.legend=T,interpolate=T,na.rm=T) +
      scale_fill_gradientn("# of taxa",colours=c("#EDB01D","#77ED1D","#1DB3ED","#391DED","#EC1DED","#ED1D3C"))+
      geom_polygon(data=usa,aes(x=long,y=lat,group=group),colour="gainsboro",fill=NA,na.rm = T, show.legend = F, size = 0.5) +
      coord_fixed(1.3)+
      geom_path(color="white") +
      coord_equal() +
     # geom_point(aes(x=lon, y=lat,col="Taxon"),size=0.3,stroke = 1,data=occ_NA_file, show.legend = TRUE,shape=17)+
     # coord_equal() + 
      ggtitle(NAME)+
      labs(x ="Longitude", y="Latitude")+
  #scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF"))+
       theme(axis.text=element_text(size=10),
            axis.title=element_text(size=10),
            legend.position="bottom",
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(colour = "grey40"),
           # panel.border = element_rect(colour = "black", fill=NA, size=1),
           panel.border = element_rect(colour = "black", fill=NA, size=1),
           legend.title=element_text(colour="black", size=16, face="bold"))
#+
           
      #scale_colour_manual(name="Taxon occurrences",values="black")+
      #guides(fill=T)
    #coord_map("albers",  at0 = 45.5, lat1 = 29.5)
    
    
    map2<-map2+coord_fixed(xlim = c(xmin,xmax),  ylim = c(ymin,ymax), ratio = 1.3)
    map2<-map2

    ggsave(paste0(out_dir,"/",as.character(taxon3),"_",Sys.Date(),"_NAREA.pdf"),units="in",width=4,height=4.8,scale=2,dpi=600)
    gc()
  }
  


##############################################################

crop<-"Industrial_sunflower_Marek"
sp_rich_files<-c(
  "Helianthus_sp_richness.asc.gz"
#   "Ipomoea_sp_richness.asc.gz",
#   "Manihot_sp_richness.asc.gz",
#   "Solanum_sp_richness.asc.gz",
#   "Xanthosoma_sp_richness.asc.gz"
)



lapply(1:length(sp_rich_files),function(i){
  sp_rich_file<-sp_rich_files[[i]]
  x<-MAPPING_CWR_FUNCTION_SP_RICH(sp_rich_file,crop)
})
