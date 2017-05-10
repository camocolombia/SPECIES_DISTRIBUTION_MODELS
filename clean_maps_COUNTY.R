

CLEAN_MAPS_FUNCTION<-function(bdir,env_dir){
  require(rgdal)
  require(raster)
  library(sp)
source(paste0(src.dir,"/","000.zipRead.R"))
source(paste0(src.dir,"/","000.zipWrite.R"))
  env_dir <- "/mnt/workspace_cluster_6/ccsosa/GAP_ANALYSIS_US_2016/gap_analysis/bio_2_5m" # !!! change accordingly !!!
  geo_data_dir<-"/mnt/workspace_cluster_6/ccsosa/GAP_ANALYSIS_US_2016/gap_analysis/geodata"
  
max_dir<-paste0(bdir,"/","maxent_modeling","/","models")
sample_dir<-paste0(bdir,"/","samples_calculations")


cat("LOADING GLOBCOVER LAYER","\n")
LAND_COVER<-raster(paste0(env_dir,"/","GLOBCOVER","/","globcover.asc"))
######################################call the land use data GIS

m<-c(0,180,1,185,290,0)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
LAND_COVER<-reclassify(LAND_COVER,rclmat)


cat("LOADING TaxaForRichness TABLE","\n")

tax_rich_table<-read.csv(paste0(bdir,"/","summary-files","/","taxaForRichness_NAREA.csv"),header = T)
tax_rich_table<-tax_rich_table[which(!is.na(tax_rich_table$IS_VALID)),]


for(i in 1:length(tax_rich_table$TAXON)){
  cat("CLEANING SDM of ",as.character(tax_rich_table$TAXON[[i]]),"\n")
    if(tax_rich_table$IS_VALID[[i]]==1){
  mod_dir<-paste0(max_dir,"/",tax_rich_table$TAXON[[i]],"/","projections")
  mod_dir2<-paste0(mod_dir,"/","NAREA_APPROACH")
  
  x<-zipRead(mod_dir2,paste0(tax_rich_table$TAXON[[i]],"_worldclim2_5_EMN_PA.asc.gz"))
  y<-x*LAND_COVER   
  zipWrite(y,mod_dir2,paste0(tax_rich_table$TAXON[[i]],"_worldclim2_5_EMN_PA.asc.gz"))
  }else{
    cat("CLEANING Sample Buffer of ",as.character(tax_rich_table$TAXON[[i]]),"\n")
    
  mod_dir<-paste0(sample_dir,"/",tax_rich_table$TAXON[[i]])
x<-zipRead(paste0(sample_dir,"/",tax_rich_table$TAXON[[i]]),"samples-buffer-na.asc.gz")
y<-x*LAND_COVER   
zipWrite(y,mod_dir2,paste0(tax_rich_table$TAXON[[i]],"_samples-buffer-na.asc.gz"))
    }
  }
gc();
}

x<-CLEAN_MAPS_FUNCTION(bdir=crop_dir,env_dir = env_dir)

