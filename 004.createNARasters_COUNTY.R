
stop("Warning: do not run the whole thing")

#crop <- "Vitis" #change according to the coded name
#crop <- "zea" #change according to the coded name
#crop<-"Fruits_nut_Preece"
#crop<-"Vegies_root_Kantar"
#crop<-"Fruits_nut_Preece"
#crop<-"Fruits_tempsmall_Hummer"
crop<-"Ornamental_ornament_Jordan"

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

#basic stuff - creating folders
biomod <- paste(crop_dir,"/biomod_modeling",sep=""); if (!file.exists(biomod)) {dir.create(biomod)}

figs <- paste(crop_dir,"/figures",sep=""); if (!file.exists(figs)) {dir.create(figs)}

msks <- paste(crop_dir,"/masks",sep=""); if (!file.exists(msks)) {dir.create(msks)}
rm(msks)

narea <- paste(biomod,"/native-areas",sep=""); if (!file.exists(narea)) {dir.create(narea)}

# polys <- paste(narea,"/polyshps",sep=""); if (!file.exists(polys)) {dir.create(polys)}

ascis <- paste(narea,"/asciigrids",sep=""); if (!file.exists(ascis)) {dir.create(ascis)}
#############################################################



















require(rgdal)
require(raster)
require(maptools)
gpclibPermit(); options(warn=-1)

source(paste(src.dir,"/000.zipWrite.R",sep=""))
source(paste(src.dir,"/000.zipRead.R",sep=""))
source(paste(src.dir,"/000.createChullBuffer.R",sep=""))

createNARaster_COUNTY <- function(spID,inDir) {
  
  cat("\n")
  cat("Taxon", spID ,"\n")
  
  inNADir <- paste(inDir, "/NAREAS/SHP", sep="")
  outNADir <- paste(inDir, "/NAREAS/ASCII", sep="")
  outFolder <- paste(outNADir, "/", spID, sep="")
  allOcc <- read.csv(paste(crop_dir, "/occurrences/",crop,".csv", sep=""))
  
  if (!file.exists(paste(outFolder, "/narea.asc.gz", sep=""))) {
    
    cat("Not processed, thus processing \n")
    
    if (!file.exists(outNADir)) {
      dir.create(outNADir)
    }
    
    if (!file.exists(outFolder)) {
      dir.create(outFolder)
    }
    
    shpName <- paste(inNADir, "/", spID, "/narea.shp", sep="")
    
    if(!file.exists(shpName)){
      #Creating convex hull instead
      cat("No shapefile available, thus creating convex hull \n")
      occFile <- paste(crop_dir, "/occurrence_files/", spID, ".csv", sep="")
      tallOcc <- allOcc[which(allOcc$Taxon == paste(spID)),]
      if(nrow(tallOcc)==0){
        cat("No coordinates for this species /n")
      }else{
        NAGrid <- chullBuffer(crop_dir, occFile, outFolder, 500000)
        zipWrite(NAGrid, outFolder, "narea.asc.gz")
      }
    }else{
      #Reading polygon shapefile and mask  
      cat("Reading and converting \n")
      pol <- readShapePoly(shpName)
      rs <- raster(paste(crop_dir, "/masks/mask.asc", sep=""))
      
      #pa <- polygonsToRaster(pol, rs)
      pa <- rasterize(pol,rs)
      
      pa[which(!is.na(pa[]))] <- 1
      pa[which(is.na(pa[]) & rs[] == 1)] <- 0
      
      cat("Writing output \n")
      pa <- zipWrite(pa, outFolder, "narea.asc.gz")
      return(pa)
    }
    
  } else {
    cat("Already processed \n")
    #pa <- zipRead(outFolder, "narea.asc.gz")
  }
}



inDir <- crop_dir
#spList <- list.files(paste(inDir,"/","biomod_modeling", "/native-areas/polyshps", sep=""))
spList <- list.dirs(paste(inDir, "/NAREAS/SHP", sep=""),full.names = F, recursive = F)

#spList <- read.csv(paste(crop_dir,"/sample_counts/sample_count_table.csv",sep=""))
#spList <- spList$TAXON

for (spp in spList) {
  cat("Processing", spp, "\n")
  ot <- createNARaster_COUNTY(spp, inDir)
}

#source(paste(src.dir,"/004.createNARasters.R",sep=""))



