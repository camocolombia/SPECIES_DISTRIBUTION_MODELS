


stop("Warning: do not run the whole thing")

#crop <- "Vitis" #change according to the coded name
#crop <- "zea" #change according to the coded name
#crop<-"Fruits_nut_Preece"
#crop<-"Vegies_root_Kantar"
#crop<-"Fruits_nut_Preece"
#crop<-"Fruits_tempsmall_Hummer"
crop<-"Vitis"

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

#== prepare taxonomic names for the analysis ==#
spList <- list.dirs(paste(crop_dir, "/maxent_modeling/models", sep=""),recursive = F)
spList<-sub(paste(crop_dir, "/maxent_modeling/models/", sep=""),"",spList)


#source(paste(src.dir,"/Ensemble_test_model_COUNTY.R",sep=""))



x <- OccFilterNArea_COUNTY(crop_dir)
