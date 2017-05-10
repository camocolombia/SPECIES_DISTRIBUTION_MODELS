

SP_RICHNESS_OPT<-function(crop,spList,out_name){
require(raster)



src.dir <- paste("X:/GAP_ANALYSIS_US_2016/gap_analysis/gap_",crop,"/_scripts",sep="") # !!! change accordingly !!!

#gap.dir <-"/curie_data/ncastaneda/gap-analysis" # !!! change accordingly !!!

gap.dir <-"X:/GAP_ANALYSIS_US_2016/gap_analysis" # !!! change accordingly !!!
#gap.dir <-"/mnt/workspace_cluster_6/ccsosa/GAP_ANALYSIS_US_2016/gap_analysis" # !!! change accordingly !!!

#crop details
crop_dir <- paste(gap.dir,"/gap_",crop,sep="")

if (!file.exists(crop_dir)) {dir.create(crop_dir)}
setwd(crop_dir)

source(paste0(src.dir,"/","000.zipRead.R"))
source(paste0(src.dir,"/","000.zipWrite.R"))


models_dir<-paste0(crop_dir,"/","maxent_modeling/models")


#path, fname
#Aphanisma_blitoides_worldclim2_5_EMN_PA.asc.gz
x_models<-lapply(1:length(spList),function(i){
  
     spID<-spList[[i]]
    # cat("Processing ",spID,"\n")
     
     mod_sp_dir<-paste0(models_dir,"/",spID,"/projections/NAREA_APPROACH")
     sp_file<-paste0(spID,"_worldclim2_5_EMN_PA.asc.gz")
     sp_file_p<-paste0(mod_sp_dir,"/",sp_file)
     
     if(!file.exists(sp_file_p)){
       
       cat("Skipping ",spID,"\n") 
     }else{
     
     cat("Loading raster for ",spID,"\n")
     
mod_sp_dir<-paste0(models_dir,"/",spID,"/projections/NAREA_APPROACH")
sp_file<-paste0(spID,"_worldclim2_5_EMN_PA.asc.gz")
  mod<-zipRead(mod_sp_dir,sp_file)
  return(mod)
}
})
gc()
x_models<-x_models[lapply(x_models,length)!=0];gc()
cat("Creating sum raster for ",out_name,"\n")
x_models_sum<-sum(stack(x_models),na.rm=T);gc()
cat("Compressing raster ",out_name,"\n")
X<-zipWrite(x_models_sum,crop_dir,paste0(out_name,".asc.gz"));gc()
}

crop<-"bean"

spList<-c(
  'Phaseolus_acutifolius',
  'Phaseolus_acutifolius_subsp._acutifolius',
  'Phaseolus_acutifolius_subsp._tenuifolius',
  'Phaseolus_albescens',
  'Phaseolus_albiflorus',
  'Phaseolus_albinervus',
  'Phaseolus_altimontanus',
  'Phaseolus_amabilis',
  'Phaseolus_amblyosepalus',
  'Phaseolus_angustissimus',
  'Phaseolus_anisophyllus',
  'Phaseolus_augusti',
  'Phaseolus_campanulatus',
  'Phaseolus_carteri',
  'Phaseolus_chiapasanus',
  'Phaseolus_coccineus',
  'Phaseolus_costaricensis',
  'Phaseolus_dasycarpus',
  'Phaseolus_dumosus',
  'Phaseolus_esperanzae',
  'Phaseolus_esquincensis',
  'Phaseolus_filiformis',
  'Phaseolus_glabellus',
  'Phaseolus_gladiolatus',
  'Phaseolus_grayanus',
  'Phaseolus_hintonii',
  'Phaseolus_jaliscanus',
  'Phaseolus_laxiflorus',
  'Phaseolus_leptostachyus',
  'Phaseolus_lunatus',
  'Phaseolus_macrolepis',
  'Phaseolus_maculatifolius',
  'Phaseolus_maculatus',
  'Phaseolus_maculatus_subsp._maculatus',
  'Phaseolus_maculatus_subsp._ritensis',
  'Phaseolus_macvaughii',
  'Phaseolus_magnilobatus',
  'Phaseolus_marechalii',
  'Phaseolus_micranthus',
  'Phaseolus_microcarpus',
  'Phaseolus_mollis',
  'Phaseolus_neglectus',
  'Phaseolus_nelsonii',
  'Phaseolus_nodosus',
  'Phaseolus_novoleonensis',
  'Phaseolus_oaxacanus',
  'Phaseolus_oligospermus',
  'Phaseolus_opacus',
  'Phaseolus_pachyrrhizoides',
  'Phaseolus_parvifolius',
  'Phaseolus_parvulus',
  'Phaseolus_pauciflorus',
  'Phaseolus_pedicellatus',
  'Phaseolus_perplexus',
  'Phaseolus_plagiocylix',
  'Phaseolus_pluriflorus',
  'Phaseolus_polymorphus',
  'Phaseolus_polystachyus',
  'Phaseolus_purpusii',
  'Phaseolus_reticulatus',
  'Phaseolus_rotundatus',
  'Phaseolus_salicifolius',
  'Phaseolus_sinuatus',
  'Phaseolus_smilacifolius',
  'Phaseolus_sonorensis',
  'Phaseolus_talamancensis',
  'Phaseolus_tenellus',
  'Phaseolus_texensis',
  'Phaseolus_trifidus',
  'Phaseolus_tuerckheimii',
  'Phaseolus_venosus',
  'Phaseolus_vulgaris',
  'Phaseolus_xanthotrichus',
  'Phaseolus_xolocotzii',
  'Phaseolus_zimapanensis'
  )


out_name<-"Phaseolus_sp_richness"



x<-SP_RICHNESS_OPT(crop,spList,out_name)

