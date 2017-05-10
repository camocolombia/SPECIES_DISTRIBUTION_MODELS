require(raster);require(dismo)
library("sp")
library("rgdal")
rs_dir<-list.files("D:/BU/AAA",pattern = ".asc$",full.names = T)
rs_dir_names<-list.files("D:/BU/AAA",pattern = ".asc$",full.names = F)
rs_dir_names<-sub("_worldclim2_5_EMN_PA.asc","",rs_dir_names)
out_dir<-"D:/BU/AAA"

lapply(1:length(rs_dir),function(i){
  cat("processing ",rs_dir_names[[i]],"\n")
  rs1<-raster(rs_dir[[i]])
  rs1[which(rs1[]==0)]<-NA
  
  rs2<-rs1
  
  rs2[]<-1:ncell(rs1)
  
  rs2[which(is.na(rs1[]))]<-NA
  rs3<-rs2[!is.na(rs2[])]
  
    points<-xyFromCell(rs1,rs3)
  
  #ch <- convHull(points)
  #ch2<-predict(ch, points)
 z<-chull(points[,1],points[,2])
 # dfHull <-cbind(points[,1][z],points[,2][z])
  
  points2<-cbind(points[,1],points[,2])
  coords <- points2[c(z, z[1]), ]  # closed polygon
  # 
  # 
  # plot(points2, pch=19)
  # lines(coords, col="red")
  
  
  sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)))
  # set coordinate reference system with SpatialPolygons(..., proj4string=CRS(...))
  # e.g. CRS("+proj=longlat +datum=WGS84")
  sp_poly_df <- SpatialPolygonsDataFrame(sp_poly, data=data.frame(ID=1))
  
  setwd(out_dir)
  cat("Saving shapefile for ",rs_dir_names[[i]],"\n")
  writeOGR(obj=sp_poly_df, dsn=rs_dir_names[[i]], layer=rs_dir_names[[i]], driver="ESRI Shapefile")
  
})
