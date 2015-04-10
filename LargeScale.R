###############################################################################
###############################################################################
#LargeScale data analysis
###############################################################################
###############################################################################

library(maptools)
library(rgdal)
library(mapdata)
library(mapplots)
library(adegenet)
library(plyr)

#loading the data, don't forget to set the right working directory
setwd("~/work/Rfichiers/Githuber/LargeScale_data")
World<-readShapePoly("CNTR_RG_03M_2010.shp",proj4string=CRS("+init=epsg:4326"))
sample_info<-read.table("sample_info3.txt", header=TRUE, sep="\t", dec=".")
patch_info<-read.table("patch_LargeScale.txt",header=TRUE,sep="\t")
patch_coord<-patch_info
coordinates(patch_coord)<-~Latitude+Longitude
#We then define the coordinates system used
proj4string(patch_coord)<-CRS("+init=epsg:4326")

plot(World,col="black")
plot(World[World$CNTR_ID,],ylim=c(43,63),xlim=c(10,20), col="grey",border="grey")



