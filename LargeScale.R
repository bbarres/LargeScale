###############################################################################
###############################################################################
#LargeScale data analysis
###############################################################################
###############################################################################

#loading the packages needed
library(maptools)
library(rgdal)
library(mapdata)
library(mapplots)
library(adegenet)
library(plyr)
library(adegenet)
library(gdata)
library(RColorBrewer)

#loading the geographic data, don't forget to set the right working directory
setwd("~/work/Rfichiers/Githuber/LargeScale_data")
World<-readShapePoly("CNTR_RG_03M_2010.shp",proj4string=CRS("+init=epsg:4326"))

#loading the population information
patch_info<-read.table("larscapatch.dat",header=TRUE,sep="\t",dec=".")
#turn the information in a spatialpoints object
patch_coord<-patch_info
coordinates(patch_coord)<-~longitude+latitude
#We then define the coordinates system used
proj4string(patch_coord)<-CRS("+init=epsg:4326")

plot(World,col="black")
plot(World[World$CNTR_ID,],ylim=c(43,63),xlim=c(10,20), 
     col="grey",border="grey")

#map of the surveyed samples
plot(World[World$CNTR_ID,],ylim=c(43,63),xlim=c(10,20), 
     col="grey",border="black")
#we add the intensity of the sampling in infeted patches
points(patch_coord,cex=sqrt(as.numeric(patch_coord@data$nb_sample)),
       bg=rgb(1,0,0,0.5),pch=21,col=rgb(1,0,0,1))
#then we add the patch without infection
points(patch_coord[as.numeric(patch_coord@data$PA)==0,],col="green3",pch=19)

#map of the relative proportion of coinfection in the different populations
plot(World[World$CNTR_ID,],ylim=c(43,63),xlim=c(10,20), 
     col="grey",border="black")
stars(cbind(patch_coord@data$nb_coinf,patch_coord@data$nb_pure),draw.segment=TRUE,
      locations=coordinates(patch_coord),add=T,len=2,radius=FALSE,full=FALSE,
      col.segments=c(6,4))


#first of all, we load the genetic dataset
larscagen<-read.table("larscagen.dat",header=T,sep="\t",dec=".")

#a summary of the different variables
summary(larscagen)
colnames(larscagen)
#number of individuals in each sampled populations
table(larscagen$patch_ID)
#total number of individuals
sum(table(larscagen$patch_ID)) #267 individuals

#here we select only a part of the samples for further analysis
#we remove the individuals with multiple allele (ie coinfection)
#and we remove the second column of the 
JDD<-MyzPeach #name of the input file
JDD<-drop.levels(JDD)
#let's define a set of color for keeping some consistency in the plots
coloor<-c("red","green","blue","yellow","orchid")


###############################################################################
###############################################################################
#DAPC on microsatellites only
###############################################################################
###############################################################################

#converting data to a genind format, first we use only the microsatellite data
JDDmicro<-df2genind(JDD[,c("MP_27","MP_39","MP_44","MP_5","MP_7","MP_23",
                           "MP_45","MP_28","MP_9","MP_13","MP_2","MP_38",
                           "MP_4","MP_46")],
                    ncode=6,ind.names=JDD$sample_ID, 
                    pop=JDD$patch_ID,missing=NA,ploidy=2)
#include the coordinates of the samples
JDDmicro@other$xy<-JDD[,c("longitude","latitude")]

#now we analyse the adegenet format dataset with dapc
JDDade<-JDDmicro
#determination of the number of clusters
clustJDDade<- find.clusters(JDDade,max.n.clust=35)
#with 40 PCs, we lost nearly no information
clustJDDade<- find.clusters(JDDade,n.pca=40,max.n.clust=35) #chose 4 clusters
#which individuals in which clusters per population
table(pop(JDDade),clustJDDade$grp)
#DAPC by itself, first we try to optimized the number of principal component (PCs) 
#to retain to perform the analysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=100)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=30)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=4,n.pca=10)
#STRUCTURE-like graphic
compoplot(dapcJDDade,lab=truenames(JDDade)$pop,legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)
#scatter plot
scatter(dapcJDDade,xax=1, yax=2,col=coloor)
#a more beautifull scatter plot
scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.3,pch=19,cex=3,scree.da=TRUE)


###############################################################################
###############################################################################
#DAPC on SNP only
###############################################################################
###############################################################################

#converting data to a genind format, first we use only the microsatellite data
JDDsnp<-df2genind(JDD[,c("MP_27","MP_39","MP_44","MP_5","MP_7","MP_23",
                           "MP_45","MP_28","MP_9","MP_13","MP_2","MP_38",
                           "MP_4","MP_46")],
                    ncode=6,ind.names=JDD$sample_ID, 
                    pop=JDD$patch_ID,missing=NA,ploidy=2)
#include the coordinates of the samples
JDDsnp@other$xy<-JDD[,c("longitude","latitude")]

#now we analyse the adegenet format dataset with dapc
JDDade<-JDDsnp
#determination of the number of clusters
clustJDDade<- find.clusters(JDDade,max.n.clust=35)
#with 40 PCs, we lost nearly no information
clustJDDade<- find.clusters(JDDade,n.pca=40,max.n.clust=35) #chose 4 clusters
#which individuals in which clusters per population
table(pop(JDDade),clustJDDade$grp)
#DAPC by itself, first we try to optimized the number of principal component (PCs) 
#to retain to perform the analysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=100)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=30)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=4,n.pca=10)
#STRUCTURE-like graphic
compoplot(dapcJDDade,lab=truenames(JDDade)$pop,legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)
#scatter plot
scatter(dapcJDDade,xax=1, yax=2,col=coloor)
#a more beautifull scatter plot
scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.3,pch=19,cex=3,scree.da=TRUE)


###############################################################################
###############################################################################
#DAPC on SNP and microsatellites
###############################################################################
###############################################################################

#converting data to a genind format, first we use only the microsatellite data
JDDall<-df2genind(JDD[,c("MP_27","MP_39","MP_44","MP_5","MP_7","MP_23",
                         "MP_45","MP_28","MP_9","MP_13","MP_2","MP_38",
                         "MP_4","MP_46","kdr","skdr","R81T","MACE")],
                  ncode=6,ind.names=JDD$sample_ID, 
                  pop=JDD$patch_ID,missing=NA,ploidy=2)
#include the coordinates of the samples
JDDall@other$xy<-JDD[,c("longitude","latitude")]

#now we analyse the adegenet format dataset with dapc
JDDade<-JDDall
#determination of the number of clusters
clustJDDade<- find.clusters(JDDade,max.n.clust=35)
#with 40 PCs, we lost nearly no information
clustJDDade<- find.clusters(JDDade,n.pca=40,max.n.clust=35) #chose 4 clusters
#which individuals in which clusters per population
table(pop(JDDade),clustJDDade$grp)
#DAPC by itself, first we try to optimized the number of principal component (PCs) 
#to retain to perform the analysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=100)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=30)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=4,n.pca=7)
#STRUCTURE-like graphic
compoplot(dapcJDDade,lab=truenames(JDDade)$pop,legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)
#scatter plot
scatter(dapcJDDade,xax=1, yax=2,col=coloor)
#a more beautifull scatter plot
scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.3,pch=19,cex=3,scree.da=TRUE)





###############################################################################
###############################################################################
#trash
###############################################################################
###############################################################################

#scatter plot with the different K groups and then plotting the population
scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.3,pch=19,cex=3,scree.da=FALSE)
#oilseed_rape
points(dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==1,1],
       dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==1,2],col="black",
       ,cex=2,bg="black",pch=21)
#peach
points(dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==3,1],
       dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==3,2],col="black",
       ,cex=2,bg="black",pch=21)
#tobacco
points(dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==4,1],
       dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==4,2],col="black",
       ,cex=2,bg="black",pch=21)
#other_crop
points(dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==2,1],
       dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==2,2],col="black",
       ,cex=2,bg="black",pch=21)


scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.0,pch=19,cex=3,scree.da=FALSE)
points(dapcJDDade$ind.coord[,1],dapcJDDade$ind.coord[,2],col=coloor[dapcJDDade$assign],
       pch=(as.numeric(as.factor(JDDade@other$host))+20),cex=2)

scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.0,pch=19,cex=3,scree.da=FALSE)
points(dapcJDDade$ind.coord[,1],dapcJDDade$ind.coord[,2],col=coloor[dapcJDDade$assign],
       pch=21,bg=coloor[(as.numeric(as.factor(JDDade@other$host)))])


plot(JDDade@other$xy,cex=3,col=dapcJDDade$assign,pch=as.numeric(as.factor(JDDade@other$host)))



# BRADEpop<-genind2genpop(BRADE,process.other=T,missing="0")
# 
# image(alt,col=brewer.pal(9,"Greys"))
# stars(table(pop(JDDade),dapcJDDade$assign),draw.segment=TRUE,
#       locations=JDDade@other$xy,
#       #locations=cbind(jitter(BRADEpop@other$xy$longitude,200),
#       #                jitter(BRADEpop@other$xy$latitude,200)),
#       add=T,len=0.5)







