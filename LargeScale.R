###############################################################################
###############################################################################
#LargeScale data analysis
###############################################################################
###############################################################################

#loading the packages needed
library(maptools)
library(rgdal)
library(scales)
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
stars(cbind(patch_coord@data$nb_coinf,patch_coord@data$nb_pure),
      draw.segment=TRUE,locations=coordinates(patch_coord),add=T,
      len=2,radius=FALSE,full=FALSE,col.segments=c(6,4))
#same map but with pie chart instead of the star-system
plot(World[World$CNTR_ID,],ylim=c(43,63),xlim=c(10,20), 
     col="grey",border="black")
draw.pie(x=patch_info$longitude,y=patch_info$latitude,
         z=cbind(patch_info$nb_pure,patch_info$nb_coinf),
         col=c(alpha("blue",0.4),alpha("red",0.4)),
         radius=(sqrt(patch_info$nb_sample)/8),labels=NA)
#another alternative map without border and smaller pie for zoom on Baltic Sea
plot(World[World$CNTR_ID,],ylim=c(43,63),xlim=c(10,20), 
     col="grey",border="transparent")
draw.pie(x=patch_info$longitude,y=patch_info$latitude,
         z=cbind(patch_info$nb_pure,patch_info$nb_coinf),
         col=c(alpha("blue",0.4),alpha("red",0.4)),
         radius=(sqrt(patch_info$nb_sample)/20),labels=NA)

#because there are still a lot of overlapping, we can merge samples by 
#geographic area
geoarea<-cbind(by(patch_info$longitude[patch_coord@data$PA==1],
                  patch_info$geo_area[patch_coord@data$PA==1],mean),
               by(patch_info$latitude[patch_coord@data$PA==1],
                  patch_info$geo_area[patch_coord@data$PA==1],mean),
               by(patch_info$nb_sample[patch_coord@data$PA==1],
                  patch_info$geo_area[patch_coord@data$PA==1],sum),
               by(patch_info$nb_coinf[patch_coord@data$PA==1],
                  patch_info$geo_area[patch_coord@data$PA==1],sum),
               by(patch_info$nb_pure[patch_coord@data$PA==1],
                  patch_info$geo_area[patch_coord@data$PA==1],sum))
colnames(geoarea)<-c("longitude","latitude","nb_sample","nb_coinf",
                     "nb_pure")
geoarea<-as.data.frame(geoarea)
plot(World[World$CNTR_ID,],ylim=c(43,63),xlim=c(10,20), 
     col="grey",border="black")
draw.pie(x=geoarea$longitude,y=geoarea$latitude,
         z=cbind(geoarea$nb_pure,geoarea$nb_coinf),
         col=c(alpha("blue",0.8),alpha("red",0.8)),
         radius=(sqrt(geoarea$nb_sample)/8),labels=NA)

#loading the genetic dataset
larscagen<-read.table("larscagen.dat",header=T,sep="\t",dec=".")
#recode the missing genotype data
larscagen[larscagen==-9]<-NA

#a summary of the different variables
summary(larscagen)
colnames(larscagen)
#number of individuals in each sampled populations
table(larscagen$patch_ID)
#total number of individuals
sum(table(larscagen$patch_ID)) #267 individuals

#here we select only a part of the samples for further analysis
#we remove the individuals with multiple allele (ie coinfection)
#and we remove the second allele of the markers because it is 
#unnecessary with no coinfection for an haploid species
JDD<-larscagen
#removing unnecessary allele2 columns
JDDhap<-subset(JDD,select=-c(uPP01,uPP05,uPP07,uPP09,uPP20,
                             h_20101214_c1217_640.5307,
                             h_20101214_c1336_788.5288,
                             h_20101214_c1421_219.5293,
                             h_20101214_c1421_455.5292,
                             h_20101214_c1720i_2036.5304,
                             h_20101214_c1892_2119.5284,
                             h_20101214_c2493_i_601.5301,
                             h_20101214_c2804_701.5309,
                             h_20101214_c3117_1457.5310,
                             h_20101214_c3926_348.5286,
                             h_20101214_c3997_i_419.5298,
                             h_20101214_c4769_1106.5305,
                             h_20101214_c5096_985.5294,
                             h_20101214_c5876_431.5299,
                             h_20101214_rep_c542_236.5287,
                             h_20101214_rep_c6068_457.5290,
                             h_20101214_rep_c664_2300.5297,
                             h_20101214_rep_c707_1118.5296,
                             h_20101214_rep_c707_1234.5303))
#removing coinfection
JDDsing<-JDDhap[JDDhap$coinf==0,]
JDDsing<-drop.levels(JDDsing)
#we also select a subset of the dataset without repeated MLG
JDDcc<-JDDsing[JDDsing$MLGrepeat==0,]
JDDcc<-drop.levels(JDDcc)

#let's define a set of color for keeping some consistency in the plots
coloor<-c("blue","red","green","violet","orange")


###############################################################################
###############################################################################
#DAPC on microsatellites only
###############################################################################
###############################################################################

#converting data to a genind format, first we use only the microsatellite data
JDDmicro<-df2genind(JDDcc[,c("uPP01_1","uPP05_1","uPP07_1","uPP09_1",
                             "uPP20_1")],
                    ncode=3,ind.names=JDDcc$ID, 
                    pop=JDDcc$geo_area,missing=NA,ploidy=1)
#include the populations names and the coordinates of the samples
JDDmicro@pop.names=c("Western Europe","Central Europe","Sweden","Finland",
                   "Estonia")
JDDmicro@other$xy<-JDDcc[,c("longitude","latitude")]

#now we analyse the adegenet format dataset with dapc
JDDade<-JDDmicro
#determination of the number of clusters
clustJDDade<- find.clusters(JDDade,max.n.clust=25)
#with 15 PCs, we lost nearly no information
clustJDDade<- find.clusters(JDDade,n.pca=15,max.n.clust=25) #chose 3/4 clusters
#which individuals in which clusters per population
table(pop(JDDade),clustJDDade$grp)
#DAPC , first we try to optimized the number of principal component (PCs) 
#to retain to perform the analysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=4,n.pca=40)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=3,n.pca=10)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=2,n.pca=3)
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
JDDsnp<-df2genind(JDDcc[,c("h_20101214_c1217_640.5307_1",
                           "h_20101214_c1336_788.5288_1",
                           "h_20101214_c1421_219.5293_1",
                           "h_20101214_c1421_455.5292_1",
                           "h_20101214_c1720i_2036.5304_1",
                           "h_20101214_c1892_2119.5284_1",
                           "h_20101214_c2493_i_601.5301_1",
                           "h_20101214_c2804_701.5309_1",
                           "h_20101214_c3117_1457.5310_1",
                           "h_20101214_c3926_348.5286_1",
                           "h_20101214_c3997_i_419.5298_1",
                           "h_20101214_c4769_1106.5305_1",
                           "h_20101214_c5096_985.5294_1",
                           "h_20101214_c5876_431.5299_1",
                           "h_20101214_rep_c542_236.5287_1",
                           "h_20101214_rep_c6068_457.5290_1",
                           "h_20101214_rep_c664_2300.5297_1",
                           "h_20101214_rep_c707_1118.5296_1",
                           "h_20101214_rep_c707_1234.5303_1")],
                    ncode=3,ind.names=JDDcc$ID, 
                    pop=JDDcc$geo_area,missing=NA,ploidy=1)
#include the population names and the coordinates of the samples
JDDsnp@pop.names=c("Western Europe","Central Europe","Sweden","Finland",
                   "Estonia")
JDDsnp@other$xy<-JDDcc[,c("longitude","latitude")]

#now we analyse the adegenet format dataset with dapc
JDDade<-JDDsnp
#determination of the number of clusters
clustJDDade<- find.clusters(JDDade,max.n.clust=25)
#with 40 PCs, we lost nearly no information
clustJDDade<- find.clusters(JDDade,n.pca=12,max.n.clust=25) #chose 4 clusters
#which individuals in which clusters per population
table(pop(JDDade),clustJDDade$grp)
#DAPC, first we try to optimized the number of principal component (PCs) 
#to retain to perform the analysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=30)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=10)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=2,n.pca=3)
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
#DAPC on both SNP and microsatellites markers
###############################################################################
###############################################################################

JDDall<-df2genind(JDDcc[,c("uPP01_1","uPP05_1","uPP07_1","uPP09_1",
                           "uPP20_1","h_20101214_c1217_640.5307_1",
                           "h_20101214_c1336_788.5288_1",
                           "h_20101214_c1421_219.5293_1",
                           "h_20101214_c1421_455.5292_1",
                           "h_20101214_c1720i_2036.5304_1",
                           "h_20101214_c1892_2119.5284_1",
                           "h_20101214_c2493_i_601.5301_1",
                           "h_20101214_c2804_701.5309_1",
                           "h_20101214_c3117_1457.5310_1",
                           "h_20101214_c3926_348.5286_1",
                           "h_20101214_c3997_i_419.5298_1",
                           "h_20101214_c4769_1106.5305_1",
                           "h_20101214_c5096_985.5294_1",
                           "h_20101214_c5876_431.5299_1",
                           "h_20101214_rep_c542_236.5287_1",
                           "h_20101214_rep_c6068_457.5290_1",
                           "h_20101214_rep_c664_2300.5297_1",
                           "h_20101214_rep_c707_1118.5296_1",
                           "h_20101214_rep_c707_1234.5303_1")],
                  ncode=3,ind.names=JDDcc$ID, 
                  pop=JDDcc$geo_area,missing=NA,ploidy=1)
#include the coordinates of the samples
JDDall@pop.names=c("Western Europe","Central Europe","Sweden","Finland",
                   "Estonia")
JDDall@other$xy<-JDDcc[,c("longitude","latitude")]

#now we analyse the adegenet format dataset with dapc
JDDade<-JDDall
#determination of the number of clusters
clustJDDade<- find.clusters(JDDade,max.n.clust=35)
#with 40 PCs, we lost nearly no information
clustJDDade<- find.clusters(JDDade,n.pca=30,max.n.clust=20) #chose 4 clusters
#which individuals in which clusters per population
table(pop(JDDade),clustJDDade$grp)
#DAPC, first we try to optimized the number of principal component (PCs) 
#to retain to perform the analysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=30)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=10)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=3,n.pca=4)
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







