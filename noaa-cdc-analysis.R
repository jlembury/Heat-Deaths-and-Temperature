#An analysis of temperature distribution and hyperthermia from 1999-2016
#Jessica Embury

#packages
if (!require(sp)){
  install.packages('sp')
  library(sp)
}
if (!require("raster")){
  install.packages('raster')
  library(raster)
}
if (!require("rgdal")){
  install.packages('rgdal')
  library(rgdal)
}
if (!require("spatstat")){
  install.packages('spatstat')
  library(spatstat)
}
if (!require("maptools")){
  install.packages('maptools')
  library(maptools)
}
if (!require(GISTools)){
  install.packages('GISTools')
  library(GISTools)
}
if (!require(spdep)){
  install.packages('spdep')
  library(spdep)
}
if (!require(gstat)){
  install.packages('gstat')
  library(gstat)
}
if (!require(splancs)){
  install.packages('splancs')
  library(splancs)
}
if (!require(lattice)){
  install.packages('lattice')
  library(lattice)
}
if (!require(pgirmess)){
  install.packages('pgirmess')
  library(pgirmess)
}
if (!require(RColorBrewer)){
  install.packages('RColorBrewer')
  library(RColorBrewer)
}
if (!require(classInt)){
  install.packages('classInt')
  library(classInt)
}
if (!require(spgwr)){
  install.packages('spgwr')
  library(spgwr)
}
if (!require(GWmodel)){
  install.packages('GWmodel')
  librGWmodel
}
if (!require(e1071)){
  install.packages('e1071')
  libre1071
}
if (!require(boot)){
  install.packages('boot')
  librboot
}
if (!require(PerformanceAnalytics)){
  install.packages('PerformanceAnalytics')
  librPerformanceAnalytics
}
if (!require(dplyr)){
  install.packages('dplyr')
  library(dplyr)
}


#path to county polygon shapefile containing temperature and death information
path <- 'output/shapefiles/heat.shp'

#read file >> creates spatialpolygondataframe
county <- rgdal::readOGR(path)

#remove La Paz County -- outlier with CDR = 5, all others CDR <=2
#county <- county[!(county$County=="La Paz"),]

#remove counties with CDR = 0
#county <- county[!(county$CDR=="0"),]

######################
#Variable Correlation#
######################

#convert spatial dataframe to dataframe
df <- as.data.frame(county)

#subset county spatialdataframe to only include desired variables
subset <- df[,c("CDR","Max_Temp","Percent90")]

#plot variable correlations
plot.new()
chart.Correlation(subset, histogram=TRUE, pch=19)

########################
#Moran I Test Statistic#
########################

#get county centroids
centroids <- coordinates(county)

#determine neighbors, create contiguity matrix
neighbors <- poly2nb(county, queen=T)

#create weighted distance matrix
wdist_matrix <- nb2listw(neighbors,style="W", zero.policy=TRUE)

#plot connections between neighbors
plot.new()
par(mar=rep(0,4))
plot(wdist_matrix,coords=centroids,pch=19, cex=0.1, col="gray")

#Calculate Global Moran's I based on deaths.
moran.test(county$CDR, listw=wdist_matrix, zero.policy=T)

###################
#Moran Scatterplot#
###################

#Create a Moran scatterplot.
plot.new()

#parameters
par(mar=c(4,4,1.5,0.5))

#scatterplot
moran.plot(county$CDR, listw=wdist_matrix, zero.policy=T, xlim=c(0,2.5),ylim=c(0,2.5), pch=16, col="black",cex=.5, quiet=F, labels=as.character(county$County),xlab="CDR", ylab="CDR (Spatial Lag)", main="Moran Scatterplot")

##############
#GWR Analysis#
##############

#GWR analysis of heat-related deaths (as CDR), maximum temperature,and percent of days above 90F, by county

#create grid: GridTopology(origin, cell size, number cells)
grid <- SpatialGrid(GridTopology(c(-125,24),c(0.1,0.1),c(585, 260)))

#plot grid
plot.new()
plot(grid)
plot(county,add=TRUE,col=adjustcolor('gray',alpha.f=0.5))

#distance matrix
#compute distances between points on grid, and points (locations)
dist_mat <- gw.dist(dp.locat=coordinates(county),rp.locat=coordinates(grid))

#basic GWR analysis
#bw = size of neighborhood, kernel = distance-decay
#dep_var~indep_var are the variables being compared (columns in car_theft dataframe)
gwr.res <- gwr.basic(CDR~Max_Temp+Percent90, data=county, regression.points=grid, bw = 10000, dMat=dist_mat,kernel='gaussian')
gwr.res

#Image of GWR result for CDR & Max_Temp
plot.new()
image(gwr.res$SDF,'Max_Temp')
plot(county,add=TRUE)

#test reliability of result (error), using bootstrap estimates for CDR & Max_Temp
plot.new()
set.seed(4676)
gwrcoef <- function(hpdf,i) gwr.basic(CDR~Max_Temp, data=county[i,], regression.points=grid, bw=10000, dMat=dist_mat[i,],kernel='gaussian')$SDF$Max_Temp
bootres <- boot(county,gwrcoef,100)
gwr.res$SDF$bseMax_Temp <- sqrt(apply(bootres$t,2,var))
image(gwr.res$SDF,'bseMax_Temp')
plot(county,add=TRUE)

#estimate bias from bootstrap sampling for CDR & Max_Temp
plot.new()
gwr.res$SDF$biasMax_Temp <- bootres$t0 - apply(bootres$t,2,mean)
image(gwr.res$SDF,'biasMax_Temp')
plot(county,add=TRUE)

#Image of GWR result for CDR & Percent90 (percentage of days above 90 degrees F)
plot.new()
image(gwr.res$SDF,'Percent90')
plot(county,add=TRUE)

#test reliability of result (error), using bootstrap estimates for CDR & Percent90
plot.new()
set.seed(4676)
gwrcoef <- function(hpdf,i) gwr.basic(CDR~Max_Temp+Percent90, data=county[i,], regression.points=grid, bw=10000, dMat=dist_mat[i,],kernel='gaussian')$SDF$Percent90
bootres <- boot(county,gwrcoef,100)
gwr.res$SDF$bsePercent90 <- sqrt(apply(bootres$t,2,var))
image(gwr.res$SDF,'bsePercent90')
plot(county,add=TRUE)

#estimate bias from bootstrap sampling for CDR & Percent90
plot.new()
gwr.res$SDF$biasPercent90 <- bootres$t0 - apply(bootres$t,2,mean)
image(gwr.res$SDF,'biasPercent90')
plot(county,add=TRUE)
