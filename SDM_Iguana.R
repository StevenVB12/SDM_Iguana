##
# The march of the Common Green Iguana (Iguana iguana): early establishment in Singapore and Thailand is facilitated by the pet trade and recreational parks
# Matthijs P. van den Burg, Steven M. Van Belleghem, Christina N. De Jesús Villanueva
##

## Check if everything is installed
# install.packages("rgdal")
# install.packages("raster")
# install.packages("dismo")
# install.packages("rJava")
# install.packages("maptools")
# install.packages("spThin")

library(dismo)
library(maptools)
library(spThin)


# read long lat tables
## Invasive record Thailand
iguana_LL <- read.table("Iguana_Thailand.txt", h = T)
iguana_LL <- iguana_LL[c("lon","lat")]

## All records
iguana_nat <- read.table('Iguana_iguana_Occurrence_Data_16Jan2020_inv_and_nat.txt', h=T)
iguana_nat <- iguana_nat[c("lon","lat")]

# remove duplicates
iguana_LL <- unique(iguana_LL)
iguana_nat <- unique(iguana_nat)


# plot Thailand data on map
data(wrld_simpl)
plot(wrld_simpl, xlim=c(100,110), ylim=c(0, 20), axes=TRUE,
       col="light yellow")
box()

points(iguana_LL$lon, iguana_LL$lat, col='red', pch=19, cex=1)


# get the bioclim data
predictors <- getData('worldclim', res=2.5, var='bio')

plot(predictors, 1, ext=c(90,120,-10,30)) # plots one bioclim layer
points(iguana_LL, col='blue', pch=19)

# Add elevation to predictors
predictors_alt <- getData('worldclim', res=2.5, var='alt')
predictors <- stack(predictors, predictors_alt)

# Remove observations above 1000m
presvals <- raster::extract(predictors, iguana_nat)
presvals <- as.data.frame(presvals)

nrow(presvals)
presvals1000 <- subset(presvals, presvals$alt >= 1000)
nrow(presvals1000)
iguana_nat <- iguana_nat[-as.numeric(rownames(presvals1000)),]
nrow(iguana_nat)

# Thin occurence data 50km (change thin parameter for higher distances)
iguana_nat$spec <- 'SPEC'

thin(iguana_nat, lat.col = "lat", long.col = "lon", spec.col = 'spec', thin.par = 50, reps = 1, out.dir = './', out.base = 'Iguana_native_thinned50km_inv')

iguana_nat_50km <- read.csv("Iguana_native_thinned50km_thin1.csv", h=T)[c("lon","lat")]

plot(predictors, 1, ext=c(-120,-20,-50,50)) # plots one bioclim layer

points(iguana_nat_50km$lon, iguana_nat_50km$lat, col='blue', pch=19, cex=0.1)


###
# maxent (machine learning model for predicting)
###

# Select bio layers of interest
predictors_bio <- subset(predictors, c(1,3,4,6,7,10,11,13,15,20)) 

# Run maxent model
xm <- maxent(predictors, iguana_nat_50km)

# See how much the bioclim variable contribute to the model
xm
plot(xm)
response(xm) # representing the probability of habitat suitability

# Run maxent and evaluate model using 80% training and 20% test data
# (Takes a long time. Only run to get AUC of final model)
AUCs <- c()
for(i in 1:10){
  
  group.k1 <- kfold(iguana_nat_50km, 5)
  iguana_nat_50km_train = iguana_nat_50km[group.k1 != 1,]
  iguana_nat_50km_test = iguana_nat_50km[group.k1 == 1,]
  
  
  xm_train <- maxent(predictors, iguana_nat_50km_train)
  
  presvals <- extract(predictors, iguana_nat_50km_test)
  
  backgr <- randomPoints(predictors, 10000)
  absvals <- raster::extract(predictors, backgr)
  
  # mod.ev <- evaluate(model = xm_train, p = presvals, a= absvals)
  
  
  
  sb <- ssb(iguana_nat_50km_test, as.data.frame(backgr), iguana_nat_50km_train)
  sb[,1] / sb[,2]
  
  # The following code takes absence values close to presence samples
  i <- pwdSample(iguana_nat_50km_test, as.data.frame(backgr), iguana_nat_50km_train, n=1, tr=0.1)
  pres_test_pwd <- iguana_nat_50km_test[!is.na(i[,1]), ]
  back_test_pwd <- as.data.frame(backgr)[na.omit(as.vector(i)), ]
  sb2 <- ssb(pres_test_pwd, back_test_pwd, iguana_nat_50km_train)
  sb2[1]/ sb2[2]
  
  mod.ev1 <- evaluate(model=xm_train, p=iguana_nat_50km_test, a=as.data.frame(backgr), x=predictors)
  mod.ev2 <- evaluate(model=xm_train, p=pres_test_pwd, a=back_test_pwd, x=predictors)
  
  tr1 <- threshold(mod.ev1, 'spec_sens')
  tr2 <- threshold(mod.ev2, 'spec_sens')
  
  
  AUCs <- rbind(AUCs, mod.ev1@auc, mod.ev2@auc, tr1, tr2)
  
  print(c(mod.ev1@auc, mod.ev2@auc, tr1, tr2))
}

colnames(AUC) <- c('AUC1','AUC2','Threshold1','Threshold2')

AUC 

## AUC1: AUC for regular model
## AUC2: AUC for model with absence locations forced to be close to presence locations
## th1: th for regular model
## th2: th for model with absence locations forced to be close to presence locations



###
# Plot predictions
###
library("rnaturalearth")
library("rnaturalearthdata")

# Make predictions Thailand
px <- predict(predictors, xm, ext=c(95,110,0,22), progress='')

colfunc <- colorRampPalette(c('white','red'))

world <- ne_countries(scale = "medium", returnclass = "sp")

par(mar=c(4,4,2,2), xpd=F)
plot(px, col=colfunc(100), xlim = c(95,110), ylim = c(0,22), zlim=c(0.25,1))
plot(world, add=T, border='gray35', col=NA, xlim = c(95,110), ylim = c(0,22), lwd=0.5)
points(iguana_LL, col='blue', pch=19)

# Click on figure to set scale bar
scalebar(d=100, xy=click(), below="km")


# Make prediction Taiwan
px <- predict(predictors, xm, ext=c(115,130,5,30), progress='')

colfunc <- colorRampPalette(c('white','red'))

world <- ne_countries(scale = "medium", returnclass = "sp")

par(mar=c(4,4,2,2), xpd=F)
plot(px, col=colfunc(100), xlim = c(115,130), ylim = c(5,30), zlim=c(0.25,1))
plot(world, add=T, border='gray35', col=NA, xlim = c(115,130), ylim = c(5,30), lwd=0.5)
points(iguana_nat, col='blue', pch=19)

# Click on figure to set scale bar
scalebar(d=100, xy=click(), below="km")


# Make prediction Americas
px <- predict(predictors, xm, ext=c(-120,-20,-40,40), progress='')

colfunc <- colorRampPalette(c('white','red'))

world <- ne_countries(scale = "medium", returnclass = "sp")

par(mar=c(4,4,2,2), xpd=F)
plot(px, col=colfunc(100), xlim = c(-120,-20), ylim = c(-40,40), zlim=c(0.25,1))
plot(world, add=T, border='gray35', col=NA, xlim = c(-120,-20), ylim = c(-40,40), lwd=0.5)
points(iguana_nat_50km$lon, iguana_nat_50km$lat, col='blue', pch=19, cex=0.1)

# Click on figure to set scale bar
scalebar(d=1000, xy=click(), below="km")


##
# Make predictions 2050
##

# Get predictors
predictors_bio <- getData('CMIP5', res=2.5, var='bio', year = 50, model = 'AC', rcp = '45')
predictors <- stack(predictors_bio, predictors_alt)
predictors <- subset(predictors, c(1,3,4,6,7,10,11,13,15,20))

# Run maxent
xm <- maxent(predictors, iguana_nat_50km)

# Make prediction Thailand 
px <- predict(predictors, xm, ext=c(95,110,0,22), progress='')

colfunc <- colorRampPalette(c('white','red'))

world <- ne_countries(scale = "medium", returnclass = "sp")

par(mar=c(4,4,2,2), xpd=F)
plot(px, col=colfunc(100), xlim = c(95,110), ylim = c(0,22), zlim=c(0.25,1))
plot(world, add=T, border='gray35', col=NA, xlim = c(95,110), ylim = c(0,22), lwd=0.5)
points(iguana_LL, col='blue', pch=19)

# Click on figure to set scale bar
scalebar(d=100, xy=click(), below="km")



###
# limiting: Shows what variables best exlain absence of Iguanas
###

library(devtools)
install_github('johnbaums/rmaxent')

library(rmaxent)
library(rasterVis)

lim <- limiting(predictors, xm)
levelplot(lim, col.regions=rainbow) +
  layer(sp.points(SpatialPoints(iguana_nat_50km), pch=20, col=1))

levelplot(lim, col.regions=rainbow, xlim=c(-120,-20), ylim=c(-50,50)) +
  layer(sp.points(SpatialPoints(iguana_nat_50km), pch=20, col=1))

levelplot(lim, col.regions=rainbow, xlim=c(75,130), ylim=c(-20,40)) +
  layer(sp.points(SpatialPoints(iguana_nat_50km), pch=20, col=1))

