####GSP salt affected soil maps
####image indices
####7/2/20 SKB
########################################################################################
###session set up

install.packages(c("raster", "sp", "rgdal", "car", "carData", "dplyr", "spacetime",
                   "gstat", "automap", "randomForest", "fitdistrplus", "e1071", "caret",
                   "soilassessment", "soiltexture", "GSIF", "aqp", "plyr", "Hmisc",
                   "corrplot", "factoextra", "spup", "purrr", "lattice", "ncf", 
                   "npsurvSS", "qrnn", "nnet", "mda", "RColorBrewer", "vcd", "readxls"
                   , "maptools", "neuralnet","psych", "tidyr"))

#package "lsei" not available

library(sp); library(foreign); library(rgdal); library(car);library(carData); library(maptools); library(spacetime); library(gstat); library(automap); library(randomForest);library(fitdistrplus); library(e1071); library(caret); library(raster); library(soilassessment); library(soiltexture); library(GSIF); library(aqp); library(plyr); library(Hmisc); library(corrplot); library(factoextra); library(spup); library(purrr); library(lattice); library(ncf); library(npsurvSS); library(nnet); library(class); library(mda); library(RColorBrewer); library(vcd); library(grid); library(neuralnet);library(readxl); library(psych); library(qrnn); library(dplyr); library (tidyr)

#library(lsei)

setwd("E:\\GSP\\1km_ covariates\\spectral")
getwd()


###read in layers and create predictors data frame
predictors <- readGDAL("LB1_blue.tif")
predictors$green <- readGDAL("LB2_green.tif")$band1
predictors$red <- readGDAL("LB3_red.tif")$band1
predictors$nir <- readGDAL("LB4_nir.tif")$band1
predictors$swir1 <- readGDAL("LB5_swir1.tif")$band1
predictors$swir2 <- readGDAL("LB6_swir2.tif")$band1
predictors$blue <- predictors$band1
predictors$band1 = NULL 
summary(predictors)

#remove NAs ##only do this if calculating PCA
#predictors@data <- na.omit(predictors@data)
#summary(predictors)

###generate salinity indices (13 total)
predictors$ndvi <- imageIndices(predictors$blue, predictors$green, predictors$red, predictors$nir, predictors$swir1, predictors$swir2, "NDVI")
summary(predictors$ndvi)
hist(predictors$ndvi)

#canned NDSI function doesn't work correctly
#predictors$ndsi <- imageIndices(predictors$blue, predictors$green, predictors$red, predictors$nir, predictors$swir1, predictors$swir2, "NDSI")
#summary(predictors$ndsi)
#hist(predictors$ndsi)

predictors$ndsi <- ((predictors$red - predictors$nir)/(predictors$red + predictors$nir))
summary(predictors$ndsi)
hist(predictors$ndsi)

predictors$si1 <- imageIndices(predictors$blue, predictors$green, predictors$red, predictors$nir, predictors$swir1, predictors$swir2, "SI1")
summary(predictors$si1)
hist(predictors$si1)

predictors$si2 <- imageIndices(predictors$blue, predictors$green, predictors$red, predictors$nir, predictors$swir1, predictors$swir2, "SI2")
summary(predictors$si2)
hist(predictors$si2)

predictors$si3 <- imageIndices(predictors$blue, predictors$green, predictors$red, predictors$nir, predictors$swir1, predictors$swir2, "SI3")
summary(predictors$si3)
hist(predictors$si3)

predictors$si4 <- imageIndices(predictors$blue, predictors$green, predictors$red, predictors$nir, predictors$swir1, predictors$swir2, "SI4")
summary(predictors$si4)
hist(predictors$si4)

predictors$si5 <- imageIndices(predictors$blue, predictors$green, predictors$red, predictors$nir, predictors$swir1, predictors$swir2, "SI5")
summary(predictors$si5)
hist(predictors$si5)

predictors$si6 <- imageIndices(predictors$blue, predictors$green, predictors$red, predictors$nir, predictors$swir1, predictors$swir2, "SI6")
summary(predictors$si6)
hist(predictors$si6)

predictors$savi <- imageIndices(predictors$blue, predictors$green, predictors$red, predictors$nir, predictors$swir1, predictors$swir2, "SAVI")
summary(predictors$savi)
hist(predictors$savi)

predictors$vssi <- imageIndices(predictors$blue, predictors$green, predictors$red, predictors$nir, predictors$swir1, predictors$swir2, "VSSI")
summary(predictors$vssi)
hist(predictors$vssi)

predictors$sr <- imageIndices(predictors$blue, predictors$green, predictors$red, predictors$nir, predictors$swir1, predictors$swir2, "SR")
summary(predictors$sr)
hist(predictors$sr)

predictors$crsi <- imageIndices(predictors$blue, predictors$green, predictors$red, predictors$nir, predictors$swir1, predictors$swir2, "CRSI")
summary(predictors$crsi)
hist(predictors$crsi)

predictors$bi <- imageIndices(predictors$blue, predictors$green, predictors$red, predictors$nir, predictors$swir1, predictors$swir2, "BI")
summary(predictors$bi)
hist(predictors$bi)

##transform skewed indices 
predictors$si1 = sqrt(predictors$si1)
hist(predictors$si1)

predictors$si2 = sqrt(predictors$si2)
hist(predictors$si2)

predictors$si3 = sqrt(predictors$si3)
hist(predictors$si3)

predictors$si5 = sqrt(predictors$si5)
hist(predictors$si5)

###PCA for image indices
#extract indices from predictors spatialgriddataframe
#predictors_c = predictors@data[ ,c("si1","si2","si3","si4","si5","si6","savi","vssi","ndsi","ndvi","sr","crsi", "bi")]
#predictors_c <- na.omit(predictors_c)
#summary(predictors_c)

#examine correlation between indices
#ind_cor = cor(predictors_c) 
#corrplot(ind_cor, method="number", number.cex = 0.8) 

#run PCA and view eigenvectors plot
#pca <- prcomp(predictors_c[], scale=TRUE)
#fviz_eig(pca)

##return the selected PCs to the predictors spatialgriddataframe
#pred_pcs <- predict(pca,predictors_c[]) 
#predictors@data$PCA1 = pred_pcs[,1]
#predictors@data$PCA2 = pred_pcs[,2]
#predictors@data$PCA3 = pred_pcs[,3]
#summary(predictors)


###convert spatialgriddataframe to raster objects
#p <- as(predictors, 'SpatialPixelsDataFrame') #only if coverting to a raster brick

ndvi <- raster(predictors, layer = 7)
ndsi <- raster(predictors, layer = 8)
si1 <- raster(predictors, layer = 9)
si2 <- raster(predictors, layer = 10)
si3 <- raster(predictors, layer = 11)
si4 <- raster(predictors, layer = 12)
si5 <- raster(predictors, layer = 13)
si6 <- raster(predictors, layer = 14)
savi <- raster(predictors, layer = 15)
vssi <- raster(predictors, layer = 16)
sr <- raster(predictors, layer = 17)
crsi <- raster(predictors, layer = 18)
bi <- raster(predictors, layer = 19)
