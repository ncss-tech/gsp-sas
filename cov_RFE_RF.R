###GSP salt affected soil maps
###Covariate preparation, stacking, extraction, reduction
###July 2020 JMP, SKB
#########################################################################

###session

#check/load packages

required.packages <- c('sp','rgdal','raster','snow','snowfall','parallel','tidyr','car','carData','dplyr','spacetime','gstat','automap','randomForest','fitdistrplus','e1071','caret','
          soilassessment','soiltexture','GSIF','aqp','plyr','Hmisc','corrplot','factoextra','spup','purrr','lattice','ncf','npsurv','lsei',
          'qrnn','nnet','mda','RColorBrewer','vcd','readxls','maptools','neuralnet','psych', 'fmsb','doParallel')

new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase active memory useable by raster package: Windows only
# memory.limit(500000)
## Raster settings: adjust based on system
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08, memfrac = 0.9)

#set working directory
####setwd("E:/DSMFocus/Salinity/Covariates_ISRIC/ISRIC/CONUS")
#setwd("E:/DSMFocus/Salinity/Covariates_ISRIC/ISRIC/testset")
setwd("K:/GSP/1km_ covariates/CONUS/landsat")
getwd()

###Create Image Indices###
library(soilassessment)

###read in layers and create predictors data frame
predictors <- readGDAL("Landsat_B1_1km_average.tif")
predictors$green <- readGDAL("Landsat_B2_1km_average.tif")$band1
predictors$red <- readGDAL("Landsat_B3_1km_average.tif")$band1
predictors$nir <- readGDAL("Landsat_B4_1km_average.tif")$band1
predictors$swir1 <- readGDAL("Landsat_B5_1km_average.tif")$band1
predictors$swir2 <- readGDAL("Landsat_B6_1km_average.tif")$band1
predictors$blue <- predictors$band1
predictors$band1 = NULL 
summary(predictors)



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

#make raster objects
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

#stack image index rasters
sinstack <- stack(ndvi, ndsi, si1, si2, si3, si4, si5, si6, savi, vssi, sr, crsi, bi)

#########################################################################
###Prepare all covariates###

#read in ISRIC CONUS raster layers and create stack
setwd("K:/GSP/1km_ covariates/CONUS")

# read in raster layers from external layers and create list
covlist <- list.files(pattern=".tif$")

#check extents
#me <- extent(raster(covlist[1]))
#mxmin <-  me@xmin
#mymin <- me@ymin
#mxmax <- me@xmax
#mymax <-  me@ymax

#for(r in covlist){
  
  #ce <-  extent(raster(r))
  #cxmin <-  ce@xmin
  #cymin <- ce@ymin
  #cxmax <- ce@xmax
  #cymax <-  ce@ymax
  
  #if(!cxmin == mxmin || !cymin == mymin || !cxmax == mxmax || !cymax == mymax){
    
    #print(paste0(r, " extent does not match ", covlist[1]))
  #}
#}

#stack ISRIC CONUS covariates
tifstack <- stack(covlist)
extent(tifstack)

###fix other rasters with extent issues###
#set extent raster
rs <- raster("B02CHE3_CONUS.tif")
extent(rs) 

#make a list of problem rasters (all rasters in list must have the same extent)
setwd ("K:/GSP/1km_ covariates/CONUS/fix")
outpath <- "K:/GSP/1km_ covariates/CONUS/fix/out/"
#summary(raster("K:/GSP/1km_ covariates/CONUS/fix/us_140evc__majority.tif"))
dir.create(outpath)
fixlist <- list.files(pattern="tif$")
outfiles <- paste0(outpath, fixlist) 
extension(outfiles) <- 'tif'
extent(raster("gsp_slopelength_rs.tif"))

#fix extent 
for(i in 1:length(fixlist)){
    r <- raster(fixlist[i])
    fr <- setExtent(r, rs, keepres = T, snap = F)
    writeRaster(fr, outfiles[i], overwrite = T)
    print(outfiles[i])
  }

#check fix
extent(raster("./out/us_140evc__majority.tif"))

#stack fixed rasters
fstack <- stack(outfiles)
extent(fstack)

### create raster stack of all covariates###
cov_stack <- stack(tifstack, sinstack, fstack)
names(cov_stack)
extent(cov_stack)

###Prepare point data###
# Load RData file with lab data
load(file = "K:/GSP/scripts/LDM-compact_20200709.RData")

#add location info to training data
site_sub <- site(spc)[ ,c('pedon_key','latitude_std_decimal_degrees','longitude_std_decimal_degrees')]
summary(site_sub)

data_site <- merge(s_mps_sf, site_sub, by = "pedon_key", all=TRUE)
str(data_site)
summary(data_site) #check for NAs
pts <- data_site %>% drop_na(latitude_std_decimal_degrees) #remove NA locations
summary(pts)

#remove imprecise locations
pts$latnchar <- nchar(abs(pts$latitude_std_decimal_degrees))
pts$longnchar <- nchar(abs(pts$longitude_std_decimal_degrees))
ptsc <- subset(pts, pts$latnchar > 5 & pts$longnchar > 6)

#turn training data into spatial file
shp.pts <- ptsc
coordinates(shp.pts) <- ~ longitude_std_decimal_degrees + latitude_std_decimal_degrees
temp.proj <- CRS("+proj=longlat +datum+WGS84")
projection(shp.pts) <- temp.proj  
  
#match points with covariates
projgrid <- raster(covlist[1])
cov.proj <- projection(projgrid)
shp.pts.proj <- spTransform(shp.pts, CRS(cov.proj))

plot(projgrid)
plot(shp.pts.proj, add=TRUE)

###extract covariates to points###
# Parallelized extract: (larger datasets)
cpus <- detectCores(all.tests = FALSE, logical = TRUE)-1
sfInit(parallel=TRUE, cpus=cpus)
sfExport("shp.pts.proj", "covlist")
sfLibrary(raster)
sfLibrary(rgdal)
ov.lst <- sfLapply(covlist, function(i){try( raster::extract(raster(i), shp.pts.proj) )})
snowfall::sfStop()
ov.lst <- as.data.frame(ov.lst)
names(ov.lst) = tools::file_path_sans_ext(basename(covlist))
ov.lst$DID <- seq.int(nrow(ov.lst))
shp.pts.proj$DID <- seq.int(nrow(shp.pts.proj))
pts.ext <- merge(as.data.frame(shp.pts.proj),ov.lst, by="DID")

#check names to ensure covariate and predictor data in same place
names(pts.ext)

#save points with extracted covariates
setwd("E:/DSMFocus/Salinity")
saveRDS(pts.ext, "pts_ext_covs.rds")

##save as shapefile if wanted
#gsppoints <- pts.ext
#gsppointscl <- gsppoints[, -c(1:11, 13:28)]
#coordinates(gsppointscl) <- ~ longitude_std_decimal_degrees + latitude_std_decimal_degrees # modify with field names in table
#projection(gsppoints) <- temp.proj

#writeOGR(gsppointscl, ".", "gsppoints", driver="ESRI Shapefile")



############ADD IN DATA SPLITTING HERE############




##clean out memory
gc()

#############################################################################

###Covariate data reduction using recursive feature elimination###
###Should be done just on training points###

##data prep

#read in points that include extracted covariate values
comp <- as.data.frame(readRDS("pts_ext_covs.rds"))
comp <- as.data.frame(pts.ext)

## testing block with smaller dataset ##
#training.pts <- readOGR("gsppointsclip.shp") #trimmed down set of points for testing
#names(training.pts)
#plot(projgrid)
#plot(training.pts, add=TRUE)
#
#comp <- as.data.frame(training.pts) #for testing
#
#
#names(comp)

#remove unneeded columns and change the first column name to Prop
comp <- comp[, -c(1:10,12:30)]
names(comp)
names(comp)[c(1)] <- c('Prop')
names(comp)

#remove NA values from Prop
compcl <- comp[complete.cases(comp), ]
summary (compcl)

#recursive feature elimination

#change subsets to match number of covariates

#check number of covariates
length(compcl) # number of covariates plus the Prop column

subsets <- c(1:(length(compcl)-1))
#subsets <- c(1:25,50)
length(subsets)

#set seeds to get reproducible results when the process in parallel
set.seed(915)
seeds <- vector(mode = "list", length=176) #length is 1 more than number of covariates CHANGE AS NEEDED
for(i in 1:15) seeds[[i]] <- sample.int(1000,16) # 1:x x=folds*repeats (for 5 folds & 3 repeats x=15) and sample.int(1000,x+1)
seeds[[16]] <- sample.int(1000, 1) #seeds = folds*repeats plus 1

#set up the rfe control
ctrl.RFE <- rfeControl(functions = rfFuncs,
                       method = "repeatedcv",
                       number = 5, #changed to from 15 folds with 5 repeats to speed up process
                       repeats = 3,
                       seeds = seeds,
                       verbose = FALSE)

#highlight and run everything from c1 to stopCluster(c1) to run RFE

c1 <- makeCluster(detectCores()-1)
registerDoParallel(c1)
set.seed(9)
rf.RFE <- rfe(x = compcl[,-1],
              y = compcl$Prop,
              sizes = subsets,
              rfeControl = ctrl.RFE,
              allowParallel = TRUE
)
stopCluster(c1)

#check results
rf.RFE