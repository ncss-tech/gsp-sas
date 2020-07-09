###GSP salt affected soil maps
###Covariate stacking and extraction
###7/7/2020 JMP
#########################################################################

###session

#check/load packages

required.packages <- c('sp','rgdal','raster','snow','snowfall','parallel','tidyr','car','carData','dplyr','spacetime','gstat','automap','randomForest','fitdistrplus','e1071','caret','
          soilassessment','soiltexture','GSIF','aqp','plyr','Hmisc','corrplot','factoextra','spup','purrr','lattice','ncf','npsurv','lsei',
          'qrnn','nnet','mda','RColorBrewer','vcd','readxls','maptools','neuralnet','psych')

new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase active memory useable by raster package: Windows only
# memory.limit(500000)
## Raster settings: adjust based on system
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08, memfrac = 0.9)
#set working directory

setwd("E:/DSMFocus/Salinity/Covariates_ISRIC/ISRIC/CONUS")
getwd()

#read in raster layers and create stack

# read in raster layers from external layers and create list
covlist <- list.files(pattern=".tif$")
covlist

# create raster stack
###cov_stack <- stack(rlist)
###names(cov_stack)

# Load RData file with lab data
load(file = "E:/DSMFocus/Salinity/LDM-compact_20200708.RData")

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

#save points with extrated covariates
setwd("E:/DSMFocus/Salinity")
saveRDS(pts.ext, "pts_ext_covs.rds")
