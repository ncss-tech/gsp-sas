####GSP salt affected soil maps
####image indices
####7/2/20 SKB
########################################################################################
###session set up

libs <- c("rgdal", "raster", "soilassessment")
lapply(libs, require, character.only = TRUE)

# setwd("D:\\GSP\\1km_ covariates\\spectral")
getwd()


# ISRIC 1km layers
vars <- c("REDL00", "NIRL00", "SW1L00", "SW2L00", "REDL14", "NIRL14", "SW1L14", "SW2L14", 
          paste0("T", sprintf("%02d", 1:12), "MSD3"),
          paste0("T", sprintf("%02d", 1:12), "MOD3"),
          paste0("N", sprintf("%02d", 1:12), "MSD3"),
          paste0("N", sprintf("%02d", 1:12), "MOD3")
          )
path <- "D:/geodata/project_data/gsp-sas/1km covariates/ISRIC/CONUS"
lf   <- list.files(path = path, pattern = ".tif$")
idx  <- sapply(vars, function(x) grep(x, lf))
il   <- stack(lapply(file.path(path, lf[idx]), raster))
names(il) <- vars

il00   <- il[[1:4]]
il14   <- il[[5:8]]
ilmsd3 <- il[[9:length(vars)]]
NAvalue(ilmsd3) <- -32768
# ilmsd3_t <- as(ilmsd3[[c(13:24, 37:48)]], "SpatRaster")


# Other 1km layers
vars <- expand.grid(source = "Landsat_B",
                    band   = c(1:6),
                    method = c("average"), #, "max"),
                    stringsAsFactors = FALSE
                    )
vars <- with(vars, paste0(source, band, "_1km_", method))
path <- "D:/geodata/project_data/gsp-sas/1km covariates/Other"
lf   <- list.files(path = path, pattern = ".tif$")
idx <- sapply(vars, function(x) grep(x, lf))

ol <- stack(lapply(file.path(path, lf[idx]), raster))
ol <- stack(ol)
names(ol) <- vars



# soilassessment indices
sa_indices <- function(blue = NULL, green =  NULL, red = NULL, nir = NULL, swir1 = NULL, swir2 = NULL) {
  
  NSI  <- imageIndices(nir   = nir, 
                       swir1 = swir1, 
                       swir2 = swir2, 
                       index = "NSI"
                       )
  SI4 <- imageIndices(nir   = nir, 
                      swir1 = swir1, 
                      index = "SI4"
                      )
  SAVI <- imageIndices(red   = red, 
                       nir   = nir, 
                       index = "SAVI"
                       )
  NDSI <- imageIndices(red   = red, 
                       nir   = nir, 
                       index = "NDSI"
                       )
  NDVI <- imageIndices(red   = red, 
                       nir   = nir, 
                       index = "NDVI"
                       )
  ROCK <- imageIndices(nir   = nir, 
                       swir1 = swir1, 
                       index = "ROCK"
                       )
  
  rs <- stack(NSI, SI4, SAVI, NDSI, NDVI, ROCK)
  names(rs) <- c("NSI", "SI4", "SAVI", "NDSI", "NDVI", "ROCK")
  
  if (!is.null(blue) & !is.null(green) & !is.null(red)) {
    
    SI1  <- sqrt(green * red)
    SI2  <- sqrt(blue * red)
    SI3  <- sqrt((green)^2 + (red)^2)
    SI5  <- blue / red
    SI6  <- red * nir / green
    VSSI <- 2 * green - 5 * (red + nir)
    SR   <- (green - red)/(blue + red)
    
    rs2  <- stack(SI1, SI2, SI3, SI5, SI6, VSSI, SR)
    names(rs2) <- c("SI1", "SI2", "SI3", "SI5", "SI6", "VSSI", "SR")
    rs   <- stack(rs, rs2)
    }
  return(rs)
}



path <- "D:/geodata/project_data/gsp-sas/1km covariates/ISRIC/CONUS"

il00_indices <- sa_indices(
  red   = il00$REDL00, 
  nir   = il00$NIRL00, 
  swir1 = il00$SW1L00, 
  swir2 = il00$SW2L00
  )
lapply(names(il00_indices), function(x) {
  writeRaster(il00_indices[[x]], 
              filename = file.path(path, paste0("L00_", x, ".tif")), 
              progress = "text"
  )
})


il14_indices <- sa_indices(
  red   = il14$REDL14, 
  nir   = il14$NIRL14, 
  swir1 = il14$SW1L14, 
  swir2 = il14$SW2L14
)
lapply(names(il14_indices), function(x) {
  writeRaster(il14_indices[[x]], 
              filename = file.path(path, paste0("L14_", x, ".tif")), 
              progress = "text"
  )
})


ilmsd3_indices <- c(
  min = min(ilmsd3_t),
  max = max(ilmsd3_t),
  avg = mean(ilmsd3_t),
  range = diff(range(ilmsd3_t)),
  cv  = stdev(ilmsd3_t) / mean(ilmsd3_t) * 100 
)
lapply(names(ilmsd3_indices), function(x) {
  writeRaster(ilmsd3_indices[[x]], 
              filename = file.path(path, paste0("MOD11A2_", x, ".tif")),
              progress = "text"
  )
})

ilmsd3_pca <- RStoolbox::rasterPCA(ilmsd3, nSamples = 2000, nComp = 10, spca = TRUE, progress = "text")
lapply(names(ilmsd3_pca$map), function(x) {
  writeRaster(ilmsd3_pca$map[[x]], 
              filename = file.path(path, paste0("MOD11A2_", x, ".tif")),
              progress = "text"
  )
})



path <- "D:/geodata/project_data/gsp-sas/1km covariates/Other"

ol_indices <- sa_indices(
  blue  = ol$Landsat_B1_1km_average,
  green = ol$Landsat_B2_1km_average,
  red   = ol$Landsat_B3_1km_average, 
  nir   = ol$Landsat_B4_1km_average, 
  swir1 = ol$Landsat_B5_1km_average, 
  swir2 = ol$Landsat_B6_1km_average
)
lapply(names(ol_indices), function(x) {
  writeRaster(ol_indices[[x]],
              filename = file.path(path, paste0("landsat_", x, "_1km_average.tif")),
              progress = "text"
  )
})

ol_pca     <- RStoolbox::rasterPCA(ol, nSamples = 2000, spca = TRUE, progress = "text")
lapply(names(ol_pca$map), function(x) {
  writeRaster(ol_pca$map[[x]],
              filename = file.path(path, paste0("landsat_", x, "_1km_average.tif")),
              progress = "text"
  )
})

# nlcd
nlcd <- "D:/geodata/land_use_land_cover/NLCD_2016_Land_Cover_L48_20190424.img"
r <- raster("D:/geodata/project_data/gsp-sas/1km covariates/SSURGO/mlra_rs.tif")
gdalUtils::gdalwarp(
  srcfile = nlcd,
  dstfile = "D:/geodata/project_data/gsp-sas/1km covariates/Other/nlcd_2016_conus.tif",
  t_srs = proj4string(r), 
  te = bbox(r),
  r = "mode",
  tr = res(r),
  ot = "Int32",
  verbose   = TRUE,
  overwrite = TRUE
)


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



# https://stackoverflow.com/questions/45428129/how-to-subset-a-raster-based-on-grid-cell-values

library(raster)

dem <- raster("D:/geodata/project_data/gsp-sas/1km covariates/ISRIC/CONUS/DEMENV5_CONUS.tif")

data(wrld_simpl, package = "maptools")

r <- raster(xmn = -180, xmx = 180, ymn = -90, ymx = 90, res = 1)
wm2 <- rasterize(wrld_simpl, r, 1)
wm3 <- mask(is.na(wm2), wm2, maskvalue = 1, updatevalue = NA, progress = "text")
wm3_dist <- distance(wm3, progress = "text")
wm3_dist <- projectRaster(wm3_dist, dem, progress = "text")

writeRaster(wm3_dist, filename = "C:/workspace2/wm3_dist.tif", overwrite = TRUE)


d <- wm3_dist * wm3 