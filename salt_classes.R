###GSP salt affected soil maps
###salt affected soil classes
###August 2020 JMP, SKB, SR
##############################################################################
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

library(soilassessment)

#############################################################################
##this section of code starts on p86 in the GSSmap_Technical_manual_3.pdf##

###predict salt classes###

# setwd("G:/GSP/predictions/predictions")
setwd("D:/geodata/project_data/gsp-sas/predictions_v2")

### Read in prediction layers###
library(raster)

f <- c(PH_030 = "ph_030_qrf_mlra_pm_NULL_final4_fill_mask.tif",
       EC_030 = "ec030bt_fill_mask.tif", 
       ESP_030 = "v2_esp030t25_coast_pm_log_nafw.tif",
       PH_100 = "ph_100_qrf_mlra_pm_NULL_final4_fill_mask.tif",
       EC_100 = "ec100bt_fill_mask.tif",
       ESP_100 = "v2_esp100t25_coast_pm_log_nafw.tif" 
       ) #change file names as needed
#f <- c("ec_030_bt_nomlra_wstatspm.tif", "cv_esp_030_nomlra_statspm.tif", "ph_030_qrf_mlra_pm_NULL_final2.tif") #change file names as needed
rs <- readAll(stack(f))
rs <- projectRaster(rs, crs = "+init=epsg:5070", progress = "text")
rs$EC_030t  <- log(rs$EC_030 + 0.1)
rs$ESP_030t <- log(rs$ESP_030 + 0.1)
rs$EC_100t  <- log(rs$EC_100 + 0.1)
rs$ESP_100t <- log(rs$ESP_100 + 0.1)
rs <- readAll(rs)


rs #inspect range in values

###back transform ec and esp if needed###

# do we still need to do this?
##rs$ec_t <- calc(rs$ec_030_t, function(x) exp(x) - 0.1) 
# rs$esp_t <- calc(rs$cv_esp_100_nomlra_statspm, function(x) exp(x) - 0.1) 
# rs #check which columns to keep in next line
# rs <- rs[[c(1,3:4)]] #only use back-transformed data
rs_030 <- as(rs[[c(1:3, 7:8)]], "SpatialGridDataFrame")
rs_100 <- as(rs[[c(4:6, 9:10)]], "SpatialGridDataFrame")



## run saltiness and salt severity for both 0-30cm and 30-100cm

# saltiness
library(soilassessment)

rs_030$salty <- saltClass(rs_030$EC_030, rs_030$PH_030, rs_030$ESP_030, "FAO") # replace ec, ph, esp variables with predicted layers for each depth interval
rs_030$saltiness <- classCode(rs_030$salty, "saltclass") # add class names
plot(rs_030$saltiness)
plot(rs_030["saltiness"])

rs_100$salty <- saltClass(rs_100$EC_100, rs_100$PH_100, rs_100$ESP_100, "FAO") # replace ec, ph, esp variables with predicted layers for each depth interval
rs_100$saltiness <- classCode(rs_100$salty, "saltclass") #add class names
plot(rs_100$saltiness)
plot(rs_100["saltiness"])



# salt severity
rs_030$salt_affected <- saltSeverity(rs_030$EC_030, rs_030$PH_030, rs_030$ESP_030, "FAO") # replace ec, ph, esp variables with predicted layers for each depth interval
rs_030$saltaffectedness <- classCode(rs_030$salt_affected, "saltseverity") # add class names
barplot(table(rs_030$saltaffectedness))
plot(rs_030$saltaffectedness)
plot(rs_030["saltaffectedness"], col = RColorBrewer::brewer.pal(11, "Spectral"))
spplot(rs_030["saltaffectedness"])

rs_100$salt_affected <- saltSeverity(rs_100$EC_100, rs_100$PH_100, rs_100$ESP_100, "FAO") #replace ec, ph, esp variables with predicted layers for each depth interval
rs_100$saltaffectedness <- classCode(rs_100$salt_affected, "saltseverity") #add class names
barplot(table(rs_100$saltaffectedness))
plot(rs_100$saltaffectedness)
plot(rs_100["saltaffectedness"], col = RColorBrewer::brewer.pal(11, "Spectral"))
spplot(rs_100["saltaffectedness"])

### export salt class maps###

#run for both 0-30cm and 30-100cm
#saltiness (not required)
#rs2$saltyclasses <- as.numeric(rs2$saltiness)
#saltiness_LUT100 <- classLUT(rs2["saltiness"], "saltclass") ##WHY WON'T THIS WORK "object 'LUT' not found" (because classLUT() is for salt affectedness)
#writeGDAL(rs2["saltyclasses"], drivername = "GTiff", "Mid100_saltiness.tif")
#write.table(saltiness_LUT100, file = "saltiness_LUT100.txt", row.names = F)

#run for both 0-30cm and 30-100cm 
rs_030$saltclasses <- as.numeric(rs_030$saltaffectedness) #convert factors to numeric
salinity_LUT30 <- classLUT(rs_030["saltaffectedness"], "saltseverity") #change object name for depth as needed
rgdal::writeGDAL(rs_030["saltclasses"], drivername = "GTiff", "Top030_saltaffected.tif")
write.table(salinity_LUT30, file = "saltaffected_LUT30.txt", row.names = FALSE)

rs_100$saltclasses <- as.numeric(rs_100$saltaffectedness) #convert factors to numeric
salinity_LUT100 <- classLUT(rs_100["saltaffectedness"], "saltseverity") #change object name for depth as needed
rgdal::writeGDAL(rs_100["saltclasses"], drivername = "GTiff", "Mid100_saltaffected.tif")
write.table(salinity_LUT100, file = "saltaffected_LUT100.txt", row.names = F)
########################################################################
###accuracy assessment of salt class maps; p88###

##import and classify validation points for salt classes; these are OBSERVED values
setwd("G:/GSP/pointdata")
getwd()
# setwd("C:/Users/stephen.roecker/Nextcloud/projects/2020_gsp-sas")

#soilv <- readOGR("filepath", "validation dataset") #change file path and validation dataset name 

#load training_data.RData

load(file = "training_data.RData")
load(file = "LDM-compact_20200709.RData")
#also load training data into environment, test is test points

soilv <- test

soilv <- soilv[c(1,10:18)] #get rid of covariate columns,keep pedon_key and observed values
head(soilv)


#### predict salt classes for each test point
soilv$saltaffected1 <- saltSeverity(EC  = soilv$ec_ptf.000_030_cm, 
                                    pH  = soilv$ph_ptf.000_030_cm, 
                                    ESP = soilv$esp.000_030_cm, 
                                    "FAO"
                                    ) #change for each depth
summary(soilv$saltaffected1)

soilv$saltaffectedness1 <- classCode(soilv$saltaffected1, "saltseverity")
summary(soilv$saltaffectedness1)

##extract the salt classes from the map using the validation samples
soilv <- subset(soilv, !is.na(soilv$saltaffectedness1)) #trim test points down to only ones with relevant data
summary(soilv)

library(sf)
s_mps_sf <- st_transform(s_mps_sf, crs = st_crs(5070))
test2 <- cbind(s_mps_sf["pedon_key"], sf::st_coordinates(s_mps_sf)) #extract coordinates from lab data
soilv <- merge(x = test2, y = soilv, by = "pedon_key", all.y = TRUE) #merge coordinates with test data
summary(soilv)

rs_030.ov <- over(as(soilv, "Spatial"), rs_030) #combine predictions with test data ####I DON'T KNOW IF THIS IS DOING ANYTHING.Nothing from soilv shows up with summary(rsxx.ov) 

saltclasses.ovv <- rs_030.ov
summary(saltclasses.ovv) #check for salt_affected and saltaffectedness
summary(soilv) #check for saltaffected1 and saltaffectedness1

soilv$salt_affected <- saltclasses.ovv$salt_affected #take predicted value for salt affected and add to the validataion points
soilv$saltaffectedness <- saltclasses.ovv$saltaffectedness #take predicted value for saltaffectedness and add to the validation points
#check summary of extracted classified pixels
summary(soilv$salt_affected)
summary(soilv$saltaffectedness)
summary(soilv$saltaffected1)
summary(soilv$saltaffectedness1)

##get rid of NAs
soilv <- subset(soilv, complete.cases(saltaffectedness1, saltaffectedness, salt_affected, saltaffected1))

##generate confusion matrix and Kappa

library(vcd); library(mda)

lv <- c(3, 6, 8:17)
soilv <- transform(soilv,
                   sa  = factor(salt_affected, levels = lv),
                   sa1 = factor(saltaffected1, levels = lv)
                   )
agreementplot(table(soilv$sa, soilv$sa1),
              main = "Accuracy assessment",xlab = "Class codes in holdout samples",
              ylab = "Class codes in map")

Kappa(table(soilv$sa, soilv$sa1))


########################################################################
###uncertainty of salt class maps with monte-carlo simulations; p90###

 ##convert input layers to raster files (ours may already be raster objects?) library(sp)
#transform EC and ESP
# rs_030$EC_030t <- log(rs_030$EC_030 + 0.1)
# rs_030$ESP_030t <- log(rs_030$ESP_030 + 0.1)
# 
# EC_030 <- raster(rs_030["EC_030t"]) #change variable name to predicted ec layer for each depth
# names(EC_030)<-c("EC") 
# EC1 <- as(EC_030, "SpatialPixelsDataFrame")
# 
#  
#  PH_030 <- raster(rs_030["PH_030"]) #change variable name to predicted ph layer for each depth
#  names(PH_030) <- c("PH")
#  PH1 <- as(PH_030, "SpatialPixelsDataFrame")
#  
#  ESP_030 <- raster(rs_030["ESP_030t"]) #change variable name to predicted esp layer for each depth
#  names(ESP_030) <- c("ESP")
#  ESP1 <- as(ESP_030, "SpatialPixelsDataFrame")
 
# #NOT SURE EXACTLY WHAT'S HAPPENING HERE; no ECte, PHt, ESPt objects created in script earlier; they might be new objects but we'll need to edit this# 
# #NEED UNCERTAINTY MAPS from the predictions


# setwd("G:/GSP/predictions/acc_unc")
setwd("D:/geodata/project_data/gsp-sas/predictions_v2/accuracy_uncert")
 
 
##bring in uncertainty layers and stack them
unc <- c(EC_unc  = "ec100_uncert_predsd.tif", 
         PH_unc  = "ph_100_uncert_predsd.tif", 
         ESP_unc = "v2_esp100t25_coast_pm_log_nafw_sd.tif"
         )
uncst <- projectRaster(readAll(stack(unc)), crs = "+init=epsg:5070", progress = "text")
uncst$EC_unct  <- log(uncst$EC_unc + 0.1)
uncst$ESP_unct <- log(uncst$ESP_unc +0.1)
uncst <- readAll(uncst)

uncst2 <- as(uncst, "SpatialGridDataFrame")
uncst2


#If ECte, PHt, ESPt are predicted values and ECsd, PHsd, and ESPsd are uncertainty standard deviations
#ECte <- raster(rs_030["EC_030"]);
# ECte <- rs_030$EC_030
# ECsd <- uncst$EC_unc
# names(ECsd) <-"ECsd"
# 
# 
# PHde <- rs$PH_030
# PHsd <- uncst$PH_unc 
# names(PHsd) <- "PHsd"
# 
# ESPt  <- rs$ESP_030
# ESPsd <- uncst$ESP_unc
# names(ESPsd) <- "ESPsd"


##obtain sample spatial autocorrelation
###############
library(automap)

EC11 <- as(rs_100["EC_100t"], "SpatialPointsDataFrame")
b1 <- nrow(EC11)
c1 <- trunc(0.01 * b1)
jj1 <- EC11[sample(b1, c1),]
ec_100_vrm <- autofitVariogram(EC_100t ~ 1, jj1)
plot(ec_100_vrm)


PH11 <- as(rs_100["PH_100"], "SpatialPointsDataFrame")
b2 <- nrow(PH11)
c2 <- trunc(0.01 * b2)
jj2 <- PH11[sample(b2, c2),]
ph_100_vrm <- autofitVariogram(PH_100 ~ 1, jj2)
plot(ph_100_vrm)


ESP11 <- as(rs_100["ESP_100t"], "SpatialPointsDataFrame")
b3 <- nrow(ESP11)
c3 <- trunc(0.01 * b3)
jj3 <- ESP11[sample(b3, c3),]
esp_100_vrm <- autofitVariogram(ESP_100t ~ 1, jj3)
plot(esp_100_vrm)



# plot autocorrelation info
library(spup)

## all crm objects have ranges = 20000 like example in manual but that is NOT the range output by the vrm; all ranges have to match for genSample()
plot(ec_030_vrm) # Note the spatial correlation model and the value of Range parameter
acf(EC11$EC_100t) ##Also note the acf0 (at lag 0)
ec_100_crm <- makeCRM(acf0 = 0.9, range = 20000, model = "Sph") #variogram model is actually "Ste" but then it won't plot...?? object 'xlim_factor' not found ; 
plot(ec_030_crm, main = "EC 30cm correlogram")

plot(ph_030_vrm)
acf(PH11$PH_100)
ph_100_crm <- makeCRM(acf0 = 0.9, range = 20000, model = "Sph") 
plot(ph_030_crm, main = "PH 30cm correlogram")

plot(esp_030_vrm)
acf(ESP11$ESP_100t)
esp_100_crm <- makeCRM(acf0 = 0.9, range = 20000, model = "Sph")
plot(esp_030_crm, main = "ESP 30cm correlogram")


# save(ec_100_vrm, ph_100_vrm, esp_100_vrm,
#      ec_100_crm, ph_100_crm, esp_100_crm,
#      file = "gsp_variograms_100.RData"
#      )
load(file = "gsp_variograms_100.RData")

###################


## Develop input marginal and joint multivariate uncertainty models for defining MC models
## ERROR Distribution parameters must be objects of the same class (change ECte to SpatialGridDataFrame?)

##create MC realizations from the distributions
MC <- 100

test_ex <- sampleRandom(rs, size = 10000, na.rm = TRUE, sp = TRUE)@data
cm <- cor(test_ex[c(9, 4, 10)])

sa <- read_sf("D:/geodata/soils/SSURGO_CONUS_FY19.gdb", layer = "SAPOLYGON")
tiles <- st_as_sf(st_make_grid(sa, n = c(6, 6)))
tiles$tile <- 1:nrow(tiles)


# compute MC simulation over tiles

lapply(1:34, function(x) {
    
    cat(as.character(Sys.time()), x, "\n")
    
    tiles_sub    <- tiles[tiles$tile == x, ]
    tiles_sub    <- st_buffer(tiles_sub, 1e4 * 10)
    rs_sub    <- crop(rs, tiles_sub)
    uncst_sub <- crop(uncst, tiles_sub)
    
    EC_UM <- defineUM(distribution = "norm", 
                      distr_param  = c(rs_sub$EC_100t, uncst_sub$EC_unct), 
                      crm          = ec_100_crm, 
                      id           = "EC"
    )
    PH_UM <- defineUM(distribution = "norm",
                      distr_param  = c(rs_sub$PH_100, uncst_sub$PH_unc),
                      crm          = ph_100_crm,
                      id = "PH"
    )
    ESP_UM <- defineUM(distribution = "norm",
                       distr_param  = c(rs_sub$ESP_100t, uncst_sub$ESP_unct),
                       crm          = esp_100_crm,
                       id           = "ESP"
    )
    MC <- 100
    salinityMUM <- defineMUM(UMlist = list(EC_UM, PH_UM, ESP_UM),
                             cormatrix = cm
    )
    input_sample <- genSample(UMobject = salinityMUM, 
                              n = MC, 
                              samplemethod = "ugs", 
                              nmax = 20, 
                              asList = FALSE, 
                              debug.level = -1
    )
    save(input_sample, file = paste0("D:/geodata/project_data/gsp-sas/predictions_v2/accuracy_uncert/input_sample_", x, ".RData"))
})
    

# back transform results

lapply(1:34, function(x){
    
    cat(as.character(Sys.time()), x, "\n")
    
    load(file = paste0("D:/geodata/project_data/gsp-sas/predictions_v2/accuracy_uncert/input_sample_", x, ".RData"))
    
    EC_sample  <- input_sample[[1:100]]
    PH_sample  <- input_sample[[101:200]]
    ESP_sample <- input_sample[[201:300]]
    
    EC_sample  <- exp(EC_sample) - 0.1
    ESP_sample <- exp(ESP_sample) - 0.1
    
    input_sample <- stack(EC_sample, PH_sample, ESP_sample)
    
    save(input_sample, file = paste0("D:/geodata/project_data/gsp-sas/predictions_v2/accuracy_uncert/input_sample_", x, "_bt.RData"))
})


# iterate over tiles and compute input sample statistics

lapply(1:34, function(x){
    
    cat(as.character(Sys.time()), x, "\n")
    
    load(file = paste0("D:/geodata/project_data/gsp-sas/predictions_v2/accuracy_uncert/input_sample_", x, "_bt.RData"))
    
    EC_sample  <- input_sample[[1:100]]
    PH_sample  <- input_sample[[101:200]]
    ESP_sample <- input_sample[[201:300]]
    
    EC_sample_mean  <- mean(EC_sample)
    PH_sample_mean  <- mean(PH_sample)
    ESP_sample_mean <- mean(ESP_sample)
    
    EC_sample_sd  <- calc(EC_sample,  fun = sd)
    PH_sample_sd  <- calc(PH_sample,   fun = sd)
    ESP_sample_sd <- calc(ESP_sample, fun = sd)

    st <- stack(EC_sample_mean, PH_sample_mean, ESP_sample_mean,
                EC_sample_sd,   PH_sample_sd,   ESP_sample_mean
                )
    
    writeRaster(st, filename = paste0("D:/geodata/project_data/gsp-sas/predictions_V2/accuracy_uncert/gsp_sample_st_", x, ".tif"), overwrite = TRUE)
    })



# load MC iterations
mc_st <- lapply(1:34, function(x) {
    stack(paste0("D:/geodata/project_data/gsp-sas/predictions_v2/accuracy_uncert/gsp_sample_st_", x, ".tif"))
    })


# mosaic and calculate averages
ec_avg_l <- lapply(mc_st, function(x) x[[1]])
ec_avg_l <- c(ec_avg_l, fun = mean, na.rm = TRUE, progress = "text")
EC_sample_mean <- do.call(mosaic, ec_avg_l)
writeRaster(EC_sample_mean, filename = "D:/geodata/project_data/gsp-sas/predictions_v2/ec_100_mc_avg.tif", overwrite = TRUE)

ph_avg_l <- lapply(mc_st, function(x) x[[2]])
ph_avg_l <- c(ph_avg_l, fun = mean, na.rm = TRUE, progress = "text")
PH_sample_mean <- do.call(mosaic, ph_avg_l)
writeRaster(PH_sample_mean, filename = "D:/geodata/project_data/gsp-sas/predictions_v2/ph_100_mc_avg.tif", overwrite = TRUE)

esp_avg_l <- lapply(mc_st, function(x) x[[3]])
esp_avg_l <- c(esp_avg_l, fun = mean, na.rm = TRUE, progress = "text")
ESP_sample_mean <- do.call(mosaic, esp_avg_l)
writeRaster(ESP_sample_mean, filename = "D:/geodata/project_data/gsp-sas/predictions_v2/esp_100_mc_avg.tif", overwrite = TRUE)


# mosaic and calculate standard deviations
ec_sd_l <- lapply(mc_st, function(x) x[[4]])
ec_sd_l <- c(ec_sd_l, fun = mean, na.rm = TRUE, progress = "text")
ec_sd_r <- do.call(mosaic, ec_sd_l)
writeRaster(ec_sd_r, filename = "D:/geodata/project_data/gsp-sas/predictions_v2/ec_100_mc_sd.tif", overwrite = TRUE)

ph_sd_l <- lapply(mc_st, function(x) x[[5]])
ph_sd_l <- c(ph_sd_l, fun = mean, na.rm = TRUE, progress = "text")
ph_sd_r <- do.call(mosaic, ph_sd_l)
writeRaster(ph_sd_r, filename = "D:/geodata/project_data/gsp-sas/predictions_v2/ph_100_mc_sd.tif", overwrite = TRUE)

esp_sd_l <- lapply(mc_st, function(x) x[[6]])
esp_sd_l <- c(esp_sd_l, fun = mean, na.rm = TRUE, progress = "text")
esp_sd_r <- do.call(mosaic, esp_sd_l)
writeRaster(esp_sd_r, filename = "D:/geodata/project_data/gsp-sas/predictions_v2/esp_100_mc_sd.tif", overwrite = TRUE)



# plot the realizations
library(tmap)

tm_shape(EC_sample_mean) + 
    tm_raster(palette = rev(viridis::viridis(5)), breaks = c(0, 2, 4, 8, 16, 300), title = "EC") + 
    tm_layout(main.title = "Mean of EC realizations", legend.outside = TRUE)
tm_shape(PH_sample_mean) + 
    tm_raster(palette = rev(viridis::viridis(11)), breaks = c(0, 3.5, 4.5, 5.1, 5.6, 6.1, 6.6, 7.4, 7.9, 8.5, 9, 14), title = "pH") + 
    tm_layout(main.title = "Mean of pH realizations", legend.outside = TRUE)
tm_shape(ESP_sample_mean) + 
    tm_raster(palette = rev(viridis::viridis(4)), breaks = round(quantile(ESP_sample_mean, p = c(0, 0.1, 0.5, 0.9, 1)), 2), title = "ESP") + 
    tm_layout(main.title = "Mean of ESP realizations", legend.outside = TRUE)



##uncertainty propagation through the classification model
lapply(1:34, function(x) {
    
    cat(as.character(Sys.time()), "propagating error for", x, "\n")
    
    load(file = paste0("D:/geodata/project_data/gsp-sas/predictions_v2/accuracy_uncert/input_sample_", x, "_bt.RData"))
    
    Salinity_model_raster <- function(EC1, PH1, ESP1){
        ww = EC1
        ww = raster(ww)
        ww$salt = saltSeverity(values(EC1),values(PH1),values(ESP1),"FAO")
        ww = ww$salt; names(ww) = c("salt")
        ww
    }
    
    v <- list()
    v[[1]] <- map(1:100, function(x){input_sample[[x]]})
    v[[2]] <- map(101:200, function(x){input_sample[[x]]})
    v[[3]] <- map(201:300, function(x){input_sample[[x]]})
    input_sample <- v
    salinity_sample <- propagate(realizations = input_sample,
                                 model        = Salinity_model_raster,
                                 n            = MC
                                 )
    
    #determine uncertainty of final classified map
    salinity_sample <- raster::stack(salinity_sample)
    salinity_freq = modal(salinity_sample, freq = TRUE)
    salinity_prop = salinity_freq / 100
    salinity_SErr = sqrt(salinity_prop * (1 - salinity_prop) / 100)
    CL = 0.95
    z_star = round(qnorm((1 - CL) / 2, lower.tail = FALSE), digits = 2)
    salinity_MErr = z_star * salinity_SErr
    
    #write final output to raster
    writeRaster(salinity_MErr, filename = paste0("D:/geodata/project_data/gsp-sas/predictions_v2/accuracy_uncert/Salinity_ME_", x, ".tif"), format = "GTiff", progress = "text")
})


mc_err <- lapply(1:34, function(x) {
    stack(paste0("D:/geodata/project_data/gsp-sas/predictions_v2/accuracy_uncert/Salinity_ME_", x, ".tif"))
})
mc_err_l <- c(mc_err, fun = mean, na.rm = TRUE, progress = "text")
mc_err_final <- do.call(mosaic, mc_err_l)
writeRaster(mc_err_final, filename = "D:/geodata/project_data/gsp-sas/predictions_v2/Salinity_ME_100.tif", overwrite = TRUE, progress = "text")

