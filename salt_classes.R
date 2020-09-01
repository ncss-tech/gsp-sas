###GSP salt affected soil maps
###salt affected soil classes
###August 2020 JMP, SKB
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

setwd("G:/GSP/predictions/predictions")

###Read in prediction layers###

f <- c("ec_100_bt_nomlra_wstatspm.tif", "cv_esp_100_nomlra_statspm.tif", "ph_100_qrf_mlra_pm_NULL_final2.tif") #change file names as needed
#f <- c("ec_030_bt_nomlra_wstatspm.tif", "cv_esp_030_nomlra_statspm.tif", "ph_030_qrf_mlra_pm_NULL_final2.tif") #change file names as needed
rs <- stack(f)

rs
###back transform ec and esp if needed###

##rs$ec_t <- calc(rs$ec_030_t, function(x) exp(x) - 0.1) 
rs$esp_t <- calc(rs$cv_esp_100_nomlra_statspm, function(x) exp(x) - 0.1)
rs #check which columns to keep in next line
rs <- rs[[c(1,3:4)]] #only use back-transformed data
rs2 <- as(rs, "SpatialGridDataFrame")
names(rs2) <- c("EC_100", "pH_100", "ESP_100")


##run saltiness and salt severity for both 0-30cm and 30-100cm
#saltiness
rs2$salty <- saltClass(rs2$EC_100, rs2$pH_100, rs2$ESP_100,"FAO") #replace ec, ph, esp variables with predicted layers for each depth interval
rs2$saltiness <- classCode(rs2$salty, "saltclass") #add class names
plot(rs2$saltiness)
spplot(rs2["saltiness"])

#salt severity
rs2$salt_affected <- saltSeverity(rs2$EC_100, rs2$pH_100, rs2$ESP_100, "FAO") #replace ec, ph, esp variables with predicted layers for each depth interval
rs2$saltaffectedness <- classCode(rs2$salt_affected, "saltseverity") #add class names
plot(rs2$saltaffectedness)
spplot(rs2["saltaffectedness"])

♥###export salt class maps###
#run for both 0-30cm and 30-100cm
#saltiness
rs2$saltyclasses <- as.numeric(rs2$saltiness)
saltiness_LUT100 <- classLUT(rs2["saltiness"], "saltclass") ##WHY WON'T THIS WORK "object 'LUT' not found" 
writeGDAL(rs2["saltyclasses"], drivername = "GTiff", "Mid100_saltiness.tif")
write.table(saltiness_LUT100, file = "saltiness_LUT100.txt", row.names = F)

#run for both 0-30cm and 30-100cm 
rs2$saltclasses <- as.numeric(rs2$saltaffectedness) #convert factors to numeric
salinity_LUT100 <- classLUT(rs2["saltaffectedness"], "saltseverity") #change object name for depth as needed
writeGDAL(rs2["saltclasses"], drivername = "GTiff", "Mid100_saltaffected.tif")
write.table(salinity_LUT100, file = "saltaffected_LUT100.txt", row.names = F)

###do we need to also export saltiness?????###


########################################################################
###accuracy assessment of salt class maps; p88###

##import and classify validation points for salt classes; these are OBSERVED values
setwd("G:/GSP/pointdata")
getwd()

#soilv <- readOGR("filepath", "validation dataset") #change file path and validation dataset name 
#load training_data.RData

load(file = "LDM-compact_20200709.RData")
#also load training data into environment, test is test points

soilv <- test

soilv$saltaffected1 <- saltSeverity(soilv$ec_ptf.030_100_cm, soilv$ph_ptf.030_100_cm, soilv$esp.030_100_cm, "FAO") 
summary(soilv$saltaffected1)

soilv$saltaffectedness1 <- classCode(soilv$saltaffected1, "saltseverity")
summary(soilv$saltaffectedness1)

##extract the salt classes from the map using the validation samples
soilv <- subset(soilv, !is.na(soilv$saltaffectedness1)) #trim test points down to only ones with relevant data
summary(soilv)

test2 <- cbind(s_mps_sf["pedon_key"], sf::st_coordinates(s_mps_sf)) #extract coordinates from lab data
soilv <- merge(x = test2, y = soilv, by = "pedon_key", all.y = TRUE) #merge coordinates with test data
summary(soilv)

rs2.ov <- over(as(soilv, "Spatial"), rs2) #combine predictions with test data

saltclasses.ovv <- rs2.ov

soilv$salt_affected <- saltclasses.ovv$salt_affected #take predicted value for salt affected and add to the validataion points
soilv$saltaffectedness <- saltclasses.ovv$saltaffectedness #take predicted value for saltaffectedness and add to the validation points
#check summary of extracted classified pixels
summary(soilv$salt_affected)
summary(soilv$saltaffectedness)



##generate confusion matrix and Kappa
#####gives error - all arguments must have the same length...but they have the same length so??????
agreementplot(confusion(soilv$salt_affected, soilv$salt_affected1),
              main = "Accuracy assessment",xlab = "Class codes in holdout samples",
              ylab = "Class codes in map")
Kappa(confusion(soilv$salt_affected, soilv$salt_affected1))


########################################################################
###uncertainty of salt class maps with monte-carlo simulations; p90###

##convert input layers to raster files (ours may already be raster objects?)

EC <- raster(rs2["EC_100"]) #change variable name to predicted ec layer for each depth
names(EC) <- c("EC")
EC1 <- as(EC, "SpatialPixelsDataFrame") 

PH <- raster(rs2["pH_100"]) #change variable name to predicted ph layer for each depth
names(PH) <- c("PH")
PH1 <- as(PH, "SpatialPixelsDataFrame")

ESP <- raster(rs2["ESP_100"]) #change variable name to predicted esp layer for each depth
names(ESP) <- c("ESP")
ESP1 <- as(ESP, "SpatialPixelsDataFrame")

#NOT SURE EXACTLY WHAT'S HAPPENING HERE; no ECte, PHt, ESPt objects created in script earlier; they might be new objects but we'll need to edit this# 
#NEED UNCERTAINTY MAPS from the predictions

##bring in uncertainty layers and stack them
unc <- c("ec_100_preduncert.tif", "ut_esp_100_t25_nomlra_prun_width.tif")
uncst <- stack(unc)
uncst2 <- as(uncst, "SpatialGridDataFrame")
names(uncst2) <- c("EC_unc", "ESP_unc")

##or bring in uncertainty layers one at a time
#ECunc <- c("ec_100_preduncert.tif")
#ECunc <- raster(ECunc)
#ECsd <- ECunc
#names(ECsd) <- c("ECsd")

#If ECte, PHt, ESPt are supposed to be observed values
#and if ECsd, PHsd, and ESPsd are values from uncertainty maps
ECte <- raster(soilv["ec_ptf.030_100_cm"])
ECsd <- uncst2$EC_unc
names(ECsd) <- c("ECsd")

PHde <- raster(soilv["ph_ptf.030_100_cm"])
PHsd <- uncst2$pH_unc
names(PHsd) <- c("PHsd")

ESPt <- raster(soilv["esp.030_100_cm"])
ESPsd <- uncst2$ESP_unc
names(ESPsd) <- c("ESPsd")

##obtain sample spatial autocorrelation
b <- nrow(EC1)
c <- trunc(0.01*b)
jj <- EC1[sample(b,c),]
vrm <- autofitVariogram(EC~1,jj) #cannot deal with non-square cells

#plot correlation info
plot(vrm)#Note the spatial correlation model and the value of Range parameter
acf((EC1$EC)) ##Also note the acf0 (at lag 0)
EC_crm <- makeCRM(acf0 = 0.85, range = 20000, model = "Sph")
plot(EC_crm, main = "EC correlogram")

##develop input marginal and joint multivariate uncertainty models for defining MC models
##takes 2.5+ hours to run each
EC_UM <- defineUM(distribution = "norm",distr_param = c(ECte,ECsd),crm = EC_crm,id = "EC")
PH_UM <- defineUM(distribution = "norm",distr_param = c(PHde,PHsd),crm = PH_crm,id = "PH")
ESP_UM <- defineUM(distribution = "norm",distr_param = c(ESPt,ESPsd),crm = ESP_crm,id = "ESP")
class(EC_UM)
class(PH_UM)
class(ESP_UM)

#get the correlation values and use them in defining the Monte Carlo Uncertainty Mode (MUM)
cor(values(ECte),values(PHde))
cor(values(ECte),values(ESPt))
cor(values(PHde),values(ESPt))

salinityMUM <- defineMUM(UMlist = list(EC_UM, PH_UM, ESP_UM), cormatrix = matrix(c(1, cor(values(ECte),values(PHde)), cor(values(ECte),values(ESPt)), cor(values(ECte), values(PHde)), 1, cor(values(PHde), values(ESPt)), cor(values(ECte), values(ESPt)), cor(values(PHde), values(ESPt)),1), nrow = 3, ncol = 3))
class(salinityMUM)

##create MC realizations from the distributions
MC <- 100
input_sample <- genSample(UMobject = salinityMUM, n = MC, samplemethod = "ugs", nmax = 20, asList = FALSE)

#compute input sample statistics
EC_sample <- input_sample[[1:MC]]
PH_sample <- input_sample[[(MC+1):(2*MC)]]
ESP_sample <- input_sample[[(2*MC+1):(3*MC)]]
EC_sample_mean <- mean(EC_sample)
PH_sample_mean <- mean(PH_sample)
ESP_sample_mean <- mean(ESP_sample)
EC_sample_sd <- calc(EC_sample, fun = sd)
PH_sample_sd <- calc(PH_sample, fun = sd)
ESP_sample_sd <- calc(ESP_sample, fun = sd)

#plot the realizations
par(mfrow=c(2,2),mar = c(1, 1, 2, 2), mgp = c(1.7, 0.5, 0), oma = c(0, 0, 0, 1), + las = 1, cex.main = 1, tcl = -0.2, cex.axis = 0.8, cex.lab = 0.8)
plot(EC_sample_mean, main = "Mean of ECt realizations", xaxt = "n", yaxt = "n")
plot(PH_sample_mean, main = "Mean of PHt realizations", xaxt = "n", yaxt = "n")
plot(ESP_sample_mean, main = "Mean of ESPt realizations", xaxt = "n", yaxt = "n")

##uncertainty propagation through the classification model
Salinity_model_raster <- function (EC1,PH1,ESP1){
    ww=EC1
    ww=raster(ww)
    ww$salt=saltSeverity(values(EC1),values(PH1),values(ESP1),"FAO")
    ww=ww$salt; names(ww)=c("salt")
    ww
   }

v <- list()
v[[1]] <- map(1:100, function(x){input_sample[[x]]})
v[[2]] <- map(101:200, function(x){input_sample[[x]]})
v[[3]] <- map(201:300, function(x){input_sample[[x]]})
input_sample <- v
salinity_sample <- propagate(realizations=input_sample,model=Salinity_model_raster,n=MC)

#determine uncertainty of final classified map
samplelist <- list()
samplelist [[1]] = map(1:100, function(x){input_sample[[x]]})
samplelist [[2]] = map(101:200, function(x){input_sample[[x]]})
samplelist [[3]] = map(201:300, function(x){input_sample[[x]]})
input_sample = samplelist
salinity_sample = propagate(realizations = input_sample, model = Salinity_model_raster, n = MC)
salinity_sample <- raster::stack(salinity_sample)
names(salinity_sample) <- paste("salt.", c(1:nlayers(salinity_sample)), sep = "")
salinity_freq = modal(salinity_sample, freq=TRUE)
salinity_prop = salinity_fre/100
salinity_SErr = sqrt(salinity_prop*(1-salinity_prop)/100)
CL = 0.95
z_star = round(qnorm((1-CL)/2,lower.tail=F),digits = 2)
salinity_MErr = z_star*salinity_SErr

#write final output to raster
writeRaster(salinity_MErr,filename="Salinity_ME.tif",format="GTiff")
