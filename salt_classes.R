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

###ADD READ IN OF EC, PH, ESP PREDICTION LAYERS###

###generate salt classes###

##run saltiness and salt severity for both 0-30cm and 30-100cm
#saltiness
salty <- saltClass(EC_depth, pH_depth, ESP_depth, "FAO") #replace ec, ph, esp variables with predicted layers for each depth interval
saltiness <- classCode(salty, "saltclass") #add class names
plot(saltiness)

#salt severity
salt_affected <- saltSeverity(EC_depth, pH_depth, ESP_depth, "FAO") #replace ec, ph, esp variables with predicted layers for each depth interval
saltaffectedness <- classCode(salt_affected, "saltseverity") #add class names
plot(saltaffectedness)

###export salt class maps###

#run for both 0-30cm and 30-100cm 
saltclasses <- as.numeric(saltaffectedness) #convert factors to numeric
salinity_LUT30 <- classLUT(saltaffectedness, "saltseverity") #change object name for depth as needed
writeGDAL(saltclasses, drivername = "GTIFF", "Top0_30saltaffected.tif")
write.table(salinity_LUT30, file = "saltaffected_LUT30.txt", row.names = F)

###do we need to also export saltiness?????###


########################################################################
###accuracy assessment of salt class maps; p88###

##import and classify validation points for salt classes; these are OBSERVED values
soilv <- readOGR("filepath", "validation dataset") #change file path and validation dataset name 
soilv <- subset(soilv, soilv$horizon==1) #change name of horizon in validation data for 0-30cm and 30-100 cm name 
soilv$saltaffected1 <- saltSeverity(soilv$ec, soilv$ph, soilv$esp, "FAO") #change ec, ph, esp to match validation data column names for properties
summary(soilv$saltaffected1)
soilv$saltaffectedness1 <- classCode(soilv$salt_affected1, "saltseverity") #add class names
summary(soilv$saltaffectedness1)

##extract predicted values to validation points; these are PREDICTED values
soilv <- subset(soilv, !is.na(soilv$saltaffectedness1))
saltclasses.ovv <- over(soilv, saltclasses) #extract class predictions to validation points
soilv$salt_affected <- saltclasses.ovv$salt_affected #take predicted value for salt affected and add to the validataion points
soilv$saltaffectedness <- saltclasses.ovv$saltaffectedness #take predicted value for saltaffectedness and add to the validation points
#check summary of extracted classified pixels
summary(soilv$salt_affected)
summary(soilv$saltaffectedness)

##generate confusion matrix and Kappa
agreementplot(confusion(soilv$salt_affected, soilv$salt_affected1),
              main = "Accuracy assessment",xlab = "Class codes in holdout samples",
              ylab = "Class codes in map")
Kappa(confusion(soilv$salt_affected, soilv$salt_affected1))


########################################################################
###uncertainty of salt class maps with monte-carlo simulations; p90###

##convert input layers to raster files (ours may already be raster objects?)
EC <- raster(EC_depth) #change variable name to predicted ec layer for each depth
names(EC) <- c("EC")
EC1 <- as(EC, "SpatialPixelsDataFrame") 

PH <- raster(pH_depth) #change variable name to predicted ph layer for each depth
names(PH) <- c("PH")
PH1 <- as(PH, "SpatialPixelsDataFrame")

ESP <- raster(esp_depth) #change variable name to predicted esp layer for each depth
names(ESP) <- c("ESP")
ESP1 <- as(ESP, "SpatialPixelsDataFrame")

#NOT SURE EXACTLY WHAT'S HAPPENING HERE; no ECte, PHt, ESPt objects created in script earlier; they might be new objects but we'll need to edit this#
ECte <- raster(predictors["ECte"])
ECsd <- pred_uncerta$pred_sd
names(ECsd) <- c("ECsd")

PHde <- raster(predictors["PHt"])
PHsd <- pred_uncertb$pred_sd
names(PHsd) <- c("PHsd")

ESPt <- raster(predictors["ESPt"])
ESPsd <- pred_uncertc$pred_sd 
names(ESPsd) <- c("ESPsd")

##obtain sample spatial autocorrelation
b <- nrow(EC1)
c <- trunc(0.01*b)
jj <- EC1[sample(b,c),]
vrm <- autofitVariogram(EC~1,jj)

#plot correlation info
plot(vrm)#Note the spatial correlation model and the value of Range parameter
acf((EC1$EC)) ##Also note the acf0 (at lag 0)
EC_crm <- makeCRM(acf0 = 0.85, range = 20000, model = "Sph")
plot(EC_crm, main = "EC correlogram")

##develop input marginal and joint multivariate uncertainty models for defining MC models
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