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
salty <- saltClass(EC_depth, pH_depth, ESP_depth, "FAO") #replace ec, ph, esp with predicted layers for each depth interval
saltiness <- classCode(salty, "saltclass") #add class names
plot(saltiness)

#salt severity
salt_affected <- saltSeverity(EC_depth, pH_depth, ESP_depth, "FAO") #replace ec, ph, esp with predicted layers for each depth interval
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
###accuracy assessment of salt class maps###

##import and classify validation points for salt classes; these are OBSERVED values
soilv <- readOGR("filepath", "validation dataset") #change file path and validation dataset name 
soilv <- subset(soilv, soilv$horizon==1) #change name of horizon in validation data for 0-30cm and 30-100 cm name 
soilv$saltaffected1 <- saltSeverity(soilv$ec, soilv$ph, soilv$esp, "FAO") #change ec, ph, esp to match validation data names for properties
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