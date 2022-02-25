load(file = "C:/Users/stephen.roecker/OneDrive - USDA/projects/gsp-sas/LDM-compact_20200709.RData")
data <- s_mps_sf[1:24]
names(data)[10:24] <- paste(rep(c("pH", "EC", "ESP", "SAR", "TSS"), each = 3), c("0.30", "30.100", "100.150"), "obs", sep = "_")


dsm_fn <- c(
  ec_030  = "USA840_CONUS_SalinityMap030.tif",
  ec_100  = "USA840_CONUS_SalinityMap30100.tif",
  ph_030  = "USA840_CONUS_pHMap030.tif",
  ph_100  = "USA840_CONUS_pHMap30100.tif",
  esp_030 = "USA840_CONUS_ESPMap030.tif",
  esp_100 = "USA840_CONUS_ESPMap30100.tif"
)
dsm_rs <- stack(file.path("D:/geodata/project_data/gsp-sas/deliverables/maps/CONUS", dsm_fn))
names(dsm_rs) <- names(dsm_fn)
gsp <- extract(dsm_rs, as(data, "Spatial"), df = TRUE, sp = TRUE)@data
names(gsp)[25:30] <- paste(rep(c("EC", "pH", "ESP"), each = 2), c("0.30", "30.100"), "pred", sep = "_")

sas030_obs <- aqp:::.rank_salts(EC = gsp$EC_0.30_obs,   pH = gsp$pH_0.30_obs,   ESP = gsp$ESP_0.30_obs)
sas100_obs <- aqp:::.rank_salts(EC = gsp$EC_30.100_obs, pH = gsp$pH_30.100_obs, ESP = gsp$ESP_30.100_obs)

sas030_pred <- aqp:::.rank_salts(EC = gsp$EC_0.30_pred,   pH = gsp$pH_0.30_pred,   ESP = gsp$ESP_0.30_pred)
sas100_pred <- aqp:::.rank_salts(EC = gsp$EC_30.100_pred, pH = gsp$pH_30.100_pred, ESP = gsp$ESP_30.100_pred)

sas <- cbind(st_coordinates(s_mps_sf), gsp, sas030_obs, sas100_obs, sas030_pred, sas100_pred)

write.csv(sas, file = "C:/workspace2/github/ncss-tech/stats_for_soil_survey/data/gsp_sas.csv", row.names = FALSE)


