
library(aqp)
library(soilDB)
library(sf)
library(terra)


# load gNATSGO ----
fp  <- file.path("D:/geodata/soils/gNATSGO_CONUS_Oct2022")
dsn <- file.path(fp, "gNATSGO_CONUS_Oct2022.gdb")

lyr <- st_layers(dsn)

f <- fetchGDB(dsn, childs = FALSE)
mu <- get_mapunit_from_GDB(dsn, stats = FALSE, stringsAsFactors = TRUE)
# save(f, mu, file = "gnatsgo_Oct22.RData")
load(file = "gnatsgo_Oct22.RData")

gnatsgo_r <- rast(file.path(fp, "gNATSGO-mukey.tif"))


h_avg <- horizons(f) |>
  within({
    hzdep_thk =  (hzdepb_r - hzdept_r) / 2
    h = 1/10^-ph1to1h2o_r
  }) |>
  segment(intervals = c(0, 30, 100), hzdepcols = c("hzdept_r", "hzdepb_r")) |>
  subset(!is.na(cokey), select = c(cokey, segment_id, ec_r, sar_r, ph1to1h2o_r, h, hzdep_thk))


co_avg <- collapse::collap(
  h_avg, 
  by   = list(cokey = h_avg$cokey, segment_id = h_avg$segment_id),
  cols = 3:6,
  keep.w = FALSE,
  FUN  = collapse::fmean, 
  w    = h_avg$hzdep_thk,
  na.rm = TRUE
  )

# co_avg2 <- h_avg %>%
#   group_by(cokey, segment_id) %>%
#   summarize(across(c(ec_r, sar_r, ph1to1h2o_r, h), ~ weighted.mean(.x, w = hzdep_thk, na.rm = TRUE))) %>%
#   ungroup() %>%
#   as.data.frame()
# 
# co_avg[is.na(co_avg)] <- 0
# co_avg1 <- co_avg
# co_avg1[3:6] <- lapply(co_avg[3:6], round, 1) 
# co_avg2[is.na(co_avg2)] <- 0
# co_avg21 <- co_avg2
# co_avg21[3:6] <- lapply(co_avg21[3:6], round, 1) 
# 
# sum(co_avg == co_avg2)
# nrow(co_avg) * ncol(co_avg)
# View(co_avg[idx])


co_avg2 <- reshape(
  data      = co_avg, 
  direction = "wide",
  idvar     = "cokey",
  timevar   = "segment_id",
  v.names   = c("ec_r", "sar_r", "ph1to1h2o_r", "h")
)


mu_avg <- merge(site(f), co_avg2, by = "cokey", all.x = TRUE)
mu_avg$statsgo <- mu_avg$mukey %in% mu$mukey[mu$areasymbol == "US"]
idx <- grep("cokey|mukey|ec_|sar_|ph1|^h.0|statsgo", names(mu_avg))
mu_avg <- mu_avg[idx]


mu_avg <- collapse::collap(
  mu_avg,
  by = list(mukey = mu_avg$mukey),
  cols = 3:10,
  keep.w = FALSE,
  FUN = collapse::fmean,
  w   = mu_avg$comppct_r,
  na.rm = TRUE
)
mu_avg <- within(mu_avg, {
  ph_030 = log10(`h.000-030`)
  ph_100 = log10(`h.030-100`)
  })

# saveRDS(mu_avg, "mu_avg_sas.rds")
mu_avg <- readRDS(file.path("C:/workspace2", "mu_avg_sas.rds"))


esp <- function(SAR) {
  (100 * (-0.0126 + 0.01475 * SAR)) / 
  (1   + (-0.0126 + 0.01475 * SAR))
}


mu_avg <- within(mu_avg, {
  ph_030[is.infinite(ph_030)] <- 0
  ph_100[is.infinite(ph_100)] <- 0
  
  ph_h2o     = ph_030
  # ph_030_sat = predict(ph_lm, data.frame(ph_h2o))
  ph_h2o     = ph_100
  # ph_100_sat = predict(ph_lm, data.frame(ph_h2o))
                        
  ESP_030 = esp(`sar_r.000-030`)
  ESP_100 = esp(`sar_r.030-100`)
  
  SAS_030 = allocate(EC = `ec_r.000-030`, pH = ph_030, ESP = ESP_030)
  SAS_100 = allocate(EC = `ec_r.030-100`, pH = ph_100, ESP = ESP_100)
  
  SAS_030x = classCode(saltSeverity(ph = ph_030, ec = `ec_r.000-030`, esp = ESP_030), "saltseverity")
  SAS_100x = classCode(saltSeverity(ph = ph_030, ec = `ec_r.030-100`, esp = ESP_100), "saltseverity")
})


idx <- !complete.cases(mu_avg$`ec_r.000-030`, mu_avg$ph_030, mu_avg$ESP_030)
test <- mu_avg[is.na(mu_avg$SAS_030), ]
table(idx)
table(!is.na(test$SAS_030))

ggplot(data.frame(ph1 = mu_avg$`ph1to1h2o_r.000-030`, ph2 = mu_avg$ph_030), aes(ph1, ph2)) + geom_hex(bins = 50) +  
  scale_fill_binned(colours = viridis::viridis(50, option = "D"), )


table(mu_avg$SAS_030, mu_avg$SAS_030x)
test <- subset(mu_avg, SAS_030x == "Saline-sodic" & SAS_030 == "slightly sodic")


dat <- with(mu_avg, data.frame(
  mukey = as.integer(mukey), 
  SAS_030 = as.integer(SAS_030), 
  SAS_100 = as.integer(SAS_100)
))



# rasterize ----
setwd(fp)

r <- gnatsgo_r
res <- sqrt(res(r)[1]^2 * ncell(r) / 50)
test <- rast(res = res, extent = ext(r), crs = crs(r))
values(test) <- 1:ncell(test)
plot(test)

makeTiles(gnatsgo_r, test, filename = "gnatsgo_.tif", overwrite = TRUE)

l <- list.files(path = fp, pattern = "gnatsgo_[0-9]{1,2}.tif")

derat <- function(l, dat, vars) {
  lapply(vars, function(var) {
    lapply(l, function(x) {
      cat("processing ", var, x, as.character(Sys.time()), "\n")
      r <- rast(x)
      val <- as.integer(values(r))
      val2 <- data.table::merge.data.table(
        data.table::data.table(mukey = val), 
        dat[c("mukey", var)], 
        by = "mukey", all.x = TRUE, 
        sort = FALSE
      )
      values(r) <- val2[[var]]
      writeRaster(r, file = gsub(".tif", paste0("_", var, ".tif"), x), overwrite = TRUE)
  })})
}
derat(l, dat, c("SAS_030", "SAS_100"))


# l2 <- list.files(path = fp, pattern = "gnatsgo_[0-9]{1,2}_cat.tif", full.names = TRUE)
# lf <- lapply(l2, rast)
# lf$filename = "test.tif"
# lf$overwrite = TRUE
# do.call("merge", lf)

lf_030 <- list.files(path = fp, pattern = "gnatsgo_[0-9]{1,2}_SAS_030.tif", full.names = TRUE)
test <- sprc(lf_030)
fn <- "gnatsgo_Oct22_SAS_030.tif"
merge(test, filename = fn, datatype = "INT1U")
# if (file.path(fp, fn) |> file.exists()) file.remove(lf_030)


lf_100 <- list.files(path = fp, pattern = "gnatsgo_[0-9]{1,2}_SAS_100.tif", full.names = TRUE)
test <- sprc(lf_100)
fn <- "gnatsgo_Oct22_SAS_100.tif"
merge(test, filename = fn, datatype = "INT1U")
# if (file.path(fp, fn) |> file.exists()) file.remove(lf_100)

fn <- "gnatsgo_Oct22_SAS_030.tif"
test <- rast(fn)
writeRaster(test, "gnatsgo_Oct22_SAS_030_INT1U.tif", datatype = "INT1U")

fn <- "gnatsgo_Oct22_SAS_100.tif"
test <- rast(fn)
writeRaster(test, "gnatsgo_Oct22_SAS_100_INT1U.tif", datatype = "INT1U")


# gdalUtilities::gdalbuildvrt(gdalfile = unlist(l2), output.vrt = "test.vrt", dryrun = TRUE)
# gdalUtilities::gdalwarp(srcfile = "test.vrt", dstfile = "test_vrt.tif")



# GSP maps ----
fp <- "D:/geodata/project_data/gsp-sas/deliverables/maps/CONUS"
l <- list.files(fp, pattern = ".tif$", full.names = TRUE)
r <- rast(c(
  EC_030  = l[[5]], EC_100  = l[[6]], 
  ph_030  = l[[3]], pH_100  = l[[4]],
  ESP_030 = l[[1]], ESP_100 = l[[2]]
  ))

test <- allocate(
  EC = values(r$USA840_CONUS_SalinityMap030), 
  pH = values(r$USA840_CONUS_pHMap030), 
  ESP = values(r$USA840_CONUS_ESPMap030), 
  to = "FAO Salt Severity"
)  

test2 <- r[[1]]
test2[!is.na(test2)] <- 0
values(test2) <- as.integer(test)
writeRaster(test2, filename = "test_dsm.tif")



# accuracy assessment ----
load(file = "C:/Users/stephen.roecker/OneDrive - USDA/projects/gsp-sas/LDM-compact_20200709.RData")

s_mps_sf$SAS <- allocate(
  EC = s_mps_sf$ec_ptf.000_030_cm, 
  pH = s_mps_sf$ph_ptf.000_030_cm, 
  ESP = s_mps_sf$esp.000_030_cm
)
s_mps_sf$SAS_int <- as.integer(s_mps_sf$SAS)
write_sf(s_mps_sf, "test.shp")

vars <- c("test2.tif", "test_dsm.tif")
l  <- file.path(fp, vars)

test2 <- rast(l[1])
test_dsm <- rast(l[2])


ex1 <- extract(test2,    st_transform(s_mps_sf, 5070))
ex2 <- extract(test_dsm, s_mps_sf)
ex <- cbind(st_drop_geometry(s_mps_sf), SAS1 = ex1[, 2], SAS2 = ex2[, 2])
vars <- c("SAS_int", "SAS1", "SAS2")
ex[vars] <- lapply(ex[vars], function(x) factor(x, levels = 1:11))


# SSURGO
caret::confusionMatrix(table(pred = ex$SAS1, obs = ex$SAS_int))$overall |> round(2) |> format(scientific = FALSE)
# DSM
caret::confusionMatrix(table(pred = ex$SAS2, obs = ex$SAS_int))$overall |> round(2) |> format(scientific = FALSE)



# spot check map units
f2 <- f[f$mukey == "2611946"]
h2 <- horizons(f2)
View(h2[grep("cokey|mukey|ec_|sar_|ph1|^h.0|statsgo", names(h2))])


