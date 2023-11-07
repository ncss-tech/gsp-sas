
library(aqp)
library(soilDB)
library(sf)
library(terra)
library(soilassessment)


# load gNATSGO ----
fp_box  <- file.path("C:/Users/stephen.roecker/Box/GSP salinity maps")
fp_geo  <- file.path("D:/geodata/soils/gNATSGO_CONUS_Oct2022")
fp_nlcd <- file.path("D:/geodata/land_use_land_cover/nlcd_2019_land_cover_l48_20210604")
fp_cdls <- file.path("D:/geodata/land_use_land_cover/2022_30m_cdls")
dsn <- file.path(fp_geo, "gNATSGO_CONUS_Oct2022.gdb")

lyr <- st_layers(dsn)

f <- fetchGDB(dsn, childs = FALSE)
mu <- get_mapunit_from_GDB(dsn, stats = FALSE, stringsAsFactors = TRUE)
mu <- subset(mu, !duplicated(mukey))
# save(f, mu, file = "gnatsgo_Oct22.RData")
load(file = "gnatsgo_Oct22.RData")

gnatsgo_r <- rast(file.path(fp_geo, "gNATSGO-mukey.tif"))


gdalUtilities::gdalwarp(
  srcfile = file.path(fp_nlcd, "nlcd_2019_land_cover_l48_20210604.img"),
  dstfile = file.path(fp_nlcd, "nlcd_2019_land_cover_l48_20210604_5070.tif"),
  t_srs = "EPSG:5070",
  r = "near",
  tr = c(30, 30),
  te = ext(gnatsgo_r)[c(1, 3, 2,4)]
)
nlcd_r <- rast(file.path(fp_nlcd, "nlcd_2019_land_cover_l48_20210604_5070.tif"))


gdalUtilities::gdalwarp(
  srcfile = file.path(fp_cdls, "2022_30m_cdls.tif"),
  dstfile = file.path(fp_cdls, "2022_30m_cdls_ext.tif"),
  r = "near",
  te = ext(gnatsgo_r)[c(1, 3, 2,4)]
)
cdls_r <- rast(file.path(fp_cdls, "2022_30m_cdls_ext.tif"))
  
  

# weighted average ----
h_seg <- horizons(f) |>
  within({
    hzdep_thk =  (hzdepb_r - hzdept_r) / 2
    h = 1/10^-ph1to1h2o_r
  }) |>
  segment(intervals = c(0, 30, 100), hzdepcols = c("hzdept_r", "hzdepb_r")) |>
  subset(!is.na(cokey), select = c(cokey, segment_id, ec_r, sar_r, ph1to1h2o_r, h, hzdep_thk))


vars <- c("mukey", "cokey", "comppct_r")
h_seg2 <- merge(h_seg, site(f)[vars], by = "cokey", all.x = TRUE)

co_avg <- collapse::collap(
  h_seg2, 
  by = ec_r + sar_r + ph1to1h2o_r + h ~ mukey + cokey + comppct_r + segment_id,
  keep.w = FALSE,
  FUN  = collapse::fmean, 
  w    =  ~ hzdep_thk,
  na.rm = TRUE
)

mu_avg <- collapse::collap(
  co_avg,
  by = ec_r + sar_r + ph1to1h2o_r +h ~ mukey + segment_id,
  keep.w = FALSE,
  FUN = collapse::fmean,
  w   = ~ comppct_r,
  na.rm = TRUE
)

mu_avg <- reshape(
  data      = mu_avg, 
  direction = "wide",
  idvar     = "mukey",
  timevar   = "segment_id",
  v.names   = c("ec_r", "sar_r", "ph1to1h2o_r", "h")
)

mu_avg <- merge(mu["mukey"], mu_avg, by = "mukey", all.x = TRUE)
mu_avg$statsgo <- mu_avg$mukey %in% mu$mukey[mu$areasymbol == "US"]

mu_avg <- within(mu_avg, {
  ph_030 = log10(`h.000-030`)
  ph_100 = log10(`h.030-100`)
  mukey  = as.integer(mukey)
})

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


s <- strsplit(tolower(mu$muname), " |, ")
idx <- sapply(s, function(x) x[2] == "water" | x[1] == "water")
mu$water <- idx

mu_avg <- merge(mu_avg, mu[c("mukey", "water")], by = "mukey", all.x = TRUE, sort = FALSE)

# saveRDS(mu_avg, file.path(fp_box, "gnatsgo_Oct2022_wavg_sas.rds"))
mu_avg <- readRDS(file.path(fp_box, "gnatsgo_Oct2022_wavg_sas.rds"))



# classify SAS ----
esp <- function(SAR) {
  (100 * (-0.0126 + 0.01475 * SAR)) / 
    (1   + (-0.0126 + 0.01475 * SAR))
}

idx <- mu_avg$water == FALSE

mu_avg <- within(mu_avg, {
  ph_030[is.infinite(ph_030) & idx] <- 0
  ph_100[is.infinite(ph_100) & idx] <- 0
  
  ph_h2o     = ph_030
  # ph_030_sat = predict(ph_lm, data.frame(ph_h2o))
  ph_h2o     = ph_100
  # ph_100_sat = predict(ph_lm, data.frame(ph_h2o))
                        
  
  EC_030 = ifelse(is.na(`ec_r.000-030`) & idx, 0, `ec_r.000-030`)
  EC_100 = ifelse(is.na(`ec_r.030-100`) & idx, 0, `ec_r.030-100`)
  
  SAR_030 = ifelse(is.na(`sar_r.000-030`) & idx, 0, `sar_r.000-030`)
  SAR_100 = ifelse(is.na(`sar_r.030-100`) & idx, 0, `sar_r.030-100`)
  
  ESP_030 = esp(SAR_030)
  ESP_100 = esp(SAR_100)
  
  SAS_030 = allocate(EC = EC_030, pH = ph_030, ESP = ESP_030)
  SAS_100 = allocate(EC = EC_100, pH = ph_100, ESP = ESP_100)
  
  SAS_030x = classCode(saltSeverity(ph = ph_030, ec = EC_030, esp = ESP_030), "saltseverity")
  SAS_100x = classCode(saltSeverity(ph = ph_030, ec = EC_100, esp = ESP_100), "saltseverity")
})



## compare with SDA ----

### get_SDA_project ----
le <- get_legend_from_SDA(WHERE = "areasymbol LIKE 'CA%'")
test <- get_SDA_property(
  "ec_r", 
  method = "Weighted Average", 
  areasymbols = le$areasymbol, 
  top_depth = 0, bottom_depth = 30, 
  include_minors = TRUE, 
  miscellaneous_areas = TRUE
)
test2 <- merge(test, mu_avg, by = "mukey", all.x = TRUE)
vars <- c("mukey", "ec_r", "ec_r.000-030")
test3 <- within(test2[vars], {
  dif = round(ec_r - `ec_r.000-030`, 2)
})


### fetchGDB ----
f2 <- f[f$mukey == "1860756"]
View(subset(horizons(f2), select = c(cokey, ec_r, sar_r, ph1to1h2o_r, hzdept_r)))
h2 <- horizons(f2)

h3 <- h2 %>% 
  select(cokey, ec_r, sar_r, ph1to1h2o_r, hzdept_r, hzdepb_r) %>% 
  mutate(thk = hzdepb_r - hzdept_r) %>% 
  left_join(site(f2), by = "cokey") %>% 
  segment(intervals = c(0, 30), hzdepcols = c("hzdept_r", "hzdepb_r")) %>% 
  group_by(mukey, cokey, comppct_r) %>% 
  summarize(ec = weighted.mean(ec_r, w = thk)) %>%
  group_by(mukey) %>%
  summarize(ec = weighted.mean(ec, w = comppct_r))

idx <- !complete.cases(mu_avg$`ec_r.000-030`, mu_avg$ph_030, mu_avg$ESP_030)
test <- mu_avg[is.na(mu_avg$SAS_030), ]
table(idx)
table(!is.na(test$SAS_030))

ggplot(data.frame(ph1 = mu_avg$`ph1to1h2o_r.000-030`, ph2 = mu_avg$ph_030), aes(ph1, ph2)) + geom_hex(bins = 50) +  
  scale_fill_binned(colours = viridis::viridis(50, option = "D"), )


table(mu_avg$SAS_030, mu_avg$SAS_030x)
test <- subset(mu_avg, SAS_030x == "Saline-sodic" & SAS_030 == "slightly sodic")


## finalize results ----
dat <- with(mu_avg, data.frame(
  mukey = as.integer(mukey), 
  SAS_030 = as.integer(SAS_030), 
  SAS_100 = as.integer(SAS_100)
))



# rasterize ----

r <- gnatsgo_r
res <- sqrt(res(r)[1]^2 * ncell(r) / 100)
test <- rast(res = res, extent = ext(r), crs = crs(r))
values(test) <- 1:ncell(test)
plot(test)

makeTiles(gnatsgo_r, test, filename = file.path(fp_geo,  "gnatsgo_.tif"), datatype = "INT4U", overwrite = TRUE)
makeTiles(nlcd_r,    test, filename = file.path(fp_nlcd, "nlcd2020_.tif"), overwrite = TRUE)

l <- list.files(path = fp_geo, pattern = "gnatsgo_[0-9]{1,2}.tif", full.names = TRUE)

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


lf_030 <- list.files(path = fp_geo, pattern = "gnatsgo_[0-9]{1,2}_SAS_030.tif", full.names = TRUE)
test <- sprc(lf_030)
fn <- file.path(fp_geo, "gnatsgo_Oct22_SAS_030.tif")
merge(test, filename = fn, datatype = "INT1U", overwrite = TRUE)
# if (file.path(fp_geo, fn) |> file.exists()) file.remove(lf_030)
gnat_030_r <- rast(fn)
# makeTiles(gnat_030_r, test, filename = file.path(fp_geo, "gnatsgo_SAS_030_.tif"), overwrite = TRUE)


lf_100 <- list.files(path = fp_geo, pattern = "gnatsgo_[0-9]{1,2}_SAS_100.tif", full.names = TRUE)
test <- sprc(lf_100)
fn <- file.path(fp_geo, "gnatsgo_Oct22_SAS_100.tif")
merge(test, filename = fn, datatype = "INT1U", overwrite = TRUE)
# if (file.path(fp_geo, fn) |> file.exists()) file.remove(lf_100)
gnat_100_r <- rast(fn)
# makeTiles(gnat_100_r, test, filename = file.path(fp_geo, "gnatsgo_SAS_100_.tif"), overwrite = TRUE)

# l2 <- list.files(path = fp_geo, pattern = "gnatsgo_[0-9]{1,2}_cat.tif", full.names = TRUE)
# lf <- lapply(l2, rast)
# lf$filename = "test.tif"
# lf$overwrite = TRUE
# do.call("merge", lf)

# gdalUtilities::gdalbuildvrt(gdalfile = unlist(l2), output.vrt = "test.vrt", dryrun = TRUE)
# gdalUtilities::gdalwarp(srcfile = "test.vrt", dstfile = "test_vrt.tif")



# GSP maps ----
fp_gsp <- "D:/geodata/project_data/gsp-sas/deliverables/maps/CONUS"
l <- list.files(fp_gsp, pattern = ".tif$", full.names = TRUE)
r <- rast(c(
  EC_030  = l[[5]], EC_100  = l[[6]], 
  ph_030  = l[[3]], pH_100  = l[[4]],
  ESP_030 = l[[1]], ESP_100 = l[[2]]
  ))

sas_030 <- allocate(
  EC = values(r$USA840_CONUS_SalinityMap030), 
  pH = values(r$USA840_CONUS_pHMap030), 
  ESP = values(r$USA840_CONUS_ESPMap030), 
  to = "FAO Salt Severity"
)  

sas_100 <- allocate(
  EC = values(r$USA840_CONUS_SalinityMap30100), 
  pH = values(r$USA840_CONUS_pHMap30100), 
  ESP = values(r$USA840_CONUS_ESPMap30100), 
  to = "FAO Salt Severity"
)  

sas_030_r <- r[[1]]
sas_030_r[!is.na(sas_030_r)] <- 0
values(sas_030_r) <- as.integer(sas_030)
# writeRaster(sas_030_r, filename = file.path(fp_box, "ssurgo_comparison/gsp_sas_030.tif"), datatype = "INT1U")
sas_030_r <- rast(file.path(fp_box, "ssurgo_comparison/gsp_sas_030.tif"))

sas_100_r <- r[[1]]
sas_100_r[!is.na(sas_100_r)] <- 0
values(sas_100_r) <- as.integer(sas_100)
# writeRaster(sas_100_r, filename = file.path(fp_box, "ssurgo_comparison/gsp_sas_100.tif"), datatype = "INT1U")
sas_100_r <- rast(file.path(fp_box, "ssurgo_comparison/gsp_sas_100.tif"))



# accuracy assessment ----
load(file = "C:/Users/stephen.roecker/OneDrive - USDA/projects/gsp-sas/LDM-compact_20200709.RData")

s_mps_sf$SAS <- allocate(
  EC = s_mps_sf$ec_ptf.000_030_cm, 
  pH = s_mps_sf$ph_ptf.000_030_cm, 
  ESP = s_mps_sf$esp.000_030_cm
)
s_mps_sf$SAS_int <- as.integer(s_mps_sf$SAS)
write_sf(s_mps_sf, "test.shp")

vars <- c(
  file.path(fp_geo, "gnatsgo_Oct22_SAS_030.tif"),
  file.path(fp_box, "ssurgo_comparison/gsp_sas_030.tif")
)
l  <- vars

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



# tally results ----

## gNATSGO ----
lf_030 <- list.files(path = fp_geo, pattern = "gnatsgo_SAS_030_[0-9]{1,3}.tif", full.names = TRUE)
lf_100 <- list.files(path = fp_geo, pattern = "gnatsgo_SAS_100_[0-9]{1,3}.tif", full.names = TRUE)
lf <- c(lf_030, lf_100)

nlcd_l <- list.files(path = fp_nlcd, pattern = "nlcd2020_[0-9]{1,3}.tif$", full.names = TRUE)
nlcd_l <- c(nlcd_l, nlcd_l)
lev <- levels(nlcd_r)[[1]]
idx <- nchar(lev$`NLCD Land Cover Class`) > 1
lev <- lev[idx, ]


# library(parallel)
# clus <- makeCluster(15)
# clusterExport(clus, list("values", "rast", "lf"))

sas_tal <- lapply(1:length(lf), function(i) {
  
  cat("processing", lf[i], as.character(Sys.time()), "\n")
  
  # x <- rast(lf[i])
  # crs(x) <- "epsg:5070"
  x <- values(rast(lf[i]))
  y <- values(rast(nlcd_l[i]))
  
  sas = x |>
    as.integer() |>
    factor(levels = 1:11)
  
  nlcd = y |>
    as.integer() |>
    factor(levels = lev$val)
  
  tb <- table(sas, nlcd, useNA = "always")
    
  return(tb)
  })
# stopCluster(clus)
# saveRDS(sas_tal, file = file.path(fp_box, "sas_tal.rds"))
sas_tal <- readRDS(file = file.path(fp_box, "sas_tal.rds"))

sas_tal2 <- lapply(sas_tal, function(x) {
  idx <- as.integer(row.names(x))
  x2 <- cbind(nm, as.matrix(x))
  return(x2)
  })
sas_tal2 <- do.call("rbind", sas_tal2) |> as.data.frame()
sas_tal2$nm[is.na(sas_tal2$nm)] <- 0L
n <- nchar(lf)
sas_tal2$dep <- sapply(lf, function(x) strsplit(x, "/|_|\\.")[[1]][9])
names(sas_tal2)[is.na(names(sas_tal2))] <- "missing"

# save(sas_tal2, file = file.path(fp_box, "sas_tally.RData"))
load(file.path(fp_box, "sas_tally.RData"))


## GSP ----
sas_030_r_tal <- sas_030_r |>
  project(y = crs("epsg:5070"), res = 1000, method = "near") |>
  values() |> 
  as.integer() |>
  cut(breaks = 0:11, labels = 1:11, right = TRUE) |>
  table(useNA = "always")

sas_100_r_tal <- sas_100_r |>
  project(y = crs("epsg:5070"), res = 1000, method = "near") |>
  values() |> 
  as.integer() |>
  cut(breaks = 0:11, labels = 1:11, right = TRUE) |>
  table(useNA = "always")


sas_gsp_tal <- rbind(sas_030_r_tal, sas_100_r_tal) |> as.data.frame()
names(sas_gsp_tal[12]) <- "missing"
sas_gsp_tal$dep <- c("030", "100")
# save(sas_gsp_tal, file = file.path(fp_box, "sas_gsp_tal.RData"))
load(file.path(fp_box, "sas_gsp_tal.RData"))



# sas_tal3 <- reshape(
#   sas_tal2,
#   direction = "wide",
#   idvar     = "id",
#   timevar   = "dep",
#   v.names   = c(as.character(1:11, "missing"))
# )


# vapply(sas_tal, function(x) {
#   {x * 30^2} |> 
#     units::set_units(value = "m2") |> 
#     units::set_units(value = "km2") |> 
#     sum()
#   },
#   FUN.VALUE = 0
# )

## CDLS ----
sas_030_r <- rast(file.path(fp_geo, "gnatsgo_Oct22_SAS_030.tif"))
sas_100_r <- rast(file.path(fp_geo, "gnatsgo_Oct22_SAS_100.tif"))
cdls_r    <- rast(file.path(fp_cdls, "2022_30m_cdls_ext.tif"))

rs <- c(sas_030_r, sas_100_r, cdls_r)
samp <- lapply(1:100, function(x) spatSample(rs, 1000))
samp <- do.call("rbind", samp)
names(samp)[1:2] <- c("sas030", "sas100")

samp2 <- data.frame(
  sas = ! samp$sas030 %in% c(5:6) & ! samp$sas100 %in% 5:6, 
  class = tolower(as.character(samp$Class_Names))
)
samp2 <- within(samp2, { 
  class2 = class
  class2 = ifelse(grepl("corn", class2),  "corn", class2)
  class2 = ifelse(grepl("wheat", class2), "wheat", class2)
})
table(samp2$class2, samp2$sas) |> as.matrix() |> as.data.frame.matrix() |> View()


## summarize ----
tal <- aggregate(.~ nm + dep, data = sas_tal2, sum)

tal2 <- by(tal, tal$dep, function(x) {
  
  idx_l <- list(
    sas = which(x$nm %in% c(1:4, 7:11)), 
    sal = which(x$nm %in% 1:4), 
    sod = which(x$nm %in% 7:11), 
    non = which(tal$nm %in% 6)
  )
  
  tal2 <- sapply(idx_l, function(i) {
    {x[i, -1 * 1:2] * 30^2} |>
      colSums(na.rm = TRUE) |>
      units::set_units(value = "m2") |>
      units::set_units(value = "ha")
  })
  tal2 <- data.frame(
    nlcd = c(lev$`NLCD Land Cover Class`, "missing"), 
    tal2
  )
})
tal2

tal_gsp2 <- sapply(idx_l, function(x) {
  {sas_gsp_tal[x] * 1000^2} |>
    rowSums(na.rm = TRUE) |>
    units::set_units(value = "m2") |>
    units::set_units(value = "ha")
}
)
cbind(dep = c("030", "100"),  tal_gsp2)


