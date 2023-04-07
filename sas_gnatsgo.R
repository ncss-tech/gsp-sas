
library(aqp)
library(soilDB)


dsn <- "D:/geodata/soils/gNATSGO_CONUS_Oct2022/gNATSGO_CONUS_Oct2022.gdb"
f <- fetchGDB(dsn, childs = FALSE)
mu <- get_mapunit_from_GDB(dsn, stats = FALSE)
# save(f, mu, file = "gnatsgo_Oct22.RData")
load(file = "gnatsgo_Oct22.RData")


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
# co_avg2[is.na(co_avg2)] <- 0
# 
# sum(co_avg == co_avg2)
# nrow(co_avg) * ncol(co_avg)
# 
# View(co_avg[idx])


co_avg2 <- reshape(
  data = co_avg, 
  direction = "wide",
  idvar = "cokey",
  timevar = "segment_id",
  v.names = c("ec_r", "sar_r", "ph1to1h2o_r", "h")
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
mu_avg <- readRDS("mu_avg_sas.rds")


mu_avg$sas030 <- allocate(pH = mu_avg$`ph1to1h2o_r.000-030`, EC = mu_avg$`ec_r.000-030`, SAR = mu_avg$`sar_r.000-030`, to = "FAO Salt Severity")
mu_avg$sas100 <- allocate(pH = mu_avg$`ph1to1h2o_r.030-100`, EC = mu_avg$`ec_r.030-100`, SAR = mu_avg$`sar_r.030-100`, to = "FAO Salt Severity")



