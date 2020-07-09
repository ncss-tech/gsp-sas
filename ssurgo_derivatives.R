


# mlra
library(dplyr)

nmo <- read.csv("nmsym_mlra_overlap.csv", stringsAsFactors = FALSE)
nmo <- nmo %>%
  group_by(nationalmusym) %>%
  summarize(muacres = sum(area_ac),
            membership = membership[which.max(membership)],
            mlra    = mlra[which.max(membership)]
            ) %>%
  ungroup() %>%
  as.data.frame()
# nmo[with(nmo, order(nationalmusym, - membership)), ]
# nmo <- nmo[! duplicated(nmo$nationalmusym), ]



# fetchSDA
library(aqp)
library(soilDB)

f <- fetchGDB(dsn ="D:/geodata/soils/gSSURGO_CONUS.gdb", WHERE = "areasymbol LIKE '%'")
mu <- get_mapunit_from_GDB(dsn = "D:/geodata/soils/gSSURGO_CONUS.gdb", WHERE = "muname LIKE '%'", stats = TRUE)
# save(f_us, mu_us, file = "C:/Users/stephen.roecker/Nextcloud/data/gnatsgo.RData")
load(file = "C:/Users/stephen.roecker/Nextcloud/data/gnatsgo.RData")


# le <- get_legend_from_SDA(WHERE = "areasymbol LIKE '%'")
# st <- unique(substr(le$areasymbol, 1, 2))
# 
# get_states <- function(x) {
#   lapply(x, function(x) {
#     cat(as.character(Sys.time()), "getting ", x)
#     tryCatch({
#       fetchSDA(WHERE = paste0("areasymbol LIKE '", x, "%'"))
#       },
#       error = function(err) {
#         print(paste("Error occured: ", err))
#         return(NULL)
#         }
#       )
#     })
# }
# test <- get_states(st[st != "MI"])
# names(test) <- st[st != "MI"]
# 
# st2 <- names(test)[unlist(lapply(test, is.null))]
# st2 <- paste0(c(st2, "MI"), collapse = "|")
# as <- le[grepl(st2, le$areasymbol), "areasymbol"]
# test2 <- get_states(as)
# names(test2) <- as
# 
# # save(test, test2, f, mu, file = "gsp_salinity_fetchSDA.RData")
# load(file = "gsp_salinity_fetchSDA.RData")
# 
# # combine fetchSDA
# idx <- unlist(lapply(test, is.null))
# test <- test[! idx]
# 
# sn <- lapply(test, function(x) siteNames(x))
# snt <- table(unlist(sn))
# sn <- dimnames(snt[snt == 55])[[1]]
# 
# s <- do.call("rbind", lapply(test, function(x) site(x)[sn]))
# h <- do.call("rbind", lapply(test, function(x) horizons(x)))
# f <- h
# depths(f) <- cokey ~ hzdept_r + hzdepb_r
# site(f) <- s[!duplicated(s$cokey), ]


site(f_us)$idx <- with(site(f_us, paste(compname, comppct_r, compkind, majcompflag, local_phase)))

mu <- aqp::horizons(aqp::slice(f_us, 0:100 ~ cokey + ec_r, strict = FALSE)) %>%
  # component
  group_by(cokey) %>%
  summarize(ec_r = weighted.mean(ec_r, w = hzdepb_r - hzdept_r)) %>%
  ungroup() %>%
  # map unit
  right_join(aqp::site(f), by = "cokey") %>%
  group_by(nationalmusym) %>%
  summarize(ec_r = weighted.mean(ec_r, w = comppct_r, na.rm = TRUE),
            landform = ifelse(all(complete.cases(landform, comppct_r)), aggregate(comppct_r ~ landform, FUN = sum, na.rm = TRUE)[1, 1], NA),
            pmkind   = ifelse(all(complete.cases(pmkind, comppct_r)),   aggregate(comppct_r ~ pmkind, FUN = sum, na.rm = TRUE)[1, 1], NA),
            pmorigin  = ifelse(all(complete.cases(pmorigin, comppct_r)), aggregate(comppct_r ~ pmorigin, FUN = sum, na.rm = TRUE)[1, 1], NA)
            ) %>%
  ungroup() %>%
  # mlra
  inner_join(nmo, by = "nationalmusym")

# save(test, f, mu, file = "gsp_salinity_fetchSDA.RData")



# aggregate by pmkind, pmorigin and landform
load(file = "gsp_salinity_fetchSDA.RData")

mu_pmorigin <- mu %>%
  # pmkind, pmorigin, and landform
  filter(muacres > 0) %>%
  group_by(mlra, pmorigin) %>%
  summarize(p = ifelse(any(ec_r > 2), sum(muacres[ec_r > 2], na.rm = TRUE) / sum(muacres, na.rm = TRUE), NA)) %>%
  ungroup() %>%
  mutate(p = round(p * 100))
  # filter(ec_r > 2)


mu_pmkind <- mu %>%
  # pmkind, pmorigin, and landform
  filter(muacres > 0) %>%
  group_by(mlra, pmkind) %>%
  summarize(p = ifelse(any(ec_r > 2), sum(muacres[ec_r > 2], na.rm = TRUE) / sum(muacres, na.rm = TRUE), NA)) %>%
  ungroup() %>%
  mutate(p = round(p * 100))
# filter(ec_r > 2)


mu_landform <- mu %>%
  # pmkind, pmorigin, and landform
  filter(muacres > 0) %>%
  group_by(mlra, landform) %>%
  summarize(p = ifelse(any(ec_r > 2), sum(muacres[ec_r > 2], na.rm = TRUE) / sum(muacres, na.rm = TRUE), NA)) %>%
  ungroup() %>%
  mutate(p = round(p * 100))
# filter(ec_r > 2)



mu$ec_c <- cut(mu$ec_r, breaks = c(0, 2, 4, 8, 16))

ec_t <- as.data.frame(xtabs(~ pmkind + ec_c, data = mu))
ec_t <- tidyr::spread(ec_t, ec_c, Freq)
          
