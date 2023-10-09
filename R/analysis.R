library(vegan)
library(rioja)
library(tidyverse)
library(geosphere)
library(broom)
library(gstat)


source("R/functions.R")


# What seasonal and vertical temperature explain most of the variance in the species assemblages ============================

# do CCA on core top and fossil assemblages
# to determine which temperature (season/depth) explains most of the variance in the assemblages
# use WOA98 temperatures from different seasons + annual mean as well as different depths down to 500 m
# perform three different CCAs
# 1) core top assemblages and atlas temperatures
# 2) core top assemblages and temperatures predicted from MAT
# 3) fossil assemblages and temperatures predicted from MAT

# the latter two follow the rationale used in randomTF (Telford & Birks 2011, QSR)
# the analyses are done for each ocean basin separately (cf. Kucera et al., 2005, QSR)
# note that there is some overlap between the regions used for the transfer model development

# where applicable ocean basins are split in tropical and extratropical regions
# since it is expected that depth will be more important in the tropics because of the
# large vertical temperature gradients

# 1) CCA core top assemblages constrained with atlas temperatures


# load species data
# note that the number of species is not the same for each region
# and add temperature data

coretopSpp <- map(readRDS("data/ForCenS_spp.RDS"), function(x){
  atlasTemp <- readRDS("data/ForCenS_atlas_temperatures.RDS")
  uID <- x$meta$forcensID
  indx <- sapply(uID, function(y) which(row.names(atlasTemp) == y))
  x$env <- atlasTemp[indx, ]
  x
})

# combine all
CCA_coretop_atlas <- bind_rows(tibble(doCCA("IND", coretopSpp$IND$meta$Latitude <= -23.5), region = "Indian extrop"),
                               tibble(doCCA("IND", coretopSpp$IND$meta$Latitude > -23.5), region = "Indian trop"),
                               tibble(doCCA("MDX", 1:nrow(coretopSpp$MDX$species)), region = "Mediterranean"),
                               tibble(doCCA("NAT", coretopSpp$NAT$meta$Latitude >= 23.5), region = "North Atlantic extrop"),
                               tibble(doCCA("NAT", coretopSpp$NAT$meta$Latitude < 23.5), region = "North Atlantic trop"),
                               tibble(doCCA("PAC", coretopSpp$PAC$meta$Latitude >= 23.5), region = "Pacific north"),
                               tibble(doCCA("PAC", coretopSpp$PAC$meta$Latitude <= -23.5), region = "Pacific south"),
                               tibble(doCCA("PAC", coretopSpp$PAC$meta$Latitude > -23.5 & coretopSpp$PAC$meta$Latitude < 23.5), region = "Pacific trop"),
                               tibble(doCCA("SAT", coretopSpp$SAT$meta$Latitude <= -23.5), region = "South Atlantic extrop"),
                               tibble(doCCA("SAT", coretopSpp$SAT$meta$Latitude > -23.5), region = "South Atlantic trop")
)
CCA_coretop_atlas$group <- "Climatology"


# 2 and 3) CCA on core top and fossil assemblages constrained with temperatures predicted from assemblages
fosSpp <- readRDS("data/LGM_spp_19082022.RDS") %>%
  map(., function(x){
  x$meta <- x$meta %>%
    group_by(lgmID) %>%
    mutate(sampleID = paste(lgmID, 1:n(), sep = "_")) %>%
    ungroup()
  x
})

# predict temperature for core tops
NAT_mod_temp <- apply(coretopSpp$NAT$env, 2, function(x) doMAT(spp = coretopSpp$NAT, env = x))
SAT_mod_temp <- apply(coretopSpp$SAT$env, 2, function(x) doMAT(spp = coretopSpp$SAT, env = x))
IND_mod_temp <- apply(coretopSpp$IND$env, 2, function(x) doMAT(spp = coretopSpp$IND, env = x))
PAC_mod_temp <- apply(coretopSpp$PAC$env, 2, function(x) doMAT(spp = coretopSpp$PAC, env = x))
MDX_mod_temp <- apply(coretopSpp$MDX$env, 2, function(x) doMAT(spp = coretopSpp$MDX, env = x))


# predict temperature for fossil data
NAT_fos_temp <- apply(coretopSpp$NAT$env, 2, function(x) doMAT(spp = coretopSpp$NAT, fos = fosSpp$NAT, env = x))
SAT_fos_temp <- apply(coretopSpp$SAT$env, 2, function(x) doMAT(spp = coretopSpp$SAT, fos = fosSpp$SAT, env = x))
IND_fos_temp <- apply(coretopSpp$IND$env, 2, function(x) doMAT(spp = coretopSpp$IND, fos = fosSpp$IND, env = x))
PAC_fos_temp <- apply(coretopSpp$PAC$env, 2, function(x) doMAT(spp = coretopSpp$PAC, fos = fosSpp$PAC, env = x))
MDX_fos_temp <- apply(coretopSpp$MDX$env, 2, function(x) doMAT(spp = coretopSpp$MDX, fos = fosSpp$MDX, env = x))

# variance explained for entire regions

# for core top assemblages
EX_NAT_mod <- map_dbl(NAT_mod_temp, function(x) varEx(spp = coretopSpp$NAT, pred = x))
EX_SAT_mod <- map_dbl(SAT_mod_temp, function(x) varEx(spp = coretopSpp$SAT, pred = x))
EX_IND_mod <- map_dbl(IND_mod_temp, function(x) varEx(spp = coretopSpp$IND, pred = x))
EX_PAC_mod <- map_dbl(PAC_mod_temp, function(x) varEx(spp = coretopSpp$PAC, pred = x))
EX_MDX_mod <- map_dbl(MDX_mod_temp, function(x) varEx(spp = coretopSpp$MDX, pred = x))

# for fossil assembalges
EX_NAT_fos <- map_dbl(NAT_fos_temp, function(x) varEx(spp = fosSpp$NAT, pred = x))
EX_SAT_fos <- map_dbl(SAT_fos_temp, function(x) varEx(spp = fosSpp$SAT, pred = x))
EX_IND_fos <- map_dbl(IND_fos_temp, function(x) varEx(spp = fosSpp$IND, pred = x))
EX_PAC_fos <- map_dbl(PAC_fos_temp, function(x) varEx(spp = fosSpp$PAC, pred = x))
EX_MDX_fos <- map_dbl(MDX_fos_temp, function(x) varEx(spp = fosSpp$MDX, pred = x))

# indices for subregions
NAT_mod_trop <- map(NAT_mod_temp, function(x) which(x$latitude <= 23.5))
NAT_mod_extrop <- map(NAT_mod_temp, function(x) which(x$latitude > 23.5))
SAT_mod_trop <- map(SAT_mod_temp, function(x) which(x$latitude >= -23.5))
SAT_mod_extrop <- map(SAT_mod_temp, function(x) which(x$latitude < -23.5))
IND_mod_trop <- map(IND_mod_temp, function(x) which(x$latitude >= -23.5))
IND_mod_extrop <- map(IND_mod_temp, function(x) which(x$latitude < -23.5))
PAC_mod_trop <- map(PAC_mod_temp, function(x) which(x$latitude <= 23.5 & x$latitude >= -23.5))
PAC_mod_north <- map(PAC_mod_temp, function(x) which(x$latitude > 23.5))
PAC_mod_south <- map(PAC_mod_temp, function(x) which(x$latitude < -23.5))

NAT_fos_trop <- map(NAT_fos_temp, function(x) which(x$latitude <= 23.5))
NAT_fos_extrop <- map(NAT_fos_temp, function(x) which(x$latitude > 23.5))
SAT_fos_trop <- map(SAT_fos_temp, function(x) which(x$latitude >= -23.5))
SAT_fos_extrop <- map(SAT_fos_temp, function(x) which(x$latitude < -23.5))
IND_fos_trop <- map(IND_fos_temp, function(x) which(x$latitude >= -23.5))
IND_fos_extrop <- map(IND_fos_temp, function(x) which(x$latitude < -23.5))
PAC_fos_trop <- map(PAC_fos_temp, function(x) which(x$latitude <= 23.5 & x$latitude >= -23.5))
PAC_fos_north <- map(PAC_fos_temp, function(x) which(x$latitude > 23.5))
PAC_fos_south <- map(PAC_fos_temp, function(x) which(x$latitude < -23.5))

# variance explained for subregions
EX_NAT_trop_mod <- map2_dbl(.f = function(x, y ) varEx(spp = coretopSpp$NAT, pred = x, indx = y), .x = NAT_mod_temp, .y = NAT_mod_trop)
EX_NAT_extrop_mod <- map2_dbl(.f = function(x, y ) varEx(spp = coretopSpp$NAT, pred = x, indx = y), .x = NAT_mod_temp, .y = NAT_mod_extrop)
EX_SAT_trop_mod <- map2_dbl(.f = function(x, y ) varEx(spp = coretopSpp$SAT, pred = x, indx = y), .x = SAT_mod_temp, .y = SAT_mod_trop)
EX_SAT_extrop_mod <- map2_dbl(.f = function(x, y ) varEx(spp = coretopSpp$SAT, pred = x, indx = y), .x = SAT_mod_temp, .y = SAT_mod_extrop)
EX_IND_trop_mod <- map2_dbl(.f = function(x, y ) varEx(spp = coretopSpp$IND, pred = x, indx = y), .x = IND_mod_temp, .y = IND_mod_trop)
EX_IND_extrop_mod <- map2_dbl(.f = function(x, y ) varEx(spp = coretopSpp$IND, pred = x, indx = y), .x = IND_mod_temp, .y = IND_mod_extrop)
EX_PAC_trop_mod <- map2_dbl(.f = function(x, y ) varEx(spp = coretopSpp$PAC, pred = x, indx = y), .x = PAC_mod_temp, .y = PAC_mod_trop)
EX_PAC_north_mod <- map2_dbl(.f = function(x, y ) varEx(spp = coretopSpp$PAC, pred = x, indx = y), .x = PAC_mod_temp, .y = PAC_mod_north)
EX_PAC_south_mod <- map2_dbl(.f = function(x, y ) varEx(spp = coretopSpp$PAC, pred = x, indx = y), .x = PAC_mod_temp, .y = PAC_mod_south)

EX_NAT_trop_fos <- map2_dbl(.f = function(x, y ) varEx(spp = fosSpp$NAT, pred = x, indx = y), .x = NAT_fos_temp, .y = NAT_fos_trop)
EX_NAT_extrop_fos <- map2_dbl(.f = function(x, y ) varEx(spp = fosSpp$NAT, pred = x, indx = y), .x = NAT_fos_temp, .y = NAT_fos_extrop)
EX_SAT_trop_fos <- map2_dbl(.f = function(x, y ) varEx(spp = fosSpp$SAT, pred = x, indx = y), .x = SAT_fos_temp, .y = SAT_fos_trop)
EX_SAT_extrop_fos <- map2_dbl(.f = function(x, y ) varEx(spp = fosSpp$SAT, pred = x, indx = y), .x = SAT_fos_temp, .y = SAT_fos_extrop)
EX_IND_trop_fos <- map2_dbl(.f = function(x, y ) varEx(spp = fosSpp$IND, pred = x, indx = y), .x = IND_fos_temp, .y = IND_fos_trop)
EX_IND_extrop_fos <- map2_dbl(.f = function(x, y ) varEx(spp = fosSpp$IND, pred = x, indx = y), .x = IND_fos_temp, .y = IND_fos_extrop)
EX_PAC_trop_fos <- map2_dbl(.f = function(x, y ) varEx(spp = fosSpp$PAC, pred = x, indx = y), .x = PAC_fos_temp, .y = PAC_fos_trop)
EX_PAC_north_fos <- map2_dbl(.f = function(x, y ) varEx(spp = fosSpp$PAC, pred = x, indx = y), .x = PAC_fos_temp, .y = PAC_fos_north)
EX_PAC_south_fos <- map2_dbl(.f = function(x, y ) varEx(spp = fosSpp$PAC, pred = x, indx = y), .x = PAC_fos_temp, .y = PAC_fos_south)


# combine all
globalEX <- bind_rows(bind_cols(EX = EX_IND_extrop_fos, region = "Indian extrop", group = "LGM"),
          bind_cols(EX = EX_IND_trop_fos, region = "Indian trop", group = "LGM"),
          bind_cols(EX = EX_IND_extrop_mod, region = "Indian extrop", group = "Core top"),
          bind_cols(EX = EX_IND_trop_mod, region = "Indian trop", group = "Core top"),
          
          bind_cols(EX = EX_MDX_fos, region = "Mediterranean", group = "LGM"),
          bind_cols(EX = EX_MDX_mod, region = "Mediterranean", group = "Core top"),
         
          bind_cols(EX = EX_NAT_extrop_fos, region = "North Atlantic extrop", group = "LGM"),
          bind_cols(EX = EX_NAT_trop_fos, region = "North Atlantic trop", group = "LGM"),
          bind_cols(EX = EX_NAT_extrop_mod, region = "North Atlantic extrop", group = "Core top"),
          bind_cols(EX = EX_NAT_trop_mod, region = "North Atlantic trop", group = "Core top"),
          
          bind_cols(EX = EX_PAC_north_fos, region = "Pacific north", group = "LGM"),
          bind_cols(EX = EX_PAC_north_mod, region = "Pacific north", group = "Core top"),
          bind_cols(EX = EX_PAC_trop_fos, region = "Pacific trop", group = "LGM"),
          bind_cols(EX = EX_PAC_trop_mod, region = "Pacific trop", group = "Core top"),
          bind_cols(EX = EX_PAC_south_fos, region = "Pacific south", group = "LGM"),
          bind_cols(EX = EX_PAC_south_mod, region = "Pacific south", group = "Core top"),
          
          bind_cols(EX = EX_SAT_extrop_fos, region = "South Atlantic extrop", group = "LGM"),
          bind_cols(EX = EX_SAT_trop_fos, region = "South Atlantic trop", group = "LGM"),
          bind_cols(EX = EX_SAT_extrop_mod, region = "South Atlantic extrop", group = "Core top"),
          bind_cols(EX = EX_SAT_trop_mod, region = "South Atlantic trop", group = "Core top")
) %>%
  mutate(SeasonDepth = rep(names(EX_IND_extrop_fos), 20)) %>%
  separate(SeasonDepth, c("season", "depth")) %>%
  bind_rows(CCA_coretop_atlas) %>%
  mutate(depth = as.numeric(depth))


# annual 50 is overall highest and also in virtually the highest in the fossil
globalEX %>%
  group_by(season, depth) %>%
  summarise(EX = mean(EX)) %>%
  ungroup() %>%
  arrange(desc(EX))

# ==============================================================================================================================

# temperature predicted from core tops ===============================================================

# annual mean at 50 m water depth 
# needed later to assess if it matters how LGM temperature anomalies are calculated
coretop_spp_temp <- list(
  IND = IND_mod_temp$annual_50 %>%
    select(samples, predT),
  MDX = MDX_mod_temp$annual_50 %>%
    select(samples, predT),
  NAT = NAT_mod_temp$annual_50 %>%
    select(samples, predT),
  PAC = PAC_mod_temp$annual_50 %>%
    select(samples, predT),
  SAT = SAT_mod_temp$annual_50 %>%
    select(samples, predT)
) %>%
  map(., ~rename(., 
                 forcensID = samples,
                 annual_50 = predT)
      )

# ==============================================================================================================================

# LGM temperatures from assemblages ===================================================================
LGM_rec_temperatures <- list(IND = IND_fos_temp$annual_50 %>%
                           select(predT),
                         MDX = MDX_fos_temp$annual_50 %>%
                           select(predT),
                         NAT = NAT_fos_temp$annual_50 %>%
                           select(predT),
                         PAC = PAC_fos_temp$annual_50 %>%
                           select(predT),
                         SAT = SAT_fos_temp$annual_50 %>%
                           select(predT)
                         )
# add meta data
LGM_rec_temperatures <- map2_df(.x = fosSpp, .y = LGM_rec_temperatures, .f = function(x, y){
  x$meta %>%
    select(Core, Latitude, Longitude, Ocean, lgmID, sampleID) %>%
    bind_cols(y)
    }) %>%
  tibble()
           
# add square chord distance to most similar core top, which can later be used for filtering
sqcd <- map2_df(coretopSpp, fosSpp, function(x, y){
  mod <- MAT(x$species, x$env$annual_0, k = 10, lean = FALSE)
  tibble(sampleID = y$meta$sampleID,
         sqcd = predict(mod, y$species)$dist.n[,1])
})

LGM_rec_temperatures <- left_join(LGM_rec_temperatures, sqcd, by = "sampleID") %>%
  arrange(sampleID)


# add analogue cutoffs based on intersample dissimilarity in core top data
LGM_rec_temperatures <- LGM_rec_temperatures %>%
   mutate(analogueCutoff = case_when(Ocean == "Indian" ~ quantile(paldist(coretopSpp$IND$species), 0.1),
                                     Ocean == "Mediterranean" ~ quantile(paldist(coretopSpp$MDX$species), 0.1),
                                     Ocean == "North Atlantic" ~ quantile(paldist(coretopSpp$NAT$species), 0.1),
                                     Ocean == "Pacific" ~ quantile(paldist(coretopSpp$PAC$species), 0.1),
                                     Ocean == "South Atlantic" ~ quantile(paldist(coretopSpp$SAT$species), 0.1))
          )

# clean up a bit
rm(list = c("CCA_coretop_atlas", "doCCA", "doMAT", "EX_IND_extrop_fos", "EX_IND_extrop_mod", "EX_IND_fos", "EX_IND_mod",
     "EX_IND_trop_fos", "EX_IND_trop_mod", "EX_MDX_fos", "EX_MDX_mod", "EX_NAT_extrop_fos", "EX_NAT_extrop_mod",
     "EX_NAT_fos", "EX_NAT_mod", "EX_NAT_trop_fos", "EX_NAT_trop_mod", "EX_PAC_fos", "EX_PAC_mod", "EX_PAC_north_fos",
     "EX_PAC_north_mod", "EX_PAC_south_fos", "EX_PAC_south_mod", "EX_PAC_trop_fos", "EX_PAC_trop_mod", "EX_SAT_extrop_fos",
     "EX_SAT_extrop_mod", "EX_SAT_fos", "EX_SAT_mod", "EX_SAT_trop_fos", "EX_SAT_trop_mod", "IND_fos_extrop", "IND_fos_temp",
     "IND_fos_trop", "IND_mod_extrop", "IND_mod_temp", "IND_mod_trop", "MDX_fos_temp", "MDX_mod_temp", "NAT_fos_extrop",
     "NAT_fos_temp", "NAT_fos_trop", "NAT_mod_extrop", "NAT_mod_temp", "NAT_mod_trop", "PAC_fos_north", "PAC_fos_south",
     "PAC_fos_temp", "PAC_fos_trop", "PAC_mod_north", "PAC_mod_south", "PAC_mod_temp", "PAC_mod_trop", "SAT_fos_extrop",
     "SAT_fos_temp", "SAT_fos_trop", "SAT_mod_extrop", "SAT_mod_temp", "SAT_mod_trop", "sqcd", "varEx")
)

# ==============================================================================================================================

# similarity decay =============================================================================================================

# prep core top assemblages because the number of species is not the same for each region
# match spp data with temperature data from World Ocean Atlas annual mean 50 m

modernSpp <- coretopSpp

modernSpp$IND$species <- replaceRuber(modernSpp$IND$species)
modernSpp$NAT$species <- replaceRuber(modernSpp$NAT$species)
modernSpp$SAT$species <- replaceRuber(modernSpp$SAT$species)
modernSpp$PAC$species <- replaceRuber(modernSpp$PAC$species)
modernSpp$IND$species <- replaceMenardii(modernSpp$IND$species)
modernSpp$NAT$species <- replaceMenardii(modernSpp$NAT$species)
modernSpp$SAT$species <- replaceMenardii(modernSpp$SAT$species)
modernSpp$MDX$species <- replaceMenardii(modernSpp$MDX$species)

# make single df
modernSppDF <- map_df(modernSpp, function(x){
  bind_cols(select(x$meta, forcensID, Latitude, Longitude), x$species, select(x$env, annual_50))
}, .id = "Ocean")

# remove duplicated samples (as regions overlap)
modernSppDF <- modernSppDF[!duplicated(modernSppDF$forcensID), ]
row.names(modernSppDF) <- modernSppDF$forcensID
rm(modernSpp)

# 5th percentile of SQCD distances between core top samples
# used for analogue filtering later on
globalNoAnalogueValue <- modernSppDF %>%
  select(Dentigloborotalia_anfracta:Globorotalia_menardii_Globorotalia_tumida) %>%
  paldist() %>%
  quantile(0.05)

# prep LGM assemblages in the same way
LGMSpp <- fosSpp
# sort taxonomy
LGMSpp$IND$species <- replaceRuber(LGMSpp$IND$species)
LGMSpp$NAT$species <- replaceRuber(LGMSpp$NAT$species)
LGMSpp$SAT$species <- replaceRuber(LGMSpp$SAT$species)
LGMSpp$PAC$species <- replaceRuber(LGMSpp$PAC$species)
LGMSpp$IND$species <- replaceMenardii(LGMSpp$IND$species)
LGMSpp$NAT$species <- replaceMenardii(LGMSpp$NAT$species)
LGMSpp$SAT$species <- replaceMenardii(LGMSpp$SAT$species)
LGMSpp$MDX$species <- replaceMenardii(LGMSpp$MDX$species)


LGMSppDF_noT <- map_df(LGMSpp, function(x){
  bind_cols(select(x$meta, lgmID, sampleID, Latitude, Longitude, ChronozoneLevel), x$species)
})
LGMSppDF_noT <- LGMSppDF_noT[order(LGMSppDF_noT$lgmID), ]

rm(LGMSpp)

# average assemblages for sites
LGMSppDF_avg_noT <- LGMSppDF_noT %>%
  select(-sampleID) %>%
  group_by(lgmID) %>%
  summarise_all(mean)

# and re-normalise
LGMSppDF_avg_noT[, 5:43] <- sweep(LGMSppDF_avg_noT[, 5:43], 1, rowSums(LGMSppDF_avg_noT[, 5:43], na.rm = TRUE), FUN = "/")

# cluster analysis on all assemblages
allSppDF <- modernSppDF %>%
  select(-Ocean, - annual_50) %>%
  tibble() %>%
  rename(uID = forcensID) %>%
  bind_rows(., LGMSppDF_avg_noT %>%
              rename(uID = lgmID))


set.seed(123)
clusterAllSpp <- allSppDF %>%
  select(-uID, - Latitude, - Longitude, - ChronozoneLevel) %>%
  kmeans(., 3, nstart = 50)

clusterOut <- tibble(allSppDF %>%
                       select(uID, Longitude, Latitude) %>%
                       mutate(set = case_when(grepl("ForCenS", uID) ~ "Modern",
                                              TRUE ~ "LGM")),
                     cluster = clusterAllSpp$cluster) %>%
  mutate(set = factor(set, levels = c("Modern", "LGM")))


# add temperatures from simulations
# note that there are more samples than temperatures because for many cores we have multiple samples for the LGM temperature interval
# calculate LGM temperature from simulated anomaly and atlas temperatures

modelsLGM <- readRDS("data/simTemp_LGM_annual50.RDS") %>%
  distinct(model)

notInEnsemble <- c("lgmID", "AWI.ESM_coast", "AWI.ESM_coast_ext", "AWI.ESM_subpolar")

LGMSppDF <- readRDS("data/LGM_atlas_temperatures.RDS") %>%
  tibble() %>%
  dplyr::select(uniqueID, annual_50) %>%
  rename(lgmID = uniqueID) %>%
  rename(T_WOA98 = annual_50) %>%
  left_join(readRDS("data/simTemp_LGM_annual50.RDS") %>%
              dplyr::select(lgmID, model, T.LGM:T.anomaly), . , by = "lgmID") %>%
  mutate(T.LGM_WOA = T_WOA98 - T.anomaly) %>%
  dplyr::select(lgmID, model, T.LGM_WOA) %>%
  pivot_wider(names_from = model, values_from = T.LGM_WOA) %>%
  mutate(ensemble = rowMeans(select(., -all_of(notInEnsemble)), na.rm = TRUE)) %>%
  left_join(LGMSppDF_noT, ., by = "lgmID") %>%
  tibble()

LGMSppDF_avg <- readRDS("data/LGM_atlas_temperatures.RDS") %>%
  dplyr::select(uniqueID, annual_50) %>%
  rename(lgmID = uniqueID) %>%
  rename(T_WOA98 = annual_50) %>%
  left_join(readRDS("data/simTemp_LGM_annual50.RDS") %>%
              dplyr::select(lgmID, model, T.LGM:T.anomaly), . , by = "lgmID") %>%
  mutate(T.LGM_WOA = T_WOA98 - T.anomaly) %>%
  select(lgmID, model, T.LGM_WOA) %>%
  pivot_wider(names_from = model, values_from = T.LGM_WOA) %>%
  mutate(ensemble = rowMeans(select(., -all_of(notInEnsemble)), na.rm = TRUE)) %>%
  left_join(LGMSppDF_avg_noT, ., by = "lgmID")


# subsetting to determine the effect of spatial sample coverage on the similarity decay pattern

# get distance between core top and nearest LGM sample 
d2LGM <- apply(as.matrix(modernSppDF[, c("Longitude", "Latitude")], ncol = 2), 1, function(x) min(distGeo(p1 = x, p2 = LGMSppDF_avg[, c("Longitude", "Latitude")])))/1000

# minimum distance between LGM sample and core top
spDistUniqueLGM <- apply(LGMSppDF_avg[, c("Longitude", "Latitude")], 1, function(z) distGeo(p1 = modernSppDF[, c("Longitude", "Latitude")], p2 = z)/1000)
d2coretopUniqueLGM <- cbind.data.frame(lgmID = LGMSppDF_avg$lgmID,
                                       minDist = apply(spDistUniqueLGM, 2, function(a) min(a)),
                                       forcensID = apply(spDistUniqueLGM, 2, function(a) modernSppDF$forcensID[which.min(a)]))

# get community and environmental distance matrices for modern

# distance for all core top samples
modDist <- makeDist(modernSppDF)

# and for core top samples with a nearby (100 km) LGM sample
modDistSub <- makeDist(modernSppDF[d2LGM <= 100,])

# and only for the core top samples that are closest to the LGM samples
modDistSub2 <- makeDist(map_df(d2coretopUniqueLGM$forcensID, function(x) modernSppDF[modernSppDF$forcensID %in% x, ]))

# get community and environmental distance matrices using simulated PI temperatures 
modSpp_PI <- readRDS("data/simTemp_Coretop_annual50.RDS") %>%
  dplyr::select(forcensID, model, T.PI) %>%
  pivot_wider(names_from = model, values_from = T.PI) %>%
  mutate(ensemble = rowMeans(select(., -forcensID), na.rm = TRUE)) %>%
  left_join(modernSppDF, ., by = "forcensID") %>%
  rename(AWI.ESM = AWI)

modelsPI <- names(modSpp_PI)[44: ncol(modSpp_PI)]

modSppDist <- as.matrix(vegdist(select(modSpp_PI, Dentigloborotalia_anfracta:Globorotalia_menardii_Globorotalia_tumida), method = "bray"))
# convert to similarity
modSppDist <- 1- modSppDist
modSppDist <- melt(modSppDist)[melt(upper.tri(modSppDist, diag = FALSE))$value,]
names(modSppDist)[3] <- "commdist"

# this file is large (1.7 GB)
PI_Dist <- modSpp_PI %>%
  select(all_of(modelsPI)) %>%
  map_dfr(., function(x){
    EnvDist <- as.matrix(dist(x))
    EnvDist <- melt(EnvDist)[melt(upper.tri(EnvDist, diag = FALSE))$value,]
    names(EnvDist)[3] <- "envdist"
    EnvDist
  }, .id = "model") %>%
  left_join(., modSppDist, by = c("Var1", "Var2")) %>%
  select(-c(Var1, Var2)) %>%
  tibble()


# get community and environmental distances matrices for LGM

# for averaged LGM assemblages
lgmSppDist_avg <- as.matrix(vegdist(select(LGMSppDF_avg, Dentigloborotalia_anfracta:Globorotalia_menardii_Globorotalia_tumida), method = "bray"))
# similarity
lgmSppDist_avg <- 1- lgmSppDist_avg
# melt to df, keep only upper triangle
lgmSppDist_avg <- melt(lgmSppDist_avg)[melt(upper.tri(lgmSppDist_avg, diag = FALSE))$value,]
names(lgmSppDist_avg)[3] <- "commdist"

# LGM species assemblages by Chronozone level
lgmSppDist_avg_chron <- LGMSppDF_avg %>%
  group_by(ChronozoneLevel) %>%
  group_split() %>%
  map_dfr(., function(x){
    comdis <- x %>%
      select(Dentigloborotalia_anfracta:Globorotalia_menardii_Globorotalia_tumida) %>%
      vegdist(method = "bray") %>%
      tidy() %>%
      mutate(commdist = 1 - distance) %>%
      select(-distance) %>%
      rename(Var1 = 1,
             Var2 = 2)
    edis <- x %>%
      reframe(across(`NCAR-CCSM4-1`:ensemble, ~as.numeric(dist(.x))))
    bind_cols(comdis, edis) %>%
      select(-(1:2)) %>%
      pivot_longer(-commdist, values_to = "envdist", names_to = "model")
  }, .id = "chronolevel") %>%
  mutate(chronolevel = factor(chronolevel))


# LGM environmental distance matrix
models <- names(LGMSppDF_avg)[44:ncol(LGMSppDF_avg)]

# get dissimilarity data for each model  
LGMDist_avg <- LGMSppDF_avg %>%
  select(all_of(models)) %>%
  map_dfr(getEnvDist_avg, .id = "model") %>%
  select(-c(Var1, Var2)) %>%
  mutate(model = as.factor(model)) %>%
  tibble()


# make grids for plotting distance decay
resEnv <- 0.5
resSim <- 0.025

modDistBin <- modDist %>%
  gridFunction(., resEnv, resSim)

modDistSubBin <- modDistSub %>%
  gridFunction(., resEnv, resSim)

modDistSub2Bin <- modDistSub2 %>%
  gridFunction(., resEnv, resSim)

LGMDist_avgBin <- LGMDist_avg %>%
  named_group_split(., model) %>%
  map_dfr(., function(x) gridFunction(x, resEnv, resSim), .id = "model")

LGMDist_avgBin_chron <- lgmSppDist_avg_chron %>%
  named_group_split(model, chronolevel) %>%
  map_dfr(., function(x) gridFunction(x, resEnv, resSim), .id = "model") %>%
  separate(model, into = c("model", "chronolevel"), sep = "_")

PI_Dist_Bin <- PI_Dist %>%
  split(., .$model) %>%
  map_dfr(., function(x) gridFunction(x, resEnv, resSim), .id = "model") %>%
  mutate(model = factor(model, levels = c("modern", modelsPI)))

# define limits of modern distance decay pattern
ecolim <- modDist %>%
  limfunction(x = ., resSim = resSim, quant = 0.99)

ecolimSub <- modDistSub %>%
  limfunction(x = ., resSim = resSim, quant = 0.99)

ecolimSub2 <- modDistSub2 %>%
  limfunction(x = ., resSim = resSim, quant = 0.99)

# which site pairs have large temperature differences in the simulations,
# but high community similarity?
extraLimitSites <- LGMDist_avg %>%
  filter(model == "ensemble") %>%
  tibble() %>%
  mutate(simBin = cut(commdist, breaks = seq(0, 1, resSim), include.lowest = TRUE, right = FALSE, labels = seq(0, 1, resSim)[-1] - 0.5*resSim)) %>%
  left_join(., ecolim, by = "simBin") %>%
  mutate(indx = envdist > Tlim) %>%
  select(indx) %>%
  bind_cols(., lgmSppDist_avg) %>%
  filter(indx == TRUE & commdist >= 0.75) %>%
  pivot_longer(-c(indx, commdist), values_to = "site", names_to = "dummy") %>%
  group_by(site) %>%
  summarise(n = n())

# modern and LGM distance decay and the difference 
diffSimDecayModLGM <- LGMDist_avgBin %>%
  split(., .$model) %>% # take this detour of splitting by model because many models have NAs for some bins
  map_dfr(., ~full_join(modDistBin, .,  by = c("envBin", "simBin")), .id = "model") %>%
  replace_na(list(lognnorm.x = 0, lognnorm.y = 0)) %>%
  mutate(diff = lognnorm.y - lognnorm.x)

# modern and LGM distance decay using ensemble mean temperature
plotDist <- bind_rows(
  modDist %>%
    mutate(model = "modern"),
  LGMDist_avg %>%
    filter(model == "ensemble")) %>%
  mutate(bin = cut(envdist, breaks = seq(0, 35, 1), include.lowest = TRUE, right = FALSE, labels = seq(0, 35, 1)[-1]-0.5*1),
         model = factor(model, levels = c("modern", "ensemble"), labels = c("Modern", "LGM modelled")))

# what does LGM similarity decay look like with modern temperatures?
LGM_atlasDist <- LGMSppDF_avg %>%
  left_join(., readRDS("data/LGM_atlas_temperatures.RDS") %>%
              tibble() %>%
              dplyr::select(uniqueID, annual_50) %>%
              rename(lgmID = uniqueID), by = "lgmID") %>%
  makeDist() %>%
  gridFunction(., resEnv, resSim)

# ==============================================================================================================================

# temporal turnover ===================================================================================

# get assemblage at nearest core top to estimate turnover

# # minimum distance to core top
spDist <- apply(LGMSppDF[, c("Longitude", "Latitude")], 1, function(z) distGeo(p1 = modernSppDF[, c("Longitude", "Latitude")], p2 = z)/1000)
d2coretop <- cbind.data.frame(lgmID = LGMSppDF$lgmID,
                              minDist = apply(spDist, 2, function(a) min(a)),
                              forcensID = apply(spDist, 2, function(a) modernSppDF$forcensID[which.min(a)])) %>%
  tibble()

# species at nearest core top
sppNearestCoretop <- map_df(d2coretop$forcensID, function(x) modernSppDF[modernSppDF$forcensID %in% x, ]) %>%
  select(Dentigloborotalia_anfracta:Globorotalia_menardii_Globorotalia_tumida)

# get LGM species only
LGMSppOnly <- select(LGMSppDF, Dentigloborotalia_anfracta:Globorotalia_menardii_Globorotalia_tumida)

# LGM-modern turnover (Bray-Curtis)
turnover <- tibble(select(LGMSppDF, lgmID, sampleID, Longitude, Latitude),
                   minDist_km = d2coretop$minDist,
                   turnover_BC = rowSums(abs((LGMSppOnly-sppNearestCoretop)))/rowSums(LGMSppOnly+sppNearestCoretop)
)

# clean up a bit again
rm(list = c("spDist", "d2coretop", "sppNearestCoretop", "LGMSppOnly"))

# ==============================================================================================================================

# assess reconstructed LGM temperatures: anomalies ================================================================

# does it matter if we use anomalies vs atlas temperatures or vs temperatures from (the nearest) core tops?

# anomaly from atlas temperature
LGM_mean <- LGM_rec_temperatures %>%
  group_by(Core, Latitude, Longitude, Ocean, lgmID) %>%
  summarise(T_LGM_foram = mean(predT), .groups = "drop") %>%
  left_join(., readRDS("data/LGM_atlas_temperatures.RDS") %>% 
              select(uniqueID, annual_50) %>%
              rename(lgmID = uniqueID, T_atlas = annual_50),
            by = "lgmID"
  ) %>%
  mutate(T_anomaly_foram_atlas = T_LGM_foram - T_atlas)

# anomaly vs temperature predicted from nearby core top
# get meta data from core top
# split again because only search for nearest core top within region
# to make sure that TF model exactly same for LGM and core top temperature estimate
LGM_mean_split <- split(LGM_mean, LGM_mean$Ocean)

# minimum distance to core top
coretop_dist <- mapply(function(x, y){
  sites <- unique(x[, c("lgmID", "Longitude", "Latitude")])
  distmatrix <- apply(sites[, c("Longitude", "Latitude")], 1, function(z) distGeo(p1 = y$meta[, c("Longitude", "Latitude")], p2 = z)/1000)
  minDist2CoreTop <- apply(distmatrix, 2, function(a) min(a))
  forcensID <- apply(distmatrix, 2, function(a) y$meta$forcensID[which.min(a)])
  cbind.data.frame(minDist2CoreTop, forcensID, lgmID = sites$lgmID)
}, x = LGM_mean_split, y = readRDS("data/ForCenS_spp.RDS"), SIMPLIFY = FALSE)

# get temperature predicted from nearest core tops
checkAnoms <- map2_df(coretop_dist, coretop_spp_temp, function(x, y){
  y %>%
    select(forcensID, annual_50) %>%
    left_join(x, ., by = "forcensID") %>%
    select(-forcensID) %>%
    rename(T_foram_coretop = annual_50)
}) %>%
  left_join(LGM_mean, ., by = "lgmID") %>%
  mutate(T_anomaly_foram_coretop = T_LGM_foram - T_foram_coretop,
         diffAnomaly = T_anomaly_foram_coretop - T_anomaly_foram_atlas)


checkAnoms %>%
  filter(minDist2CoreTop <= 100) %>%
  summarise(mean(diffAnomaly), sd(diffAnomaly))

checkAnoms %>%
  filter(minDist2CoreTop <= 100) %>%
  ggplot(aes(diffAnomaly)) +
  geom_histogram()

# the difference is small, random and has no spatial pattern --> we can use the entire data set and use anomalies vs atlas

rm(list = c("checkAnoms", "coretop_dist", "LGM_mean_split", "LGM_mean"))

# ==============================================================================================================================

# assess reconstructed LGM temperatures: analogues ================================================================


# fraction of no analogue samples using global threshold
LGM_rec_temperatures %>%
  mutate(goodAnalogue = sqcd <= globalNoAnalogueValue) %>%
  summarise(sum(goodAnalogue)/n())

# median dissimilarity to most similar core tops
LGM_rec_temperatures %>%
  summarise(median(sqcd))

# assess effect of analogue quality and save LGM temperatures
LGM_rec_temperatures %>%
  ggplot(aes(sqcd)) +
  geom_histogram()

LGM_mean_NoFiltering <- LGM_rec_temperatures %>%
  left_join(.,
            turnover %>%
              select(sampleID, turnover_BC, minDist_km),
            by = "sampleID"
            ) %>%
  group_by(lgmID, Core, Latitude, Longitude, Ocean) %>%
  summarise(T_LGM_foram = mean(predT),
            turnover = mean(turnover_BC),
            minDist_km = mean(minDist_km),
            n = n(),
            .groups = "drop") %>%
  left_join(., readRDS("data/LGM_atlas_temperatures.RDS") %>% 
              select(uniqueID, annual_50) %>%
              rename(lgmID = uniqueID, T_atlas = annual_50),
            by = "lgmID"
  ) %>%
  mutate(T_anomaly_foram_atlas = T_LGM_foram - T_atlas)

LGM_mean_GlobalFiltering <- LGM_rec_temperatures %>%
  left_join(.,
            turnover %>%
              select(sampleID, turnover_BC, minDist_km),
            by = "sampleID"
  ) %>%
  filter(sqcd <= globalNoAnalogueValue) %>% 
  group_by(lgmID, Core, Latitude, Longitude, Ocean) %>%
  summarise(T_LGM_foram = mean(predT),
            turnover = mean(turnover_BC),
            minDist_km = mean(minDist_km),
            n = n(),
            .groups = "drop") %>%
  left_join(., readRDS("data/LGM_atlas_temperatures.RDS") %>% 
              select(uniqueID, annual_50) %>%
              rename(lgmID = uniqueID, T_atlas = annual_50),
            by = "lgmID"
  ) %>%
  mutate(T_anomaly_foram_atlas = T_LGM_foram - T_atlas) 

LGM_mean_RegionalFiltering <- LGM_rec_temperatures %>%
  left_join(.,
            turnover %>%
              select(sampleID, turnover_BC, minDist_km),
            by = "sampleID"
  ) %>%
  filter(sqcd <= analogueCutoff) %>%
  group_by(lgmID, Core, Latitude, Longitude, Ocean) %>%
  summarise(T_LGM_foram = mean(predT),
            turnover = mean(turnover_BC),
            minDist_km = mean(minDist_km),
            n = n(),
            .groups ="drop") %>%
  left_join(., readRDS("data/LGM_atlas_temperatures.RDS") %>% 
              select(uniqueID, annual_50) %>%
              rename(lgmID = uniqueID, T_atlas = annual_50),
            by = "lgmID"
  ) %>%
  mutate(T_anomaly_foram_atlas = T_LGM_foram - T_atlas) 


# how many samples are lost when filtering?
LGM_mean_NoFiltering %>%
  summarise(sum(n), n())

# how many samples and sites are filtered out with global threshold?
LGM_mean_NoFiltering %>%
  summarise(sum(n), n()) -
  LGM_mean_GlobalFiltering %>%
  summarise(sum(n), n())

# with regional threshold?
LGM_mean_NoFiltering %>%
  summarise(sum(n), n()) -
  LGM_mean_RegionalFiltering %>%
  summarise(sum(n), n())

# what is the effect on the temperature anomalies?
LGM_mean_NoFiltering %>%
  reframe(avg = mean(T_anomaly_foram_atlas),
            q = quantile(T_anomaly_foram_atlas, c(0.025, 0.975)))


LGM_mean_GlobalFiltering %>%
  reframe(avg = mean(T_anomaly_foram_atlas),
            q = quantile(T_anomaly_foram_atlas, c(0.025, 0.975)))

LGM_mean_RegionalFiltering %>%
  reframe(avg = mean(T_anomaly_foram_atlas),
            q = quantile(T_anomaly_foram_atlas, c(0.025, 0.975)))

# ==============================================================================================================================

# uncertainty of the reconstructions ==============================================================

# make TF model to assess spatial autocorrelation in the residuals
TFmodel_an50 <- map(coretopSpp, function(x){
  MAT(x$species, x$env$annual_50, k = 10, lean = FALSE)
})

# get residuals
TFM_residuals <- map(TFmodel_an50, function(x) cbind.data.frame(res = x$fitted.values[, ncol(x$fitted.values)] - x$x, x =  x$x))

# add metadata
residuals <- map2(coretopSpp, TFM_residuals,  function(x, y) cbind(select(x$meta, c(forcensID, Longitude, Latitude)), y))

# detrend residuals
residuals <- map(residuals, function(x){
  detRes <- resid(loess(x$res ~ x$x, span = 0.1))
  cbind(x, detRes) %>%
    st_as_sf(coords = c("Longitude", "Latitude"),
             crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    
})

# get variograms
variograms <- residuals %>%
  map(., function(x){
    x %>%
      variogram(detRes ~ 1, ., cutoff = 5000, width = 100)
    }
  )

# fit using Matern model
fit_variograms <- variograms %>%
  map(~fit.variogram(., model = vgm("Mat"), fit.kappa = TRUE))

# make spatial models
spModels <- map2(.x = residuals, .y = fit_variograms, function(x, y){
  gstat(formula = detRes ~ 1, data = x, beta = mean(x$detRes), dummy = TRUE, model = y, nmax = 100)
})

names(spModels) <- c("Indian", "Mediterranean", "North Atlantic", "Pacific", "South Atlantic")

# predict spatially correlated noise
set.seed(123)
spResNoise <- names(spModels) %>%
  map_df(., function(x){
    mod <- spModels[[x]]
    newdata <- LGM_rec_temperatures %>%
          filter(Ocean == x) %>%
          distinct(Latitude, Longitude, Ocean, lgmID) %>%
          st_as_sf(coords = c("Longitude", "Latitude"),
                   crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
                   )
  predict(
    object = mod,
    newdata = newdata,
    nsim = 1000
  ) %>%
    mutate(Ocean = newdata$Ocean,
           lgmID = newdata$lgmID)
}) %>%
  tibble() %>%
  select(-geometry, -Ocean)

# join with reconstructions
recWithNoise <- left_join(LGM_mean_NoFiltering, spResNoise, by = c("lgmID")) %>%
  pivot_longer(contains("sim"), names_to = "ens", values_to = "uncertainty") %>%
  mutate(T_anomaly = T_anomaly_foram_atlas + uncertainty) %>%
  select(lgmID, Longitude, Latitude, ens, T_anomaly)

recWithNoise %>%
  group_by(ens) %>%
  reframe(avgT = mean(T_anomaly)) %>%
  summarise(avg =  mean(avgT),
            q2.5 = quantile(avgT, 0.025),
            q97.5 = quantile(avgT, 0.975))

# export reconstruction with spatially coherent noise
recWithNoise %>%
  mutate(T_anomaly = round(T_anomaly, 3)) %>%
  rename(T_anomaly_degC = T_anomaly) %>%
  write_csv("out/LGM_foram_reconstruction_ensembles.csv")
  
LGM_mean_NoFiltering <- left_join(LGM_mean_NoFiltering, spResNoise, by = c("lgmID")) %>%
  pivot_longer(contains("sim"), names_to = "sim", values_to = "uncertainty") %>%
  group_by(lgmID, Latitude, Longitude, Ocean, T_LGM_foram, turnover, n, minDist_km, T_atlas, T_anomaly_foram_atlas) %>%
  summarise(avgUncertainty = sd(uncertainty)*2, .groups = "drop") %>%
  mutate(uncertainty2sd = avgUncertainty/sqrt(n))

# export mean reconstruction
fosSpp %>%
  map_df(function(x){
    x$meta %>%
      select(lgmID, Core)
  }) %>%
  distinct() %>%
  left_join(., LGM_mean_NoFiltering %>%
              select(lgmID, Longitude, Latitude, T_LGM_foram, T_anomaly_foram_atlas, turnover, n, uncertainty2sd) %>%
              mutate(across(c(T_LGM_foram, T_anomaly_foram_atlas, turnover, uncertainty2sd), ~round(., 2))) %>%
              rename(T_degC = T_LGM_foram,
                     T_anomaly_degC = T_anomaly_foram_atlas),
            by = "lgmID"
            ) %>%
  write_csv("out/LGM_foram_recontruction.csv")

rm(list = c("TFmodel_an50", "TFM_residuals", "residuals", "variograms", "fit_variograms", "spModels", "spResNoise"))
 
# ==============================================================================================================================


# compare simulated and reconstructed temperatures ===================================================
modelsLGM <- readRDS("data/simTemp_LGM_annual50.RDS") %>%
  distinct(model)

notInEnsemble <- c("lgmID", "AWI.ESM_coast_ext", "AWI.ESM_coast", "AWI.ESM_subpolar")

LGM_sim <- readRDS("data/simTemp_LGM_annual50.RDS") %>%
  select(-c(Longitude, Latitude, T.LGM, T.PI, Depth, Season)) %>%
  pivot_wider(names_from = model, values_from = T.anomaly) %>%
  mutate(ensemble = rowMeans(select(., -all_of(notInEnsemble)), na.rm = TRUE)) %>%
  pivot_longer(-lgmID, names_to = "model", values_to = "T_anomaly_sim") %>%
  mutate(T_anomaly_sim = -T_anomaly_sim) # ensure that this is in the right place

# mean temperature anomaly per simulation
LGM_sim %>%
  group_by(model) %>%
  filter(!model %in% c("ensemble", "AWI.coast", "AWI.subpolar", "AWI.MKZ", "AWI.MSSP")) %>%
  summarise(avg = mean(T_anomaly_sim, na.rm = TRUE))

# join reconstructed and simulated data
# no analogue filtering
LGM_noAF <- LGM_sim %>%
  inner_join(., LGM_mean_NoFiltering, by = "lgmID") %>%
  mutate(difDataModel = T_anomaly_sim - T_anomaly_foram_atlas)

# analogue filtering using global threshold
LGM_globalAF <- LGM_sim %>%
  inner_join(., LGM_mean_GlobalFiltering, by = "lgmID") %>%
  mutate(difDataModel = T_anomaly_sim - T_anomaly_foram_atlas)

# analogue filtering using regional threshold
LGM_regionalAF <- LGM_sim %>%
  inner_join(., LGM_mean_RegionalFiltering, by = "lgmID") %>%
  mutate(difDataModel = T_anomaly_sim - T_anomaly_foram_atlas)

# model data diff in NAT
LGM_noAF %>%
  filter(Ocean == "North Atlantic") %>%
  group_by(model) %>%
  summarise(min(T_anomaly_foram_atlas), max(difDataModel, na.rm = TRUE))

# averaged
LGM_noAF %>%
  filter(Ocean == "North Atlantic") %>%
  filter(Latitude >= 40 & Latitude < 60) %>%
  group_by(model) %>%
  summarise(round(mean(T_anomaly_foram_atlas), 1), round(mean(difDataModel, na.rm = TRUE), 1))

# prep for figures
latRangeData <- LGM_noAF %>%
  slice(c(which.min(Latitude), which.max(Latitude))) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform(crs = "+proj=robin") %>%
  mutate(y = purrr::map_dbl(geometry, 2),
         x = 0)

latRangeMap <- tibble(y = c(-90, 90), x = c(0, 0)) %>%
  st_as_sf(coords = c("x", "y"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform(crs = "+proj=robin") %>%
  mutate(y = purrr::map_dbl(geometry, 2))

lat_noAF <- LGM_noAF %>%
  filter(model == "ensemble") %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform(crs = "+proj=robin") %>%
  mutate(y = purrr::map_dbl(geometry, 2))

# analyse hosing exps
RMSEreduction <- left_join(LGM_noAF %>%
                             filter(model == "AWI.ESM_coast_ext" | model =="AWI.ESM_coast" | model == "AWI.ESM_subpolar" | model == "AWI.ESM") %>%
                             group_by(model) %>%
                             summarise(RMSE = sqrt(sum(difDataModel^2)/n())) %>%
                             pivot_wider(values_from = RMSE, names_from = model) %>%
                             mutate(across(everything(), ~round((1-./AWI.ESM)*100))) %>%
                             pivot_longer(-AWI.ESM, names_to = "model", values_to = "RMSEred"),
                           LGM_noAF %>%
                             filter(model == "AWI.ESM_coast_ext" | model =="AWI.ESM_coast" | model == "AWI.ESM_subpolar" | model == "AWI.ESM") %>%
                             filter(Ocean == "North Atlantic") %>%
                             filter(Latitude >= 40 & Latitude < 60) %>%
                             group_by(model) %>%
                             summarise(RMSE = sqrt(sum(difDataModel^2)/n())) %>%
                             pivot_wider(values_from = RMSE, names_from = model) %>%
                             mutate(across(everything(), ~round((1-./AWI.ESM)*100))) %>%
                             pivot_longer(-AWI.ESM, names_to = "model", values_to = "RMSEredNA")
) %>%
  mutate(lab = paste0(RMSEredNA, "%"),
         x = 1,
         y = 50) %>%
  mutate(model = factor(model, labels = c("Newfoundland", "Newfoundland ext.", "Subpolar")))
