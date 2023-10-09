library(tidyverse)
library(RColorBrewer)
library(sf)
# devtools::install_github("eliocamp/tagger")
library(tagger) 
library(patchwork)
library(rnaturalearth)

# theme and map data =================================================================================
theme_set(theme_bw() + 
            theme(line = element_line(colour = "#464646"),
                  rect = element_rect(fill = NA, colour = NA),
                  text = element_text(colour = "#000000", size = 7),
                  axis.ticks = element_line(colour = "#464646"), 
                  legend.key = element_rect(colour = NA, fill = NA), 
                  panel.border = element_rect(colour = "#464646", fill = NA),
                  panel.grid = element_blank(),
                  plot.background = element_blank(), 
                  panel.background = element_blank(),
                  strip.background = element_blank(),
                  plot.title = element_text(hjust = 0.5, size = 7),
                  plot.subtitle = element_text(hjust = 0.5))
)

# edits for maps
theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.background = element_blank(),
                         panel.border = element_blank(),
                         axis.line = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         legend.key = element_blank(),
                         strip.background = element_rect(fill = "transparent"),
                         legend.position = "bottom",
                         legend.key.height = unit(0.2, "cm"),
                         legend.key.width = unit(0.3, "cm"),
                         legend.margin = margin(0, 0, 0, 0),
                         legend.box.margin = margin(t = -1, r = 10, b = 0, l = 10, unit = "mm")
                         )
)


# data for maps
wmap <- ne_download(scale = 110, type = "land", category = "physical", returnclass = "sf")

bbox <- st_sf(geometry = st_sfc(
  st_polygon(x = list(cbind(c(-180, rep(180, 100), rep(-180, 100)),
                            c(-90, seq(-90, 90, length = 100), 
                              seq(90, -90, length = 100))))),
  crs = "WGS84"))

grat <- list(
  st_sf(geometry = st_sfc(
    st_multilinestring(x = map(seq(-180, 180, 30), 
                               function(x) cbind(x, seq(-90, 90, 1)))), crs = "WGS84")),
  st_sf(geometry = st_sfc(
    st_multilinestring(x = map(seq(-90, 90, 30), function(x) {
      cbind(seq(-180, 180, 1), x)
    })), crs = "WGS84"))) %>% 
  bind_rows()

# labels for maps
ylabs <- tibble(x = -180,
                y = seq(-60, 60, 30), 
                label = expression(
                  60 * degree * S,
                  30 * degree * S,
                  0 * degree, 
                  30 * degree * N,
                  60 * degree * N
                )
) %>%
  st_as_sf(
    coords = c("x", "y"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    )

xlabstop <- tibble(y = 90,
                x = seq(-180, 180, 90), 
                label = expression(
                  180 * degree * W,
                  90 * degree * W,
                  0 * degree, 
                  90 * degree * E,
                  180 * degree * E
                )
) %>%
  st_as_sf(
    coords = c("x", "y"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    )

xlabsbot <- tibble(y = -90,
                   x = seq(-180, 180, 90), 
                   label = expression(
                     180 * degree * W,
                     90 * degree * W,
                     0 * degree, 
                     90 * degree * E,
                     180 * degree * E
                   )
) %>%
  st_as_sf(
    coords = c("x", "y"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  )
# ===============================================================================================

# Figure 1 clusters =============================================================================
clustPlot <- clusterOut %>%
  left_join(., 
            map_df(fosSpp, function(x) x$meta %>%
                     select(lgmID, inMARGO)) %>%
              rename(uID = lgmID),
            by = "uID"
            ) %>%
  mutate(inMARGO = !inMARGO) %>%
  mutate(cluster = case_when(cluster == 1 ~ "Transitional",
                             cluster == 2 ~ "Warm",
                             cluster == 3 ~ "Cold")) %>%
  replace_na(list(inMARGO = FALSE)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  ggplot() +
  geom_sf(data = bbox, fill = "grey80", size = 0.3, colour = "#464646") +
  geom_sf(data = grat, colour = "grey20", linewidth = 0.1, linetype = "dotted") +
  geom_sf(aes(fill = cluster, colour = inMARGO), size = 1, shape = 21, stroke = 0.15) +
  geom_sf(data = wmap, fill = "white", colour = "black", size = 0.1) +
  geom_sf_text(data = ylabs, aes(label = label),
               nudge_x = -1e6,
               size = 1.5,
               parse = TRUE,
               hjust = 1) +
  geom_sf_text(data = xlabsbot, 
               aes(label = label),
               nudge_y = -5e5,
               size = 1.5,
               parse = TRUE,
               vjust = 1) +
  coord_sf(crs = "+proj=robin") +
  scale_fill_manual(values = c("dodgerblue4", "darkgoldenrod2", "red4")) +
  scale_colour_manual(values = c("transparent", "white"), guide = "none") +
  facet_wrap(~set) +
  labs(fill = NULL) +
  theme_opts

ggsave("out/Fig1.pdf", plot = clustPlot,  units = "mm", width = 180, height = 70)

# ===============================================================================================

# Figure 2 distance decay in core top and ensemble mean, difference ===============================
Fig2A <- modDistBin %>%
  mutate(model = "modern") %>%
  ggplot() +
  geom_tile(aes(envBin, simBin, fill = lognnorm)) +
  scale_fill_viridis_c(option = "magma", guide = "none") +
  geom_boxplot(data = plotDist %>%
                 mutate(bin = cut(envdist, breaks = seq(0, 35, 1), include.lowest = TRUE, right = FALSE, labels = seq(0, 35, 1)[-1]-0.5*1))
               %>% # only to speed things up, remove later!!
                 slice_sample(n = 1e5)
               ,
               aes(envdist, commdist, group = bin), outlier.colour = NA, fill = NA, colour = "black", size = 0.25) +
  xlim(c(0, 32)) +
  labs(title = "Modern",
       x = NULL,
       y = "Bray-Curtis similarity",
       fill = NULL) +
  theme(legend.box.margin = margin(t = 0, r = 1, b = 0, l = -5, unit = "mm"))


Fig2B <- LGMDist_avgBin %>%
  filter(model == "ensemble") %>%
  ggplot() +
  geom_tile(aes(envBin, simBin, fill = lognnorm)) +
  scale_fill_viridis_c(option = "magma", breaks = c(0, 1), labels = c("min", "max")) +
  geom_boxplot(data = LGMDist_avg %>%
                 filter(model == "ensemble") %>%
                 mutate(bin = cut(envdist, breaks = seq(0, 35, 1), include.lowest = TRUE, right = FALSE, labels = seq(0, 35, 1)[-1]-0.5*1)),
               aes(envdist, commdist, group = bin), outlier.colour = NA, fill = NA, colour = "black", size = 0.25) +
  xlim(c(0, 32)) +
  labs(title = "LGM",
       x = expression("Temperature difference between sites ["*degree*"C]"),
       y = NULL,
       fill = "Number\nof pairs") +
  guides(fill = guide_colourbar(label.theme = element_text(size = 5))) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.height = unit(1.5, "mm"),
        legend.key.width = unit(7, "points"),
        legend.box.margin = margin(t = -2, r = 2, b = 0, l = 2, unit = "mm"),
        legend.position = "bottom")

Fig2C <- diffSimDecayModLGM %>%
  filter(model == "ensemble") %>%
  ggplot(aes(simBin, envBin)) +
  geom_tile(aes(fill = diff)) +
  geom_step(data = ecolim, aes(simBinPlot, Tlim), colour = "black", linewidth = 0.25) +
  scale_fill_distiller(palette = "RdGy", limits = c(-.8, .8), breaks = c(-0.8, 0.8), labels = c("less", "more")) +
  ylim(c(0, 32)) +
  labs(title = "Difference",
       fill = NULL,
       y = NULL,
       x = NULL) +
  coord_flip() +
  guides(fill = guide_colourbar(label.theme = element_text(size = 5))) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.height = unit(1.5, "mm"),
        legend.key.width = unit(7, "points"),
        legend.box.margin = margin(t = -2, r = 5, b = 0, l = 5, unit = "mm"),
        legend.position = "bottom")

Fig2 <- Fig2A + Fig2B + Fig2C + plot_annotation(tag_levels = "a")

ggsave("out/Fig2.pdf", plot = Fig2,  units = "mm", width = 180, height = 80)

# ===============================================================================================

# Figure 3 turnover, reconstructions ======================================================================================

# turnover
mapTurnover <- LGM_noAF %>%
  filter(minDist_km <= 100) %>%
  filter(model == "ensemble") %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  ggplot() +
  geom_sf(data = bbox, fill = "grey80", size = 0.3, colour = "#464646") +
  geom_sf(data = grat, colour = "grey20", linewidth = 0.1, linetype = "dotted") +
  geom_sf(aes(colour =  turnover), size = 0.3) +
  geom_sf(data = wmap, fill = "white", colour = "black", size = 0.1) +
  geom_sf_text(data = ylabs, aes(label = label),
               nudge_x = -1e6,
               size = 1.5,
               parse = TRUE,
               hjust = 1) +
  geom_sf_text(data = xlabstop, aes(label = label),
               size = 1.5,
               parse = TRUE,
               vjust = -1) +
  coord_sf(crs = "+proj=robin",
           clip = "off") +
  scale_colour_viridis_c(option = "magma", limits = c(0, 1), breaks = c(1, .5, 0)) +
  labs(colour = expression("Turnover   ")) +
  theme_opts +
  theme(legend.key.height = unit(0.1, "cm"),
        legend.box.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "mm"))

latTurnover <- lat_noAF %>%
  filter(minDist_km <= 100) %>%
  ggplot(aes(turnover, y)) +
  geom_smooth(orientation = "y", method = "gam", colour = "black", fill = "grey60", alpha = .5, method.args = list(family = binomial), linewidth = 0.7, lineend = "round") + 
  coord_cartesian(xlim = c(1, 0), ylim = latRangeMap$y) +
  scale_x_continuous(breaks = c(1, 0.5, 0)) +
  labs(x = NULL) +
  theme(axis.line.x = element_line(colour = "#464646", linewidth = 0.3),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank())

plotTurnover <- mapTurnover + latTurnover + plot_layout(nrow = 1, widths = unit(c(3, -1), c("null","null")))

# reconstructed temperature anomaly
xlimRec <- c(-7, 2)
xlimDif <- c(-4, 5)

# cannot get geom_tile to work with x labels
mapTrec <- LGM_noAF %>%
  filter(model == "ensemble") %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  ggplot() +
  geom_sf(data = bbox, fill = "grey80", size = 0.3, colour = "#464646") +
  geom_sf_text(data = xlabstop, aes(label = label),
               size = 1.5,
               parse = TRUE,
               vjust = -1) +
  geom_tile(data = readRDS("data/gridded_LGM_temp_rec.RDS") %>%
              drop_na() #%>%
              # slice_sample(n = 300) # remove for final figure
            ,
            aes(x, y, width = 1, height = 1, fill = sst_anomaly)) +
  geom_sf(data = grat, colour = "grey20", linewidth = 0.1, linetype = "dotted") +
  geom_sf(data = wmap, fill = "white", colour = "black", size = 0.1) +
  geom_sf(aes(fill =  T_anomaly_foram_atlas), size = 0.7, shape = 21, stroke = 0.2, colour = "#464646") +
  coord_sf(crs = "+proj=robin",
           default_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
           clip = "off"
           ) +
  scale_fill_gradient2(low = "dodgerblue1", high = "red2", limits = c(-14, 4), na.value = "grey90") + # can tinker with na.value to visualise areas that have not been gridded
  labs(fill = expression("Reconstructed\ntemperature anomaly ["*degree*"C]")) +
  theme_opts +
  theme(legend.key.height = unit(0.1, "cm"),
        legend.box.margin = margin(t = -7, r = 0, b = 0, l = 0, unit = "mm"))


latTrec <- lat_noAF %>%
  ggplot(aes(T_anomaly_foram_atlas, y)) +
  geom_line(data = latRangeData, aes(x, y), linetype = "dashed") +
  geom_smooth(orientation = "y", method = "gam", colour = "black", fill = "grey60", alpha = .5, lineend = "round", linewidth = 0.7) +
  coord_cartesian(xlim = xlimRec, ylim = latRangeMap$y) +
  labs(x = NULL) +
  theme(axis.line.x = element_line(colour = "#464646", linewidth = 0.3),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank())

plotReconstruction <- mapTrec + latTrec + plot_layout(nrow = 1, widths = unit(c(3, -1), c("null","null")))

# modelled temperature anomaly
mapTmod <-  LGM_noAF %>%
  filter(model == "ensemble") %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  ggplot() +
  geom_sf(data = bbox, fill = "grey80", size = 0.3, colour = "#464646") +
  geom_sf(data = grat, colour = "grey20", linewidth = 0.1, linetype = "dotted") +
  geom_sf(aes(colour =  T_anomaly_sim), size = 0.3) +
  geom_sf(data = wmap, fill = "white", colour = "black", size = 0.1) +
  geom_sf_text(data = ylabs, aes(label = label),
               nudge_x = -1e6,
               size = 1.5,
               parse = TRUE,
               hjust = 1) +
  geom_sf_text(data = xlabstop, aes(label = label),
               size = 1.5,
               parse = TRUE,
               vjust = -1) +
  coord_sf(crs = "+proj=robin",
           clip = "off") +
  scale_colour_gradient2(low = "dodgerblue1", high = "red2", limits = c(-14, 4), na.value = "grey90") +
  labs(colour = expression("Simulated\ntemperature anomaly ["*degree*"C]"),
       y = "No analogue filtering") +
  theme_opts  +
  theme(legend.key.height = unit(0.1, "cm"),
        legend.box.margin = margin(t = -7, r = 0, b = 0, l = 0, unit = "mm"))

latTmod <- lat_noAF %>%
  ggplot(aes(T_anomaly_sim, y)) +
  geom_line(data = latRangeData, aes(x, y), linetype = "dashed") +
  geom_smooth(orientation = "y", method = "gam", colour = "black", fill = "grey60", alpha = .5, linewidth = 0.7, lineend = "round") +
  coord_cartesian(xlim = xlimRec, ylim = latRangeMap$y) +
  labs(x = NULL) +
  theme(axis.line.x = element_line(colour = "#464646", linewidth = 0.3),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank())

plotSimul <- mapTmod + latTmod  + plot_layout(nrow = 1, widths = unit(c(3, -1), c("null","null")))

# difference between reconstruction and simulation
mapDifMod <- LGM_noAF %>%
  filter(model == "ensemble") %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  ggplot() +
  geom_sf(data = bbox, fill = "grey80", size = 0.3, colour = "#464646") +
  geom_sf(data = grat, colour = "grey20", linewidth = 0.1, linetype = "dotted") +
  geom_sf(aes(colour =  difDataModel), size = 0.3) +
  geom_sf(data = wmap, fill = "white", colour = "black", size = 0.1) +
  geom_sf_text(data = ylabs, aes(label = label),
               nudge_x = -1e6,
               size = 1.5,
               parse = TRUE,
               hjust = 1) +
  geom_sf_text(data = xlabstop, aes(label = label),
               size = 1.5,
               parse = TRUE,
               vjust = -1) +
  coord_sf(crs = "+proj=robin",
           clip = "off") +
  scale_colour_gradient2(low = "dodgerblue1", high = "red2", limits = c(-7, 11), na.value = "grey90") + # can tinker with na.value to visualise areas that have not been gridded
  labs(colour = expression("Simulated - reconstructed\ntemperature anomaly ["*degree*"C]")) +
  theme_opts +
  theme(legend.key.height = unit(0.1, "cm"),
        legend.box.margin = margin(t = -7, r = 0, b = 0, l = 0, unit = "mm"),
        plot.margin = margin(t = -40, r = 5.5, b = 5.5, l = 5.5, unit = "pt"))

latDifMod <-  lat_noAF %>%
  ggplot(aes(difDataModel, y)) +
  geom_line(data = latRangeData, aes(x, y), linetype = "dashed") +
  geom_smooth(orientation = "y", method = "gam", colour = "black", fill = "grey60", alpha = .5, linewidth = 0.7, lineend = "round") +
  coord_cartesian(xlim = xlimDif, ylim = latRangeMap$y) +
  labs(x = NULL) +
  theme(axis.line.x = element_line(colour = "#464646", linewidth = 0.3),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank())

plotDifDM <- mapDifMod + latDifMod + plot_layout(nrow = 1, widths = unit(c(3, -1), c("null","null")))

Fig3 <- plotTurnover/plotReconstruction/plotSimul/plotDifDM + 
  plot_annotation(tag_levels = "a") + 
  plot_layout(widths = unit(c(3, -1), c("null","null")))

ggsave("out/Fig3.pdf", plot = Fig3,  units = "mm", width = 88, height = 180)

# ===============================================================================================

# Figure 4 comparison with hosing simulations ===================================================

# distance decay in hsing simulations
hosing_DD <- LGMDist_avgBin %>%
  filter(model %in% c("AWI.ESM_coast", "AWI.ESM_subpolar", "AWI.ESM_coast_ext")) %>%
  ggplot(., aes(envBin, simBin)) +
  geom_tile(aes(fill = lognnorm)) +
  scale_fill_viridis_c(option = "magma", breaks = c(0, 1), labels = c("min", "max")) +
  geom_boxplot(data = LGMDist_avg %>%
                 filter(model %in% c("AWI.ESM_coast", "AWI.ESM_subpolar", "AWI.ESM_coast_ext")) %>%
                 mutate(bin = cut(envdist, breaks = seq(0, 35, 1), include.lowest = TRUE, right = FALSE, labels = seq(0, 35, 1)[-1]-0.5*1)),
               aes(envdist, commdist, group = bin), outlier.colour = NA, fill = NA, colour = "black", size = 0.25) +
  labs(title = "LGM",
       x = expression("Temperature\ndifference between sites ["*degree*"C]"),
       y = "Bray-Curtis similarity",
       fill = "Number\nof pairs   ") +
  guides(fill = guide_colourbar(label.theme = element_text(size = 5))) +
  facet_wrap(.~model, ncol = 1) +
  theme(axis.title.x.bottom = element_text(margin = margin(t = 10, r = 2, b = 0, l = 2, unit = "points")),
        legend.key.height = unit(1.5, "mm"),
        legend.key.width = unit(10, "points"),
        legend.position = "bottom",
        strip.text = element_blank())

# draw rectangle around tiles that show excess sim when using LGM simulation
# the idea is to draw segments where the first derivative changes
# define block using threshold
block <- LGMDist_avgBin %>%
  filter(grepl("AWI", model)) %>%
  split(., .$model) %>% # take this detour of splitting by model because many models have NAs for some bins
  map_dfr(., ~full_join(modDistBin, .,  by = c("envBin", "simBin")), .id = "model") %>%
  replace_na(list(lognnorm.x = 0, lognnorm.y = 0)) %>%
  mutate(diff = lognnorm.y - lognnorm.x) %>%
  filter(model == "AWI.ESM") %>%
  mutate(block = diff > 0.1) %>%
  select(envBin, simBin, block) %>%
  pivot_wider(names_from = envBin, values_from = block) %>%
  mutate(`31.75` = FALSE)

# convert to matrix
blockM <- as.matrix(block[, -1])
blockM[is.na(blockM)] <- FALSE

# vertical walls
BMv <- cbind(FALSE, blockM)
vw <- map_df(seq(nrow(BMv)), function(i){
  tibble(y = block$simBin[i] - resSim/2,
         yend = block$simBin[i] + resSim/2,
         x = as.numeric(colnames(blockM)[which(diff(BMv[i,]) != 0)]) - resEnv/2,
         xend = x
  )
})

# horizontal walls
BMh <- rbind(FALSE, blockM, FALSE)
hw <- map_df(seq(ncol(BMh)), function(i){
  indx <- which(diff(BMh[,i]) != 0)
  EnvBreaks <- c(0, as.numeric(colnames(blockM)) + resEnv/2)
  simBreaks <- c(0, block$simBin + resSim/2)
  tibble(x = EnvBreaks[i],
         xend = EnvBreaks[i + 1],
         y = simBreaks[indx],
         yend = y
  )
})


hosing_DDdiff <- LGMDist_avgBin %>%
  filter(grepl("AWI", model)) %>%
  split(., .$model) %>% # take this detour of splitting by model because many models have NAs for some bins
  map_dfr(., ~full_join(modDistBin, .,  by = c("envBin", "simBin")), .id = "model") %>%
  replace_na(list(lognnorm.x = 0, lognnorm.y = 0)) %>%
  mutate(diff = lognnorm.y - lognnorm.x) %>%
  filter(model %in% c("AWI.ESM_coast", "AWI.ESM_subpolar", "AWI.ESM_coast_ext")) %>%
  ggplot(., aes(envBin, simBin)) +
  geom_tile(aes(fill = diff)) +
  geom_segment(data = vw, aes(x = x, xend = xend, y = y, yend = yend), colour = "black", linewidth = 0.2, lineend = "round") +
  geom_segment(data = hw, aes(x = x, xend = xend, y = y, yend = yend), colour = "black", linewidth = 0.2, lineend = "round") +
  scale_fill_distiller(palette = "RdGy", limits = c(-.8, .8), breaks = c(-0.8, 0.8), labels = c("less", "more")) +
  labs(title = "LGM - modern",
       x = expression("Temperature\ndifference between sites ["*degree*"C]"),
       y = NULL,
       fill = "Difference   ") +
  guides(fill = guide_colourbar(label.theme = element_text(size = 5))) +
  theme(axis.title.x.bottom = element_text(margin = margin(t = 10, r = 2, b = 0, l = 2, unit = "points")),
        legend.key.height = unit(1.5, "mm"),
        legend.key.width = unit(10, "points"),
        legend.position = "bottom",
        strip.text = element_blank()
  ) +
  facet_wrap(.~model, ncol = 1)

# labels
a1 <- tibble(x = 2,
             y = -30,
             model = factor("Subpolar"),
             text = "equilibrium")

a2 <- tibble(x = -2,
             y = -30,
             model = factor("Subpolar"),
             text = "hosing")

latHosing <- LGM_noAF %>%
  filter(model == "AWI.ESM_coast_ext" | model =="AWI.ESM_coast" | model == "AWI.ESM_subpolar" ) %>%
  mutate(model = factor(model, labels = c("Newfoundland", "Newfoundland ext.", "Subpolar"))) %>%
  ggplot(aes(difDataModel, Latitude)) +
  geom_line(data = tibble( x = 0, y = range(LGM_noAF$Latitude)),
            aes(x, y),
            linetype = "dashed", inherit.aes = FALSE) +
  geom_smooth(data = LGM_noAF %>%
                filter(model == "AWI.ESM") %>%
                select(-model),
              aes(difDataModel, Latitude),
              orientation = "y", colour = "#F4A582", fill = "#F4A582", se = TRUE, method = "gam", lineend = "round") +
  geom_smooth(orientation = "y", method = "gam", colour = "black", fill = "grey60", alpha = .5, lineend = "round") +
  labs(title = "Simulated -\nreconstructed",
       x = expression("Simulated - reconstructed\ntemperature anomaly ["*degree*"C]")) +
  facet_wrap(~model, ncol = 1, strip.position = "right") +
  geom_text(data = a1, aes(x, y, label = text), colour = "#F4A582", size = 1.8, hjust = 0) +
  geom_text(data = a2, aes(x, y, label = text), size = 1.8, hjust = 0) +
  geom_text(data = RMSEreduction, aes(x, y, label = lab), size = 1.8, hjust = 0) +
  theme(axis.title.x.bottom = element_text(margin = margin(t = 10, r = 2, b = 0, l = 2, unit = "points")),
        strip.text = element_text(face = "bold", size = 6))


Fig4 <- hosing_DD + hosing_DDdiff +  latHosing + plot_annotation(tag_levels = "a") + plot_layout(widths = c(3, 3, 2))

ggsave("out/Fig4.pdf", plot = Fig4,  units = "mm", width = 180, height = 150)

# ===============================================================================================

# EDF 1 LGM similarity decay in PMIP simulations ==============================================
LGM_DD_models <- LGMDist_avgBin %>%
  filter(!model %in% c("ensemble", "AWI.ESM_coast", "AWI.ESM_subpolar", "AWI.ESM_coast_ext")) %>%
  ggplot(., aes(envBin, simBin)) +
  geom_tile(aes(fill = lognnorm)) +
  scale_fill_viridis_c(option = "magma", breaks = c(0, 1), labels = c("min", "max")) +
  geom_boxplot(data = LGMDist_avg %>%
                 filter(!model %in% c("ensemble", "AWI.ESM_coast", "AWI.ESM_subpolar", "AWI.ESM_coast_ext")) %>%
                 mutate(bin = cut(envdist, breaks = seq(0, 35, 1), include.lowest = TRUE, right = FALSE, labels = seq(0, 35, 1)[-1]-0.5*1)),
               aes(envdist, commdist, group = bin), outlier.colour = NA, fill = NA, colour = "black", size = 0.25) +
  xlim(c(0, 32)) +
  labs(x = expression("Temperature difference between sites ["*degree*"C]"),
       y = "Bray-Curtis similarity",
       fill = "Number of pairs") +
  facet_wrap(.~model, nrow = 2) +
  theme(legend.key.width = unit(1, "mm"),
        legend.key.height = unit(0.25, "cm"))

LGM_DD_models_diff <- diffSimDecayModLGM %>%
  filter(!model %in% c("ensemble", "AWI.ESM_coast", "AWI.ESM_subpolar", "AWI.ESM_coast_ext")) %>%
  ggplot(., aes(envBin, simBin)) +
  geom_tile(aes(fill = diff)) +
  scale_fill_distiller(palette = "RdGy", limits = c(-.8, .8), breaks = c(-0.8, 0.8), labels = c("less", "more")) +
  xlim(c(0, 32)) +
  labs(x = expression("Temperature difference between sites ["*degree*"C]"),
       y = "Bray-Curtis similarity",
       fill = NULL) +
  facet_wrap(.~model, nrow = 2) +
  theme(legend.key.width = unit(1, "mm"),
        legend.key.height = unit(0.25, "cm"))



ExtF1 <- LGM_DD_models / LGM_DD_models_diff + plot_annotation(tag_levels = "a")

ggsave("out/Ext_Fig1.eps", plot = ExtF1,  units = "mm", width = 184, height = 180)

# ===============================================================================================

# EDF 2 sites outside modern limits =============================================================
extraLimitMap <- bind_cols(LGMSppDF_avg[extraLimitSites$site, 1:3], n = extraLimitSites$n) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  ggplot() +
  geom_sf(data = bbox, fill = "grey80", size = 0.3, colour = "#464646") +
  geom_sf(data = grat, colour = "grey20", linewidth = 0.1, linetype = "dotted") +
  geom_sf(aes(alpha = n), size = 0.7, colour = "#491078", shape = 16) +
  geom_sf(data = wmap, fill = "white", colour = "black", size = 0.1) +
  geom_sf_text(data = ylabs, aes(label = label),
               nudge_x = -1e6,
               size = 1.5,
               parse = TRUE,
               hjust = 1) +
  geom_sf_text(data = xlabstop, 
               aes(label = label),
               # nudge_y = -5e5,
               size = 1.5,
               parse = TRUE,
               vjust = -1) +
  coord_sf(crs = "+proj=robin",
           clip = "off") +
  labs(colour = NULL) +
  theme_opts +
  theme(legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"),
        legend.key.width = unit(1, "mm")
  )


ggsave("out/Ext_Fig2.tiff", plot = extraLimitMap,  units = "mm", width = 88, height = 60)

# ===============================================================================================

# EDF 3 distance decay pattern using simulated PI temperature ====================================
DD_PI <- PI_Dist_Bin %>%
  filter(model != "annual_50" & model != "ensemble") %>%
  ggplot(., aes(envBin, simBin)) +
  geom_tile(aes(fill = lognnorm)) +
  geom_boxplot(data = PI_Dist %>%
                 mutate(model = factor(model)) %>%
                 filter(model != "annual_50" & model != "ensemble") %>%
                 mutate(bin = cut(envdist, breaks = seq(0, 35, 1), include.lowest = TRUE, right = FALSE, labels = seq(0, 35, 1)[-1]-0.5*1)),
               aes(envdist, commdist, group = bin), outlier.colour = NA, fill = NA, colour = "black", size = 0.25) +
  scale_fill_viridis_c(option = "magma", breaks = c(0, 1), labels = c("min", "max")) +
  xlim(c(0, 32)) +
  labs(x = expression("Temperature difference between sites ["*degree*"C]"),
       y = "Bray-Curtis similarity",
       fill = "Number of pairs") +
  facet_wrap(.~model, nrow = 2) +
  theme(legend.key.width = unit(1, "mm"),
        legend.key.height = unit(0.25, "cm")) +
  tag_facets(tag_pool = letters[1:8], position = "tr")

# difference from pattern using atlas temperatures
DD_PI_dif <- PI_Dist_Bin %>%
  select(model, envBin, simBin, lognnorm) %>%
  pivot_wider(., names_from = model, values_from = lognnorm) %>%
  replace(is.na(.), 0) %>%
  mutate(across(-c(1, 2), ~ . - annual_50)) %>%
  select(-annual_50) %>%
  pivot_longer(-c(envBin, simBin), values_to = "diff", names_to = "model") %>%
  mutate(model = factor(model, levels = c("modern", modelsPI))) %>%
  filter(model != "Tierney2020" & model != "ensemble") %>%
  ggplot(., aes(envBin, simBin)) +
  geom_tile(aes(fill = diff)) +
  scale_fill_distiller(palette = "RdGy", limits = c(-.8, .8), breaks = c(-0.8, 0.8), labels = c("less", "more")) +
  xlim(c(0, 32)) +
  labs(x = expression("Temperature difference between sites ["*degree*"C]"),
       y = "Bray-Curtis similarity",
       fill = NULL) +
  facet_wrap(.~model, nrow = 2) +
  theme(legend.key.width = unit(1, "mm"),
        legend.key.height = unit(0.25, "cm")) +
  tag_facets(tag_pool = letters[9:16], position = "tr")

DD_PI_plot <- DD_PI / DD_PI_dif

ggsave("out/Ext_Fig3.eps", plot = DD_PI_plot,  units = "mm", width = 184, height = 150)

# ===============================================================================================

# EDF 4 LGM assemblages and modern temp =========================================================
LGM_atlasDD <- LGM_atlasDist %>%
  ggplot(., aes(envBin, simBin)) +
  geom_tile(aes(fill = lognnorm)) +
  scale_fill_viridis_c(option = "magma", breaks = c(0, 1), labels = c("min", "max")) +
  geom_boxplot(data = LGMSppDF_avg %>%
                 left_join(., readRDS("data/LGM_atlas_temperatures.RDS") %>%
                             tibble() %>%
                             dplyr::select(uniqueID, annual_50) %>%
                             rename(lgmID = uniqueID), by = "lgmID") %>%
                 makeDist() %>%
                 mutate(bin = cut(envdist, breaks = seq(0, 35, 1), include.lowest = TRUE, right = FALSE, labels = seq(0, 35, 1)[-1]-0.5*1)),
               aes(envdist, commdist, group = bin), outlier.colour = NA, fill = NA, colour = "black", size = 0.25) +
  xlim(c(0, 32)) +
  labs(x = NULL,
       y = "Bray-Curtis similarity", 
       fill = "Number\nof pairs") +
  guides(fill = guide_colourbar(label.theme = element_text(size = 5))) +
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(0.1, "cm"),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm")
  )

# difference from modern  
PLotDiffDD_atlasModern <- full_join(LGM_atlasDist, modDistBin, by = c("envBin", "simBin")) %>%
  replace_na(list(lognnorm.x = 0, lognnorm.y = 0)) %>%
  mutate(diff = lognnorm.x - lognnorm.y)  %>%
  ggplot(., aes(simBin, envBin)) +
  geom_tile(aes(fill = diff)) +
  geom_step(data = ecolim, aes(simBinPlot, Tlim), colour = "black", size = 0.3) +
  scale_fill_distiller(palette = "RdGy", limits = c(-.8, .8), guide = "none") +
  ylim(c(0, 32)) +
  labs(y = expression("Temperature difference between sites ["*degree*"C]"),
       x = NULL
  ) +
  coord_flip() 

# difference from LGM ensemble
PLotDiffDD_atlasLGM <- full_join(LGM_atlasDist,
                                 LGMDist_avgBin %>%
                                   filter(model == "ensemble") %>%
                                   select(-model),
                                 by = c("envBin", "simBin")) %>%
  replace_na(list(lognnorm.x = 0, lognnorm.y = 0)) %>%
  mutate(diff = lognnorm.x - lognnorm.y)  %>%
  ggplot(., aes(simBin, envBin)) +
  geom_tile(aes(fill = diff)) +
  geom_step(data = ecolim, aes(simBinPlot, Tlim), colour = "black", size = 0.3) +
  scale_fill_distiller(palette = "RdGy", limits = c(-.8, .8), breaks = c(-0.8, 0.8), labels = c("less", "more")) +
  guides(fill = guide_colourbar(label.theme = element_text(size = 5))) +
  ylim(c(0, 32)) +
  labs(fill = NULL,
       y = NULL,
       x = NULL
  ) +
  coord_flip() +
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(0.1, "cm"),
        axis.title.x = element_text(vjust = -3),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"))

ExtF4 <- LGM_atlasDD + PLotDiffDD_atlasModern + PLotDiffDD_atlasLGM + plot_layout(nrow = 1) + plot_annotation(tag_levels = "a")

ggsave("out/Ext_Fig4.eps", plot = ExtF4,  units = "mm", width = 180, height = 80)

# ===============================================================================================

# EDF 5 LGM temperature gradients and difference with models ====================================
tgradients <- LGM_noAF %>%
  filter(!model %in% c("ensemble", "AWI.ESM_coast", "AWI.ESM_subpolar", "AWI.ESM_coast_ext")) %>%
  mutate(T_sim_LGM = T_atlas + T_anomaly_sim,
         model = factor(model)) %>%
  pivot_longer(c(T_LGM_foram, T_atlas, T_sim_LGM), names_to = "group", values_to = "Temperature") %>%
  mutate(latBin = cut(Latitude, breaks = seq(-90, 90, 10), include.lowest = TRUE, right = FALSE, labels = seq(-90, 90, 10)[-1] - 0.5*10)) %>%
  group_by(group, latBin, model) %>%
  summarise(AvgT = mean(Temperature, na.rm = TRUE),
            n = n(),
            seAvgT = sd(Temperature, na.rm= TRUE)/sqrt(n),
            .groups = "drop") %>%
  mutate(latBin = as.numeric(levels(latBin))[latBin],
         Group = factor(group, levels = c("T_atlas", "T_LGM_foram", "T_sim_LGM"), labels = c("Modern", "LGM reconstruction", "LGM simulation"))
         ) %>%
  ggplot(aes(AvgT, latBin, colour = Group)) +
  geom_path() +
  geom_point(size = 1) +
  geom_errorbarh(aes(xmin = AvgT - seAvgT, xmax = AvgT + seAvgT), height = 2) +
  scale_colour_manual(values = c("red2", "grey80", "dodgerblue1")) +
  labs(x = expression("Temperature ["*degree*"C]"),
       y = "Latitude",
       colour = NULL) +
  facet_wrap(~model, nrow = 2) +
  theme(legend.position = "top")

mapDiffDM <- LGM_noAF %>%
  filter(! model %in% c("ensemble", "AWI.ESM_subpolar", "AWI.ESM_coast", "AWI.ESM_coast_ext")) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  ggplot() +
  geom_sf(data = bbox, fill = "grey80", size = 0.3, colour = "#464646") +
  geom_sf(data = grat, colour = "grey20", linewidth = 0.1, linetype = "dotted") +
  geom_sf(aes(colour =  difDataModel), size = 0.3) +
  geom_sf(data = wmap, fill = "white", colour = "black", size = 0.1) +
  coord_sf(crs = "+proj=robin",
           clip = "off") +
  scale_colour_gradient2(low = "dodgerblue1", high = "red2") +
  labs(colour = expression("Simulated - reconstructed temperature anomaly ["*degree*"C]")) +
  facet_wrap(.~model, nrow = 2) +
  theme_opts +
  theme(legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.15, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(t = -2, r = 10, b = 0, l = 10, unit = "mm"))

ExtF5 <- tgradients / mapDiffDM + plot_annotation(tag_levels = "a")

ggsave("out/Ext_Fig5.eps", plot = ExtF5,  units = "mm", width = 180, height = 180)

# ===============================================================================================

# EDF 6 LGM temperature anomalies in hosing simulations ========================================
Ext_Fig6 <- LGM_noAF %>%
  filter(model == "AWI.ESM_coast_ext" | model =="AWI.ESM_coast" | model == "AWI.ESM_subpolar" ) %>%
  mutate(model = factor(model, labels = c("Newfoundland", "Newfoundland ext.", "Subpolar"))) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  ggplot() +
  geom_sf(data = bbox, fill = "grey80", size = 0.3, colour = "#464646") +
  geom_sf(data = grat, colour = "grey20", linewidth = 0.1, linetype = "dotted") +
  geom_sf(aes(colour =  difDataModel), size = 0.3) +
  geom_sf(data = wmap, fill = "white", colour = "black", size = 0.1) +
  geom_sf_text(data = ylabs, aes(label = label),
               nudge_x = -1e6,
               size = 1.5,
               parse = TRUE,
               hjust = 1) +
  geom_sf_text(data = xlabsbot, 
               aes(label = label),
               nudge_y = -5e5,
               size = 1.5,
               parse = TRUE,
               vjust = 1) +
  coord_sf(crs = "+proj=robin") +
  scale_colour_gradient2(low = "dodgerblue1", high = "red2") + 
  facet_wrap(.~model, ncol = 3) +
  labs(colour = expression("Simulated - reconstructed temperature anomaly ["*degree*"C]    ")) +
  theme_opts +
  theme(legend.key.height = unit(1.5, "mm"),
        legend.key.width = unit(10, "points"),
        legend.box.margin = margin(t = 0, r = 2, b = 0, l = 2, unit = "mm"),
        legend.position = "bottom")


ggsave("out/Ext_Fig6.eps", plot = Ext_Fig6,  units = "mm", width = 180, height = 100)

# ===============================================================================================

# EDF 7 how robust is the modern pattern against effect of subsetting ===========================
modernDD <- bind_rows(
  modDistSubBin %>%
    mutate(group = "LGM <100km"),
  modDistSub2Bin %>%
    mutate(group = "Nearest to LGM")
) %>%
  ggplot() +
  geom_tile(aes(envBin, simBin, fill = lognnorm)) +
  geom_boxplot(data = bind_rows(
    modDistSub %>%
      mutate(group = "LGM <100km"),
    modDistSub2 %>%
      mutate(group = "Nearest to LGM")) %>%
      mutate(bin = cut(envdist, breaks = seq(0, 35, 1), include.lowest = TRUE, right = FALSE, labels = seq(0, 35, 1)[-1]-0.5*1)),
    aes(envdist, commdist, group = bin), outlier.colour = NA, fill = NA, colour = "black", size = 0.25) +
  scale_fill_viridis_c(option = "magma", breaks = c(0, 1), labels = c("min", "max")) +
  xlim(c(0, 32)) +
  labs(x = expression("Temperature difference between sites ["*degree*"C]"),
       y = "Bray-Curtis similarity",
       fill = NULL) +
  facet_wrap(.~group) +
  theme(legend.key.width = unit(1, "mm"),
        legend.key.height = unit(2.5, "mm")) +
  tag_facets(tag_pool = letters[1:2], position = "tr")


difModernDD <- bind_rows(
  full_join(modDistBin, modDistSubBin, by = c("envBin", "simBin")) %>%
    replace(is.na(.), 0) %>%
    mutate(diff = lognnorm.y - lognnorm.x,
           group = "All vs <100 km"),
  full_join(modDistBin, modDistSub2Bin, by = c("envBin", "simBin")) %>%
    replace(is.na(.), 0) %>%
    mutate(diff = lognnorm.y - lognnorm.x,
           group = "All vs nearest to LGM")
) %>%
  ggplot() +
  geom_tile(aes(simBin, envBin, fill = diff)) +
  geom_step(data = ecolim,
            aes(simBinPlot, Tlim),
            colour = "black",
            linewidth = 0.25
            ) +
  geom_step(data = bind_rows(ecolimSub %>%
                               mutate(group = "All vs <100 km"),
                             ecolimSub2 %>%
                               mutate(group = "All vs nearest to LGM")
                             ),
            aes(simBinPlot, Tlim),
            colour = "orange",
            linewidth = 0.25
            ) +
  scale_fill_distiller(palette = "RdGy", limits = c(-.8, .8), breaks = c(-0.8, 0.8), labels = c("less", "more")) +
  ylim(c(0, 32)) +
  coord_flip() +
  labs(y = expression("Temperature difference between sites ["*degree*"C]"),
       x = "Bray-Curtis similarity",
       fill = NULL) +
  facet_wrap(.~group) +
  theme(
    legend.key.width = unit(1, "mm"),
    legend.key.height = unit(2.5, "mm")) +
  tag_facets(tag_pool = letters[3:4], position = "tr")

modDDsens <- modernDD / difModernDD

ggsave("out/Ext_Fig7.eps", plot = modDDsens,  units = "mm", width = 121, height = 100)

# ===============================================================================================

# EDF 8 robustness of LGM DD pattern against chronological confidence ===========================

chronMap <- LGMSppDF_avg %>%
  select(Longitude, Latitude, ChronozoneLevel) %>%
  mutate(ChronozoneLevel = paste0("Chrono zone level ", ChronozoneLevel)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  ggplot() +
  geom_sf(data = bbox, fill = "grey80", size = 0.3, colour = "#464646") +
  geom_sf(data = grat, colour = "grey20", linewidth = 0.1, linetype = "dotted") +
  geom_sf(colour = "#7e4e90b2", size = 0.03) +
  geom_sf(data = wmap, fill = "white", colour = "black", size = 0.1) +
  geom_sf_text(data = ylabs, aes(label = label),
               nudge_x = -1e6,
               size = 1.5,
               parse = TRUE,
               hjust = 1) +
  geom_sf_text(data = xlabsbot, 
               aes(label = label),
               nudge_y = -5e5,
               size = 1.5,
               parse = TRUE,
               vjust = 1) +
  coord_sf(crs = "+proj=robin",
           clip = "off") +
  facet_wrap(~ChronozoneLevel, ncol = 1) +
  theme_opts


chronDD <- LGMDist_avgBin_chron %>%
  filter(model == "ensemble") %>%
  mutate(chronolevel = paste0("Chrono zone level ", chronolevel)) %>%
  ggplot(., aes(simBin, envBin)) +
  geom_tile(aes(fill = lognnorm)) +
  geom_step(data = ecolim, aes(simBinPlot, Tlim), colour = "grey", linewidth = 0.3) +
  coord_flip() +
  scale_fill_viridis_c(option = "magma", breaks = c(0, 1), labels = c("min", "max")) +
  facet_wrap(~chronolevel, ncol = 1) +
  labs(y = expression("Temperature difference between sites ["*degree*"C]"),
       x = "Bray-Curtis similarity",
       fill = "Number of pairs") +
  theme(legend.position = "none") +
  theme(strip.text = element_blank())


chronSensitivity <- chronDD + chronMap + plot_layout(widths = c(2, 3))

ggsave("out/Ext_Fig8.tiff", plot = chronSensitivity,  units = "mm", width = 121, height = 120)

# ===============================================================================================

# EDF 9  attribution ==============================================================================
Cols <- brewer.pal(4, "RdBu")
cols <- c("black", Cols[c(3, 2, 1, 4)])

Ext_Fig9 <- globalEX %>%
  mutate(group = factor(group, levels = c("Climatology", "Core top", "LGM"))) %>%
  ggplot() +
  geom_path(aes(EX, depth, colour = season)) +
  scale_color_manual(values = cols) +
  ylim(c(500, 0)) +
  ylab("Depth [m]") +
  xlab("Fraction variance explained") +
  facet_grid(vars(region), vars(group), labeller = label_wrap_gen(width = 10)) +
  theme(legend.position = "bottom")

ggsave("out/Ext_Fig9.eps", plot = Ext_Fig9,  units = "mm", width = 89, height = 150)

# ===============================================================================================

# EDF 10 analogue quality =======================================================================
mapTrecNoGrid <- LGM_noAF %>%
  filter(model == "ensemble") %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  ggplot() +
  geom_sf(data = bbox, fill = "grey80", size = 0.3, colour = "#464646") +
  geom_sf(data = grat, colour = "grey20", linewidth = 0.1, linetype = "dotted") +
  geom_sf(aes(colour =  T_anomaly_foram_atlas), size = 0.3) +
  geom_sf(data = wmap, fill = "white", colour = "black", size = 0.1) +
  geom_sf_text(data = ylabs, aes(label = label),
               nudge_x = -1e6,
               size = 1.5,
               parse = TRUE,
               hjust = 1) +
  geom_sf_text(data = xlabstop, aes(label = label),
               # nudge_y = 5e5,
               size = 1.5,
               parse = TRUE,
               vjust = -1) +
  coord_sf(crs = "+proj=robin",
           clip = "off") +
  scale_colour_gradient2(low = "dodgerblue1", high = "red2", limits = c(-14, 4), na.value = "grey90") +
  labs(colour = expression("Reconstructed\ntemperature anomaly ["*degree*"C]"),
       y = "No filtering") +
  theme_opts +
  theme(axis.title.y = element_text(face = "bold",
                                    margin = margin(r = 5)))


plotRecAll <- mapTrecNoGrid + theme(legend.position = "none") + latTrec + mapDifMod + theme(legend.position = "none") + latDifMod +
  plot_layout(nrow = 1, widths = unit(c(3, -1, 3, -1), c("null","null", "null", "null")))


# for global analogue filtering 
mapTrec_globalAF <- LGM_globalAF %>%
  filter(model == "ensemble") %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  ggplot() +
  geom_sf(data = bbox, fill = "grey80", size = 0.3, colour = "#464646") +
  geom_sf(data = grat, colour = "grey20", linewidth = 0.1, linetype = "dotted") +
  geom_sf(aes(colour =  T_anomaly_foram_atlas), size = 0.5) +
  geom_sf(data = wmap, fill = "white", colour = "black", size = 0.1) +
  geom_sf_text(data = ylabs, aes(label = label),
               nudge_x = -1e6,
               size = 1.5,
               parse = TRUE,
               hjust = 1) +
  geom_sf_text(data = xlabstop, aes(label = label),
               # nudge_y = 5e5,
               size = 1.5,
               parse = TRUE,
               vjust = -1) +
  coord_sf(crs = "+proj=robin",
           clip = "off") +
  scale_colour_gradient2(low = "dodgerblue1", high = "red2", limits = c(-14, 4)) + 
  labs(colour = expression("Reconstructed\ntemperature anomaly ["*degree*"C]"),
       y = "Global filtering") +
  theme_opts +
  theme(axis.title.y = element_text(face = "bold",
                                    margin = margin(r = 5)))


latTrec_globalAF <- LGM_globalAF %>%
  filter(model == "ensemble") %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform(crs = "+proj=robin") %>%
  mutate(y = purrr::map_dbl(geometry, 2)) %>%
  ggplot(aes(T_anomaly_foram_atlas, y)) +
  geom_line(data = latRangeData, aes(x, y), linetype = "dashed") +
  geom_smooth(data = lat_noAF,
              orientation = "y", colour = "red", fill = "pink", se = TRUE, method = "gam") +
  geom_smooth(orientation = "y", colour = "black", fill = "grey60", method = "gam") +
  labs(x = expression("Temperature anomaly ["*degree*"C]")) +
  labs(x = NULL) +
  coord_cartesian(xlim = xlimRec, ylim = latRangeMap$y) +
  theme(axis.line.x = element_line(colour = "#464646", linewidth = 0.3),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank())



mapDifMod_globalAF <- LGM_globalAF %>%
  filter(model == "ensemble") %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  ggplot() +
  geom_sf(data = bbox, fill = "grey80", size = 0.3, colour = "#464646") +
  geom_sf(data = grat, colour = "grey20", linewidth = 0.1, linetype = "dotted") +
  geom_sf(aes(colour =  difDataModel), size = 0.5) +
  geom_sf(data = wmap, fill = "white", colour = "black", size = 0.1) +
  geom_sf_text(data = ylabs, aes(label = label),
               nudge_x = -1e6,
               size = 1.5,
               parse = TRUE,
               hjust = 1) +
  geom_sf_text(data = xlabstop, aes(label = label),
               # nudge_y = 5e5,
               size = 1.5,
               parse = TRUE,
               vjust = -1) +
  coord_sf(crs = "+proj=robin",
           clip = "off") +
  scale_colour_gradient2(low = "dodgerblue1", high = "red2", limits = c(-7, 11)) + 
  labs(colour = expression("Simulated - reconstructed\ntemperature anomaly ["*degree*"C]"),
       y = "Global filtering") +
  theme_opts


latDifMod_globalAF <- LGM_globalAF %>%
  filter(model == "ensemble") %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform(crs = "+proj=robin") %>%
  mutate(y = purrr::map_dbl(geometry, 2)) %>%
  ggplot(aes(difDataModel, y)) +
  geom_line(data = latRangeData, aes(x, y), linetype = "dashed") +
  geom_smooth(data = lat_noAF,
              orientation = "y", colour = "red", fill = "pink", se = TRUE, method = "gam") +
  geom_smooth(orientation = "y", colour = "black", fill = "grey60", method = "gam") +
  labs(x = expression("Temperature anomaly ["*degree*"C]")) +
  labs(x = NULL) +
  coord_cartesian(xlim = xlimDif, ylim = latRangeMap$y) +
  theme(axis.line.x = element_line(colour = "#464646", linewidth = 0.3),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank())


plotRec_globalAF <- mapTrec_globalAF + theme(legend.position = "none") + latTrec_globalAF +
  mapDifMod_globalAF + theme(legend.position = "none") + latDifMod_globalAF +
  plot_layout(nrow = 1, widths = unit(c(3, -1, 3, -1), c("null","null", "null", "null")))

# for regional analogue filtering
mapTrec_regionalAF <- LGM_regionalAF %>%
  filter(model == "ensemble") %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  ggplot() +
  geom_sf(data = bbox, fill = "grey80", size = 0.3, colour = "#464646") +
  geom_sf(data = grat, colour = "grey20", linewidth = 0.1, linetype = "dotted") +
  geom_sf(aes(colour =  T_anomaly_foram_atlas), size = 0.5) +
  geom_sf(data = wmap, fill = "white", colour = "black", size = 0.1) +
  geom_sf_text(data = ylabs, aes(label = label),
               nudge_x = -1e6,
               size = 1.5,
               parse = TRUE,
               hjust = 1) +
  geom_sf_text(data = xlabstop, aes(label = label),
               # nudge_y = 5e5,
               size = 1.5,
               parse = TRUE,
               vjust = -1) +
  coord_sf(crs = "+proj=robin",
           clip = "off") +
  scale_colour_gradient2(low = "dodgerblue1", high = "red2", limits = c(-14, 4)) + 
  labs(colour = expression("Reconstructed\ntemperature anomaly ["*degree*"C]"),
       y = "Regional filtering") +
  theme_opts +
  theme(axis.title.y = element_text(face = "bold",
                                    margin = margin(r = 5)))


latTrec_regionalAF <- LGM_regionalAF %>%
  filter(model == "ensemble") %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform(crs = "+proj=robin") %>%
  mutate(y = purrr::map_dbl(geometry, 2)) %>%
  ggplot(aes(T_anomaly_foram_atlas, y)) +
  geom_line(data = latRangeData, aes(x, y), linetype = "dashed") +
  geom_smooth(data = lat_noAF,
              orientation = "y", colour = "red", fill = "pink", se = TRUE, method = "gam") +
  geom_smooth(orientation = "y", colour = "black", fill = "grey60", method = "gam") +
  labs(x = expression("Temperature anomaly ["*degree*"C]")) +
  labs(x = NULL) +
  coord_cartesian(xlim = xlimRec, ylim = latRangeMap$y) +
  theme(axis.line.x = element_line(colour = "#464646", linewidth = 0.3),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank())

mapDifMod_regionalAF <- LGM_regionalAF %>%
  filter(model == "ensemble") %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  ggplot() +
  geom_sf(data = bbox, fill = "grey80", size = 0.3, colour = "#464646") +
  geom_sf(data = grat, colour = "grey20", linewidth = 0.1, linetype = "dotted") +
  geom_sf(aes(colour =  difDataModel), size = 0.5) +
  geom_sf(data = wmap, fill = "white", colour = "black", size = 0.1) +
  geom_sf_text(data = ylabs, aes(label = label),
               nudge_x = -1e6,
               size = 1.5,
               parse = TRUE,
               hjust = 1) +
  geom_sf_text(data = xlabstop, aes(label = label),
               # nudge_y = 5e5,
               size = 1.5,
               parse = TRUE,
               vjust = -1) +
  coord_sf(crs = "+proj=robin",
           clip = "off") +
  scale_colour_gradient2(low = "dodgerblue1", high = "red2", limits = c(-7, 11)) + 
  labs(colour = expression("Simulated - reconstructed\ntemperature anomaly ["*degree*"C]"),
       y = "Regional filtering") +
  theme_opts

latDifMod_regionalAF <- LGM_regionalAF %>%
  filter(model == "ensemble") %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform(crs = "+proj=robin") %>%
  mutate(y = purrr::map_dbl(geometry, 2)) %>%
  ggplot(aes(difDataModel, y)) +
  geom_line(data = latRangeData, aes(x, y), linetype = "dashed") +
  geom_smooth(data = lat_noAF,
              orientation = "y", colour = "red", fill = "pink", se = TRUE, method = "gam") +
  geom_smooth(orientation = "y", colour = "black", fill = "grey60", method = "gam") +
  labs(x = expression("Temperature anomaly ["*degree*"C]")) +
  labs(x = NULL) +
  coord_cartesian(xlim = xlimDif, ylim = latRangeMap$y) +
  theme(axis.line.x = element_line(colour = "#464646", linewidth = 0.3),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank())


plotRec_regionalAF <- mapTrec_regionalAF + latTrec_regionalAF +
  mapDifMod_regionalAF + latDifMod_regionalAF +
  plot_layout(nrow = 1, widths = unit(c(3, -1, 3, -1), c("null","null", "null", "null")))

# plot of affect of analogue filtering
ExtF10 <- plotRecAll / plotRec_globalAF / plotRec_regionalAF

ggsave("out/Ext_Fig10.eps", plot = ExtF10,  units = "mm", width = 183, height = 120)

# ===============================================================================================
