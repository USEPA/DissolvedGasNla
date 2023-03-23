## This script generates data to be referenced in inline-code in the manuscript.
## This makes the inline code more succint and the .Rmd easier to read.


# SETUP----------------------
# load libraries-----------------
library(sf) # spatial data
library(tidyverse) # dplyr, ggplot
library(janitor) # clean names
library(USAboundaries) # state boundaries
library(tictoc) # processing time

# Identify local path for each user
localPath <- Sys.getenv("USERPROFILE")


# DATA---------------
#### load sample data (dg object)
load(paste0(localPath,  # object name dg
            "/Environmental Protection Agency (EPA)/",
            "ORD NLA17 Dissolved Gas - Documents/",
            "inputData/dg.2021-02-01.RData"))

# To enable spatial analysis of the data, the dataframe will be converted to a 'simple features' (sf) object.
# Define coordinates
coords <- data.frame(longitude = dg$map.lon.dd, latitude = dg$map.lat.dd)

dg.sf <- st_as_sf(dg, coords = c("map.lon.dd", "map.lat.dd"), 
                  crs = 4269) %>% # standard for lat/lon
  st_transform(5070) # project to CONUS ALBERS for plotting


##### load ecoregions
# read in ecoregion polygons
ecoR <- st_read(dsn = paste0(localPath, 
                             "/Environmental Protection Agency (EPA)/",
                             "ORD NLA17 Dissolved Gas - Documents/inputData"),
                layer = "aggr_ecoregions_simple")

# Check CRS
st_crs(ecoR) # 3857
ecoR <- st_transform(ecoR, 5070) # convert to CONUS Albers
st_crs(ecoR) # 5070

# SET UP CUSTOM COLORS FOR ECOREGIONS---------
# Custom color pallette for ecoregion polygons. Attempted to mirror
# https://www.epa.gov/national-aquatic-resource-surveys/
# ecoregional-results-national-lakes-assessment-2012
cols <- c("Coastal Plains" = "orange1",
          "Northern Appalachians" = "lightpink1",
          "Northern Plains" = "darksalmon",
          "Southern Appalachians" = "mediumturquoise",
          "Southern Plains" = "khaki4",
          "Temperate Plains" = "forestgreen", 
          "Upper Midwest" = "deepskyblue4",
          "Western Mountains" = "saddlebrown",
          "Xeric" = "lightskyblue4")


#### Load state boundaries
states <- USAboundaries::us_states() %>%
  dplyr::filter(!state_name %in% c("Alaska", "District of Columbia", "Hawaii", "Puerto Rico")) %>%
  st_transform(5070) # convert to CONUS Albers


### load population data
tic()
load(paste0(localPath, "\\Environmental Protection Agency (EPA)\\ORD NLA17 Dissolved Gas - Documents\\inputData\\all_predictions.rda"))
toc()


# FIGURES-----------------


# # Figure 1
# # Commented to expedite script
 ggplot() +
   geom_sf(data = ecoR, color = NA, aes(fill = WSA9_NAME)) +
   geom_sf(data = dg.sf %>% filter(!is.na(n2o.src.snk)), 
           aes(size = dissolved.n2o.nmol, color = n2o.src.snk),
           show.legend = "point") +
   geom_sf(data = states, fill = NA, color = "cornsilk3", size = 0.1) +
   # guide argument below removes points from the boxes in the ecoregion legend
   # https://aosmith.rbind.io/2020/07/09/ggplot2-override-aes/
   scale_fill_manual("Ecoregion", values = cols, 
                     guide = guide_legend(override.aes = list(shape = NA))) +
   scale_color_manual(values = c("white", "black"), name = "source/sink") +
   scale_size(name = expression(N[2]*O~(nM)),
              range = c(0.1, 10), # custom size range
              breaks = c(1, 10, 25, 50, 100)) + # custom breaks
   theme(legend.key.size = unit(0.4, "cm"), # size of boxes in legend
         legend.title = element_text(size = 8))
 ggsave("manuscript/manuscript_figures/figure1.png", width = 8, height = 4, units = "in")



# MANUSCRIPT DATA-----

# paragraph 1
undersatN2oPercent <- round(((dg %>% dplyr::filter(n2o.src.snk == "sink") %>% 
                         {nrow(.)} / dg %>% distinct(site.id) %>% 
                         {nrow(.)}) * 100), 1)

propN2oSinkMean <- all_predictions %>% 
  group_by(.draw) %>% 
  summarise(prop_sat = sum(n2osat < 1) / length(.row)) %>% 
  summarise(estimate = round(median(prop_sat), 3) * 100) %>% 
  pull()

propN2oSinkMean <- all_predictions %>% 
  group_by(.draw) %>% 
  summarise(prop_sat = sum(n2osat < 1) / length(.row)) %>% 
  summarise(estimate = round(median(prop_sat), 3) * 100) %>% 
  pull()

propN2oSinkLower <- all_predictions %>% 
  group_by(.draw) %>% 
  summarise(prop_sat = sum(n2osat < 1) / length(.row)) %>% 
  summarise(LCL = round(quantile(prop_sat, probs = 0.025), 3) * 100) %>% 
  pull()

propN2oSinkUpper <- all_predictions %>% 
  group_by(.draw) %>% 
  summarise(prop_sat = sum(n2osat < 1) / length(.row)) %>% 
  summarise(UCL = round(quantile(prop_sat, probs = 0.975), 3) * 100) %>% 
  pull()

propN2oSinkMeanWM <- all_predictions %>% 
  group_by(WSA9, .draw) %>% 
  summarise(prop_sat = sum(n2osat < 1) / length(.row)) %>% 
  summarise(estimate = round(median(prop_sat), 3) * 100) %>% 
  filter(WSA9 == "WMT") %>% 
  select(estimate) %>% 
  pull()


propN2oSinkLowerWM <- all_predictions %>% 
  group_by(WSA9, .draw) %>% 
  summarise(prop_sat = sum(n2osat < 1) / length(.row)) %>% 
  summarise(LCL = round(quantile(prop_sat, probs = 0.025), 3) * 100) %>% 
  filter(WSA9 == "WMT") %>% 
  pull()


propN2oSinkUpperWM <- all_predictions %>% 
  group_by(WSA9, .draw) %>% 
  summarise(prop_sat = sum(n2osat < 1) / length(.row)) %>% 
  summarise(UCL = round(quantile(prop_sat, probs = 0.975), 3) * 100) %>% 
  filter(WSA9 == "WMT") %>% 
  pull()

propN2oSinkMeanNPL <- all_predictions %>% 
  group_by(WSA9, .draw) %>% 
  summarise(prop_sat = sum(n2osat < 1) / length(.row)) %>% 
  summarise(estimate = round(median(prop_sat), 3) * 100) %>% 
  filter(WSA9 == "NPL") %>% 
  select(estimate) %>% 
  pull()


propN2oSinkLowerNPL <- all_predictions %>% 
  group_by(WSA9, .draw) %>% 
  summarise(prop_sat = sum(n2osat < 1) / length(.row)) %>% 
  summarise(LCL = round(quantile(prop_sat, probs = 0.025), 3) * 100) %>% 
  filter(WSA9 == "NPL") %>% 
  pull()


propN2oSinkUpperNPL <- all_predictions %>% 
  group_by(WSA9, .draw) %>% 
  summarise(prop_sat = sum(n2osat < 1) / length(.row)) %>% 
  summarise(UCL = round(quantile(prop_sat, probs = 0.975), 3) * 100) %>% 
  filter(WSA9 == "NPL") %>% 
  pull()
