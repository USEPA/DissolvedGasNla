## This script generates data to be referenced in inline-code in the manuscript.
## This makes the inline code more succint and the .Rmd easier to read.


# SETUP----------------------
## load libraries
library(sf) # spatial data
library(tidyverse) # dplyr, ggplot
library(janitor) # clean names
library(USAboundaries) # state boundaries
library(tictoc) # processing time

## Identify local path for each user
localPath <- Sys.getenv("USERPROFILE")


# DATA---------------
## load sample data (dg object)----
load(paste0(localPath,  # object name dg
            "/Environmental Protection Agency (EPA)/",
            "ORD NLA17 Dissolved Gas - Documents/",
            "inputData/dg.2021-02-01.RData"))

# To enable spatial analysis of the data, the dataframe will be converted to a 
# 'simple features' (sf) object.
# Define coordinates
coords <- data.frame(longitude = dg$map.lon.dd, latitude = dg$map.lat.dd)

dg.sf <- st_as_sf(dg, coords = c("map.lon.dd", "map.lat.dd"), 
                  crs = 4269) %>% # standard for lat/lon
  st_transform(5070) # project to CONUS ALBERS for plotting


## load ecoregions----
# read in ecoregion polygons
ecoR <- st_read(dsn = paste0(localPath, 
                             "/Environmental Protection Agency (EPA)/",
                             "ORD NLA17 Dissolved Gas - Documents/inputData"),
                layer = "aggr_ecoregions_simple",
                quiet = TRUE)

# Check CRS
st_crs(ecoR) # 3857
ecoR <- st_transform(ecoR, 5070) # convert to CONUS Albers
st_crs(ecoR) # 5070

# SET UP CUSTOM COLORS FOR ECOREGIONS
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


## Load state boundaries----
states <- USAboundaries::us_states() %>%
  dplyr::filter(!state_name %in% c("Alaska", "District of Columbia", "Hawaii", "Puerto Rico")) %>%
  st_transform(5070) # convert to CONUS Albers


## load population data----
# this object is created in dgIndicatorAnalysis.Rmd.  It is Roy's
# all_predictions object with emission rates and total flux added.
#tic()
load(paste0(localPath, "\\Environmental Protection Agency (EPA)\\",
            "ORD NLA17 Dissolved Gas - Documents\\",
            "inputData\\all_predictions.rda"))
#toc()

# add full ecoregion names to all_predictions
wsa9names <- structure(list(WSA9 = c("CPL", "NAP", "NPL", "SAP", "SPL", 
                                     "TPL","UMW", "WMT", "XER"), 
                            WSA9_NAME = c("Coastal Plains", "Northern Appalachians", 
                                          "Northern Plains", "Southern Appalachians", 
                                          "Southern Plains", "Temperate Plains", 
                                          "Upper Midwest", "Western Mountains", "Xeric")), 
                       row.names = c(NA, 9L), class = "data.frame")

all_predictions <- full_join(all_predictions, wsa9names)


## emission rates----
# code copied from dgIndicatorAnalysis.Rmd
### read wind data----
# Mean index period wind speed was calculated on VM via 
# nlaSampleFrameWind.R.  The resulting data object was written out and read in below.
wind <- read.table(paste0(Sys.getenv("USERPROFILE"), 
                          "/Environmental Protection Agency (EPA)/",
                          "ORD NLA17 Dissolved Gas - Documents/",
                          "inputData/windSampleFrame.txt"))
# print first 6 lines
wind %>% head()

# print first 6 .row values of all_predictions
all_predictions %>%
  filter(.row %in% 1:6) %>% # each row is a unique waterbody
  group_by(.row) %>%
  slice(1) # grabs first record

# Row order in wind and all_predictions is inherited from sample frame.  Thus row 1
# of wind corresponds to .row == 1 in all_predictions.  This is a sketchy way to merge
# things, but will do for now.

# add .id to wind which is simply row number
wind <- wind %>% mutate(.id = row_number())

# Merge wind and all_predictions
#tic() # 21 seconds on VM
all_predictions <- full_join(all_predictions, wind, by = c(".row" = ".id")) # 40 seconds on laptop
#toc()
dim(wind) # 465897
dim(all_predictions) #232948500, 500 observations for each lake
465897*500 #232948500



### estimate gas exchange rate for population----

# Using gas exchange model from:
# Vachon, D., and Y. T. Prairie (2013), The ecosystem size and shape 
# dependence of gas transfer velocity versus wind speed relationships in 
# lakes, Can. J. Fish. Aquat. Sci., 70(12), 1757-1764, doi:10.1139/cjfas-2013-0241.
# Table 2, model B.

# define k600 prediction function
pred.k600 <- function(ws, area) { 
  # predict k600 in cm h-1
  # ws = wind speed (m s-1)
  # area = lake area (km2)
  2.51 + (1.48 * ws) + (0.39 * ws * log10(area))
}

### calculate emission and flux----
all_predictions <- all_predictions %>%
  mutate(k600.cm.h = pred.k600(ws = ws, area = area_ha/100), # 100ha = 1km2
         # 1000 L to m3. 44ng to 1nmol. ng->ug. ug->mg
         e.n2o.mg.m2.d = (((n2o - n2oeq) * 1000 * 44) / (1000 * 1000)) * 
           (k600.cm.h * (24/100)), #h->day. cm->m
         f.n2o.Mg.y = (e.n2o.mg.m2.d * 365 * (area_ha*10000)) / # day to year. ha to m2.
           (1000*1000*1000)) #mg-g, g-kg, kg-Mg


## ecoregion and national means----
ecoreg_stats <- all_predictions %>%
  group_by(WSA9, WSA9_NAME, .draw) %>%
  summarise( mean_n2o = mean(n2o),
             mean_n2oeq = mean(n2oeq),
             mean_sat = mean(n2osat),
             prop_sat = sum(n2osat < 1) / length(unique(.row)),
             mean_en2o = mean(e.n2o.mg.m2.d),
             fn2o_Mg = sum(f.n2o.Mg.y)) %>%
  summarise( post_med_n2o = round(mean(mean_n2o), 1),
             LCI_n2o = round(quantile(mean_n2o, probs = 0.025), 1),
             UCI_n2o = round(quantile(mean_n2o, probs = 0.975), 1),
             post_med_n2oeq = round(mean(mean_n2oeq), 2),
             LCI_n2oeq = round(quantile(mean_n2oeq, probs = 0.025), 2),
             UCI_n2oeq = round(quantile(mean_n2oeq, probs = 0.975), 2),
             post_med_sat = round(mean(mean_sat), 2),
             LCI_sat = round(quantile(mean_sat, probs = 0.025), 2),
             UCI_sat = round(quantile(mean_sat, probs = 0.975), 2),
             estimate = round(mean(prop_sat), 3),
             LCL = round(quantile(prop_sat, probs = 0.025), 3),
             UCL = round(quantile(prop_sat, probs = 0.975), 3),
             post_med_en2o = round(mean(mean_en2o), 3),
             LCL_en2o = round(quantile(mean_en2o, probs = 0.025), 3),
             UCL_en2o = round(quantile(mean_en2o, probs = 0.975), 3),
             post_med_fn2o = round(mean(fn2o_Mg), 3), 
             LCL_fn2o = round(quantile(fn2o_Mg, probs = 0.025), 3),
             UCL_fn2o = round(quantile(fn2o_Mg, probs = 0.975), 3),
             n = dplyr::n())

national_stats <- all_predictions %>%
  group_by(.draw) %>%
  summarise( mean_n2o = mean(n2o),
             mean_n2oeq = mean(n2oeq),
             mean_sat = mean(n2osat),
             prop_sat = sum(n2osat < 1) / length(unique(.row)),
             mean_en2o = mean(e.n2o.mg.m2.d),
             fn2o_Mg = sum(f.n2o.Mg.y)) %>%
  summarise( post_med_n2o = round(mean(mean_n2o), 1),
             LCI_n2o = round(quantile(mean_n2o, probs = 0.025), 1),
             UCI_n2o = round(quantile(mean_n2o, probs = 0.975), 1),
             post_med_n2oeq = round(mean(mean_n2oeq), 2),
             LCI_n2oeq = round(quantile(mean_n2oeq, probs = 0.025), 2),
             UCI_n2oeq = round(quantile(mean_n2oeq, probs = 0.975), 2),
             post_med_sat = round(mean(mean_sat), 2),
             LCI_sat = round(quantile(mean_sat, probs = 0.025), 2),
             UCI_sat = round(quantile(mean_sat, probs = 0.975), 2),
             estimate = round(mean(prop_sat), 3),
             LCL = round(quantile(prop_sat, probs = 0.025), 3),
             UCL = round(quantile(prop_sat, probs = 0.975), 3),
             post_med_en2o = round(mean(mean_en2o), 3),
             LCL_en2o = round(quantile(mean_en2o, probs = 0.025), 3),
             UCL_en2o = round(quantile(mean_en2o, probs = 0.975), 3),
             post_med_fn2o = round(mean(fn2o_Mg), 3), 
             LCL_fn2o = round(quantile(fn2o_Mg, probs = 0.025), 3),
             UCL_fn2o = round(quantile(fn2o_Mg, probs = 0.975), 3),
             n = dplyr::n()) %>%
  mutate(WSA9 = "national", WSA9_NAME = "national")


# FIGURES-----------------

# The images are written to manuscript/manuscript_figures/...
# the first time these scripts are run.  Unless the code has changed,
# there is no need to run the scripts again.  Commenting out the script
# will expedite knitting of manuscript_file.Rmd

## Figure 1: site map----
# if the image is already on computer, then nothing, else create image
if(!("figure1.png" %in% list.files("manuscript/manuscript_figures"))) {
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
}


## Figure 2: N2O saturation ratio density plot----
# code taken from dgIndicatorAnalysis.Rmd, "N2O saturation ratio: continuous variable
# if the image is already on computer, then nothing, else create image
if(!("figure2.tiff" %in% list.files("manuscript/manuscript_figures"))) {
  
  # density plots of predicted distribution, mean, and median
  # using xlim to zoom in on values close to 1
  
  dummy <- all_predictions %>%
    group_by(WSA9) %>%
    summarize(mean = mean(n2osat),
              median = median(n2osat)) %>%
    mutate(dens.med = c(2, 2.5, 2.3, 2.05, 1.75, 2.15, 2.2, 1.9, 1.9),
           dens.mean = c(.4, 1.9, 2, .6, .65, .45, 1.5, 1.45, 0.9),
           .draw = 1) # needed to match grouping aesthetic in ggplot call
  
  # tic() # 40 seconds for 100 draws
  all_predictions %>%
    filter(.draw %in% 1:100) %>% # subset for practice
    #mutate(ecoregion = factor(WSA9)) %>%
    #mutate(ecoregion = fct_reorder(ecoregion, preds_sat)) %>%
    ggplot(aes(x = n2osat, group = .draw, color = .draw ) ) +
    geom_density(show.legend = FALSE) +
    geom_vline(xintercept = 1, color = "red") +
    geom_segment(data = dummy,
                 aes(x = mean, xend = mean, y = 0,
                     yend = dens.mean),
                 linetype = "solid") +
    geom_segment(data = dummy,
                 aes(x = median, xend = median, y = 0,
                     yend = dens.med),
                 linetype = "dashed") +
    ylab ("density") + 
    xlab(expression(N[2]*O~saturation~ratio)) +
    xlim(0,3) + # data range is 0-600, but sample max is 30
    facet_wrap(~WSA9) +
    theme_bw() +
    theme(legend.position = "none")
  # toc()
  ggsave("manuscript/manuscript_figures/figure2.tiff", width = 8.5, height = 5)
  
}

## Figure 3: delta N2O by waterbody size----
# See 'N2O saturation ratio: continuous variable' section of dgIndicatorAnalysis.Rmd
# if the image is already on computer, then nothing, else create image
if(!("figure3.tiff" %in% list.files("manuscript/manuscript_figures"))) {
  


# PLOT N2O SAT RATIO BY WSA9 AND SIZE
# point and linerange
p1.data <- all_predictions %>%
  group_by(WSA9, size_cat, .draw) %>% # group by iteration
  summarise(mean_sat = mean(n2osat)) %>% # 500 means for each WSA9
  # now summarize to 1 statistic per WSA9
  summarise( estimate = round(median(mean_sat), 3), 
             LCL = round(quantile(mean_sat, probs = 0.025), 3),
             UCL = round(quantile(mean_sat, probs = 0.975), 3)) %>% 
  mutate(ecoregion = factor(WSA9)) %>%
  mutate(ecoregion = fct_reorder(ecoregion, estimate)) 

# data fpr arrow segment.  Only smallest size category
p1.data.arrow <- p1.data %>%
  filter(size_cat == "min_4")

p1 <- p1.data %>%
  ggplot(aes(x=estimate, y=size_cat)) +
  geom_point() +
  geom_linerange( aes( xmin = LCL, xmax = UCL)) +
  geom_vline(xintercept = 1, color='blue') +
  geom_segment(data = p1.data.arrow, 
               aes(x=1, y = "min_4", xend = estimate, yend = "min_4"),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "red") +
  xlab(expression(mean~N[2]*O~saturation~ratio)) +
  ylab("waterbody size category (ha)") +
  theme_bw() +
  facet_wrap(~ecoregion)


# PLOT N2O* BY SIZE
# point and linerange
p2 <- all_predictions %>%
  mutate(n2ostar = abs(n2o - n2oeq)) %>%
  group_by(size_cat, .draw) %>% # group by iteration
  #group_by(.draw) %>%
  summarise(mean_star = mean(n2ostar)) %>% # 500 means for each WSA9
  # now summarize to 1 statistic per WSA9
  summarise( estimate = round(median(mean_star), 3), 
             LCL = round(quantile(mean_star, probs = 0.025), 3),
             UCL = round(quantile(mean_star, probs = 0.975), 3)) %>% 
  ggplot( aes( x = estimate, y = size_cat ) ) +
  geom_point() +
  geom_linerange( aes( xmin = LCL, xmax = UCL)) +
  xlab(expression("["*Delta~N[2]*O*"]"~"("*nmol~L^{-1}*")")) +
  #coord_flip() + 
  theme_bw() +
  theme(axis.title.y = element_blank())
#ggsave("output/figures/n2oStarBySize.tiff")

ggpubr::ggarrange(p1, p2, ncol=2, nrow = 1, widths = c(0.7, 0.3), labels = c("A", "B"))

ggsave("manuscript/manuscript_figures/figure3.tiff", width = 8.5, height = 5)
}


# MANUSCRIPT DATA-----

## paragraph 1----
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
