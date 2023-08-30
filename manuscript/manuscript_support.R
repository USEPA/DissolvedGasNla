## This script generates data to be referenced in inline-code in the manuscript.
## This makes the inline code more succint and the .Rmd easier to read.


# SETUP----------------------
## load libraries
library(sf) # spatial data
library(tidyverse) # dplyr, ggplot
library(janitor) # clean names
library(USAboundaries) # state boundaries
library(tictoc) # processing time
library(ggpubr) # multiple plots
library(ggallin) # psuedlolog transformation for negative values

## Identify local path for each user
localPath <- Sys.getenv("USERPROFILE")


# DATA---------------
## load sample data (dg object)----
if(localPath == "") { # if DMAP, then
  load("inputData/dg.2021-02-01.RData")
} else { # if not DMAP, then...
load(paste0(localPath,  # object name dg
            "/Environmental Protection Agency (EPA)/",
            "ORD NLA17 Dissolved Gas - Documents/",
            "inputData/dg.2021-02-01.RData"))
}

# To enable spatial analysis of the data, the dataframe will be converted to a 
# 'simple features' (sf) object.
# Define coordinates
coords <- data.frame(longitude = dg$map.lon.dd, latitude = dg$map.lat.dd)

dg.sf <- st_as_sf(dg, coords = c("map.lon.dd", "map.lat.dd"), 
                  crs = 4269) %>% # standard for lat/lon
  st_transform(5070) # project to CONUS ALBERS for plotting


## load ecoregions----
# read in ecoregion polygons
if(localPath == ""){ # if DMAP, then
  ecoR <- st_read(dsn = "inputData",
                  layer = "aggr_ecoregions_simple",
                  quiet = TRUE)
  } else { # if not DMAP, then...
ecoR <- st_read(dsn = paste0(localPath, 
                             "/Environmental Protection Agency (EPA)/",
                             "ORD NLA17 Dissolved Gas - Documents/inputData"),
                layer = "aggr_ecoregions_simple",
                quiet = TRUE)
}
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
#tic() # 7.6min on DMAP
if(localPath == ""){ # if DMAP, then
  load("../../shared/jbeaulie/all_predictions.rda")
} else { # if not DMAP
load(paste0(localPath, "\\Environmental Protection Agency (EPA)\\",
            "ORD NLA17 Dissolved Gas - Documents\\",
            "inputData\\all_predictions.rda"))
}
#toc()

# remove Great Salt Lake
dim(all_predictions) #  931794000        13
all_predictions <- all_predictions %>%
  filter(.row != 258807) # this row id is for Great Salt lake, 431066 ha
dim(all_predictions) # 931792000 13, good, removed 2000 instances of salt lake

format(object.size(all_predictions), "Gb") # 72.9 GB

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
if(localPath == ""){ # if DMAP, then
  wind <- read.table("inputData/windSampleFrame.txt")
} else { # if not DMAP, then...
wind <- read.table(paste0(Sys.getenv("USERPROFILE"), 
                          "/Environmental Protection Agency (EPA)/",
                          "ORD NLA17 Dissolved Gas - Documents/",
                          "inputData/windSampleFrame.txt"))
}

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
# left join needed to exclude great salt lake in wind data.  see line 88-93
all_predictions <- left_join(all_predictions, wind, by = c(".row" = ".id")) # 40 seconds on laptop
#toc()
dim(wind) # 465897
dim(all_predictions) #931792000, 2000 observations for each lake
465896*2000 #931792000

summary(wind$ws)

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
         ScN2o = 2055.6 - 137.11*surftemp + 4.3173*surftemp^2 - 0.05435*surftemp^3, # schmidt number (Wanninkhof et al 1992)
         #k600 to kn2o
         kn2o.cm.h = case_when(surftemp <= 30 ~ k600.cm.h * (ScN2o/600)^-(2/3), # conversion only good up to 30C. see paper!
                               TRUE ~ k600.cm.h),
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
             median_en2o = median(e.n2o.mg.m2.d),
             fn2o_Mg = sum(f.n2o.Mg.y),
             mean_no3 = mean(no3_cat),
             total_sa = sum(area_ha),
             prop_small = sum(size_cat == "min_4") / length(unique(.row))) %>%
  summarise( post_med_n2o = round(median(mean_n2o), 1),
             LCI_n2o = round(quantile(mean_n2o, probs = 0.025), 1),
             UCI_n2o = round(quantile(mean_n2o, probs = 0.975), 1),
             post_med_n2oeq = round(median(mean_n2oeq), 2),
             LCI_n2oeq = round(quantile(mean_n2oeq, probs = 0.025), 2),
             UCI_n2oeq = round(quantile(mean_n2oeq, probs = 0.975), 2),
             post_med_sat = round(median(mean_sat), 2),
             LCI_sat = round(quantile(mean_sat, probs = 0.025), 2),
             UCI_sat = round(quantile(mean_sat, probs = 0.975), 2),
             estimate = round(median(prop_sat), 3),
             LCL = round(quantile(prop_sat, probs = 0.025), 3),
             UCL = round(quantile(prop_sat, probs = 0.975), 3),
             post_med_en2o = round(median(mean_en2o), 3),
             post_median_en2o = round(median(median_en2o), 3),
             LCL_en2o = round(quantile(mean_en2o, probs = 0.025), 3),
             UCL_en2o = round(quantile(mean_en2o, probs = 0.975), 3),
             post_med_fn2o = round(median(fn2o_Mg), 3), 
             LCL_fn2o = round(quantile(fn2o_Mg, probs = 0.025), 3),
             UCL_fn2o = round(quantile(fn2o_Mg, probs = 0.975), 3),
             post_med_no3 = round(median(mean_no3), 4),
             LCL_no3 = round(quantile(mean_no3, probs = 0.025), 4),
             UCL_no3 = round(quantile(mean_no3, probs = 0.975), 4),
             post_prop_small = round(median(prop_small), 1),
             LCL_prop_small = round(quantile(prop_small, probs = 0.025), 1),
             UCL_prop_small = round(quantile(prop_small, probs = 0.975), 1),
             total_sa = round(mean(total_sa), digits = 0),
             n = dplyr::n()) %>%
  ungroup()

national_stats <- all_predictions %>%
  group_by(.draw) %>%
  summarise( mean_n2o = mean(n2o),
             mean_n2oeq = mean(n2oeq),
             mean_sat = mean(n2osat),
             prop_sat = sum(n2osat < 1) / length(unique(.row)),
             mean_en2o = mean(e.n2o.mg.m2.d),
             median_en2o = median(e.n2o.mg.m2.d),
             fn2o_Mg = sum(f.n2o.Mg.y),
             total_sa = sum(area_ha)) %>%
  summarise( post_med_n2o = round(median(mean_n2o), 1),
             LCI_n2o = round(quantile(mean_n2o, probs = 0.025), 1),
             UCI_n2o = round(quantile(mean_n2o, probs = 0.975), 1),
             post_med_n2oeq = round(median(mean_n2oeq), 2),
             LCI_n2oeq = round(quantile(mean_n2oeq, probs = 0.025), 2),
             UCI_n2oeq = round(quantile(mean_n2oeq, probs = 0.975), 2),
             post_med_sat = round(median(mean_sat), 2),
             LCI_sat = round(quantile(mean_sat, probs = 0.025), 2),
             UCI_sat = round(quantile(mean_sat, probs = 0.975), 2),
             estimate = round(mean(prop_sat), 3),
             LCL = round(quantile(prop_sat, probs = 0.025), 3),
             UCL = round(quantile(prop_sat, probs = 0.975), 3),
             post_median_en2o = round(median(median_en2o), 3),
             post_med_en2o = round(median(mean_en2o), 3),
             LCL_en2o = round(quantile(mean_en2o, probs = 0.025), 3),
             UCL_en2o = round(quantile(mean_en2o, probs = 0.975), 3),
             post_med_fn2o = round(median(fn2o_Mg), 3), 
             LCL_fn2o = round(quantile(fn2o_Mg, probs = 0.025), 3),
             UCL_fn2o = round(quantile(fn2o_Mg, probs = 0.975), 3),
             total_sa = round(mean(total_sa), digits = 0),
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

## Figure X: delta N2O by waterbody size----
# See 'N2O saturation ratio: continuous variable' section of dgIndicatorAnalysis.Rmd
# if the image is already on computer, then nothing, else create image
if(!("deltaN2ObySize.tiff" %in% list.files("manuscript/manuscript_figures"))) {
  

  
# PLOT N2O SAT RATIO BY WSA9 AND SIZE
# point and linerange
p1.data <- all_predictions %>%
  filter(WSA9 %in% c("CPL", "NPL")) %>% # only show two examples.  comment out to include all ecoregions
  group_by(WSA9_NAME, size_cat, .draw) %>% # group by iteration
  summarise(mean_sat = mean(n2osat)) %>% # 500 means for each WSA9
  # now summarize to 1 statistic per WSA9
  summarise( estimate = round(median(mean_sat), 3),
             LCL = round(quantile(mean_sat, probs = 0.025), 3),
             UCL = round(quantile(mean_sat, probs = 0.975), 3)) 
  # mutate(ecoregion = factor(WSA9_NAME)) %>%
  # mutate(ecoregion = fct_reorder(ecoregion, estimate))

# data for arrow segment.  Only smallest size category.  Only CPL and NPL.
p1.data.arrow <- p1.data %>%
  filter(size_cat == "min_4")

# pushing CPL and NPL to separate ggplot images to allow ggpubr
# to label each plot A, B, and C.  When pushing to two panels
# via faceting, ggpubr couldn't separately label the panels.

# CPL plot
p1.cpl <- p1.data %>%
  filter(WSA9_NAME == "Coastal Plains") %>%
  ggplot(aes(x=estimate, y=size_cat)) +
  geom_point() +
  geom_linerange( aes( xmin = LCL, xmax = UCL)) +
  geom_vline(xintercept = 1, color='blue') +
  # geom_segment(data = p1.data.arrow %>% filter(WSA9_NAME == "Coastal Plains"),
  #              aes(x=1, y = size_cat, xend = estimate, yend = size_cat),
  #              arrow = arrow(length = unit(0.2, "cm")),
  #              color = "red") +
  xlab(expression(mean~N[2]*O~saturation~ratio)) +
  xlim(0.84, 1.5) +
  ylab("waterbody size category (ha)") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  facet_grid(rows = vars(WSA9_NAME))

# NPL plot
p1.npl <- p1.data %>%
  filter(WSA9_NAME == "Northern Plains") %>%
  ggplot(aes(x=estimate, y=size_cat)) +
  geom_point() +
  geom_linerange( aes( xmin = LCL, xmax = UCL)) +
  geom_vline(xintercept = 1, color='blue') +
  # geom_segment(data = p1.data.arrow %>% filter(WSA9_NAME == "Northern Plains"),
  #              aes(x=1, y = size_cat, xend = estimate, yend = size_cat),
  #              arrow = arrow(length = unit(0.2, "cm")),
  #              color = "red") +
  xlab(expression(mean~N[2]*O~saturation~ratio)) +
  xlim(0.84, 1.5) +
  ylab("waterbody size category (ha)") +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  facet_grid(rows = vars(WSA9_NAME))

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
  ylab("waterbody size category (ha)") +
  #coord_flip() + 
  theme_bw() 

#ggsave("output/figures/n2oStarBySize.tiff")

ggpubr::ggarrange(p2, # first column
                  ggpubr::ggarrange(p1.cpl, p1.npl, ncol=1, # second column with plots in 2 rows
                                    nrow = 2, labels = c("B", "C")),
                  ncol = 2, 
                  labels = "A") # label for first plot
ggsave("manuscript/manuscript_figures/deltaN2ObySize.tiff", width = 6, height = 5)


}

## Figure X:  Flux by lake, emission rate by lake, and nation flux.  VS continuous size-------
if(!("n2oFluxAndEmissionRateVsContinuousArea.tiff" %in% list.files("manuscript/manuscript_figures"))) {

  #  flux vs waterbody size continuous
b1 <- all_predictions %>%
  #filter(.draw %in% 1:200) %>% # subset for practice
  group_by(.row) %>%
  # mean and CI of all realizations for each lake
  summarise(mean.flux = mean(f.n2o.Mg.y),
            LCI_flux = quantile(f.n2o.Mg.y, probs = 0.025),
            UCI_flux = quantile(f.n2o.Mg.y, probs = 0.975)) %>% 
  full_join(., # add other data back in
            all_predictions %>% 
              filter(.draw == 1) %>%# just 1 realization
              select(.row, WSA9_NAME, area_ha)) %>%
  ggplot(., aes(area_ha, mean.flux)) +
  geom_point(size = 0.1) +
  geom_errorbar(aes(xmin = area_ha, ymin = LCI_flux,
                    xmax = area_ha, ymax = UCI_flux),
                linewidth = 0.1) +
  scale_x_log10(labels=scales::comma) + 
  scale_y_continuous(trans = ggallin::pseudolog10_trans, 
                     breaks =c(-1000, -500, -50, -5, 0, 5, 50, 500, 1000)) +
  xlab("waterbody size (Ha)") +
  ylab(expression(N[2]*O~flux~"("*metric~tons~year^-1*")"))


# N2O emission rate vs waterbody size continuous
b2 <- all_predictions %>%
  #filter(.draw %in% 1:10) %>% # subset for practice
  group_by(.row) %>% # goup by lake
  # mean and CI of all realizations for each lake
  summarise(mean.e = mean(e.n2o.mg.m2.d),
            LCI_e = quantile(e.n2o.mg.m2.d, probs = 0.025),
            UCI_e = quantile(e.n2o.mg.m2.d, probs = 0.975)) %>% 
  full_join(., # add other data back in
            all_predictions %>% 
              filter(.draw == 1) %>%# just 1 realization
              select(.row, WSA9_NAME, area_ha)) %>%
  ggplot(., aes(area_ha, mean.e)) +
  geom_point(size = 0.1) +
  geom_errorbar(aes(xmin = area_ha, ymin = LCI_e,
                    xmax = area_ha, ymax = UCI_e),
                linewidth = 0.1) +
  scale_x_log10(labels=scales::comma) +
  scale_y_continuous(trans = ggallin::pseudolog10_trans) +
  xlab("waterbody size (Ha)") +
  ylab(expression(N[2]*O~emission~rate~"("*mg~N[2]*O~m^-2~day^-1*")"))


# National flux by lake size class: point and line range
b3 <- all_predictions %>%
  group_by(size_cat, .draw) %>% # group by iteration
  summarise(mean_f.n2o.Mg.y = sum(f.n2o.Mg.y)) %>% # 2000 means for size cat
  # now summarize to 1 statistic per size_cat
  summarise( estimate = round(median(mean_f.n2o.Mg.y), 3), 
             LCL = round(quantile(mean_f.n2o.Mg.y, probs = 0.025), 3),
             UCL = round(quantile(mean_f.n2o.Mg.y, probs = 0.975), 3)) %>% 
  ggplot(., aes(x = estimate, y = size_cat)) +
  geom_point() +
  geom_linerange(aes(xmin = LCL, xmax = UCL)) +
  geom_vline(xintercept = 0, color='blue') +
  xlab(expression(N[2]*O~flux~(metric~tons~year^{-1}))) + # 1Mg = 1 metric ton
  ylab("Lake size class (ha)") +
  theme_bw() 

# ggpubr::ggarrange(b1, b2, # first two rows
#                   ggarrange(b3, ncol = 2, nrow=1, labels = "AUTO"), # single plot in 3rd row, only one column
#                   nrow=3, labels = "AUTO")

ggpubr::ggarrange(b2, b1, b3, ncol = 3, nrow=1, labels = "AUTO")
ggsave("manuscript/manuscript_figures/n2oFluxAndEmissionRateVsContinuousArea.tiff", width = 8, height = 4)   

}

all_predictions %>%
  group_by(WSA9_NAME, .draw) %>% # group by iteration
  summarise(mean_f.n2o.Mg.y = sum(f.n2o.Mg.y)) %>% # 2000 means for size cat
  # now summarize to 1 statistic per WSA9
  summarise( estimate = round(median(mean_f.n2o.Mg.y), 3), 
             LCL = round(quantile(mean_f.n2o.Mg.y, probs = 0.025), 3),
             UCL = round(quantile(mean_f.n2o.Mg.y, probs = 0.975), 3)) %>% 
  ggplot(., aes(x = estimate, y = WSA9_NAME)) +
  geom_point() +
  geom_linerange(aes(xmin = LCL, xmax = UCL)) +
  geom_vline(xintercept = 0, color='blue') +
  xlab(expression(N[2]*O~flux~(metric~tons~year^{-1}))) + # 1Mg = 1 metric ton
  ylab("Lake size class (ha)") +
  theme_bw() 



## SI Figure 1: N2O emission rate distribution----
if(!("SIfigure1.tiff" %in% list.files("manuscript/manuscript_figures"))) {
  
# density plot by ecoregion
all_predictions %>%
  #filter(.draw %in% 1:10) %>% # subset for practice
  #mutate(ecoregion = factor(WSA9)) %>%
  #mutate(ecoregion = fct_reorder(ecoregion, preds_sat)) %>%
  ggplot(aes(x = e.n2o.mg.m2.d, group = .draw, color = .draw )) +
  geom_density(aes(y=after_stat(scaled)), show.legend = FALSE) +
  #geom_rug() + 
  xlab(expression("N"[2]*"O"~ "emission"~ "rate" ~ "(mg" ~ N[2]*O ~ m^-2 ~ d^-1 ~ ")")) + 
  ylab("density") +
  facet_wrap(~WSA9_NAME, scales = "free_x") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("manuscript/manuscript_figures/SIfigure1.tiff", width = 8.5, height = 5)

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


## Lake size distribution----
# proportion of total SA attributable to lakes <=50 ha
proportionSmallSA <- all_predictions %>%
  filter(.draw == 1) %>% # just grab one realization
  select(area_ha) %>%
  arrange(area_ha) %>% 
  mutate(cumulative = cumsum(area_ha), # cumulative surface area from small to large
         distribution = cumulative / sum(area_ha), # prop of total SA
         delta50 = abs(50-area_ha)) %>% # how close to 50 Ha
  filter(delta50 == min(delta50)) %>% # select closest to 50 Ha
  pull(distribution) # 50Ha lake
  
# proportion of population <=50ha
proportionSmallNumber <- all_predictions %>%
  filter(.draw == 1) %>% # just grab one realization
  select(area_ha) %>%
  arrange(area_ha) %>% 
  mutate(number = 1:nrow(.), # number each lake from smallest to largest
         proportion = number / nrow(.)) %>% # proportion of all lakes 
  filter(abs(50-area_ha) == min(abs(50-area_ha))) %>% # grab 50Ha lake
  pull(proportion) # extract proportion of all lakes represented by this lake
  
# flux from largest lake size class------
FluxBySize <- all_predictions %>%
  group_by(size_cat, .draw) %>% # group by iteration
  summarise(mean_f.n2o.Mg.y = sum(f.n2o.Mg.y)) %>% # 2000 means for size cat
  # now summarize to 1 statistic per size_cat
  summarise( estimate = round(median(mean_f.n2o.Mg.y), 3), 
             LCL = round(quantile(mean_f.n2o.Mg.y, probs = 0.025), 3),
             UCL = round(quantile(mean_f.n2o.Mg.y, probs = 0.975), 3))
  
  
# Mean and max emission rates-----
min.max <- all_predictions %>%
  group_by(.row) %>%
  summarize(e.n2o.mg.m2.d = mean(e.n2o.mg.m2.d)) %>%
  ungroup() %>%
  summarize(min = min(e.n2o.mg.m2.d),
            max = max(e.n2o.mg.m2.d))



## IPCC Indirect N2O Leaching and Runoff----
# 1990 - 2021 inventory reports 15.2 + 2.9 MMT CO2-Eq per year for CONUS (Table 5-19) in 2021.
# https://www.epa.gov/system/files/documents/2023-04/US-GHG-Inventory-2023-Chapter-5-Agriculture.pdf
# 15.2 and 2.9 via tables A-181 and A-182 in Annex 3.12, 
# https://www.epa.gov/system/files/documents/2023-04/US-GHG-Inventory-2023-Annex-3-Additional-Source-or-Sink-Categories-Part-B.pdf

# 1 mass of N2O = 265 mass CO2 www.epa.gov/system/files/documents/2023-04/US-GHG-Inventory-2023-Chapter-1-Introduction.pdf
# 1 MMT (million metric ton) = 1000 kt (kilo ton) www.epa.gov/system/files/documents/2023-04/US-GHG-Inventory-2023-Chapter-1-Introduction.pdf
# 1 ton = 1000 kg
# 1000 kg = 1 Mg = 1 metric ton
# (15.2+2.9) = 18 MMT CO2 eq * (1MMT N2O / 265 MMT CO2 Eq) * (1000kt N2O/1MMt N2O) * (1000t N2O/1kt N2O) * (1000kg N2O/1t N2O) * (1Mg/1000kg)
ipcc.n2o.indirect <- (15.2+2.9)*(1/265)*(1000/1)*(1000/1)*(1000/1)*(1/1000) #68,302 Mg N2O

# this is indirect emissions from groundwater (EF5-g), streams/rivers (EF5-r), and estuaries (EF5-e).  The default EF is
# 0.0025 for each, for an EF5 of 0.0075.  Therefore, only 1/3 of indirect N2O emissions reported above
# should be prescribed to streams/rivers/lakes/reservoirs.
ipcc.n2o.surface = round(ipcc.n2o.indirect * (1/3), digits = 0) # 22,767
