## This script generates data to be referenced in inline-code in the manuscript.
## This makes the inline code more succinct and the .Rmd easier to read.

# 32 minutes to run on DMAP with images already stored locally

# SETUP----------------------
## Identify local path for each user
#localPath <- Sys.getenv("USERPROFILE")

## load libraries
# DMAP not under renv control.  Should probably specify same
# versions as used for renv library, but causing issues.  Just 
# install newest version.
# if(localPath == "") { # if DMAP, then
#   # require returns FALSE and gives a warning 
#   # (rather than an error as library() does by default) 
#   # if the package does not exist.
#   
#   # older version of dplyr already installed on instance.  The only way I could
#   # figure out how to update was to use devtools::install_version
#   if(!require("devtools")) install.packages("devtools")
#   devtools::install_version("dplyr", "1.1.4", 
#                             repos = "https://cran.r-project.org") # get newest version
#   if(!require("ggallin")) install.packages("ggallin")
#   if(!require("sf")) install.packages("sf")
#   if(!require("tidyverse")) install.packages("tidyverse")
#   if(!require("janitor")) install.packages("janitor")
#   if(!require("tictoc")) install.packages("tictoc")
#   if(!require("pbkrtest")) install.packages("pbkrtest") 
#   if(!require("ggpubr")) install.packages("ggpubr")
#   if(!require("ggallin")) install.packages("ggallin")
#   if(!require("tidybayes")) install.packages("tidybayes")
# }

library(sf) # spatial data
library(tidyverse) # dplyr, ggplot
library(janitor) # clean names
library(usmap) # state boundaries 
library(tictoc) # processing time
library(ggpubr) # multiple plots
library(ggallin) # psuedlolog transformation for negative values
library(tidybayes) # prop flux by size figure


if(!("figure1.png" %in% list.files("manuscript/manuscript_figures"))) {

# DATA---------------
## load sample data (dg object)----
load("./inputData/dg.RData")
save(dg, file = "./manuscript/manuscript_files/dg.rda")  # save dg file to ms folder for faster knitting
  
# To enable spatial analysis of the data, the dataframe will be converted to a 
# 'simple features' (sf) object.
# Define coordinates
coords <- data.frame(longitude = dg$map.lon.dd, latitude = dg$map.lat.dd)

dg.sf <- st_as_sf(dg, coords = c("map.lon.dd", "map.lat.dd"), 
                  crs = 4269) %>% # standard for lat/lon
  st_transform(5070) # project to CONUS ALBERS for plotting

## load ecoregions----
# read in ecoregion polygons
ecoR <- st_read("/vsizip/inputData/aggr_ecoregions_simple.zip", quiet = TRUE)
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


states <- usmap::us_map() %>%
  dplyr::filter(!full %in% c("Alaska", "District of Columbia", "Hawaii", "Puerto Rico")) %>%
  sf::st_transform(5070) # convert to CONUS Albers
} # fig 1

if(!("ecoreg_stats.rda" %in% list.files("manuscript/manuscript_files"))) {
  if(!("national_stats.rda" %in% list.files("manuscript/manuscript_files"))) {
    if(!("figure2.tiff" %in% list.files("manuscript/manuscript_figures"))) {
      if(!("n2oFluxAndEmissionRateVsContinuousArea.tiff" %in% list.files("manuscript/manuscript_figures"))) {
        if(!("n2oStarBySize.tiff" %in% list.files("manuscript/manuscript_figures"))) {
          if(!("SIfigure1.tiff" %in% list.files("manuscript/manuscript_figures"))) {
            if(!("FluxBySize.rda" %in% list.files("manuscript/manuscript_files"))) {
              if(!("em_min_max.rda" %in% list.files("manuscript/manuscript_files"))) {
      
## load posterior predictions from modeling----
#tic() # 7.6min on DMAP.  5 minutes on Dell Precision workstation
if(!("all_predictions.rda" %in% list.files("inputData"))){
  inUrl1  <- "https://zenodo.org/records/15159394/files/all_predictions.rda?download=1"
  infile1 <- tempfile()
  try(download.file(inUrl1,infile1,method="curl")) 
  if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")
  
  save(all_predictions, file = "inputData/all_predictions.rda")
} else {
  load("inputData/all_predictions.rda")
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


all_predictions_ms <- full_join(all_predictions, wsa9names)

rm(all_predictions) # free up RAM

## emission rates----
# code copied from dgIndicatorAnalysis.Rmd
### read wind data----
# Mean index period wind speed was calculated on VM via 
# nlaSampleFrameWind.R.  The resulting data object was written out and read in below.
wind <- read.table("./inputData/windSampleFrame.txt")

# print first 6 lines
wind %>% head()

# print first 6 .row values of all_predictions
#all_predictions_ms %>%
#  filter(.row %in% 1:6) %>% # each row is a unique waterbody
#  group_by(.row) %>%
#  slice(1) # grabs first record

# Row order in wind and all_predictions is inherited from sample frame.  Thus row 1
# of wind corresponds to .row == 1 in all_predictions.  This is a sketchy way to merge
# things, but will do for now.

# add .id to wind which is simply row number
wind <- wind %>% mutate(.id = row_number())

# Merge wind and all_predictions
#tic() # 21 seconds on VM
# left join needed to exclude great salt lake in wind data.  see line 88-93
all_predictions_ms <- left_join(all_predictions_ms, wind, by = c(".row" = ".id")) # 40 seconds on laptop
#toc()
#dim(wind) # 465897
#dim(all_predictions_ms) #931792000, 2000 observations for each lake
#465896*2000 #931792000

#summary(wind$ws)

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
all_predictions_ms <- all_predictions_ms %>%
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

all_predictions_ms  <- all_predictions_ms %>%
  select(.row, 
         .draw,
         WSA9,
         WSA9_NAME,
         state,
         size_cat, 
         area_ha, 
         no3_cat, 
         n2o, 
         n2oeq, 
         n2osat,
         e.n2o.mg.m2.d,
         f.n2o.Mg.y)
              } # ecoreg_stats
            } # national_stats
          } #fig2
        } # flux
      } # n2o star
    } # SI figure
  } # FluxBySize
} # emissions min/max

if(!("ecoreg_stats.rda" %in% list.files("manuscript/manuscript_files"))) {
## ecoregion and national means----
ecoreg_stats <- all_predictions_ms %>%
  group_by(WSA9, WSA9_NAME, .draw) %>%
  summarise( mean_n2o = mean(n2o),
             mean_n2oeq = mean(n2oeq),
             mean_sat = mean(n2osat),
             median_sat = median(n2osat),
             prop_sat = sum(n2osat < 1) / length(unique(.row)),
             mean_en2o = mean(e.n2o.mg.m2.d),
             median_en2o = median(e.n2o.mg.m2.d),
             fn2o_Mg = sum(f.n2o.Mg.y),
             median_no3 = median(no3_cat),
             total_sa = sum(area_ha),
             prop_small = sum(size_cat == "min_4") / length(unique(.row))) %>%
  summarise( post_m_mean_n2o = round(mean(mean_n2o), 1),
             LCL_mean_n2o = round(quantile(mean_n2o, probs = 0.025), 1),
             UCL_mean_n2o = round(quantile(mean_n2o, probs = 0.975), 1),
             post_m_mean_n2oeq = round(mean(mean_n2oeq), 2),
             LCL_mean_n2oeq = round(quantile(mean_n2oeq, probs = 0.025), 2),
             UCL_mean_n2oeq = round(quantile(mean_n2oeq, probs = 0.975), 2),
             post_m_mean_sat = round(mean(mean_sat), 2),
             post_m_median_sat = round(mean(median_sat), 2),
             LCL_mean_sat = round(quantile(mean_sat, probs = 0.025), 2),
             UCL_mean_sat = round(quantile(mean_sat, probs = 0.975), 2),
             post_m_prop_sat = round(mean(prop_sat), 3),
             LCL_prop_sat = round(quantile(prop_sat, probs = 0.025), 3),
             UCL_prop_sat = round(quantile(prop_sat, probs = 0.975), 3),
             post_m_mean_en2o = round(mean(mean_en2o), 3),
             post_m_median_en2o = round(mean(median_en2o), 3),
             LCL_median_en2o = round(quantile(median_en2o, probs = 0.025), 3),
             UCL_median_en2o = round(quantile(median_en2o, probs = 0.975), 3),
             LCL_mean_en2o = round(quantile(mean_en2o, probs = 0.025), 3),
             UCL_mean_en2o = round(quantile(mean_en2o, probs = 0.975), 3),
             post_m_fn2o = round(mean(fn2o_Mg), 3), 
             LCL_fn2o = round(quantile(fn2o_Mg, probs = 0.025), 3),
             UCL_fn2o = round(quantile(fn2o_Mg, probs = 0.975), 3),
             post_m_median_no3 = round(mean(median_no3), 4),
             LCL_median_no3 = round(quantile(median_no3, probs = 0.025), 4),
             UCL_median_no3 = round(quantile(median_no3, probs = 0.975), 4),
             post_m_prop_small = round(mean(prop_small), 1),
             LCL_prop_small = round(quantile(prop_small, probs = 0.025), 1),
             UCL_prop_small = round(quantile(prop_small, probs = 0.975), 1),
             total_sa = round(mean(total_sa), digits = 0),
             n = dplyr::n()) %>%
  ungroup()

save(ecoreg_stats, file = "manuscript/manuscript_files/ecoreg_stats.rda")
} # ecoreg_stats

if(!("national_stats.rda" %in% list.files("manuscript/manuscript_files"))) {
national_stats <- all_predictions_ms %>%
  group_by(.draw) %>%
  summarise( mean_n2o = mean(n2o),
             mean_n2oeq = mean(n2oeq),
             mean_sat = mean(n2osat),
             prop_sat = sum(n2osat < 1) / length(unique(.row)),
             mean_en2o = mean(e.n2o.mg.m2.d),
             median_en2o = median(e.n2o.mg.m2.d),
             fn2o_Mg = sum(f.n2o.Mg.y),
             total_sa = sum(area_ha)) %>%
  summarise( post_m_mean_n2o = round(mean(mean_n2o), 1),
             LCL_mean_n2o = round(quantile(mean_n2o, probs = 0.025), 1),
             UCL_mean_n2o = round(quantile(mean_n2o, probs = 0.975), 1),
             post_m_mean_n2oeq = round(mean(mean_n2oeq), 2),
             LCL_mean_n2oeq = round(quantile(mean_n2oeq, probs = 0.025), 2),
             UCL_mean_n2oeq = round(quantile(mean_n2oeq, probs = 0.975), 2),
             post_m_mean_sat = round(mean(mean_sat), 2),
             LCL_mean_sat = round(quantile(mean_sat, probs = 0.025), 2),
             UCL_mean_sat = round(quantile(mean_sat, probs = 0.975), 2),
             post_m_prop_sat = round(mean(prop_sat), 3),
             LCL_prop_sat = round(quantile(prop_sat, probs = 0.025), 3),
             UCL_prop_sat = round(quantile(prop_sat, probs = 0.975), 3),
             post_m_median_en2o = round(mean(median_en2o), 3),
             LCL_median_en2o = round(quantile(median_en2o, probs = 0.025), 3),
             UCL_median_en2o = round(quantile(median_en2o, probs = 0.975), 3),
             post_m_mean_en2o = round(mean(mean_en2o), 3),
             LCL_mean_en2o = round(quantile(mean_en2o, probs = 0.025), 3),
             UCL_mean_en2o = round(quantile(mean_en2o, probs = 0.975), 3),
             post_m_fn2o = round(mean(fn2o_Mg), 3), 
             LCL_fn2o = round(quantile(fn2o_Mg, probs = 0.025), 3),
             UCL_fn2o = round(quantile(fn2o_Mg, probs = 0.975), 3),
             total_sa = round(mean(total_sa), digits = 0),
             n = dplyr::n()) %>%
  mutate(WSA9 = "National", WSA9_NAME = "National")

save(national_stats, file = "manuscript/manuscript_files/national_stats.rda")
} # national_stats


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
    geom_sf(data = dg.sf %>% 
              filter(!is.na(n2o.src.snk.error), sitetype == "PROB", visit.no == 1) %>%
              arrange(desc(n2o.src.snk.error)), # grey on top of black
            aes(size = dissolved.n2o.nmol, color = n2o.src.snk.error),
            show.legend = "point") +
    geom_sf(data = states, fill = NA, color = "cornsilk3", size = 0.1) +
    # guide argument below removes points from the boxes in the ecoregion legend
    # https://aosmith.rbind.io/2020/07/09/ggplot2-override-aes/
    scale_fill_manual("Ecoregion", values = cols,
                      guide = guide_legend(override.aes = list(shape = NA))) +
    scale_color_manual(values = c("red4", "black", "grey"), name = "source/sink") +
    scale_size(name = expression(N[2]*O~(nM)),
               range = c(0.1, 10), # custom size range
               breaks = c(1, 10, 25, 50, 100)) + # custom breaks
    theme(legend.key.size = unit(0.4, "cm"), # size of boxes in legend
          legend.title = element_text(size = 8),
          plot.margin = unit(c(0.3, 0, 0.2, 0), 
                             "inches"))
  ggsave("manuscript/manuscript_figures/figure1.png", width = 8, height = 4, units = "in")
} # fig 1


## Figure 2: N2O saturation ratio density plot----
# code taken from dgIndicatorAnalysis.Rmd, "N2O saturation ratio: continuous variable
# if the image is already on computer, then nothing, else create image
if(!("figure2.tiff" %in% list.files("manuscript/manuscript_figures"))) {
  
  # density plots of predicted distribution, mean, and median
  # using xlim to zoom in on values close to 1
  
  dummy <- ecoreg_stats %>%
                    #CPL   NAP   NPL    SAP    SPL   TPL  UWM    WMT    XER
    mutate(dens.med = c(2.5,  4,    3.1,   3.5,   2.5,  3,   3,      2,     2.5),
           dens.mean = c(.5,  3,    2.5,   1.5,   1,    1,   2,      1.5,   1.1),
           # dens.equ = rep(1, 9),
           # `median legend` = "median",
           # `mean legend` = "mean",
           # `eq legend` = "equilibrium",
           .draw = 1) # needed to match grouping aesthetic in ggplot call
  
  tic() # 60 seconds for 100 draws, and 16 minutes for 2000 draws, on memory intensive DMAP, 2797 secs for 2000 draws on Precision workstation
  # After spending many hours trying to get the legend right, I gave up.  Write
  # image to disk, then add legend in powerpoint.
  all_predictions_ms %>%
    #filter(.draw %in% 1:100) %>% # subset for practice
    ggplot(aes(x = n2osat, group = .draw, color = .draw )) +
    geom_density(show.legend = FALSE) + # this argument doesn't work when faceted, hence guides call below
    geom_vline(xintercept = 1, color = "red") +
    geom_segment(data = dummy,
                 aes(x = post_m_mean_sat, xend = post_m_mean_sat,
                     y = 0, yend = dens.mean),
                 color = "black", linetype = "solid") +
    geom_segment(data = dummy,
                 aes(x = post_m_median_sat, xend = post_m_median_sat,
                     y = 0, yend = dens.med),
                     color = "black", linetype = "dotted") +
    # geom_segment(#data = dummy,
    #              aes(x = dens.equ, xend = dens.equ,
    #                  y = 0, yend = 4, linetype = `eq legend`),
    #              color = "red") +
    ylab ("density") +
    xlab(expression(N[2]*O~saturation~ratio)) +
    xlim(0,3) + # data range is 0-600, but sample max is 30
    facet_wrap(~WSA9_NAME) +
    theme_bw()
  toc()
  ggsave("manuscript/manuscript_figures/figure2.tiff", width = 8.5, height = 5)
  
} # fig 2

## fig 3 Conditional effects of NO3 on sat ratio----
if(!("cond_effect_no3.tiff" %in% list.files("manuscript/manuscript_figures"))) {
  load("./modelFiles/n2o_mod6.rda")
  load("./modelFiles/df_model.rda")
  effect_no3_surftemp <- expand.grid(no3_cat = levels(df_model$no3_cat),
                                     surftemp = round(quantile(df_model$surftemp, 
                                                               probs = seq(0.10, 0.90, 0.2)), 0),
                                     log_area = mean(df_model$log_area),
                                     log_elev = mean(df_model$log_elev)) %>%
    add_epred_draws(n2o_mod6, resp = c("n2o", "n2oeq"), re_formula = ~ 1) %>%
    pivot_wider(names_from = .category, values_from = .epred) %>%
    mutate(n2o_pred = n2o,
           n2oeq_pred = n2oeq,
           sat_pred = n2o / n2oeq,
           surftemp = factor(surftemp))
  
  effect_surftemp_no3 <- expand.grid(no3_cat = levels(df_model$no3_cat),
                                     surftemp = seq(range(df_model$surftemp)[1],
                                                    range(df_model$surftemp)[2],
                                                    length.out = 101),
                                     log_area = mean(df_model$log_area),
                                     log_elev = mean(df_model$log_elev)) %>%
    add_epred_draws(n2o_mod6, resp = c("n2o", "n2oeq"), re_formula = ~ 1) %>%
    pivot_wider(names_from = .category, values_from = .epred) %>%
    mutate(n2o_pred = n2o,
           n2oeq_pred = n2oeq,
           sat_pred = n2o / n2oeq,
           surftemp = surftemp)
  
  
  effect_no3_area <- expand.grid(no3_cat = levels(df_model$no3_cat),
                                 log_area = round(quantile(df_model$log_area, 
                                                           probs = seq(0.10, 0.90, 0.2)), 0),
                                 surftemp = mean(df_model$surftemp),
                                 log_elev = mean(df_model$log_elev)) %>%
    add_epred_draws(n2o_mod6, resp = c("n2o", "n2oeq"), re_formula = ~ 1) %>%
    pivot_wider(names_from = .category, values_from = .epred) %>%
    mutate(n2o_pred = n2o,
           n2oeq_pred = n2oeq,
           sat_pred = n2o / n2oeq,
           log_area = factor(log_area))
  
  effect_area_no3 <- expand.grid(no3_cat = levels(df_model$no3_cat),
                                 log_area = seq(range(df_model$log_area)[1],
                                                range(df_model$log_area)[2],
                                                length.out = 101),
                                 surftemp = mean(df_model$surftemp),
                                 log_elev = mean(df_model$log_elev)) %>%
    add_epred_draws(n2o_mod6, resp = c("n2o", "n2oeq"), re_formula = ~ 1) %>%
    pivot_wider(names_from = .category, values_from = .epred) %>%
    mutate(n2o_pred = n2o,
           n2oeq_pred = n2oeq,
           sat_pred = n2o / n2oeq,
           log_area = log_area)
  
  plt_sat_no3_temp <-
    effect_no3_surftemp %>%
    ggplot(aes(x = no3_cat, 
               y = sat_pred, 
               group = surftemp,
               color = surftemp)) +
    stat_pointinterval(shape = 95,
                       point_interval = "median_qi",
                       position = position_dodge(width = 1.4, preserve = "single")) +
    scale_color_brewer(palette = "PuBu") +
    xlab(expression(paste("NO"[3], " category"))) +
    ylab("") +
    theme_pubr() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 7), breaks = seq(0, 7, 1)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(legend.position = c(0.4, 0.6),
          axis.text.y = element_blank()) +
    labs(color = "Surface temperature (\u00B0C)")
  
  plt_sat_temp_no3 <-
    effect_surftemp_no3 %>%
    ggplot() +
    stat_lineribbon(aes(x = surftemp,
                        y = sat_pred, 
                        group = no3_cat,
                        color = no3_cat),
                    point_interval = "median_qi",
                    show.legend = TRUE,
                    size = 0.7,
                    .width = 0.95) +
    scale_color_brewer(palette = "OrRd") +
    xlab("Surface temperature (\u00B0C)") +
    ylab("Saturation ratio") +
    theme_pubr() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 7), breaks = seq(0, 7, 1)) +
    scale_x_continuous(expand = c(0, 0)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_fill_brewer(palette = "Greys", guide = "none") +
    theme(legend.position = c(0.4, 0.6)) + 
    labs(color = expression(paste("NO"[3], " category")))
  
  plt_sat_no3_area <-
    effect_no3_area %>%
    ggplot(aes(x = no3_cat, 
               y = sat_pred, 
               group = log_area,
               color = log_area)) +
    stat_pointinterval(shape = 95,
                       point_interval = "median_qi",
                       position = position_dodge(width = 1.4, preserve = "single")) +
    scale_color_brewer(palette = "Greys") +
    xlab(expression(paste("NO"[3], " category"))) +
    ylab("") +
    theme_pubr() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 7), breaks = seq(0, 7, 1)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(legend.position = c(0.4, 0.6),
          axis.text.y = element_blank()) +
    labs(color = "log(surface area (ha))")
  
  plt_sat_area_no3 <-
    effect_area_no3 %>%
    ggplot() +
    stat_lineribbon(aes(x = log_area,
                        y = sat_pred, 
                        group = no3_cat,
                        color = no3_cat),
                    point_interval = "median_qi",
                    show.legend = TRUE,
                    size = 0.7,
                    .width = 0.95) +
    scale_color_brewer(palette = "OrRd") +
    xlab("log(surface area (ha))") +
    ylab("Saturation ratio") +
    theme_pubr() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 7), breaks = seq(0, 7, 1)) +
    scale_x_continuous(expand = c(0, 0)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_fill_brewer(palette = "Greys", guide = "none") +
    theme(legend.position = c(0.6, 0.6)) + 
    labs(color = expression(paste("NO"[3], " category")))
  
  plt_temp_area_effects <- ggarrange(plt_sat_temp_no3,
                                     plt_sat_no3_temp,
                                     plt_sat_area_no3,
                                     plt_sat_no3_area,
                                     ncol = 2,
                                     nrow = 2,
                                     widths = c(0.9, 0.9, 0.9, 0.9),
                                     heights = c(0.9, 0.9, 0.9, 0.9),
                                     labels = c("A", "B", "C", "D"),
                                     hjust = -7,
                                     vjust = 3
  )
  
#  plt_temp_area_effects <- plt_temp_area_effects %>%
#    annotate_figure(plt_temp_area_effects,
#                    left = text_grob("Derived posterior mean", 
#                                     size = 18, 
#                                     rot = 90)
#    )
  
  # ggsave only works in R console (i.e., copy + paste to save)
  ggsave(filename = "./manuscript/manuscript_figures/cond_effect_no3.tiff",
         plot = plt_temp_area_effects,
         device = "tiff",
         width = 210 / 25.4,
         height = 210 / 25.4,
         units = "in",
         bg = "white")
  
} # fig 3

## Figure X: delta N2O by waterbody size----
# See 'N2O saturation ratio: continuous variable' section of dgIndicatorAnalysis.Rmd
# if the image is already on computer, then nothing, else create image
if(!("n2oStarBySize.tiff" %in% list.files("manuscript/manuscript_figures"))) {
  
  
  # originally drafted as three panel figure, but decided to scale down to one panel.
  # PLOT N2O SAT RATIO BY WSA9 AND SIZE
  # # point and linerange
  # p1.data <- all_predictions %>%
  #   filter(WSA9 %in% c("CPL", "NPL")) %>% # only show two examples.  comment out to include all ecoregions
  #   group_by(WSA9_NAME, size_cat, .draw) %>% # group by iteration
  #   summarise(mean_sat = mean(n2osat)) %>% # 500 means for each WSA9
  #   # now summarize to 1 statistic per WSA9
  #   summarise( estimate = round(median(mean_sat), 3),
  #              LCL = round(quantile(mean_sat, probs = 0.025), 3),
  #              UCL = round(quantile(mean_sat, probs = 0.975), 3)) 
  #   # mutate(ecoregion = factor(WSA9_NAME)) %>%
  #   # mutate(ecoregion = fct_reorder(ecoregion, estimate))
  # 
  # # data for arrow segment.  Only smallest size category.  Only CPL and NPL.
  # p1.data.arrow <- p1.data %>%
  #   filter(size_cat == "min_4")
  # 
  # # pushing CPL and NPL to separate ggplot images to allow ggpubr
  # # to label each plot A, B, and C.  When pushing to two panels
  # # via faceting, ggpubr couldn't separately label the panels.
  # 
  # # CPL plot
  # p1.cpl <- p1.data %>%
  #   filter(WSA9_NAME == "Coastal Plains") %>%
  #   ggplot(aes(x=estimate, y=size_cat)) +
  #   geom_point() +
  #   geom_linerange( aes( xmin = LCL, xmax = UCL)) +
  #   geom_vline(xintercept = 1, color='blue') +
  #   # geom_segment(data = p1.data.arrow %>% filter(WSA9_NAME == "Coastal Plains"),
  #   #              aes(x=1, y = size_cat, xend = estimate, yend = size_cat),
  #   #              arrow = arrow(length = unit(0.2, "cm")),
  #   #              color = "red") +
  #   xlab(expression(mean~N[2]*O~saturation~ratio)) +
  #   xlim(0.84, 1.5) +
  #   ylab("waterbody size category (ha)") +
  #   theme_bw() +
  #   theme(axis.title.y = element_blank(),
  #         axis.title.x = element_blank()) +
  #   facet_grid(rows = vars(WSA9_NAME))
  # 
  # # NPL plot
  # p1.npl <- p1.data %>%
  #   filter(WSA9_NAME == "Northern Plains") %>%
  #   ggplot(aes(x=estimate, y=size_cat)) +
  #   geom_point() +
  #   geom_linerange( aes( xmin = LCL, xmax = UCL)) +
  #   geom_vline(xintercept = 1, color='blue') +
  #   # geom_segment(data = p1.data.arrow %>% filter(WSA9_NAME == "Northern Plains"),
  #   #              aes(x=1, y = size_cat, xend = estimate, yend = size_cat),
  #   #              arrow = arrow(length = unit(0.2, "cm")),
  #   #              color = "red") +
  #   xlab(expression(mean~N[2]*O~saturation~ratio)) +
  #   xlim(0.84, 1.5) +
  #   ylab("waterbody size category (ha)") +
  #   theme_bw() +
  #   theme(axis.title.y = element_blank()) +
  #   facet_grid(rows = vars(WSA9_NAME))
  
  # PLOT N2O* BY SIZE
  # point and linerange
  p2 <- all_predictions_ms  %>%
    mutate(
      size_cat = fct_recode(as.factor(size_cat),
                            "< 4" = "min_4",
                            "4 to < 10" = "4_10",
                            "10 to < 20" = "10_20",
                            "20 to < 50" = "20_50",
                            "> 50" = "50_max")) %>%
    mutate(n2ostar = abs(n2o - n2oeq)) %>%
    group_by(size_cat, .draw) %>% # group by iteration
    #group_by(.draw) %>%
    summarise(mean_star = mean(n2ostar)) %>% # 2000 means for each WSA9
    # now summarize to 1 statistic per WSA9
    summarise( post_m_mean = round(mean(mean_star), 3),
               LCL_mean = round(quantile(mean_star, probs = 0.025), 3),
               UCL_mean = round(quantile(mean_star, probs = 0.975), 3)) %>% 
    ggplot( aes( x = post_m_mean, y = size_cat ) ) +
    geom_point() +
    geom_linerange( aes( xmin = LCL_mean, xmax = UCL_mean)) +
    xlab(expression("["*Delta~N[2]*O*"]"~"("*nmol~L^{-1}*")")) +
    ylab("waterbody size category (ha)") +
    #coord_flip() + 
    theme_bw() 
  
  ggsave(plot = p2, filename = "manuscript/manuscript_figures/n2oStarBySize.tiff", width = 3, height = 3)
  
  # ggpubr::ggarrange(p2, # first column
  #                   ggpubr::ggarrange(p1.cpl, p1.npl, ncol=1, # second column with plots in 2 rows
  #                                     nrow = 2, labels = c("B", "C")),
  #                   ncol = 2, 
  #                   labels = "A") # label for first plot
  # ggsave("manuscript/manuscript_figures/deltaN2ObySize.tiff", width = 6, height = 5)
  
  
} # fig 4: n2o star

## Figure X:  Flux by lake, emission rate by lake, and national flux VS continuous size-------
if(!("n2oFluxAndEmissionRateVsContinuousArea.tiff" %in% list.files("manuscript/manuscript_figures"))) {
  
  # Predicted (mean, L95CI, U95CI of PPD) waterbody N2O emission rates vs waterbody size continuous 
  b1.dat <- all_predictions_ms %>%
    #filter(.draw %in% 1:10) %>% # subset for practice
    group_by(.row) %>% # goup by waterbody
    # median and CI of all realizations for each waterbody
    summarise(post_m_pred = mean(e.n2o.mg.m2.d),
              LCL_pred = quantile(e.n2o.mg.m2.d, probs = 0.025),
              UCL_pred = quantile(e.n2o.mg.m2.d, probs = 0.975)) %>%
    full_join(., # add other data back in
              all_predictions_ms %>% 
                filter(.draw == 1) %>%# just 1 realization
                select(.row, WSA9_NAME, area_ha))
  
  b1 <- ggplot(b1.dat, aes(x = area_ha, y = post_m_pred)) +
    geom_errorbar(aes(xmin = area_ha, ymin = LCL_pred,
                      xmax = area_ha, ymax = UCL_pred),
                  linewidth = 0.1) +
    geom_point(size = 0.1, color = "red") +
    scale_x_log10(labels=scales::comma) +
    scale_y_continuous(trans = ggallin::pseudolog10_trans) +
    ylab(expression(atop(Areal~Emission~rate, "("*mg~N[2]*O~m^-2~day^-1*")"))) +
    theme_bw()  +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank())
  
  #  predicted flux vs waterbody size continuous
  b2.dat <- all_predictions_ms %>%
    #filter(.draw %in% 1:200) %>% # subset for practice
    group_by(.row) %>%
    # mean and CI of all realizations for each lake
    summarise(post_m_pred = mean(f.n2o.Mg.y),
              LCL_pred = quantile(f.n2o.Mg.y, probs = 0.025),
              UCL_pred = quantile(f.n2o.Mg.y, probs = 0.975)) %>% 
    full_join(., # add other data back in
              all_predictions_ms %>% 
                filter(.draw == 1) %>%# just 1 realization
                select(.row, WSA9_NAME, area_ha))
  
  b2 <- ggplot(b2.dat, aes(x = area_ha, y = post_m_pred)) +
    geom_errorbar(aes(xmin = area_ha, ymin = LCL_pred,
                      xmax = area_ha, ymax = UCL_pred),
                  linewidth = 0.1) + # size with ggplot 3.3.3.  linewidth with other versions?
    geom_point(size = 0.1, color = "red") +
    scale_x_log10(labels=scales::comma) + 
    scale_y_continuous(trans = ggallin::pseudolog10_trans, 
                       breaks =c(-1000, -500, -50, -5, 0, 5, 50, 500, 1000)) +
    #xlab("waterbody size (Ha)") +
    ylab(expression(atop(Annual~Flux,"("*metric~tons~N[2]*O~year^-1*")"))) +
    theme_bw()  +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank())
  
  
  # Proportion of total CONUS flux by lake size class------
  PropFluxBySize <- all_predictions_ms %>%
    group_by(.draw) %>% # group by iteration
    mutate(sum_f.n2o.Mg.y = sum(f.n2o.Mg.y)) %>% # 2000 estimates for total flux across CONUS
    #left_join(all_predictions_ms, by = ".draw") %>% # rejoin to full predictions (note duplicate sum flux value for same draw)
    # now sum flux by size cat (and draw) then calc proportion of total CONUS flux for each size class for each draw
    group_by(size_cat, .draw) %>%
    mutate(size_f = sum(f.n2o.Mg.y),
           size_p_f = (size_f / sum_f.n2o.Mg.y)) %>% # This step leaves duplicates for .draw
    #distinct(.draw, .keep_all = TRUE) %>% # remove duplicate data
    group_by(size_cat) %>%
    # Now summarize (over draws) the posterior distributions of proportions by mean and 95% CI
    summarise(post_m_prop = round(mean(size_p_f), 3), 
              LCL = round(quantile(size_p_f, probs = 0.025), 3),
              UCL = round(quantile(size_p_f, probs = 0.975), 3))
  
  save(PropFluxBySize, file = "manuscript/manuscript_files/PropFluxBySize.rda")
  
  ## Figure on proportion of flux (ratio size class / CONUS)
  ############## ROY
  ##############
  # OBJECT PropFluxBySize DOES NOT HAVE size_f and size_p_f variables.
  # they were lost during the final summarize above (551-553)
  #b3 <- PropFluxBySize %>% 
  #  mutate(total_f = size_f / size_p_f) %>% # CODE BREAKS HERE
  #  select(size_cat, .draw, total_f, size_f, size_p_f) %>%
  #  mutate(
  #    size_cat = fct_recode(as.factor(size_cat),
  #                          "< 4" = "min_4",
  #                          "4 to < 10" = "4_10",
  #                          "10 to < 20" = "10_20",
  #                          "20 to < 50" = "20_50",
  #                          "> 50" = "50_max")) %>%
  #  ggplot(aes(x = factor(size_cat), y = size_p_f)) +
    #ylim(-1, 1) +
    #geom_line(aes(x = size_cat, y = size_p_f, group = .draw), color = "grey80") +
  #  tidybayes::stat_lineribbon(color = 'red', .width = 0.95, alpha = 0.7) +  
    #tidybayes::stat_pointinterval(color = 'black', linewidth = 2) +
  #  scale_fill_manual(values = 'black') +
  #  guides(fill = "none") +
  #  coord_cartesian(ylim = c(-3, 3)) +
  #  xlab("Size category (ha)") +
  #  ylab(expression(atop(Flux~ratio,"(Size category flux / CONUS flux)"))) +
  #  theme_bw()
  
  b3.dat <- all_predictions_ms %>%
    group_by(.draw) %>% # group by iteration
    mutate(sum_f.n2o.Mg.y = sum(f.n2o.Mg.y)) %>% # 2000 estimates for total flux across CONUS
    ungroup() %>%
    mutate(p_total_f = f.n2o.Mg.y / sum_f.n2o.Mg.y) %>%
    select(.row, .draw, area_ha, p_total_f) %>%
    group_by(.row) %>% 
    mean_qi()

  b3 <-  ggplot(b3.dat, aes(x = area_ha, y = p_total_f)) +
    geom_errorbar(aes(xmin = area_ha, ymin = p_total_f.lower,
                      xmax = area_ha, ymax = p_total_f.upper),
                  linewidth = 0.1) + 
    geom_point(size = 0.1, color = "red") +
    scale_x_log10(labels = scales::comma) +
    coord_cartesian(ylim = c(-3, 3)) +
    xlab("Waterbody size (ha)") +
    ylab(expression(atop(Proportion, of~CONUS~flux))) +
    theme_bw()  +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  # # National flux by lake size class: point and line range
  # b3 <- all_predictions %>%
  #   group_by(size_cat, .draw) %>% # group by iteration
  #   summarise(mean_f.n2o.Mg.y = sum(f.n2o.Mg.y)) %>% # 2000 means for size cat
  #   # now summarize to 1 statistic per size_cat
  #   summarise( estimate = round(median(mean_f.n2o.Mg.y), 3), 
  #              LCL = round(quantile(mean_f.n2o.Mg.y, probs = 0.025), 3),
  #              UCL = round(quantile(mean_f.n2o.Mg.y, probs = 0.975), 3)) %>% 
  #   ggplot(., aes(x = estimate, y = size_cat)) +
  #   geom_point() +
  #   geom_linerange(aes(xmin = LCL, xmax = UCL)) +
  #   geom_vline(xintercept = 0, color='blue') +
  #   xlab(expression(N[2]*O~flux~(metric~tons~year^{-1}))) + # 1Mg = 1 metric ton
  #   ylab("Lake size class (ha)") +
  #   theme_bw() 
  
  # ggpubr::ggarrange(b1, b2, # first two rows
  #                   ggarrange(b3, ncol = 2, nrow=1, labels = "AUTO"), # single plot in 3rd row, only one column
  #                   nrow=3, labels = "AUTO")
  
  # ggpubr::ggarrange(b2, b1, b3, ncol = 3, nrow=1, labels = "AUTO")
  ggpubr::ggarrange(b1, b2, b3, ncol = 1, nrow=3, labels = "AUTO", align = "v", hjust = -7, vjust = 2) #, hjust = -5, vjust = 2
  ggsave("manuscript/manuscript_figures/n2oFluxAndEmissionRateVsContinuousArea.tiff", width = 8, height = 5)   
  
} # fig 5: Flux and Flux by size


## SI Figure 1: N2O emission rate distribution----
if(!("SIfigure1.tiff" %in% list.files("manuscript/manuscript_figures"))) {
  
  # density plot by ecoregion
  all_predictions_ms %>%
    #filter(.draw %in% 1:10) %>% # subset for practice
    #mutate(ecoregion = factor(WSA9)) %>%
    #mutate(ecoregion = fct_reorder(ecoregion, preds_sat)) %>%
    ggplot(aes(x = e.n2o.mg.m2.d, group = .draw, color = .draw)) +
    geom_density(aes(y=after_stat(scaled)), show.legend = FALSE) +
    #geom_rug() + 
    xlab(expression("N"[2]*"O"~ "emission"~ "rate" ~ "(mg" ~ N[2]*O ~ m^-2 ~ d^-1 ~ ")")) + 
    ylab("density") +
    facet_wrap(~WSA9_NAME, scales = "free_x") +
    theme_bw() +
    theme(legend.position = "none") 
  
  
  ggsave("manuscript/manuscript_figures/SIfigure1.tiff", width = 8.5, height = 5)
  
} # SI figure


# MANUSCRIPT DATA-----

## paragraph 1----
#undersatN2oPercent <- round(((dg %>% dplyr::filter(n2o.src.snk == "sink") %>% 
#                         {nrow(.)} / dg %>% distinct(site.id) %>% 
#                         {nrow(.)}) * 100), 1)

#propN2Osummary <- all_predictions_ms %>% 
#  group_by(.draw) %>% 
#  summarise(propSink = sum(n2osat < 1) / length(.row)) %>% 
#  summarise(post_m_prop = round(mean(propSink), 3) * 100,
#            LCL = round(quantile(propSink, probs = 0.025), 3) * 100,
#            UCL = round(quantile(propSink, probs = 0.975), 3) * 100)

#save(propN2Osummary, file = "manuscript/manuscript_files/propN2Osummary.rda")

#propN2OsummaryWSA9 <- all_predictions_ms %>% 
#  group_by(WSA9, .draw) %>% 
#  summarise(propSink = sum(n2osat < 1) / length(.row)) %>% 
#  summarise(post_m_prop = round(mean(propSink), 3) * 100,
#            LCL = round(quantile(propSink, probs = 0.025), 3) * 100,
#            UCL = round(quantile(propSink, probs = 0.975), 3) * 100)

#save(propN2OsummaryWSA9, file = "manuscript/manuscript_files/propN2OsummaryWSA9.rda")

## Lake size distribution----
if(!("FluxBySize.rda" %in% list.files("manuscript/manuscript_files"))) {
  if(!("em_min_max.rda" %in% list.files("manuscript/manuscript_files"))) {
# Total flux by lake size class------
FluxBySize <- all_predictions_ms %>%
  group_by(size_cat, .draw) %>% # group by iteration
  summarise(sum_f.n2o.Mg.y = sum(f.n2o.Mg.y)) %>% # 2000 sums of flux by size cat
  # now summarize to 1 statistic per size_cat
  summarise( post_m_flux = round(mean(sum_f.n2o.Mg.y), 3), 
             LCL = round(quantile(sum_f.n2o.Mg.y, probs = 0.025), 3),
             UCL = round(quantile(sum_f.n2o.Mg.y, probs = 0.975), 3))

save(FluxBySize, file = "manuscript/manuscript_files/FluxBySize.rda")

# Min and max estimates for distribution of CONUS emission rates--
em_min_max <- all_predictions_ms %>%
  group_by(.draw) %>%
  summarise(min = min(e.n2o.mg.m2.d),
            max = max(e.n2o.mg.m2.d)) %>%
  summarise(min_LCL = round(quantile(min, probs = 0.025), 3),
            min_UCL = round(quantile(min, probs = 0.975), 3),
            max_LCL = round(quantile(max, probs = 0.025), 3),
            max_UCL = round(quantile(max, probs = 0.975), 3))

save(em_min_max, file = "manuscript/manuscript_files/em_min_max.rda")
  } # Flux by size
} # emissions min/max


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
# ipcc.n2o.indirect <- (15.2+2.9)*(1/265)*(1000/1)*(1000/1)*(1000/1)*(1/1000) #68,302 Mg N2O

# this is indirect emissions from groundwater (EF5-g), streams/rivers (EF5-r), and estuaries (EF5-e).  The default EF is
# 0.0025 for each, for an EF5 of 0.0075.  Therefore, only 1/3 of indirect N2O emissions reported above
# should be prescribed to streams/rivers/lakes/reservoirs.
# ipcc.n2o.surface = round(ipcc.n2o.indirect * (1/3), digits = 0) # 22,767

#ipcc.n2o <- tibble(indirect = (15.2+2.9)*(1/265)*(1000/1)*(1000/1)*(1000/1)*(1/1000),
#                   surface = round(ipcc.n2o.indirect * (1/3), digits = 0))

#save(ipcc.n2o, file = "manuscript/manuscript_files/ipcc.n2o.rda")
