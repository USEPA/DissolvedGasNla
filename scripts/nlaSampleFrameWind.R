# THIS SCRIPT IS TO CALCULATE WIND DATA FOR ALL ~600,000 WATERBODIES IN THE 
# NLA SAMPLE FRAME.  THESE DATA WILL BE USED TO CALCULATE A K VALUE FOR
# EACH SITE.  WHEN COMBINED WITH ROY'S MODEL PREDICTIONS OF DISSOLVED N2O 
# WILL ENABLE EMISSION RATE PREDICTION FOR EACH SITE.  THESE PREDICTIONS CAN
# BE USED TO ESTIMATE POPULATION LEVEL MEAN AND VARIANCE.  
# THIS CODE SHOULD BE RUN ON VM.  SCRIPT IS NOT YET COMPLETE.
# SEE 'readWind.R' IN THE NLA PROJECT FOR A VERSION TAILORED TO THE 1200
# NLA OBSERVATIONS.  THE NLA VERSION UTILIZES LARGE DATA FRAMES, BUT CAN
# BE RUN ON LAPTOP.  VERSION BELOW HAS A STRONGER RELIANCE ON LISTS, BUT
# THE DATAFRAME APPROACH MIGHT WORK ON VM.  NEED TO DECIDE HOW WE ARE GOING
# TO AGGREGATE DATA.  OPTIONS INCLUDE 1) INDEX PERIOD MEAN PER SITE, 2) 
# INDEX PERIOD MEAN OVER DAYLIGHT HOURS, AND 3) MEAN DAYLIGHT HOURS PER 
# DAY.  OPTION 3 WOULD ENABLE CALCULATING A K AND EMISSION RATE PER DAY
# FOR EACH SITE.


#R 4.0.3
library(tidyverse)
library(sf)
library(ncdf4)
library(ncdf4.helpers)
library(suncalc)
library(lubridate)
library(pbapply)

# ENFORCE EPA FORMAT ON NLA NAMES
toEPA <- function(X1){
  names(X1) = tolower(names(X1))
  names(X1) = gsub(pattern = "_", replacement = ".", x = names(X1))
  X1
}


# NLA17 sample frame--------------------------
sampFrm <- read.csv(file = paste0(Sys.getenv("USERPROFILE"), 
                                    "/Environmental Protection Agency (EPA)/",
                                    "ORD NLA17 Dissolved Gas - Documents/",
                                    "inputData/NLA_Sample_Frame.csv"), header = TRUE,
                    na.strings = "") %>%
  toEPA()

dim(sampFrm) # 614730, 71

sampFrm <- sampFrm %>%
  filter(nla17.sf != "Exclude2017") %>%
  filter(nla17.sf != "Exclude2017_Include2017NH") %>%
  filter(state != "DC") %>%
  filter(state != "HI")

dim(sampFrm) # 465,897, 71 

# netCDF file------------------------------------------
# awesome package vignette: 
# https://cran.r-project.org/web/packages/futureheatwaves/vignettes/starting_from_netcdf.html#:~:text=You%20can%20read%20netCDF%20data,connection%20to%20a%20netCDF%20file.
# ncdf file downloaded  from Climate Data Center contains hourly wind and lake temperature
# for US. (https://confluence.ecmwf.int/display/CKB/ERA5-Land?src=contextnavpagetreemode)
wind <- nc_open(paste0(Sys.getenv("USERPROFILE"), 
                       "/Environmental Protection Agency (EPA)/",
                       "ORD NLA17 Dissolved Gas - Documents/",
                       "inputData/10u10v24hour.nc")) #

wind # good!  two variables with three dimensions!

lon <- ncvar_get(wind, varid = "longitude")
lat <- ncvar_get(wind, varid = "latitude")
summary(lon);summary(lat)


# inspect time dimension
wind$dim$time
wind$dim$time$units # hours since 1900-01-01 00:00:00.0
wind$dim$time$calendar # gregorian.  Good, consistent with base R as.Date

wind.time <- as.POSIXct(wind$dim$time$vals * 3600, # convert hours to seconds, which are basis for POSIX class (60*60=3600)
                        origin = as.POSIXct("1900-01-01 00:00:00"), tz = "UTC") # UTC per CDS manual
wind.time # looks good
range(wind.time) # May 1 - Nov. 1, range of NLA sample dates

windU <- ncvar_get(wind, varid = "u10") # pull u10 data. in m/s

# This variable is in a 3 dimensional array ordered as longitude, then latitude, then time
dim(windU) # 701 461 4416
dim(lon) # 701
dim(lat) # 461
dim(wind.time) # 4416

{# demo
# # demo:  pull the modeled U10 at a certain location and time step.
# lon.index <- which.min(abs(lon - sampFrm$lon.dd83[1])) # minimal distance between grid and desired location
# lat.index <- which.min(abs(lat - sampFrm$lat.dd83[1])) # minimal distance between grid and desired location
# time.index <- which(wind.time == as.POSIXct("2017-07-13 12:00:00")) # hourly data available.
# 
# windU[lon.index, lat.index, time.index] # -0.9, cool!
# 
# # demo: get and plot full time series at grid point closes to lake:
# windU.ts <- nc.get.var.subset.by.axes(wind, "u10",
#                                       axis.indices = list(X = lon.index,
#                                                           Y = lat.index))
# dim(windU.ts) # 1 1 2208  1 lon, 1 lat, 2208 times
# 
# data_frame(time = wind.time,
#            u10 = as.vector(windU.ts)) %>%
#   ggplot(aes(time, u10)) +
#   geom_line()
}

# pull entire time series for each waterbody in sample frame
# create miniature to test with
#sampFrm.min <- sampFrm[sample(x = 1:nrow(sampFrm), size = 10), ] %>% as.data.frame() # sample() yields a weird object class

# calculate lon and lat indices for each row
sampFrm <- sampFrm %>% 
  rowwise() %>% # mutate by row!!!!
  mutate(lon.index = which.min(abs(lon - lon.dd83)), # index for minimum distance between target and grid longitude
         lat.index = which.min(abs(lat - lat.dd83)), # index for minimum distance between target and grid latitude
         lon.grid = lon[lon.index], # extract lon of closest grid, not needed for calcs, just to double check map matches observed
         lat.grid = lat[lat.index]) %>% # extract lat of closest grid, not needed for calcs, just to double check map matches observed
  ungroup()

# double check that extracted grids are spatially matched to lake targets
sampFrm %>%
  select(lon.dd83, lat.dd83, lon.grid, lat.grid) # yes, perfect match


# split dataframe into list
sampFrm.l <- split(sampFrm, seq(nrow(sampFrm)))

# extract wind speed time series for each waterbody and place in list element.
Sys.time()
sampFrm.l.w <- pblapply(sampFrm.l, function(x) {
  data.frame(u10 = nc.get.var.subset.by.axes(wind, "u10", # extract u10 component
                                             axis.indices = list(X = x$lon.index,
                                                                 Y = x$lat.index)) %>%
               as.vector(),
             v10 = nc.get.var.subset.by.axes(wind, "v10", # extract v10 component
                                             axis.indices = list(X = x$lon.index,
                                                                 Y = x$lat.index)) %>%
               as.vector(),
             time = wind.time, # add time stamp
             date = as.Date(wind.time), # collapse time stamp to date
             lat.dd83 = x$lat.dd83, # add latitude
             lon.dd83 = x$lon.dd83) %>% # add longitude
    mutate(ws = sqrt(u10^2 + v10^2)) %>% # calculate wind speed from north and east components (http://colaweb.gmu.edu/dev/clim301/lectures/wind/wind-uv)
    select(-v10, -u10)
})
Sys.time() # xx min on laptop

# calculate average wind speed per site
sampFrm.l.w.mean <- lapply(sampFrm.l.w, function(x){
  poo <- x %>%
    summarize(lat.unique = unique(lat.dd83), 
              lon.unique = unique(lon.dd83),
              ws = mean(ws))
  poo
})

# coerce to dataframe and write to disk
wind.sampFr <- do.call("rbind", sampleFrm.l.w.mean)









