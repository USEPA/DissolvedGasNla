# READ IN DISSOLVED GAS DATA
# This file contains the aggregated (by duplicates) dissolved gas data
# from NLA17.  See data dictionary in documents library for details.


dg <- read.table(file = paste0(localPath, 
                               "/Environmental Protection Agency (EPA)/",
                               "ORD NLA17 Dissolved Gas - Documents/",
                               "inputData/nla17gasDataAggregated_2019-11-12.txt"), 
                 header = TRUE, sep = "\t", as.is = TRUE)

# need to drop WSA9 = AK, WSA9_NAME = Alaska

# sp
# Define coordinates
coords <- data.frame(longitude = dg$map.lon.dd, latitude = dg$map.lat.dd)

# to spatial points df
dg.sp <- SpatialPointsDataFrame(coords, dg)
class(dg.sp)
summary(dg.sp)
names(dg.sp@data)
# projection
dg.sp@proj4string <- CRS('+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs')
dg.sp@bbox

dg.sf <- st_as_sf(dg.sp)
class(dg.sf)
# remove air data, 
# get error message as trying to sf-like function, filter, on an sp object
# convert dg.sp to sf
dg.sf <- filter(dg.sf, sample.source == "DG")
dim(dg.sf)
# sf
adg.sf <- st_as_sf(dg.sp, coords = c("map.lon.dd", "map.lat.dd"), crs = 4269)
