# READ IN DISSOLVED GAS DATA
# This file contains the aggregated (by duplicates) dissolved gas data
# from NLA17.  See data dictionary in documents library for details.


dg <- read.table(file = paste0(localPath, 
                               "/Environmental Protection Agency (EPA)/",
                               "ORD NLA17 Dissolved Gas - Documents/",
                               "inputData/nla17gasDataAggregated_2019-11-12.txt"), 
                 header = TRUE, sep = "\t", as.is = TRUE)


# sp
# Define coordinates
coords <- data.frame(longitude = dg$map.lon.dd, latitude = dg$map.lat.dd)

# to spatial points df
dg.sp <- SpatialPointsDataFrame(coords, dg)

# projection
dg.sp@proj4string <- CRS('+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs')

# remove air data
dg.sp <- filter(dg.sp, sample.source == "DG")

# sf
adg.sf <- st_as_sf(dg.sp, coords = c("map.lon.dd", "map.lat.dd"), crs = 4269)
