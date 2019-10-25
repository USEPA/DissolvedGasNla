# READ IN DISSOLVED GAS DATA
# This file contains the aggregated (by duplicates) dissolved gas data
# from NLA17.  See data dictionary in documents library for details.

dg <- read.table(file = paste0("C:/Users/JBEAULIE/Environmental Protection Agency (EPA)/",
                               "ORD NLA17 Dissolved Gas - Documents/",
                               "inputData/nla17gasDataAggregated_2019-10-25.txt"), 
                 header = TRUE, sep = "\t", as.is = TRUE)

# sp
# Define coordinates
coords <- data.frame(longitude = dg$map.lon.dd, latitude = dg$map.lat.dd)
# spatial points df
dg.sp <- SpatialPointsDataFrame(coords, dg)
# projection
dg.sp@proj4string <- CRS('+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs')



# sf
dg.sf <- st_as_sf(dg, coords = c("map.lon.dd", "map.lat.dd"), crs = 4269)
