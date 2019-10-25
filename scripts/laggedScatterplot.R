# Exploratory
bubble(dg.sf['dissolved.ch4']) # fails with sf

bubble(dg.sp['dissolved.ch4']) # fails due to missing values

bubble(dg.sp[!is.na(dg.sp$dissolved.ch4), "dissolved.ch4"]) # works

# lagged scatterplots
# What are distance between points
gdis <- spDists(dg.sp)
max(gdis) # 5617, this must be km.  5617 km = 3490 miles, AK to FL

hscat(dissolved.ch4 ~ 1,  
      dg.sp[!is.na(dg.sp$dissolved.ch4), ], # remove missing
      # 5617275m = 3490 miles. AK to FL
      c(0, 1000, 2000, 3000, 4000, 5000, 6000)) 

# Do this by ecoregion?