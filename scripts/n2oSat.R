n2oSat <- function(inputFile, waterTemp, airN2O, baro) {

# Input variables
# inputFile: Name of the data fram containing the information needed to calculate the dissolved gas concentrations. 
# waterTemp: Temperature of the waterbody when sampled [celsius]
# airN2O: Concentration of nitrous oxide in atmosphere [ppmv]
# baro: Barometric pressure at the time of equilibration [kPa]

##### Constants #####
cGas<-8.3144598 #universal gas constant (J K-1 mol-1)
cKelvin <- 273.15 #Conversion factor from Kelvin to Celsius
cPresConv <- 0.000001 # Constant to convert mixing ratio from umol/mol (ppmv) to mol/mol. Unit conversions from kPa to Pa, m^3 to L, cancel out.
cT0 <- 298.15#Henry's law constant T0

#Henry's law constants and temperature dependence from Sander (2015) DOI: 10.5194/acp-15-4399-2015
ckHN2O <- 0.00024 #mol m-3 Pa, range: 0.00018 - 0.00025
cdHdTN2O <- 2700 #K, range: 2600 - 3600


##### Calculate dissolved gas concentration at 100% saturation ##### 

# 100% saturation occurs when the dissolved gas concentration is in equilibrium
# with the atmosphere.
satN2O <- (ckHN2O * exp(cdHdTN2O*(1/(waterTemp + cKelvin) - 1/cT0))) * airN2O * baro * cPresConv

return(satN2O)
}

