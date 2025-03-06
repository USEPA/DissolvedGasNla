Data Munging
================
2025-03-05

- [1 Purpose](#1-purpose)
- [2 Data](#2-data)
  - [2.1 Read data](#21-read-data)
  - [2.2 Manipulate](#22-manipulate)
    - [2.2.1 Implement unit conversion and create source/sink
      column:](#221-implement-unit-conversion-and-create-sourcesink-column)
    - [2.2.2 Calculate emission rate.](#222-calculate-emission-rate)
  - [2.3 Write out data](#23-write-out-data)

# 1 Purpose

The purpose of this .rmd was to prepare the data set for subsequent
modeling. The data were written as .RData object that is stored in the
repository associated with this project:
<https://github.com/USEPA/DissolvedGasNla>

The dissolved gas data were originally curated in the ‘NLA’ dissolved
gas project in RStudio. The code can be found at a private github
repository (<https://github.com/USEPA/NLA>). After aggregating across
duplicate samples, the data were written to
nla17gasDataAggregated_2021-01-13.txt, which is also avaialble in the
GitHub repository for this project:
<https://github.com/USEPA/DissolvedGasNla>

The data object (‘dg.Rdata’) munged in this document was subsequently
used to estimate measurement uncertainty in the N2O observations, which
is documented in the ‘DG_sensitivity_to_measurement_error.Rmd’
(Saturation ratio/Standard GC sections) file that is also available in
the repository. From there, an updated ‘dg.RData’ object including this
new information on uncertainty was imported into additional RStudio
document workflows for (1) making design-based estimates from the survey
data (Design_Based_Estimates.Rmd) and (2) generating population
estimates using a Bayesian hierarchical model (NLA_N2O_models.Rmd). The
object is also imported for the files used to generate the manuscript
(manuscript_file.Rmd and manuscript_support.Rmd). All of these files
(and more) are included in the GitHub repository for this project in
order to ensure the results were reproducible starting from this base
file of the observed data.

# 2 Data

## 2.1 Read data

Below we read in the data files.

``` r
# Read dissolved gas data file
dg <- read.table(file = "./../inputData/nla17gasDataAggregated_2021-02-01.txt",
                 header = TRUE, sep = "\t", as.is = TRUE) %>%
  filter(!grepl("AK", site.id)) # omit Alaska sites

# Read simulated saturation ratio values (See DG_sensitivity_to_measurement_error.Rmd, Saturation ratio/Standard GC)
#S_a_me_abs <- readRDS("inputData/S_a_me_abs.RDS")
```

## 2.2 Manipulate

### 2.2.1 Implement unit conversion and create source/sink column:

``` r
dg <- dg %>%
  # unit conversion
  mutate(dissolved.ch4.umol = dissolved.ch4 * 1000000, # mol/L -> umol/L
         dissolved.co2.umol = dissolved.co2 * 1000000, # mol/L -> umol/L,
         dissolved.n2o.nmol = dissolved.n2o * 1000000000, # mol/L -> nmol/L
         sat.ch4.umol = sat.ch4 * 1000000, # mol/L -> umol/L
         sat.co2.umol = sat.co2 * 1000000, # mol/L -> umol/L
         sat.n2o.nmol =  sat.n2o * 1000000000, # mol/L -> nmol/L
         # add source/sink column
         co2.src.snk = ifelse(co2.sat.ratio > 1, "source", "sink"),
         ch4.src.snk = ifelse(ch4.sat.ratio > 1, "source", "sink"),
         n2o.src.snk = ifelse(n2o.sat.ratio > 1, "source", "sink")
         ) %>%
  # remove fields no longer needed
  select(-dissolved.ch4, -dissolved.co2, -dissolved.n2o,
         -sat.ch4, -sat.co2, - sat.n2o)
```

### 2.2.2 Calculate emission rate.

The emission rate (E<sub>gas</sub>) is the rate at which a gas crosses
the air water interface and is expressed on an areal basis (i.e. mass
m<sup>-2</sup> day<sup>-1</sup>). The emission rate can be converted to
a flux (F<sub>gas</sub>; mass per unit time; mass d<sup>-1</sup>). In
the context of emissions from waterbodies, the flux (F<sub>gas</sub>) is
calculated as the product of the emission rate (E<sub>gas</sub>) and
waterbody area.

E<sub>gas</sub> can be calculated as the product of gas<sup>\*</sup> and
the gas transfer velocity (k):

E<sub>gas</sub> = gas<sup>\*</sup> \* k

where gas<sup>\*</sup> is the difference between the observed
(gas<sub>obs</sub>) and equilibrium (gas<sub>eq</sub>) dissolved gas
concentration:

gas^\* = (gas<sub>obs</sub>) - (gas<sub>eq</sub>)

where a positive value indicates the waterbody is a source of the gas
and a negative value indicates the waterbody is a sink.

The gas transfer velocity (k; cm h<sup>-1</sup>) is a measure of the
physical interaction between the waterbody and the atmosphere. Wavy and
turbulent waterbodies have high gas transfer velocities whereas calm
waterbodies have low values. k was not measured in the 2017 National
Lakes Assessment, but was estimated from wind speed and lake area
(Vachon, D., and Y. T. Prairie (2013), The ecosystem size and shape
dependence of gas transfer velocity versus wind speed relationships in
lakes, Can. J. Fish. Aquat. Sci., 70(12), 1757-1764,
<doi:10.1139/cjfas-2013-0241>):

k<sub>600</sub> = 2.51 + 1.48 \* U<sub>10</sub> + 0.39 \* U<sub>10</sub>
\* log<sub>10</sub>Lake area

where k<sub>600</sub> is the gas transfer velocity normalized to
CO<sub>2</sub> at 20 <sup>o</sup>C, U10 is wind speed 10m above the
water surface. Wind speed was obtained on a 7.5km grid from the ERA-5
Land database and overlaid on the NLA sampling points. U<sub>10</sub>
was averaged between sunrise and sunset for each site, assuming that
daytime wind conditions are best matched with the dissolved gas
measurements which were conducted during the day.

k<sub>600</sub> is a standardized value that must be corrected for the
differences in diffusivity among gases and water temperature at the
sampling sites. These corrections will be implemented at a later date
and emission rates will be calculated using k<sub>600</sub> for now.

``` r
dg <- dg %>% 
  mutate(
    # e.ch4 umol CH4/m2/d 
    e.ch4.umol.d = ((dissolved.ch4.umol - sat.ch4.umol) * 1000) * # 1000 L to m3
           (k600.day * (24/100)), # 24 hour to day.  100 cm to m
    # e.co2 umol CO2/m2/d
    e.co2.umol.d = ((dissolved.co2.umol - sat.co2.umol) * 1000) * # 1000 L to m3
           (k600.day * (24/100)), # 24 hour to day.  100 cm to m
        # e.n2o nmol N2O/m2/d
    e.n2o.nmol.d = ((dissolved.n2o.nmol - sat.n2o.nmol) * 1000) * # 1000 L to m3
           (k600.day * (24/100)), # 24 hour to day.  100 cm to m
    # total flux of CH4 per day kmol CH4 day-1
    f.ch4.km.d = e.ch4.umol.d * area.ha * 10000 * (1/10^9), # 1ha = 10,000 m2. 10^9 umol->mmol->mol->kmol
    # total flux of CO2 per day kmol CO2 day-1
    f.co2.km.d = e.co2.umol.d * area.ha * 10000 * (1/10^9), # 1ha = 10,000 m2. 10^9 umol->mmol->mol->kmol
    # total flux of N2O per day mol N2O day-1
    f.n2o.m.d = e.n2o.nmol.d * area.ha * 10000 * (1/10^9))  # 1ha = 10,000 m2. 10^9 nmol->umol->mmol->mol
```

## 2.3 Write out data

Write out final .RData object.

``` r
save(dg, file = "./../inputData/dg.RData")
```
