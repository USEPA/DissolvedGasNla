Modeling workflow: Lake N2O survey data
================
Roy Martin, Jake Beaulieu, Michael McManus
2023-03-17

# 1 Background and Objectives

This document details the modeling workflow implemented for estimating
dissolved and equilibrium N2O gas concentrations and saturation ratios
for all freshwater lakes and reservoirs in the lower 48 US states larger
than 4 ha. Data for these estimates came from the 2017 Nation Lakes
Assessment (NLA) survey, wherein waterbodies were sampled according to a
spatially balanced, stratified, and unequal probability design.
Stratification was among categories of lake size (surface area in
hectares), WSA9 ecoregion, and US state (excluding AK and HI); and
larger lakes were intentionally over-sampled relative to smaller lakes.
In order for this sample to be useful for making inferences relevant to
the population of interest, any estimates derived from it needed to be
adjusted for potential biases arising from these design factors. It can
be useful to consider the complete data likelihood when imagining how
such biases may arise in sampling designs([Link and Barker
2010](#ref-Link_Barker_2010); [Zachmann et al.
2022](#ref-Zachmann_etal_2022)). For the 2017 NLA survey, the complete
data was considered to be all US lakes and reservoirs in the lower 48
states larger than 4 hectares. The survey samples, by comparison, are
then considered a subset of the complete data with a known pattern of
“missingess” assigned by the survey. For example, some lakes in the
target population were missing from the sample not at random, but
conditional on the design parameters. Their “missingness” is not random
and, therefore, not ignorable when making inferences from the sample to
the population of interest ([Gelman et al. 2014](#ref-Gelman_etal_2014),
Ch. 8; [Zachmann et al. 2022](#ref-Zachmann_etal_2022)).

One model-based strategy to account for this selection bias is to
include the survey design variables as predictors in a regression model
and then adjust the resulting estimates based on a known target
population distribution; a post-fitting adjustment referred to as
poststratification ([Gelman, Hill, and Vehtari
2020](#ref-Gelman_etal_2020), Ch. 17). A particularly popular version of
this approach uses multilevel regression models and is aptly referred to
as multilevel regression and poststratification, or MRP ([Park, Gelman,
and Bafumi 2004](#ref-Park_etal_2004); [Gelman and Little
1997](#ref-Gelman_Little_1997); [Gelman, Hill, and Vehtari
2020](#ref-Gelman_etal_2020), Ch. 17; [Kennedy and Gelman
2021](#ref-Kennedy_Gelman_2021)). Multilevel models are often
recommended because they can provide regularized estimates along the
design groupings, which can improve out-of-sample inferences (e.g.,
[Gelman and Little 1997](#ref-Gelman_Little_1997); [Kennedy and Gelman
2021](#ref-Kennedy_Gelman_2021)). Estimates for group levels that may be
missing from the sample, but are part of the population of interest, are
also straightforward using the multilevel approach ([Gelman, Hill, and
Vehtari 2020](#ref-Gelman_etal_2020) Ch. 17; [McElreath
2020](#ref-McElreath_2020)).

With a stratified, unequal probability design, such as the one used for
the 2017 NLA, the selection design prescribes a number of random draws
for each of a number of pre-specified groupings (i.e., sampling unit
types). Within a particular stratification level or grouping, the draws
may be considered a simple random sample (SRS) from that group type in
the target population. A regression model conditioning the sample
observations on the design grouping structure could, therefore, be
parameterized to provide group-specific estimates that are unbiased with
respect to the target population. To make inferences to the broader
target population, however, those regression estimates would need to be
re-weighted by the known proportions of units across groups in the
target population. If the estimates are not adjusted in this manner, the
estimates may be biased, particularly if the response of interest varies
importantly across the grouping structure. For example, for the 2017 NLA
survey, the proportion of samples prescribed to small lakes was much
lower, relative to their proportion in the target population, compared
to larger lakes, which were sampled in greater proportion relative to
the target population. Small lakes are so numerous in the US that, if
this unequal weighting by lake size weren’t prescribed, nearly all the
samples from the survey would be from small lakes and inferences
regarding larger lakes could be highly uncertain for any reasonably
achievable number of samples. On the other hand, if lake size had no
real influence on the N2O responses, the stratification across size
class would make little difference with regard to inferences. Sampling
any sizes of lakes would result in the same marginal estimates and
standard errors.

For a recent applied example and helpful tutorial employing these
concepts in a model-based approach using MRP, see Kennedy and Gelman
([2021](#ref-Kennedy_Gelman_2021)). For an example of a recent
application in the context of national surveys of environmental
resources, see Zachmann et al. ([2022](#ref-Zachmann_etal_2022)). In
this study, a similar model-based approach was used to make
population-level estimates from the 2017 NLA dissolved gas data. The
specific workflow, data, models, and code used to make these estimates
are documented in the remainder of this document. First, Bayesian
multilevel regression models were used to fit the sample data to the
design variable structure. Next, because eventual flux estimates needed
to be estimated from lake-level measures (i.e., surface area), instead
of predicting to a typical postratification table, predictions were made
to each individual lake in the population of interest. This meant
predicting to 465,897 individual natural and man made US lakes larger
than 4 hectares in the lower 48 states. These predictions were assumed
relevant to average conditions during the biological index period for
each lake in 2017. The specific objective of the modeling effort was to
provide population estimates for (1) dissolved and equilibrium N2O
concentrations; (2) the N2O saturation ratio (i.e., dissolved
N2O/equilibrium N2O); and (3) the proportion of under-saturated water
bodies (i.e., saturation ratio \< 1). As previously indicated, the gas
estimates would also be used to estimate the total flux of N2O gas
attributable to the target population of lakes over the index period.
Saturation ratio was calculated as the ratio of dissolved to equilibrium
N2O. Because dissolved and equilibrium N2O were observed on the same
sample units (lake sites), models were fit to their joint distribution
and the ratio estimates were assembled as a derived quantity. The
response variable was, therefore, multivariate in order to account for
any potential statistical dependencies between dissolved and equilibrium
N2O due to, for example, common dependencies on geography. Although
predictions of the mean marginal probabilities from separate models may
have provided comparable estimates, a joint model allowing correlated
errors was expected to better capture uncertainty and potentially
improve out-of-sample predictions, should the variables be conditionally
correlated ([Warton et al. 2015](#ref-Warton_etal_2015); [Poggiato et
al. 2021](#ref-Poggiato_etal_2021)). All of the models fit were
constructed using the `brms` package ([Bürkner 2017](#ref-Burkner_2017))
in `R` ([R Core Team 2021](#ref-R_Core_Team_2021)) as an interface to
Stan, a software package for fitting fully Bayesian models via
Hamiltonian Monte Carlo (HMC, [Stan Development Team
2018b](#ref-Stan_Development_Team_2018_a),
[2018c](#ref-Stan_Development_Team_2018_b),
[2018a](#ref-Stan_Development_Team_2018_c)). More specific details on
the model structures is provided throughout the “Model fitting” section
below.

# 2 Data

As explained in a previous data munging document document
(<https://github.com/USEPA/DissolvedGasNla/blob/master/scripts/dgIndicatorAnalysis.html>),
duplicate dissolved gas samples were collected at a depth of \~0.1m at
designated index sites distributed across 1091 lakes nationwide, of
which 95 were sampled twice as repeat visits. This subset of revisit
sites was used as a test set for assessing model fit and out-of-sample
performance.

Gas samples were analyzed via gas chromotography and concentrations were
recorded to the nearest 0.001 nmol/L. As part of the survey design, each
gas observation was indexed to an individual lake selected with unequal
probability from 5 different lake size categories,
![j \\in j=1,...,J = 5](https://latex.codecogs.com/svg.image?j%20%5Cin%20j%3D1%2C...%2CJ%20%3D%205 "j \in j=1,...,J = 5"),
according to surface area (ha), and from within a state,
![k \\in k=1,...,K = 48](https://latex.codecogs.com/svg.image?k%20%5Cin%20k%3D1%2C...%2CK%20%3D%2048 "k \in k=1,...,K = 48"),
situated within an aggregated, WSA9 or Omernik ecoregion,
![l \\in l=1,...,L = 9](https://latex.codecogs.com/svg.image?l%20%5Cin%20l%3D1%2C...%2CL%20%3D%209 "l \in l=1,...,L = 9").
All 9 WSA9 ecoregions were represented in the sample, including Xeric
(XER), Western Mountain (WMT), Northern Plains (NPL), Southern Plains
(SPL), Temperate Plains (TPL), Coastal Plains (CPL), Upper Midwest
(UMW), Northern Appalachian (NAP), and Southern Appalachian (SAP)
regions. As shown below, the data from the initial and revisit samples
were separately compiled into data frame objects in
![\\textbf{R}](https://latex.codecogs.com/svg.image?%5Ctextbf%7BR%7D "\textbf{R}"),
with ![n=984](https://latex.codecogs.com/svg.image?n%3D984 "n=984") and
![n=95](https://latex.codecogs.com/svg.image?n%3D95 "n=95") rows,
respectively, of gas observations indexed to the survey design variables
and several potentially relevant covariates.

## 2.1 Import

The gas data and covariates were previously described and munged at
<https://github.com/USEPA/DissolvedGasNla/blob/master/scripts/dataMunge.html>.
That dataset was imported below.

``` r
load( file = paste0( localPath,
              "/Environmental Protection Agency (EPA)/",
              "ORD NLA17 Dissolved Gas - Documents/",
              "inputData/dg.2021-02-01.RData")
      )

save(dg, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/dg.rda") 
```

A new data frame for modeling was constructed from the original file
including only the variables of interest: (1) the N2O gas observations;
(2) the survey design variables indexed to those observations; and (3)
additional covariates considered potentially useful for improving the
fit of the model. The modeling data frame below excluded the
second-visit observations, which would later be used for model checking.
Some variables from the imported data were renamed for convenience. In
addition, the NO3 covariate was rounded according to the documented
measurement precision and an alternative version was also created by
log-transforming and re-coding the variable as an ordered factor with
five levels. The left-most cut point separated observations below the
detection limit from the completely observed samples. The remaining cut
points were drawn at approximately equal distances in the positive
direction along the log scale. Finally, one lake had no N2O gas
information and was removed from the data frame.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/dg.rda")

dg %>%
  filter(!is.na(dissolved.n2o.nmol)) %>% # 1 obs with missing measurement
  nrow() # number of observations before filtering
```

    ## [1] 1185

``` r
df_model <- dg %>%
  filter(!is.na(dissolved.n2o.nmol)) %>%
  filter(sitetype == "PROB") %>% # probability samples only
  filter(visit.no == 1) %>%
  mutate(n2o = round(dissolved.n2o.nmol, 2),
         n2o_eq = round(sat.n2o.nmol, 2),
         n2o_sat = n2o.sat.ratio,
         n2o_em = e.n2o.nmol.d,
         n2o_flux = f.n2o.m.d,
         WSA9 = factor(ag.eco9),
         state = factor(state.abb[match(state.nm, state.name)]),
         area_ha = area.ha,
         log_area = log(area_ha),
         chla = chla.result,
         log_chla = log(chla),
         elev = elevation,
         log_elev = log(elev + 1),
         do_surf = o2.surf,
         log_do = log(do_surf),
         bf_max = max.bf,
         sqrt_bf = sqrt(bf_max),
         size_cat = recode(area.cat6, 
                           "(1,4]" = "min_4" ,
                           "(10,20]" = "10_20",
                           "(20,50]" = "20_50",
                           "(4,10]" = "4_10",
                           ">50" = "50_max")) %>%
  mutate(size_cat = factor(size_cat,
                           levels = c("min_4", "4_10", "10_20", "20_50", "50_max"),
                           ordered = TRUE)) %>%
  mutate(no3 = ifelse(nitrate.n.result <= 0.0005, 0.0005, round(nitrate.n.result, 4))) %>%# 1/2 mdl 0.01
  mutate(no3_cat = cut(log(no3), # convert no3 to ordered factor with 5 levels
                       breaks = c(-Inf, -7.5, -5.5, -3.5, -1.5, Inf),
                       labels =seq(1, 5, 1))) %>%
  mutate(no3_cat = factor(no3_cat,
                          levels = seq(1, 5, 1),
                          ordered = TRUE)) %>%
  mutate(date = as.Date(date.col)) %>%
  mutate(jdate = as.numeric(format(date, "%j"))) %>% 
  mutate(lat = map.lat.dd,
         lon = map.lon.dd) %>% # longitude
  mutate(surftemp = surftemp,
         log_surftemp = log(surftemp)) %>% 
  select(WSA9,
         state,
         size_cat,
         site.id,
         lat,
         lon,
         date,
         jdate,
         surftemp,
         log_surftemp,
         area_ha,
         log_area,
         elev,
         log_elev,
         chla,
         log_chla,
         do_surf,
         log_do,
         bf_max,
         sqrt_bf,
         n2o,
         n2o_eq,
         no3,
         no3_cat
         )

save(df_model, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/df_model.rda") 

nrow(df_model) # number of obs after filtering
```

    ## [1] 984

``` r
df_model %>%
  head(10)
```

    ## # A tibble: 10 x 24
    ##    WSA9  state size_cat site.id       lat   lon date       jdate surft~1 log_s~2
    ##    <fct> <fct> <ord>    <chr>       <dbl> <dbl> <date>     <dbl>   <dbl>   <dbl>
    ##  1 SAP   AL    50_max   NLA17_AL-1~  33.3 -87.4 2017-07-10   191    23.9    3.17
    ##  2 SAP   AL    50_max   NLA17_AL-1~  33.7 -86.2 2017-08-14   226    28.1    3.34
    ##  3 CPL   AL    50_max   NLA17_AL-1~  32.5 -87.8 2017-07-12   193    27.2    3.30
    ##  4 CPL   AL    20_50    NLA17_AL-1~  31.6 -88.4 2017-07-14   195    29      3.37
    ##  5 CPL   AL    4_10     NLA17_AL-1~  33.4 -88.2 2017-07-11   192    28      3.33
    ##  6 CPL   AL    10_20    NLA17_AL-1~  32.2 -87.8 2017-07-13   194    29      3.37
    ##  7 SAP   AL    10_20    NLA17_AL-1~  33.2 -87.2 2017-08-15   227    27.1    3.30
    ##  8 CPL   AL    min_4    NLA17_AL-1~  31.2 -85.6 2017-08-17   229    29.5    3.38
    ##  9 CPL   AR    50_max   NLA17_AR-1~  34.5 -92.3 2017-06-07   158    28      3.33
    ## 10 SAP   AR    50_max   NLA17_AR-1~  34.4 -93.1 2017-06-12   163    25.2    3.23
    ## # ... with 14 more variables: area_ha <dbl>, log_area <dbl>, elev <dbl>,
    ## #   log_elev <dbl>, chla <dbl>, log_chla <dbl>, do_surf <dbl>, log_do <dbl>,
    ## #   bf_max <dbl>, sqrt_bf <dbl>, n2o <dbl>, n2o_eq <dbl>, no3 <dbl>,
    ## #   no3_cat <ord>, and abbreviated variable names 1: surftemp, 2: log_surftemp

A second data frame including only the second visit observations was
constructed below. These data were later used as a “test set” to assess
the out-of-sample performance of the model developed on the first-visit
or “training set”.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/dg.rda")

# number of observations before filtering probability samples
dg %>%
  filter(!is.na(dissolved.n2o.nmol)) %>% # remove obs with missing response measurements
  nrow()
```

    ## [1] 1185

``` r
df_test <- dg %>%
  filter(!is.na(dissolved.n2o.nmol)) %>%
  filter(sitetype == "PROB") %>% # probability samples only
  filter(visit.no == 2) %>%
  mutate(n2o = round(dissolved.n2o.nmol, 2),
         n2o_eq = round(sat.n2o.nmol, 2),
         n2o_sat = n2o.sat.ratio,
         n2o_em = e.n2o.nmol.d,
         n2o_flux = f.n2o.m.d,
         WSA9 = factor(ag.eco9),
         state = factor(state.abb[match(state.nm, state.name)]),
         area_ha = area.ha,
         log_area = log(area_ha),
         chla = chla.result,
         log_chla = log(chla),
         elev = elevation,
         log_elev = log(elev + 1),
         do_surf = o2.surf,
         log_do = log(do_surf),
         bf_max = max.bf,
         sqrt_bf = sqrt(bf_max),
         size_cat = recode(area.cat6, 
                           "(1,4]" = "min_4" ,
                           "(10,20]" = "10_20",
                           "(20,50]" = "20_50",
                           "(4,10]" = "4_10",
                           ">50" = "50_max")) %>%
  mutate(size_cat = factor(size_cat,
                           levels = c("min_4", "4_10", "10_20", "20_50", "50_max"),
                           ordered = TRUE)) %>%
  mutate(no3 = ifelse(nitrate.n.result <= 0.0005, 0.0005, round(nitrate.n.result, 4))) %>%# 1/2 mdl 0.01
  mutate(no3_cat = cut(log(no3), # convert no3 to ordered factor with 5 levels
                       breaks = c(-Inf, -7.5, -5.5, -3.5, -1.5, Inf),
                       labels =seq(1, 5, 1))) %>%
  mutate(no3_cat = factor(no3_cat,
                          levels = seq(1, 5, 1),
                          ordered = TRUE)) %>%
  mutate(date = as.Date(date.col)) %>%
  mutate(jdate = as.numeric(format(date, "%j"))) %>% 
  mutate(lat = map.lat.dd,
         lon = map.lon.dd) %>% # longitude
  mutate(surftemp = surftemp,
         log_surftemp = log(surftemp)) %>% 
  select(WSA9,
         state,
         size_cat,
         site.id,
         lat,
         lon,
         date,
         jdate,
         surftemp,
         log_surftemp,
         area_ha,
         log_area,
         elev,
         log_elev,
         chla,
         log_chla,
         do_surf,
         log_do,
         bf_max,
         sqrt_bf,
         n2o,
         n2o_eq,
         no3,
         no3_cat
         )

save(df_test, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/df_test.rda") 

nrow(df_test) # number of obs after filtering for probability samples, first visits, and removing one site missing ecoregion (WSA9) info.
```

    ## [1] 95

``` r
df_test %>%
  head(10)
```

    ## # A tibble: 10 x 24
    ##    WSA9  state size_cat site.id      lat    lon date       jdate surft~1 log_s~2
    ##    <fct> <fct> <ord>    <chr>      <dbl>  <dbl> <date>     <dbl>   <dbl>   <dbl>
    ##  1 SAP   AL    50_max   NLA17_AL-~  33.3  -87.4 2017-08-14   226    28.2    3.34
    ##  2 CPL   AL    20_50    NLA17_AL-~  31.6  -88.4 2017-08-16   228    29      3.37
    ##  3 CPL   AR    50_max   NLA17_AR-~  34.5  -92.3 2017-07-20   201    28.8    3.36
    ##  4 CPL   AR    4_10     NLA17_AR-~  33.1  -92.7 2017-09-06   249    24.1    3.18
    ##  5 XER   AZ    50_max   NLA17_AZ-~  33.6 -112.  2017-08-30   242    28      3.33
    ##  6 WMT   AZ    20_50    NLA17_AZ-~  33.8 -109.  2017-08-08   220    20.2    3.01
    ##  7 XER   CA    50_max   NLA17_CA-~  38.1 -123.  2017-08-07   219    21.3    3.06
    ##  8 XER   CA    50_max   NLA17_CA-~  33.8 -118.  2017-08-23   235    26.7    3.28
    ##  9 WMT   CO    4_10     NLA17_CO-~  39.0 -108.  2017-07-31   212    17.6    2.87
    ## 10 SPL   CO    50_max   NLA17_CO-~  40.3 -105.  2017-07-25   206    23.6    3.16
    ## # ... with 14 more variables: area_ha <dbl>, log_area <dbl>, elev <dbl>,
    ## #   log_elev <dbl>, chla <dbl>, log_chla <dbl>, do_surf <dbl>, log_do <dbl>,
    ## #   bf_max <dbl>, sqrt_bf <dbl>, n2o <dbl>, n2o_eq <dbl>, no3 <dbl>,
    ## #   no3_cat <ord>, and abbreviated variable names 1: surftemp, 2: log_surftemp

## 2.2 Target population

Below. the NLA sampling frame was imported and filtered to include all
lakes the target population. The resulting target population above
included a total of 465,897 waterbodies.

``` r
df_pop <- read.csv(file = paste0(localPath,
              "/Environmental Protection Agency (EPA)/",
              "ORD NLA17 Dissolved Gas - Documents/",
              "inputData/NLA_Sample_Frame.csv"), header = T)

sframe <- df_pop %>%
  filter(nla17_sf != "Exclude2017") %>%
  filter(nla17_sf != "Exclude2017_Include2017NH") %>%
  filter(state != "DC") %>%
  filter(state != "HI") %>%
  droplevels() %>%
  mutate(WSA9 = factor(ag_eco9),
         WSA9 = forcats::fct_drop(WSA9), # remove NA level
         state = factor(state),
         size_cat = factor(area_cat6),
         lat = lat_dd83,
         lon = lon_dd83,
         log_area = log(area_ha),
         elev = elevation,
         log_elev = ifelse(elev <= 0, 0, elev), # assumed elev < 0 to be elev = 0
         log_elev = log(log_elev + 1)
         ) %>% 
  mutate(size_cat = recode(size_cat, 
                           "(1,4]" = "min_4" ,
                           "(10,20]" = "10_20",
                           "(20,50]" = "20_50",
                           "(4,10]" = "4_10",
                           ">50" = "50_max")) %>%
  mutate(size_cat = factor(size_cat, 
                           levels = c("min_4", "4_10", "10_20", "20_50", "50_max"),
                           ordered = TRUE)) %>%
  select(WSA9, state, size_cat, lat, lon, area_ha, log_area, elev, log_elev)

rm(df_pop)

save(sframe, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/sframe.rda") 

sframe %>%
  head(10)
```

    ## # A tibble: 10 x 9
    ##    WSA9  state size_cat   lat   lon area_ha log_area  elev log_elev
    ##    <fct> <fct> <ord>    <dbl> <dbl>   <dbl>    <dbl> <int>    <dbl>
    ##  1 NAP   MA    50_max    42.4 -72.3   9544.     9.16   159     5.08
    ##  2 NAP   VT    50_max    42.8 -72.9    812.     6.70   454     6.12
    ##  3 NAP   ME    50_max    44.9 -69.2   1888.     7.54    60     4.11
    ##  4 NAP   NY    50_max    44.3 -74.1    525.     6.26   465     6.14
    ##  5 NAP   NY    50_max    43.2 -74.1   9493.     9.16   234     5.46
    ##  6 CPL   VA    50_max    37.4 -76.8    361.     5.89    17     2.89
    ##  7 SAP   VA    50_max    37.3 -77.7   1317.     7.18    45     3.83
    ##  8 CPL   VA    50_max    36.8 -76.6    339.     5.83     8     2.20
    ##  9 SAP   NC    50_max    36.5 -78.9   1102.     7.01   131     4.88
    ## 10 CPL   NC    50_max    35.9 -77.9    287.     5.66    36     3.61

Cross tabulations below describe the structure of the target population
with respect to the survey design variables. The cross-tabulation shows
that each ecoregion does not contain each state. Therefore, in the
statistical sense, states were nested in ecoregions.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/sframe.rda")

sframe %>%
  group_by(WSA9, state) %>%
  summarise(n = n(), .groups = "drop") %>%
  spread(state, n) %>%
  head(10)
```

    ## # A tibble: 9 x 49
    ##   WSA9     AL    AR    AZ    CA    CO    CT    DE    FL    GA    IA    ID    IL
    ##   <fct> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int>
    ## 1 CPL    7326  6395    NA    NA    NA    NA   529 37888 23761    NA    NA    44
    ## 2 NAP      NA    NA    NA    NA    NA  2143    NA    NA    NA    NA    NA    NA
    ## 3 NPL      NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## 4 SAP    3877  2362    NA    NA    NA    NA    10    NA  9014    NA    NA   275
    ## 5 SPL      NA    NA    NA    NA  2006    NA    NA    NA    NA    NA    NA    NA
    ## 6 TPL      NA    NA    NA    NA    NA    NA    NA    NA    NA  5015    NA 10308
    ## 7 UMW      NA    NA    NA    NA    NA    NA    NA    NA    NA   205    NA    68
    ## 8 WMT      NA    NA   332  4261  2958    NA    NA    NA    NA    NA  1998    NA
    ## 9 XER      NA    NA   578  5043   606    NA    NA    NA    NA    NA   923    NA
    ## # ... with 36 more variables: IN <int>, KS <int>, KY <int>, LA <int>, MA <int>,
    ## #   MD <int>, ME <int>, MI <int>, MN <int>, MO <int>, MS <int>, MT <int>,
    ## #   NC <int>, ND <int>, NE <int>, NH <int>, NJ <int>, NM <int>, NV <int>,
    ## #   NY <int>, OH <int>, OK <int>, OR <int>, PA <int>, RI <int>, SC <int>,
    ## #   SD <int>, TN <int>, TX <int>, UT <int>, VA <int>, VT <int>, WA <int>,
    ## #   WI <int>, WV <int>, WY <int>

The cross-tablulation below indicates that lake size category was nested
in state (which was nested in ecoregion). That is, not every
ecoregion:state contained every size category.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/sframe.rda")

sframe %>%
  group_by(WSA9, state, size_cat) %>%
  summarise(n = n(), .groups = "drop") %>%
  spread(size_cat, n) %>%
  head(10)
```

    ## # A tibble: 10 x 7
    ##    WSA9  state min_4 `4_10` `10_20` `20_50` `50_max`
    ##    <fct> <fct> <int>  <int>   <int>   <int>    <int>
    ##  1 CPL   AL     5812   1078     253     119       64
    ##  2 CPL   AR     4178   1246     473     314      184
    ##  3 CPL   DE      380     72      42      28        7
    ##  4 CPL   FL    27613   6138    1918    1294      925
    ##  5 CPL   GA    20135   2814     481     238       93
    ##  6 CPL   IL       25     14      NA       2        3
    ##  7 CPL   KY      398     45      20       9       NA
    ##  8 CPL   LA    14429   3670    1163     699      636
    ##  9 CPL   MA      422    155      87      44       44
    ## 10 CPL   MD     1111    234      75      41       11

Below, the sampling frame was munged to create a post-stratification
table. There were 536 “types” or groupings of lakes in the population of
interest with respect to the sampling design. The total counts of those
lake types (n_lakes) and their proportions (prop_cell) relative to the
counts in the target population were tabulated.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/sframe.rda")

pframe <- sframe %>%
  mutate(obs = 1) %>%
  group_by(WSA9, state, size_cat) %>%
  summarise(n_lakes = sum(obs), .groups = "drop") %>%
  ungroup() %>%
  mutate(prop_cell = n_lakes/sum(n_lakes)) %>%
  mutate(type = "population") 

save(pframe, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/pframe.rda")

pframe %>%
  head(10)
```

    ## # A tibble: 10 x 6
    ##    WSA9  state size_cat n_lakes prop_cell type      
    ##    <fct> <fct> <ord>      <dbl>     <dbl> <chr>     
    ##  1 CPL   AL    min_4       5812  0.0125   population
    ##  2 CPL   AL    4_10        1078  0.00231  population
    ##  3 CPL   AL    10_20        253  0.000543 population
    ##  4 CPL   AL    20_50        119  0.000255 population
    ##  5 CPL   AL    50_max        64  0.000137 population
    ##  6 CPL   AR    min_4       4178  0.00897  population
    ##  7 CPL   AR    4_10        1246  0.00267  population
    ##  8 CPL   AR    10_20        473  0.00102  population
    ##  9 CPL   AR    20_50        314  0.000674 population
    ## 10 CPL   AR    50_max       184  0.000395 population

## 2.3 Sample vs. population

Below, the lake distributions in the population of interest were
compared to the proportions in the observed sample. There were 352 lake
types in the sample compared to the 536 in the population of of
interest. There were 984 observations distributed across these 352 lake
types in the sample; and the number of samples was not distributed
evenly across the types. Some cells were represented by as few as 1
lake. In total, 536-352 = 184 lake types in the population of interest
were not represented in the sample.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/df_model.rda")

samp_props <- df_model %>%
  mutate(obs = 1) %>%
  group_by(WSA9, state, size_cat) %>%
  summarize(n_lakes = sum(obs), .groups = "drop") %>%
  ungroup() %>%
  mutate(prop_cell = round(n_lakes / sum(n_lakes), 7)) %>%
  mutate(type = "sample") 

save(samp_props, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/samp_props.rda")

samp_props %>%
  head(10)
```

    ## # A tibble: 10 x 6
    ##    WSA9  state size_cat n_lakes prop_cell type  
    ##    <fct> <fct> <ord>      <dbl>     <dbl> <chr> 
    ##  1 CPL   AL    min_4          1   0.00102 sample
    ##  2 CPL   AL    4_10           1   0.00102 sample
    ##  3 CPL   AL    10_20          1   0.00102 sample
    ##  4 CPL   AL    20_50          1   0.00102 sample
    ##  5 CPL   AL    50_max         1   0.00102 sample
    ##  6 CPL   AR    min_4          1   0.00102 sample
    ##  7 CPL   AR    4_10           1   0.00102 sample
    ##  8 CPL   AR    10_20          3   0.00305 sample
    ##  9 CPL   AR    20_50          1   0.00102 sample
    ## 10 CPL   AR    50_max         1   0.00102 sample

Below, a graphical comparison was constructed to depict the distribution
of cells in the population of interest *vs.* those in the sample.

<img src="NLA_N2O_models_files/figure-gfm/compare_sample_pop_cells-1.png" style="display: block; margin: auto;" />

Another comparison between population and sample was constructed by
ecoregion. The samples were not balanced across ecoregions. Lakes in the
Coastal Plains (CPL) ecoregion, for example, were clearly undersampled
relative to their proportion of the population.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/pframe.rda")

pframe_eco <- pframe %>%
  group_by(WSA9) %>%
  summarise(n_lakes = sum(n_lakes)) %>%
  ungroup() %>%
  mutate(prop_cell = round(n_lakes/sum(n_lakes), 7)) %>%
  ungroup() %>%
  mutate(type = 'population') 

save(pframe_eco, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/pframe_eco.rda")
```

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/samp_props.rda")

samp_props_eco <- samp_props %>%
  group_by(WSA9) %>%
  summarise(n_lakes = sum(n_lakes)) %>%
  ungroup() %>%
  mutate(prop_cell = round(n_lakes/sum(n_lakes), 7)) %>%
  ungroup() %>%
  mutate(type = 'sample')

save(samp_props_eco, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/samp_props_eco.rda")
```

<img src="NLA_N2O_models_files/figure-gfm/compare_eco_sample_pop_cells-1.png" style="display: block; margin: auto;" />

A similar comparison by state was constructed below.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/pframe.rda")

pframe_state <- pframe %>%
  group_by(state) %>%
  summarise(n_lakes = sum(n_lakes)) %>%
  ungroup() %>%
  mutate(prop_cell = round(n_lakes/sum(n_lakes), 7)) %>%
  ungroup() %>%
  mutate(type = 'population')

save(pframe_state, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/pframe_state.rda")
```

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/samp_props.rda")
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/samp_props.rda")

samp_props_state <- samp_props %>%
  group_by(state) %>%
  summarise(n_lakes = sum(n_lakes)) %>%
  ungroup() %>%
  mutate(prop_cell = round(n_lakes/sum(n_lakes), 7)) %>%
  ungroup() %>%
  mutate(type = 'sample')

save(samp_props_state, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/samp_props_state.rda")
```

<img src="NLA_N2O_models_files/figure-gfm/compare_state_sample_pop_cells-1.png" style="display: block; margin: auto;" />

Finally, a comparison by lake size category is shown below. Small lakes
were under-sampled relative to larger lakes.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/pframe.rda")

pframe_size <- pframe %>%
  group_by(size_cat) %>%
  summarise(n_lakes = sum(n_lakes)) %>%
  ungroup() %>%
  mutate(prop_cell = round(n_lakes/sum(n_lakes), 7)) %>%
  ungroup() %>%
  mutate(type = 'population')

save(pframe_size, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/pframe_size.rda")
```

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/samp_props.rda")

samp_props_size <- samp_props %>%
  group_by(size_cat) %>%
  summarise(n_lakes = sum(n_lakes)) %>%
  ungroup() %>%
  mutate(prop_cell = round(n_lakes/sum(n_lakes), 7)) %>%
  ungroup() %>%
  mutate(type = 'sample')

save(samp_props_size, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/samp_props_size.rda")
```

<img src="NLA_N2O_models_files/figure-gfm/compare_size_sample_pop_cells-1.png" style="display: block; margin: auto;" />

## 2.4 Sample-based estimates

Below are *naive* estimates, based only on the sample, for national
means for dissolved and equilibrium N2O and the saturation ratio.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/df_model.rda")

df_model %>%
  summarise(mean = mean(n2o),
             sd = sd(n2o)) %>%
  print()
```

    ##       mean      sd
    ## 1 8.720661 9.52093

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/df_model.rda")

df_model %>%
  summarise(mean = mean(n2o_eq),
             sd = sd(n2o_eq)) %>%
  print()
```

    ##       mean        sd
    ## 1 7.483567 0.8453779

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/df_model.rda")

df_model %>%
  summarise(mean = mean(n2o / n2o_eq),
             sd = sd(n2o / n2o_eq)) %>%
  print()
```

    ##       mean       sd
    ## 1 1.170095 1.302608

Roughly 67% of lakes in the sample were under-saturated (i.e.,
saturation ratio \< 1):

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/df_model.rda")

df_model %>%
  summarise(prop_undersat = sum((n2o / n2o_eq) < 1) / 984) %>%
  print()
```

    ##   prop_undersat
    ## 1     0.6656504

Using only the sample observations again, a plot was constructed of the
overall mean (dashed line) along with the ecoregion-specific means
(black circles). The shaded areas indicate +/- 1 standard deviation.
Neither dissolved N2O nor the saturation ratio were clearly structured
by ecoregion, but there did appear to be some structure along this
variable in the equilibrium N2O observations.

<img src="NLA_N2O_models_files/figure-gfm/sample_summary_eco-1.png" style="display: block; margin: auto;" />

The same summary by state is below.

<img src="NLA_N2O_models_files/figure-gfm/sample_summary_state-1.png" style="display: block; margin: auto;" />

Finally, the same summary by size category.

<img src="NLA_N2O_models_files/figure-gfm/sample_summary_size-1.png" style="display: block; margin: auto;" />

## 2.5 Sample data exploration

Below, the distribution N2O concentrations for the sample was summarized
using a density and rug plot. Note the natural log scale of the x-axis.
Both the dissolved and equilibrium N2O data had considerable right skew
even after the log transformation. This was not unexpected and has been
noted in other studies ([Webb et al. 2019](#ref-Webb_etal_2019)). The
saturation ratio was also skewed since it was derived from the other two
observed variables (i.e., sat_ratio = n2o / n2o_eq).

<img src="NLA_N2O_models_files/figure-gfm/summary_N2O-1.png" style="display: block; margin: auto;" />

Below are plots of dissolved N2O vs. NO3. The first plot shows log(N2O)
vs. log(NO3), as well as the ordinal categories assigned to NO3
(vertical lines). The leftmost vertical line is dashed and separates the
NO3 observations below the detection limit. The trend is increasing and
nonlinear on the log scale, with increasing variance in N2O as NO3
increased.

<img src="NLA_N2O_models_files/figure-gfm/summary_N2O_vs_NO3-1.png" style="display: block; margin: auto;" />

Below are plots of dissolved N2O vs. NO3 for 6 quantiles of the surface
temperature measurements (quantiles increasing from 1 to 6). This plot
suggested that the NO3 effect on N2O may have been stronger in lakes
with higher observed temperatures.

<img src="NLA_N2O_models_files/figure-gfm/summary_N2O_vs_NO3_surftemp-1.png" style="display: block; margin: auto;" />

The next plot below shows the relationship between dissolved N2O and NO3
at 6 different quantiles (increasing 1 to 6) of the log-scaled lake
surface area estimates.

<img src="NLA_N2O_models_files/figure-gfm/summary_N2O_vs_NO3_logarea-1.png" style="display: block; margin: auto;" />

Similar plots are below, but with NO3 expressed as an ordered
categorical variable with 5 levels. The positive and monotonic trends
area similar to the previous plots where NO3 was treated as continuous.
Note the large number of observations in the first NO3 category (no3_cat
= 1). This category represented all of the observations for NO3 that
were below the detection limit, which was most of the data.

<img src="NLA_N2O_models_files/figure-gfm/summary_N2O_vs_NO3cat-1.png" style="display: block; margin: auto;" />

<img src="NLA_N2O_models_files/figure-gfm/summary_N2O_vs_NO3cat_surftemp-1.png" style="display: block; margin: auto;" />

<img src="NLA_N2O_models_files/figure-gfm/summary_N2O_vs_NO3cat_logarea-1.png" style="display: block; margin: auto;" />

Below is a plot of log(dissolved N2O) vs. log(NO3) by ecoregion, which
suggested that the NO3 effect on N2O may have varied by ecoregion.

<img src="NLA_N2O_models_files/figure-gfm/summary_N2O_vs_NO3_ecoregion-1.png" style="display: block; margin: auto;" />

Below is the same plot as above but for the ordered categorical version
of NO3.

<img src="NLA_N2O_models_files/figure-gfm/summary_N2O_vs_NO3cat_ecoregion-1.png" style="display: block; margin: auto;" />

A plot below shows trends by state within just the Temperate Plains
(TPL) ecoregion. Within states, the number of observations were
relatively small, but the trends appeared closer to linear.

<img src="NLA_N2O_models_files/figure-gfm/summary_N2O_vs_NO3_wsa9state3-1.png" style="display: block; margin: auto;" />

# 3 Model fitting

The first regression model was constructed to estimate the joint
distribution of log-transformed dissolved and equilibrium N2O
conditional on the the design factors. Each log-transformed observation,
![i \\in 1,..,N=984](https://latex.codecogs.com/svg.image?i%20%5Cin%201%2C..%2CN%3D984 "i \in 1,..,N=984"),
for each response,
![p \\in 1:P=2](https://latex.codecogs.com/svg.image?p%20%5Cin%201%3AP%3D2 "p \in 1:P=2"),
was assumed to be drawn from a multivariate normal distribution with the
parameters ![\\nu](https://latex.codecogs.com/svg.image?%5Cnu "\nu") and
![\\Sigma](https://latex.codecogs.com/svg.image?%5CSigma "\Sigma"),
where ![\\nu](https://latex.codecogs.com/svg.image?%5Cnu "\nu") is the
multivariate mean estimated conditional on the design effects and
![\\Sigma](https://latex.codecogs.com/svg.image?%5CSigma "\Sigma") is a
covariance matrix containing the observation-level variances and
residual correlation:

![Y \\sim MVN(\\nu, \\Sigma)](https://latex.codecogs.com/svg.image?Y%20%5Csim%20MVN%28%5Cnu%2C%20%5CSigma%29 "Y \sim MVN(\nu, \Sigma)")

The multivariate mean is a vector of location parameters,
![\\nu:\[\\mu\_{p=1}, \\mu\_{p=2}\]](https://latex.codecogs.com/svg.image?%5Cnu%3A%5B%5Cmu_%7Bp%3D1%7D%2C%20%5Cmu_%7Bp%3D2%7D%5D "\nu:[\mu_{p=1}, \mu_{p=2}]"),
for each response. Each location parameter was further defined by a
linear combination of parameters where, for each response
![p](https://latex.codecogs.com/svg.image?p "p") and observation
![i](https://latex.codecogs.com/svg.image?i "i"):

![\\mu\_{pi} = \\alpha\_{0(pi)} + \\alpha\_{1(pij)} + \\alpha\_{2(pijk)} + \\alpha\_{3(pijkl)} \\\\
\\alpha_1 \\sim MVN(0, \\Lambda_1) \\\\
\\alpha_2 \\sim MVN(0, \\Lambda_2) \\\\
\\alpha_3 \\sim MVN(0, \\Lambda_3)](https://latex.codecogs.com/svg.image?%5Cmu_%7Bpi%7D%20%3D%20%5Calpha_%7B0%28pi%29%7D%20%2B%20%5Calpha_%7B1%28pij%29%7D%20%2B%20%5Calpha_%7B2%28pijk%29%7D%20%2B%20%5Calpha_%7B3%28pijkl%29%7D%20%5C%5C%0A%5Calpha_1%20%5Csim%20MVN%280%2C%20%5CLambda_1%29%20%5C%5C%0A%5Calpha_2%20%5Csim%20MVN%280%2C%20%5CLambda_2%29%20%5C%5C%0A%5Calpha_3%20%5Csim%20MVN%280%2C%20%5CLambda_3%29 "\mu_{pi} = \alpha_{0(pi)} + \alpha_{1(pij)} + \alpha_{2(pijk)} + \alpha_{3(pijkl)} \\
\alpha_1 \sim MVN(0, \Lambda_1) \\
\alpha_2 \sim MVN(0, \Lambda_2) \\
\alpha_3 \sim MVN(0, \Lambda_3)")

The linear combination included a fixed global intercept,
![a_0](https://latex.codecogs.com/svg.image?a_0 "a_0"), estimated
directly from the data, and three separate, latent group-level effects
matrices,
![\\alpha_1, \\alpha_2, \\alpha_3](https://latex.codecogs.com/svg.image?%5Calpha_1%2C%20%5Calpha_2%2C%20%5Calpha_3 "\alpha_1, \alpha_2, \alpha_3").
The group effects were assumed to be multivariate normal and were
centered on zero in the two-dimensional multivariate space. The spread
of the effects around zero were determined by a covariance matrix,
![\\Lambda_1, \\Lambda_2, \\text{or } \\Lambda_3](https://latex.codecogs.com/svg.image?%5CLambda_1%2C%20%5CLambda_2%2C%20%5Ctext%7Bor%20%7D%20%5CLambda_3 "\Lambda_1, \Lambda_2, \text{or } \Lambda_3"),
which were estimated directly from the data. The covariance terms were
further defined where:

![\\Lambda = \\begin{pmatrix} 1 & \\tau^2\_{p=1} \\\\ \\tau^2\_{p=2} & 1 \\end{pmatrix} \\chi \\begin{pmatrix} 1 & \\tau^2\_{p=1} \\\\ \\tau^2\_{p=2} & 1 \\end{pmatrix}](https://latex.codecogs.com/svg.image?%5CLambda%20%3D%20%5Cbegin%7Bpmatrix%7D%201%20%26%20%5Ctau%5E2_%7Bp%3D1%7D%20%5C%5C%20%5Ctau%5E2_%7Bp%3D2%7D%20%26%201%20%5Cend%7Bpmatrix%7D%20%5Cchi%20%5Cbegin%7Bpmatrix%7D%201%20%26%20%5Ctau%5E2_%7Bp%3D1%7D%20%5C%5C%20%5Ctau%5E2_%7Bp%3D2%7D%20%26%201%20%5Cend%7Bpmatrix%7D "\Lambda = \begin{pmatrix} 1 & \tau^2_{p=1} \\ \tau^2_{p=2} & 1 \end{pmatrix} \chi \begin{pmatrix} 1 & \tau^2_{p=1} \\ \tau^2_{p=2} & 1 \end{pmatrix}")

The ![\\tau](https://latex.codecogs.com/svg.image?%5Ctau "\tau")
parameters captured the group-level standard deviations, which constrain
the spread of group-level effects for each response, and
![\\chi](https://latex.codecogs.com/svg.image?%5Cchi "\chi") was the
group-level residual correlation matrix:

![\\chi = \\begin{pmatrix} 1 & \\varrho \\\\ \\varrho & 1 \\end{pmatrix}](https://latex.codecogs.com/svg.image?%5Cchi%20%3D%20%5Cbegin%7Bpmatrix%7D%201%20%26%20%5Cvarrho%20%5C%5C%20%5Cvarrho%20%26%201%20%5Cend%7Bpmatrix%7D "\chi = \begin{pmatrix} 1 & \varrho \\ \varrho & 1 \end{pmatrix}")

wherein
![\\varrho](https://latex.codecogs.com/svg.image?%5Cvarrho "\varrho")
captured the group-level residual correlation between effects.

The explicit indexing in the notation above conveys the relationship
between the parameters and each observation,
![i](https://latex.codecogs.com/svg.image?i "i"), and emphasizes the
nested structure of the observations within the group effects.
Specifically, each observation,
![i](https://latex.codecogs.com/svg.image?i "i"), was nested in a lake
size category, ![l](https://latex.codecogs.com/svg.image?l "l"), which
was nested in a state, ![k](https://latex.codecogs.com/svg.image?k "k"),
and ecoregion, ![j](https://latex.codecogs.com/svg.image?j "j"). The
parameter
![\\alpha_1](https://latex.codecogs.com/svg.image?%5Calpha_1 "\alpha_1"),
therefore, accounted for ecoregion-scale group effects or deviations
from the global mean;
![\\alpha_2](https://latex.codecogs.com/svg.image?%5Calpha_2 "\alpha_2")
accounted for state-level group effects nested in ecoregions; and
![\\alpha_3](https://latex.codecogs.com/svg.image?%5Calpha_3 "\alpha_3")
accounted for lake size group effects within states and ecoregions.

Finally, the observation-level covariance term,
![\\Sigma](https://latex.codecogs.com/svg.image?%5CSigma "\Sigma"), was
parameterized as:

![\\Sigma = \\begin{pmatrix} 1 & \\sigma^2\_{p=1} \\\\ \\sigma^2\_{p=2} & 1 \\end{pmatrix} \\Omega \\begin{pmatrix} 1 & \\sigma^2\_{p=1} \\\\ \\sigma^2\_{p=2} & 1 \\end{pmatrix}](https://latex.codecogs.com/svg.image?%5CSigma%20%3D%20%5Cbegin%7Bpmatrix%7D%201%20%26%20%5Csigma%5E2_%7Bp%3D1%7D%20%5C%5C%20%5Csigma%5E2_%7Bp%3D2%7D%20%26%201%20%5Cend%7Bpmatrix%7D%20%5COmega%20%5Cbegin%7Bpmatrix%7D%201%20%26%20%5Csigma%5E2_%7Bp%3D1%7D%20%5C%5C%20%5Csigma%5E2_%7Bp%3D2%7D%20%26%201%20%5Cend%7Bpmatrix%7D "\Sigma = \begin{pmatrix} 1 & \sigma^2_{p=1} \\ \sigma^2_{p=2} & 1 \end{pmatrix} \Omega \begin{pmatrix} 1 & \sigma^2_{p=1} \\ \sigma^2_{p=2} & 1 \end{pmatrix}")

wherein the
![\\sigma](https://latex.codecogs.com/svg.image?%5Csigma "\sigma")
parameters were the observation-level standard deviations for each
response and
![\\Omega](https://latex.codecogs.com/svg.image?%5COmega "\Omega") was
the observation-level residual correlation matrix:

![\\Omega = \\begin{pmatrix} 1 & \\rho \\\\ \\rho & 1 \\end{pmatrix}](https://latex.codecogs.com/svg.image?%5COmega%20%3D%20%5Cbegin%7Bpmatrix%7D%201%20%26%20%5Crho%20%5C%5C%20%5Crho%20%26%201%20%5Cend%7Bpmatrix%7D "\Omega = \begin{pmatrix} 1 & \rho \\ \rho & 1 \end{pmatrix}")

wherein ![\\rho](https://latex.codecogs.com/svg.image?%5Crho "\rho")
captured the residual correlation between responses.

For model fitting, priors were needed for all parameters conditioned
directly on the data, which included the global intercept, the scale
parameters, and the correlation matrices. A normal or Gaussian prior,
![N(\\mu = 2, \\sigma = 1)](https://latex.codecogs.com/svg.image?N%28%5Cmu%20%3D%202%2C%20%5Csigma%20%3D%201%29 "N(\mu = 2, \sigma = 1)")
centered near the (log-scale) data means, was used for the global mean
parameter for each response. This prior was considered minimally
informative as it placed most (\~80%) of the prior mass over values
between approximately 2 and 27 ng/L for median N2O or N2O equilibrium
concentration and included support in the tails for values approaching 0
ng/L on the lower end and 80 ng/L on the high end. We placed
![Exp(2)](https://latex.codecogs.com/svg.image?Exp%282%29 "Exp(2)")
priors over all scale parameters, which placed most of the support
between values very close to 0 and values near 1 (central 80% density
interval from approximately 0.005 to 1.15). Finally, for the correlation
matrices, an
![LKJ(\\eta =2)](https://latex.codecogs.com/svg.image?LKJ%28%5Ceta%20%3D2%29 "LKJ(\eta =2)")
prior was used, which, for a 2-dimensional response, placed most support
for correlations between approximately -0.9 and 0.9. This prior seemed
reasonable as there were no apparant causal mechanisms to ensure a very
strong correlation between the N2O measures. Any potential residual
dependence was expected to be indirect due to, for example, a common
correlate (e.g., elevation, temperature). For more information on prior
choice recommendations in Stan, see:
<https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations>

The
![\\textbf{brms}](https://latex.codecogs.com/svg.image?%5Ctextbf%7Bbrms%7D "\textbf{brms}")
package ([Bürkner 2017](#ref-Burkner_2017)) for
![\\textbf{R}](https://latex.codecogs.com/svg.image?%5Ctextbf%7BR%7D "\textbf{R}")
([R Core Team 2021](#ref-R_Core_Team_2021)) was used to fit all of the
models in a fully Bayesian setting. The formula syntax of the
![\\textbf{brms}](https://latex.codecogs.com/svg.image?%5Ctextbf%7Bbrms%7D "\textbf{brms}")
package is similar to the syntax used in the
![\\textbf{lme4}](https://latex.codecogs.com/svg.image?%5Ctextbf%7Blme4%7D "\textbf{lme4}")
package that is widely used to fit mixed effects models in frequentist
settings. In either package, the linear predictor for
![\\mu](https://latex.codecogs.com/svg.image?%5Cmu "\mu") described
above could be expressed as:

![\\sim 1 + (1\|WSA9) + (1\|WSA9:state) + (1\|WSA9:state:size)](https://latex.codecogs.com/svg.image?%5Csim%201%20%2B%20%281%7CWSA9%29%20%2B%20%281%7CWSA9%3Astate%29%20%2B%20%281%7CWSA9%3Astate%3Asize%29 "\sim 1 + (1|WSA9) + (1|WSA9:state) + (1|WSA9:state:size)")

In the
![\\textbf{brms}](https://latex.codecogs.com/svg.image?%5Ctextbf%7Bbrms%7D "\textbf{brms}")
package, there is additional functionality and syntax for multivariate
responses and for allowing the varying intercepts in a multivariate
model to be correlated, e.g.,:

![
\\begin{aligned} 
  N_2O\_{diss} \\sim 1 + (1\|a\|WSA9) + (1\|b\|WSA9:state) + (1\|c\|WSA9:state:size) \\\\
  N_2O\_{equi} \\sim 1 + (1\|a\|WSA9) + (1\|b\|WSA9:state) + (1\|c\|WSA9:state:size) 
\\end{aligned}
](https://latex.codecogs.com/svg.image?%0A%5Cbegin%7Baligned%7D%20%0A%20%20N_2O_%7Bdiss%7D%20%5Csim%201%20%2B%20%281%7Ca%7CWSA9%29%20%2B%20%281%7Cb%7CWSA9%3Astate%29%20%2B%20%281%7Cc%7CWSA9%3Astate%3Asize%29%20%5C%5C%0A%20%20N_2O_%7Bequi%7D%20%5Csim%201%20%2B%20%281%7Ca%7CWSA9%29%20%2B%20%281%7Cb%7CWSA9%3Astate%29%20%2B%20%281%7Cc%7CWSA9%3Astate%3Asize%29%20%0A%5Cend%7Baligned%7D%0A "
\begin{aligned} 
  N_2O_{diss} \sim 1 + (1|a|WSA9) + (1|b|WSA9:state) + (1|c|WSA9:state:size) \\
  N_2O_{equi} \sim 1 + (1|a|WSA9) + (1|b|WSA9:state) + (1|c|WSA9:state:size) 
\end{aligned}
")

The above syntax would indicate that the linear predictor for both
responses in the multivariate model have the same group-level varying
effects, and that those effects may be correlated between responses.

For the remainder of this document, only this simplified syntax is
presented to describe the model structure. For more information on
![\\textbf{brms}](https://latex.codecogs.com/svg.image?%5Ctextbf%7Bbrms%7D "\textbf{brms}")
functionality and syntax with multivariate response models, the package
vignette may be helpful, and can be found at:
<https://cran.r-project.org/web/packages/brms/vignettes/brms_multivariate.html>.

## 3.1 Model 1

The first model fit was as explained above.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/df_model.rda")

bf_n2o <- bf(log(n2o) ~ 1 + 
               (1 | a | WSA9) + 
               (1 | b | WSA9:state) + 
               (1 | c | WSA9:state:size_cat),
             family = gaussian())

bf_n2oeq <- bf(log(n2o_eq) ~ 1 + 
               (1 | a | WSA9) + 
               (1 | b | WSA9:state) +
               (1 | c | WSA9:state:size_cat),
             family = gaussian())

priors <- c(
  prior(normal(2, 1), class = "Intercept", resp = "logn2o"), # centered near data mean
  prior(exponential(2), class = "sd", resp = "logn2o"),
  prior(exponential(2), class = "sigma", resp = "logn2o"),
  prior(normal(2, 1), class = "Intercept", resp = "logn2oeq"), # centered near data mean
  prior(exponential(2), class = "sd", resp = "logn2oeq"),
  prior(exponential(2), class = "sigma", resp = "logn2oeq"),
  prior(lkj(2), class = "rescor"),
  prior(lkj(2), class = "cor")
  )

n2o_mod1 <- brm(bf_n2o + bf_n2oeq + set_rescor(rescor = TRUE),
                data = df_model,
                prior = priors,
                control = list(adapt_delta = 0.99, max_treedepth = 14),
                #sample_prior = "only",
                save_pars = save_pars(all = TRUE),
                seed = 145,
                chains=4, 
                iter=5000, 
                cores=4)

save(n2o_mod1, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/n2o_mod1.rda")
```

### 3.1.1 Summarize fit

The summaries of the estimated parameters and key HMC convergence
diagnostics for the fitted model are printed below. There were no
obvious issues with the HMC sampling. All
![\\hat{R}](https://latex.codecogs.com/svg.image?%5Chat%7BR%7D "\hat{R}")
values were less than 1.01 and effective sample size
(![ESS](https://latex.codecogs.com/svg.image?ESS "ESS")) calculations
suggested that the posterior contained a sufficient number of effective
samples for conducting inference.

    ##  Family: MV(gaussian, gaussian) 
    ##   Links: mu = identity; sigma = identity
    ##          mu = identity; sigma = identity 
    ## Formula: log(n2o) ~ 1 + (1 | a | WSA9) + (1 | b | WSA9:state) + (1 | c | WSA9:state:size_cat) 
    ##          log(n2o_eq) ~ 1 + (1 | a | WSA9) + (1 | b | WSA9:state) + (1 | c | WSA9:state:size_cat) 
    ##    Data: df_model (Number of observations: 984) 
    ##   Draws: 4 chains, each with iter = 5000; warmup = 2500; thin = 1;
    ##          total post-warmup draws = 10000
    ## 
    ## Priors: 
    ## Intercept_logn2o ~ normal(2, 1)
    ## Intercept_logn2oeq ~ normal(2, 1)
    ## L ~ lkj_corr_cholesky(2)
    ## Lrescor ~ lkj_corr_cholesky(2)
    ## <lower=0> sd_logn2o ~ exponential(2)
    ## <lower=0> sd_logn2oeq ~ exponential(2)
    ## <lower=0> sigma_logn2o ~ exponential(2)
    ## <lower=0> sigma_logn2oeq ~ exponential(2)
    ## 
    ## Group-Level Effects: 
    ## ~WSA9 (Number of levels: 9) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(logn2o_Intercept)                         0.04      0.03     0.00     0.13
    ## sd(logn2oeq_Intercept)                       0.06      0.02     0.03     0.11
    ## cor(logn2o_Intercept,logn2oeq_Intercept)     0.05      0.44    -0.80     0.83
    ##                                          Rhat Bulk_ESS Tail_ESS
    ## sd(logn2o_Intercept)                     1.00     2475     4230
    ## sd(logn2oeq_Intercept)                   1.00     3326     5332
    ## cor(logn2o_Intercept,logn2oeq_Intercept) 1.01     1141     3129
    ## 
    ## ~WSA9:state (Number of levels: 96) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(logn2o_Intercept)                         0.26      0.03     0.21     0.32
    ## sd(logn2oeq_Intercept)                       0.04      0.01     0.03     0.05
    ## cor(logn2o_Intercept,logn2oeq_Intercept)     0.24      0.15    -0.07     0.52
    ##                                          Rhat Bulk_ESS Tail_ESS
    ## sd(logn2o_Intercept)                     1.00     3199     4775
    ## sd(logn2oeq_Intercept)                   1.00     3549     5476
    ## cor(logn2o_Intercept,logn2oeq_Intercept) 1.00     3332     5553
    ## 
    ## ~WSA9:state:size_cat (Number of levels: 352) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(logn2o_Intercept)                         0.09      0.04     0.01     0.16
    ## sd(logn2oeq_Intercept)                       0.02      0.01     0.00     0.03
    ## cor(logn2o_Intercept,logn2oeq_Intercept)    -0.22      0.37    -0.85     0.58
    ##                                          Rhat Bulk_ESS Tail_ESS
    ## sd(logn2o_Intercept)                     1.00      863     1397
    ## sd(logn2oeq_Intercept)                   1.00      928     1365
    ## cor(logn2o_Intercept,logn2oeq_Intercept) 1.00     1661     2807
    ## 
    ## Population-Level Effects: 
    ##                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## logn2o_Intercept       2.02      0.04     1.95     2.09 1.00     2641     3963
    ## logn2oeq_Intercept     2.00      0.02     1.96     2.04 1.00     3446     4550
    ## 
    ## Family Specific Parameters: 
    ##                Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma_logn2o       0.40      0.01     0.37     0.42 1.00     2749     5523
    ## sigma_logn2oeq     0.08      0.00     0.08     0.09 1.00     4790     6698
    ## 
    ## Residual Correlations: 
    ##                         Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
    ## rescor(logn2o,logn2oeq)     0.21      0.04     0.14     0.27 1.00     6172
    ##                         Tail_ESS
    ## rescor(logn2o,logn2oeq)     7243
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

In the summary above, the estimated standard deviations for the varying
group effects on the mean behavior of the dissolved N2O response
suggested fairly low, but non-zero variability across each of the three
levels. The standard deviations estimated for the same varying effects
for equilibrium N2O were also relatively small. Finally, note the
relatively small, but positive residual correlation between the two N2O
responses.

Before investing too much into the inferences from this model, however,
the model fit was evaluated below using a series of graphical posterior
predictive checks (PPC, [Gelman et al. 2014](#ref-Gelman_etal_2014);
[Gelman, Hill, and Vehtari 2020](#ref-Gelman_etal_2020), Ch. 11).

### 3.1.2 Model checks

#### 3.1.2.1 Dissolved N2O

Below are a series of panels illustrating graphical PPCs for the
dissolved N2O component of the multivariate model. The top left panel
compares a density plot of the observed data (black line) to density
estimates drawn for 200 samples from the posterior predictive
distribution (PPD; blue lines) of the fitted model. The top right panel
compares the cumulative density distributions in the same manner. The
left middle panel compares means *vs.* standard deviations for 1000
draws from the PPD (blue dots) to the sample mean *vs.* standard
deviation (black dot). The right middle panel compares skewness *vs.*
kurtosis for 1000 draws from the PPD to the skewness *vs.* kurtosis
estimates for the observed data. The bottom left panel compares max
*vs.* min values for 1000 draws from the PPD to the max *vs.* min values
of the sample data. Finally, the bottom right panel shows the observed
*vs.* average predicted values for each observation in the sample. The
average predicted values were calculated as the mean prediction for each
observation in the PPD based on 1000 draws.

<img src="NLA_N2O_models_files/figure-gfm/ppc_n2o1-1.png" style="display: block; margin: auto;" />

The general takeaway from the PPCs above was that the model replicated
the central tendency of the observed data fairly well, but failed to
sufficiently replicate other important aspects of the distribution, such
as skewness and kurtosis. The observed *vs.* average predictions
scatterplot also suggested substantial heteroscedasticity in the errors.

The same checks were run below, but for the test set of 95 held-out,
second-visit data points. The patterns in misfit were similar to the
patterns indicated in the PPCs with the training data above.

<img src="NLA_N2O_models_files/figure-gfm/ppc_n2o1_test-1.png" style="display: block; margin: auto;" />

#### 3.1.2.2 Equilibrium N2O

Below are PPCs for the equilibrium N2O component of the model. As with
the dissolved N2O PPCs, the model seemed to do an OK job at replicating
the central tendency, but performed less well at replicating other
important aspects of the overall distribution.

<img src="NLA_N2O_models_files/figure-gfm/ppc_n2oeq1-1.png" style="display: block; margin: auto;" />

Below are the same PPCs for equilibrium N2O in the re-visit sites.

<img src="NLA_N2O_models_files/figure-gfm/ppc_n2oeq1_test-1.png" style="display: block; margin: auto;" />

#### 3.1.2.3 Bivariate

The graphical check below compared bivariate density contours estimated
from the observed data (black lines) to density contours estimated for
each of 20 draws from the PPD. The model appeared to do a good job of
replicating the bivariate mean, but was poor at representing the overall
joint distribution.

<img src="NLA_N2O_models_files/figure-gfm/ppc_biv1-1.png" style="display: block; margin: auto;" />

The same bivariate check is shown below for the re-visit data.

<img src="NLA_N2O_models_files/figure-gfm/ppc_biv1_test-1.png" style="display: block; margin: auto;" />

#### 3.1.2.4 Saturation

The graphical PPCs below were aimed at evaluating how well the
multivariate model did at representing the observed saturation ratio:

![\\dfrac{N_2O\_{diss}} {N_2O\_{equi}}](https://latex.codecogs.com/svg.image?%5Cdfrac%7BN_2O_%7Bdiss%7D%7D%20%7BN_2O_%7Bequi%7D%7D "\dfrac{N_2O_{diss}} {N_2O_{equi}}")

This quantity was estimated as a derived variable by dividing the
dissolved N2O PPD by the equilibrium N2O PPD. The proportion of
under-saturated lakes in the sample was estimated by summing the number
of lakes from each posterior predictive draw wherein the ratio was \< 1
and dividing that number by the total number of lakes in the sample,
which was 984. Overall, these checks indicated that properly
representing the tails of the dissolved N2O and N2O-eq observations
would likely be necessary in order to better replicate the observed
saturation metrics. The observed proportion of under-saturated lakes was
underestimated by more than 10 percentage points, on average.

The top left panel, below, is a density plot of the observed saturation
ratio (black line) compared to an estimate using 50 draws from the
derived PPD (blue lines). The top right panel shows the observed
proportion of under-saturated lakes compared to a model estimate based
on 1000 draws from the PPD. The left middle panel shows the mean *vs.*
standard deviation of the saturation ratio for the observed data
compared to the same estimates for 500 draws from the PPD. The right
middle panel shows the max *vs.* min for the sample compared to 500
draws from the PPD. Finally, the bottom left panel shows the observed
*vs.* average predicted saturation ratio for all 984 lakes sampled in
the dataset.

<img src="NLA_N2O_models_files/figure-gfm/ppc_sat1-1.png" style="display: block; margin: auto;" />

The same PPCs are shown below for the revisit data. These checks
indicated that the model did a similarly underwhelming job of
replicating the re-visit data.

<img src="NLA_N2O_models_files/figure-gfm/ppc_sat1_test-1.png" style="display: block; margin: auto;" />

#### 3.1.2.5 R-square

Below, the Bayesian
![R^2](https://latex.codecogs.com/svg.image?R%5E2 "R^2") values are
reported for each response in the model.

    ##          Estimate Est.Error  Q2.5 Q97.5
    ## R2logn2o    0.247     0.031 0.187 0.309

    ##            Estimate Est.Error  Q2.5 Q97.5
    ## R2logn2oeq    0.377     0.025 0.328 0.425

The ![R^2](https://latex.codecogs.com/svg.image?R%5E2 "R^2") were also
estimated for the re-visit data.

    ##          Estimate Est.Error  Q2.5 Q97.5
    ## R2logn2o    0.413      0.04 0.331 0.489

    ##            Estimate Est.Error  Q2.5 Q97.5
    ## R2logn2oeq    0.322     0.026 0.271 0.374

## 3.2 Model 2

In an attempt to better fit the observed data, the next model included a
distributional model for each sub-model that allowed for heterogeneous
variances. The distributional terms were each fit as a function of the
survey design structure. The same structure as for the models for the
mean components.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/df_model.rda")

bf_n2o <- bf(log(n2o) ~ 1 +
               (1 | a | WSA9) + 
               (1 | b | WSA9:state) + 
               (1 | c | WSA9:state:size_cat),
             sigma ~ 1 +
               (1 | WSA9) + 
               (1 | WSA9:state) + 
               (1 | WSA9:state:size_cat), 
             family = gaussian())

bf_n2oeq <- bf(log(n2o_eq) ~ 1 +
                 (1 | a | WSA9) + 
                 (1 | b | WSA9:state) +
                 (1 | c | WSA9:state:size_cat),
             sigma ~ 1 +
               (1 | WSA9) + 
               (1 | WSA9:state) + 
               (1 | WSA9:state:size_cat),
             family = gaussian())

priors <- c(
  prior(normal(2, 1), class = "Intercept", resp = "logn2o"),
  prior(exponential(2), class = "sd", resp = "logn2o"),
  prior(normal(-1, 2), class = "Intercept", dpar = "sigma", resp = "logn2o"),
  prior(exponential(2), class = "sd", dpar = "sigma", resp = "logn2o"),
  prior(normal(2, 1), class = "Intercept", resp = "logn2oeq"), 
  prior(exponential(2), class = "sd", resp = "logn2oeq"),
  prior(normal(-1, 2), class = "Intercept", dpar = "sigma", resp = "logn2oeq"),
  prior(exponential(2), class = "sd", dpar = "sigma", resp = "logn2oeq"),
  
  prior(lkj(2), class = "rescor"),
  prior(lkj(2), class = "cor")
  )

n2o_mod2 <- brm(bf_n2o + bf_n2oeq + set_rescor(rescor = TRUE),
                data = df_model, 
                prior = priors,
  control = list(adapt_delta = 0.975, max_treedepth = 12),
  #sample_prior = "only",
  save_pars = save_pars(all = TRUE),
  seed = 84512,
  chains=4, 
  iter=5000, 
  cores=4)

save(n2o_mod2, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/n2o_mod2.rda")
```

### 3.2.1 Summarize fit

The summaries of the estimated parameters and key HMC convergence
diagnostics for the fitted model are printed below.

    ##  Family: MV(gaussian, gaussian) 
    ##   Links: mu = identity; sigma = log
    ##          mu = identity; sigma = log 
    ## Formula: log(n2o) ~ 1 + (1 | a | WSA9) + (1 | b | WSA9:state) + (1 | c | WSA9:state:size_cat) 
    ##          sigma ~ 1 + (1 | WSA9) + (1 | WSA9:state) + (1 | WSA9:state:size_cat)
    ##          log(n2o_eq) ~ 1 + (1 | a | WSA9) + (1 | b | WSA9:state) + (1 | c | WSA9:state:size_cat) 
    ##          sigma ~ 1 + (1 | WSA9) + (1 | WSA9:state) + (1 | WSA9:state:size_cat)
    ##    Data: df_model (Number of observations: 984) 
    ##   Draws: 4 chains, each with iter = 5000; warmup = 2500; thin = 1;
    ##          total post-warmup draws = 10000
    ## 
    ## Priors: 
    ## Intercept_logn2o ~ normal(2, 1)
    ## Intercept_logn2o_sigma ~ normal(-1, 2)
    ## Intercept_logn2oeq ~ normal(2, 1)
    ## Intercept_logn2oeq_sigma ~ normal(-1, 2)
    ## L ~ lkj_corr_cholesky(2)
    ## Lrescor ~ lkj_corr_cholesky(2)
    ## <lower=0> sd_logn2o ~ exponential(2)
    ## <lower=0> sd_logn2o_sigma ~ exponential(2)
    ## <lower=0> sd_logn2oeq ~ exponential(2)
    ## <lower=0> sd_logn2oeq_sigma ~ exponential(2)
    ## 
    ## Group-Level Effects: 
    ## ~WSA9 (Number of levels: 9) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(logn2o_Intercept)                         0.06      0.03     0.02     0.12
    ## sd(logn2oeq_Intercept)                       0.05      0.02     0.03     0.09
    ## sd(sigma_logn2o_Intercept)                   0.22      0.12     0.02     0.50
    ## sd(sigma_logn2oeq_Intercept)                 0.16      0.08     0.04     0.34
    ## cor(logn2o_Intercept,logn2oeq_Intercept)     0.57      0.30    -0.18     0.95
    ##                                          Rhat Bulk_ESS Tail_ESS
    ## sd(logn2o_Intercept)                     1.00     2247     1772
    ## sd(logn2oeq_Intercept)                   1.00     3423     4761
    ## sd(sigma_logn2o_Intercept)               1.00     1892     3094
    ## sd(sigma_logn2oeq_Intercept)             1.00     2413     1991
    ## cor(logn2o_Intercept,logn2oeq_Intercept) 1.00     2696     3538
    ## 
    ## ~WSA9:state (Number of levels: 96) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(logn2o_Intercept)                         0.08      0.02     0.04     0.11
    ## sd(logn2oeq_Intercept)                       0.04      0.01     0.03     0.05
    ## sd(sigma_logn2o_Intercept)                   0.56      0.08     0.40     0.74
    ## sd(sigma_logn2oeq_Intercept)                 0.22      0.05     0.11     0.32
    ## cor(logn2o_Intercept,logn2oeq_Intercept)     0.42      0.19    -0.00     0.74
    ##                                          Rhat Bulk_ESS Tail_ESS
    ## sd(logn2o_Intercept)                     1.00     1276     1733
    ## sd(logn2oeq_Intercept)                   1.00     2695     5090
    ## sd(sigma_logn2o_Intercept)               1.00     2007     3462
    ## sd(sigma_logn2oeq_Intercept)             1.00     1851     1881
    ## cor(logn2o_Intercept,logn2oeq_Intercept) 1.00     1076     1960
    ## 
    ## ~WSA9:state:size_cat (Number of levels: 352) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(logn2o_Intercept)                         0.06      0.01     0.03     0.09
    ## sd(logn2oeq_Intercept)                       0.02      0.01     0.00     0.03
    ## sd(sigma_logn2o_Intercept)                   0.61      0.05     0.51     0.72
    ## sd(sigma_logn2oeq_Intercept)                 0.19      0.06     0.07     0.29
    ## cor(logn2o_Intercept,logn2oeq_Intercept)    -0.00      0.34    -0.70     0.62
    ##                                          Rhat Bulk_ESS Tail_ESS
    ## sd(logn2o_Intercept)                     1.00     1329     2182
    ## sd(logn2oeq_Intercept)                   1.00     1052     1452
    ## sd(sigma_logn2o_Intercept)               1.00     2622     4939
    ## sd(sigma_logn2oeq_Intercept)             1.00     1603     1725
    ## cor(logn2o_Intercept,logn2oeq_Intercept) 1.00      996     2063
    ## 
    ## Population-Level Effects: 
    ##                          Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
    ## logn2o_Intercept             1.93      0.03     1.88     1.99 1.00     3338
    ## sigma_logn2o_Intercept      -1.40      0.12    -1.63    -1.17 1.00     3511
    ## logn2oeq_Intercept           2.00      0.02     1.96     2.03 1.00     3384
    ## sigma_logn2oeq_Intercept    -2.60      0.07    -2.75    -2.45 1.00     4004
    ##                          Tail_ESS
    ## logn2o_Intercept             4462
    ## sigma_logn2o_Intercept       5077
    ## logn2oeq_Intercept           4145
    ## sigma_logn2oeq_Intercept     4848
    ## 
    ## Residual Correlations: 
    ##                         Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
    ## rescor(logn2o,logn2oeq)     0.36      0.03     0.29     0.43 1.00     6474
    ##                         Tail_ESS
    ## rescor(logn2o,logn2oeq)     7733
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

### 3.2.2 Model checks

Below the same PPCs were performed as with the initial model (see above
for more details on each panel). Though the checks below suggested some
improvement in replicating the tails of the observed data, the overall
fit again suggested room for improvement.

##### 3.2.2.0.1 Dissolved N2O

<img src="NLA_N2O_models_files/figure-gfm/ppc_n2o2-1.png" style="display: block; margin: auto;" />

#### 3.2.2.1 Equilibrium N2O

The checks below suggested that this model offered little no improvement
upon the initial model for equilibrium N2O.

<img src="NLA_N2O_models_files/figure-gfm/ppc_n2oeq2-1.png" style="display: block; margin: auto;" />

#### 3.2.2.2 Bivariate

This check perhaps suggested an improvement with regard to replicating
the joint density. However, the predictions were still clearly
over-dispersed relative to the observations.

<img src="NLA_N2O_models_files/figure-gfm/ppc_biv2-1.png" style="display: block; margin: auto;" />

#### 3.2.2.3 Saturation

The PPCs for the saturation metrics below indicated that including the
distributional models was perhaps an improvement on the initial model in
some aspects; in particular, the bias in the predicted proportion of
under-saturated lakes was substantially decreased. However, there
appeared to still be issues in replicating the tails as well as issues
with central tendency.

<img src="NLA_N2O_models_files/figure-gfm/ppc_sat2-1.png" style="display: block; margin: auto;" />

#### 3.2.2.4 R-square

Relative to model 1, there was a substantial decrease in the
![R^2](https://latex.codecogs.com/svg.image?R%5E2 "R^2") estimate for
the dissolved N2O component of this model. The estimate for the
equilibrium N2O-eq component was similar to the model 1.

    ##          Estimate Est.Error  Q2.5 Q97.5
    ## R2logn2o    0.056     0.009 0.038 0.075

    ##            Estimate Est.Error  Q2.5 Q97.5
    ## R2logn2oeq    0.379     0.022 0.335 0.421

## 3.3 Model 3

In the next model, the categorical version of the NO3 covariate and an
surface temperature covariate were included to try to improve the fit.
The ordinal NO3 variable was used as a monotonic, ordinal effect and
only in the dissolved N2O component of the model. For the equlibrium N2O
component, surface temperature and log-transformed elevation were used,
along with their interaction. The model also retained the same
distributional specifications included in model 2 above.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/df_model.rda")

bf_n2o <- bf(log(n2o) ~ mo(no3_cat) +
               surftemp +
               (mo(no3_cat) | a | WSA9) + 
               (mo(no3_cat) | b | WSA9:state) + 
               (1 | c | WSA9:state:size_cat),
             sigma ~ 1 +
               (1 | WSA9) + 
               (1 | WSA9:state) + 
               (1 | WSA9:state:size_cat), 
             family = gaussian())

bf_n2oeq <- bf(log(n2o_eq) ~ surftemp +
                 log_elev +
                 surftemp:log_elev +
                 (1 | a | WSA9) + 
                 (1 | b | WSA9:state) +
                 (1 | c | WSA9:state:size_cat),
             sigma ~ 1 +
               (1 | WSA9) + 
               (1 | WSA9:state) + 
               (1 | WSA9:state:size_cat),
             family = gaussian())

priors <- c(
  prior(normal(2, 1), class = "Intercept", resp = "logn2o"),
  prior(normal(0, 1), class = "b", resp = "logn2o"),
  prior(exponential(2), class = "sd", resp = "logn2o"),
  prior(normal(-1, 2), class = "Intercept", dpar = "sigma", resp = "logn2o"),
  prior(exponential(2), class = "sd", dpar = "sigma", resp = "logn2o"),
  prior(normal(2, 1), class = "Intercept", resp = "logn2oeq"), 
  prior(normal(0, 1), class = "b", resp = "logn2oeq"), 
  prior(exponential(2), class = "sd", resp = "logn2oeq"), 
  prior(normal(-1, 2), class = "Intercept", dpar = "sigma", resp = "logn2oeq"),
  prior(exponential(2), class = "sd", dpar = "sigma", resp = "logn2oeq"),
  
  prior(lkj(2), class = "rescor"),
  prior(lkj(2), class = "cor")
  )

n2o_mod3 <- brm(bf_n2o + bf_n2oeq + set_rescor(rescor = TRUE),
                data = df_model, 
                prior = priors,
  control = list(adapt_delta = 0.975, max_treedepth = 12),
  #sample_prior = "only",
  save_pars = save_pars(all = TRUE),
  seed = 98456,
  chains=4, 
  iter=5000, 
  cores=4)

save(n2o_mod3, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/n2o_mod3.rda")
```

### 3.3.1 Summarize fit

The fitted parameters and MCMC diagnostics are below.

    ##  Family: MV(gaussian, gaussian) 
    ##   Links: mu = identity; sigma = log
    ##          mu = identity; sigma = log 
    ## Formula: log(n2o) ~ mo(no3_cat) + surftemp + (mo(no3_cat) | a | WSA9) + (mo(no3_cat) | b | WSA9:state) + (1 | c | WSA9:state:size_cat) 
    ##          sigma ~ 1 + (1 | WSA9) + (1 | WSA9:state) + (1 | WSA9:state:size_cat)
    ##          log(n2o_eq) ~ surftemp + log_elev + surftemp:log_elev + (1 | a | WSA9) + (1 | b | WSA9:state) + (1 | c | WSA9:state:size_cat) 
    ##          sigma ~ 1 + (1 | WSA9) + (1 | WSA9:state) + (1 | WSA9:state:size_cat)
    ##    Data: df_model (Number of observations: 984) 
    ##   Draws: 4 chains, each with iter = 5000; warmup = 2500; thin = 1;
    ##          total post-warmup draws = 10000
    ## 
    ## Priors: 
    ## b_logn2o ~ normal(0, 1)
    ## b_logn2oeq ~ normal(0, 1)
    ## Intercept_logn2o ~ normal(2, 1)
    ## Intercept_logn2o_sigma ~ normal(-1, 2)
    ## Intercept_logn2oeq ~ normal(2, 1)
    ## Intercept_logn2oeq_sigma ~ normal(-1, 2)
    ## L ~ lkj_corr_cholesky(2)
    ## Lrescor ~ lkj_corr_cholesky(2)
    ## <lower=0> sd_logn2o ~ exponential(2)
    ## <lower=0> sd_logn2o_sigma ~ exponential(2)
    ## <lower=0> sd_logn2oeq ~ exponential(2)
    ## <lower=0> sd_logn2oeq_sigma ~ exponential(2)
    ## simo_logn2o_mono3_cat1 ~ dirichlet(1)
    ## 
    ## Group-Level Effects: 
    ## ~WSA9 (Number of levels: 9) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(logn2o_Intercept)                         0.05      0.02     0.02     0.10
    ## sd(logn2o_mono3_cat)                         0.14      0.05     0.06     0.26
    ## sd(logn2oeq_Intercept)                       0.04      0.01     0.02     0.07
    ## sd(sigma_logn2o_Intercept)                   0.12      0.08     0.01     0.31
    ## sd(sigma_logn2oeq_Intercept)                 0.36      0.11     0.19     0.64
    ## cor(logn2o_Intercept,logn2o_mono3_cat)      -0.14      0.33    -0.72     0.53
    ## cor(logn2o_Intercept,logn2oeq_Intercept)     0.36      0.29    -0.28     0.83
    ## cor(logn2o_mono3_cat,logn2oeq_Intercept)     0.37      0.29    -0.26     0.82
    ##                                          Rhat Bulk_ESS Tail_ESS
    ## sd(logn2o_Intercept)                     1.00     4265     3917
    ## sd(logn2o_mono3_cat)                     1.00     4747     5482
    ## sd(logn2oeq_Intercept)                   1.00     4798     5704
    ## sd(sigma_logn2o_Intercept)               1.00     3419     4599
    ## sd(sigma_logn2oeq_Intercept)             1.00     5027     6638
    ## cor(logn2o_Intercept,logn2o_mono3_cat)   1.00     4970     6018
    ## cor(logn2o_Intercept,logn2oeq_Intercept) 1.00     5073     6636
    ## cor(logn2o_mono3_cat,logn2oeq_Intercept) 1.00     7256     7538
    ## 
    ## ~WSA9:state (Number of levels: 96) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(logn2o_Intercept)                         0.03      0.02     0.00     0.07
    ## sd(logn2o_mono3_cat)                         0.14      0.02     0.10     0.18
    ## sd(logn2oeq_Intercept)                       0.03      0.00     0.03     0.04
    ## sd(sigma_logn2o_Intercept)                   0.30      0.09     0.09     0.47
    ## sd(sigma_logn2oeq_Intercept)                 0.28      0.06     0.17     0.40
    ## cor(logn2o_Intercept,logn2o_mono3_cat)      -0.30      0.32    -0.82     0.45
    ## cor(logn2o_Intercept,logn2oeq_Intercept)     0.00      0.27    -0.53     0.55
    ## cor(logn2o_mono3_cat,logn2oeq_Intercept)     0.45      0.13     0.16     0.69
    ##                                          Rhat Bulk_ESS Tail_ESS
    ## sd(logn2o_Intercept)                     1.00      894     2100
    ## sd(logn2o_mono3_cat)                     1.00     3795     5396
    ## sd(logn2oeq_Intercept)                   1.00     4164     5957
    ## sd(sigma_logn2o_Intercept)               1.00     1157     1279
    ## sd(sigma_logn2oeq_Intercept)             1.00     2747     4313
    ## cor(logn2o_Intercept,logn2o_mono3_cat)   1.00      550      693
    ## cor(logn2o_Intercept,logn2oeq_Intercept) 1.01      596      971
    ## cor(logn2o_mono3_cat,logn2oeq_Intercept) 1.00     2536     4431
    ## 
    ## ~WSA9:state:size_cat (Number of levels: 352) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(logn2o_Intercept)                         0.06      0.01     0.04     0.08
    ## sd(logn2oeq_Intercept)                       0.00      0.00     0.00     0.01
    ## sd(sigma_logn2o_Intercept)                   0.58      0.06     0.47     0.70
    ## sd(sigma_logn2oeq_Intercept)                 0.28      0.05     0.18     0.39
    ## cor(logn2o_Intercept,logn2oeq_Intercept)     0.26      0.39    -0.60     0.87
    ##                                          Rhat Bulk_ESS Tail_ESS
    ## sd(logn2o_Intercept)                     1.00     1429     2386
    ## sd(logn2oeq_Intercept)                   1.00     1779     3895
    ## sd(sigma_logn2o_Intercept)               1.00     2003     4421
    ## sd(sigma_logn2oeq_Intercept)             1.00     1812     3745
    ## cor(logn2o_Intercept,logn2oeq_Intercept) 1.00     4523     5554
    ## 
    ## Population-Level Effects: 
    ##                            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
    ## logn2o_Intercept               2.40      0.05     2.29     2.50 1.00     5046
    ## sigma_logn2o_Intercept        -1.70      0.08    -1.87    -1.55 1.00     5162
    ## logn2oeq_Intercept             3.10      0.05     3.01     3.19 1.00     8903
    ## sigma_logn2oeq_Intercept      -3.54      0.14    -3.81    -3.27 1.00     4065
    ## logn2o_surftemp               -0.02      0.00    -0.03    -0.02 1.00     5529
    ## logn2oeq_surftemp             -0.04      0.00    -0.04    -0.04 1.00     9722
    ## logn2oeq_log_elev             -0.07      0.01    -0.09    -0.06 1.00     9200
    ## logn2oeq_surftemp:log_elev     0.00      0.00     0.00     0.00 1.00     9527
    ## logn2o_mono3_cat               0.23      0.05     0.12     0.34 1.00     3948
    ##                            Tail_ESS
    ## logn2o_Intercept               6296
    ## sigma_logn2o_Intercept         5452
    ## logn2oeq_Intercept             7808
    ## sigma_logn2oeq_Intercept       5338
    ## logn2o_surftemp                6684
    ## logn2oeq_surftemp              7658
    ## logn2oeq_log_elev              7505
    ## logn2oeq_surftemp:log_elev     8028
    ## logn2o_mono3_cat               4777
    ## 
    ## Simplex Parameters: 
    ##                      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
    ## logn2o_mono3_cat1[1]     0.02      0.01     0.00     0.04 1.00     5535
    ## logn2o_mono3_cat1[2]     0.09      0.02     0.04     0.13 1.00     4281
    ## logn2o_mono3_cat1[3]     0.21      0.05     0.13     0.32 1.00     3683
    ## logn2o_mono3_cat1[4]     0.69      0.05     0.58     0.77 1.00     3498
    ##                      Tail_ESS
    ## logn2o_mono3_cat1[1]     5409
    ## logn2o_mono3_cat1[2]     5251
    ## logn2o_mono3_cat1[3]     4649
    ## logn2o_mono3_cat1[4]     4251
    ## 
    ## Residual Correlations: 
    ##                         Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
    ## rescor(logn2o,logn2oeq)     0.15      0.04     0.07     0.23 1.00    11381
    ##                         Tail_ESS
    ## rescor(logn2o,logn2oeq)     7610
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

### 3.3.2 Model checks

#### 3.3.2.1 Dissolved N2O

The PPCs below indicated a better fit compared to the previous models.
The central tendency and tail behavior looked to be reasonably
replicated by comparison. However, the observed *vs.* predicted plot
suggested that larger overserved values were likely being systematically
underestimated.

<img src="NLA_N2O_models_files/figure-gfm/ppc_full_checks_mod_n2o3-1.png" style="display: block; margin: auto;" />

#### 3.3.2.2 Equilibrium N2O

The PPCs below indicated that this model appeared to be an improvement
for equilibrium N2O as well. However, some checks (e.g., skewness)
suggested some room for additional improvement.

<img src="NLA_N2O_models_files/figure-gfm/ppc_full_checks_mod_n2oeq3-1.png" style="display: block; margin: auto;" />

#### 3.3.2.3 Bivariate

The check for the joint distribution below also suggested an improvement
upon the previous models.

<img src="NLA_N2O_models_files/figure-gfm/ppc_bv_check_mod_n2o3-1.png" style="display: block; margin: auto;" />

#### 3.3.2.4 Saturation

This model looked to be an improvement with regard to the PPCs for the
saturation metrics. However, the proportion of under-saturated lakes
remained biased low and other checks indicated that further improvements
would be ideal.

<img src="NLA_N2O_models_files/figure-gfm/ppc_sat_check_mod_n2o3-1.png" style="display: block; margin: auto;" />

#### 3.3.2.5 R-square

The ![R^2](https://latex.codecogs.com/svg.image?R%5E2 "R^2") estimates
for this model are below and suggested substantial improvements on the
previous models.

    ##          Estimate Est.Error  Q2.5 Q97.5
    ## R2logn2o    0.626     0.017 0.591  0.66

    ##            Estimate Est.Error Q2.5 Q97.5
    ## R2logn2oeq    0.879     0.004 0.87 0.886

### 3.3.3 Covariate effects

Below are plots illustrating the modeled effects of covariates on both
N2O and equilibrium N2O.

#### 3.3.3.1 Dissolved N2O

The conditional effects plots below for dissolved N2O illustrated a
positive, monotonic, and non-linear relationship with NO3; and a
negative, linear relationship with surface temperature.

<img src="NLA_N2O_models_files/figure-gfm/conditional_effects_mod_n2o3-1.png" style="display: block; margin: auto;" />

#### 3.3.3.2 Equilibrium N2O

The modeled effects for the equilibrium N2O component illustrated a
negative relationship with both predictors and an interaction such that
the surface temperature effect became slightly steeper at lower
elevations.

<img src="NLA_N2O_models_files/figure-gfm/conditional_effects_mod_n2oeq3-1.png" style="display: block; margin: auto;" />

## 3.4 Model 4

In the next model, covariate terms were also included in the
![\\sigma](https://latex.codecogs.com/svg.image?%5Csigma "\sigma")
components of both models in order to try to better capture remaining
heterogeneity in the variances of both N2O and N2O-eq.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/df_model.rda")

bf_n2o <- bf(log(n2o) ~ mo(no3_cat) +
               surftemp +
               (mo(no3_cat) | a | WSA9) + 
               (mo(no3_cat) | b | WSA9:state) + 
               (1 | c | WSA9:state:size_cat),
             sigma ~ mo(no3_cat) +
               surftemp +
               (1 | WSA9) + 
               (1 | WSA9:state) + 
               (1 | WSA9:state:size_cat), 
             family = gaussian())

bf_n2oeq <- bf(log(n2o_eq) ~ surftemp +
                 log_elev +
                 surftemp:log_elev +
                 (1 | a | WSA9) + 
                 (1 | b | WSA9:state) +
                 (1 | c | WSA9:state:size_cat),
             sigma ~ surftemp +
               log_elev +
               (1 | WSA9) + 
               (1 | WSA9:state) + 
               (1 | WSA9:state:size_cat),
             family = gaussian())

priors <- c(
  prior(normal(2, 1), class = "Intercept", resp = "logn2o"),
  prior(normal(0, 1), class = "b", resp = "logn2o"),
  prior(exponential(2), class = "sd", resp = "logn2o"),
  prior(normal(-1, 2), class = "Intercept", dpar = "sigma", resp = "logn2o"),
  prior(normal(0, 1), class = "b", dpar = "sigma", resp = "logn2o"),
  prior(exponential(2), class = "sd", dpar = "sigma", resp = "logn2o"),
  prior(normal(2, 1), class = "Intercept", resp = "logn2oeq"), 
  prior(normal(0, 1), class = "b", resp = "logn2oeq"), 
  prior(exponential(2), class = "sd", resp = "logn2oeq"), 
  prior(normal(-1, 2), class = "Intercept", dpar = "sigma", resp = "logn2oeq"),
  prior(normal(0, 1), class = "b", dpar = "sigma", resp = "logn2oeq"),
  prior(exponential(2), class = "sd", dpar = "sigma", resp = "logn2oeq"),
  
  prior(lkj(2), class = "rescor"),
  prior(lkj(2), class = "cor")
  )

n2o_mod4 <- brm(bf_n2o + bf_n2oeq + set_rescor(rescor = TRUE),
                data = df_model, 
                prior = priors,
  control = list(adapt_delta = 0.975, max_treedepth = 12),
  #sample_prior = "only",
  save_pars = save_pars(all = TRUE),
  seed = 15851,
  chains=4, 
  iter=5000, 
  cores=4)

save(n2o_mod4, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/n2o_mod4.rda")
```

### 3.4.1 Summarize fit

Below is a summary of the fitted parameters along with some convergence
diagnostics.

    ##  Family: MV(gaussian, gaussian) 
    ##   Links: mu = identity; sigma = log
    ##          mu = identity; sigma = log 
    ## Formula: log(n2o) ~ mo(no3_cat) + surftemp + (mo(no3_cat) | a | WSA9) + (mo(no3_cat) | b | WSA9:state) + (1 | c | WSA9:state:size_cat) 
    ##          sigma ~ mo(no3_cat) + surftemp + (1 | WSA9) + (1 | WSA9:state) + (1 | WSA9:state:size_cat)
    ##          log(n2o_eq) ~ surftemp + log_elev + surftemp:log_elev + (1 | a | WSA9) + (1 | b | WSA9:state) + (1 | c | WSA9:state:size_cat) 
    ##          sigma ~ surftemp + log_elev + (1 | WSA9) + (1 | WSA9:state) + (1 | WSA9:state:size_cat)
    ##    Data: df_model (Number of observations: 984) 
    ##   Draws: 4 chains, each with iter = 5000; warmup = 2500; thin = 1;
    ##          total post-warmup draws = 10000
    ## 
    ## Priors: 
    ## b_logn2o ~ normal(0, 1)
    ## b_logn2o_sigma ~ normal(0, 1)
    ## b_logn2oeq ~ normal(0, 1)
    ## b_logn2oeq_sigma ~ normal(0, 1)
    ## Intercept_logn2o ~ normal(2, 1)
    ## Intercept_logn2o_sigma ~ normal(-1, 2)
    ## Intercept_logn2oeq ~ normal(2, 1)
    ## Intercept_logn2oeq_sigma ~ normal(-1, 2)
    ## L ~ lkj_corr_cholesky(2)
    ## Lrescor ~ lkj_corr_cholesky(2)
    ## <lower=0> sd_logn2o ~ exponential(2)
    ## <lower=0> sd_logn2o_sigma ~ exponential(2)
    ## <lower=0> sd_logn2oeq ~ exponential(2)
    ## <lower=0> sd_logn2oeq_sigma ~ exponential(2)
    ## simo_logn2o_mono3_cat1 ~ dirichlet(1)
    ## simo_logn2o_sigma_mono3_cat1 ~ dirichlet(1)
    ## 
    ## Group-Level Effects: 
    ## ~WSA9 (Number of levels: 9) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(logn2o_Intercept)                        0.050     0.020    0.020    0.100
    ## sd(logn2o_mono3_cat)                        0.145     0.054    0.065    0.272
    ## sd(logn2oeq_Intercept)                      0.036     0.012    0.019    0.065
    ## sd(sigma_logn2o_Intercept)                  0.113     0.080    0.005    0.303
    ## sd(sigma_logn2oeq_Intercept)                0.208     0.098    0.039    0.435
    ## cor(logn2o_Intercept,logn2o_mono3_cat)     -0.186     0.327   -0.755    0.489
    ## cor(logn2o_Intercept,logn2oeq_Intercept)    0.344     0.290   -0.277    0.817
    ## cor(logn2o_mono3_cat,logn2oeq_Intercept)    0.364     0.286   -0.269    0.819
    ##                                           Rhat Bulk_ESS Tail_ESS
    ## sd(logn2o_Intercept)                     1.000     4733     4348
    ## sd(logn2o_mono3_cat)                     1.000     4386     6208
    ## sd(logn2oeq_Intercept)                   1.001     3977     5215
    ## sd(sigma_logn2o_Intercept)               1.001     2543     4693
    ## sd(sigma_logn2oeq_Intercept)             1.001     2777     2177
    ## cor(logn2o_Intercept,logn2o_mono3_cat)   1.000     4616     6092
    ## cor(logn2o_Intercept,logn2oeq_Intercept) 1.000     5667     6074
    ## cor(logn2o_mono3_cat,logn2oeq_Intercept) 1.001     7634     7624
    ## 
    ## ~WSA9:state (Number of levels: 96) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(logn2o_Intercept)                        0.035     0.017    0.003    0.068
    ## sd(logn2o_mono3_cat)                        0.117     0.022    0.076    0.162
    ## sd(logn2oeq_Intercept)                      0.033     0.003    0.027    0.040
    ## sd(sigma_logn2o_Intercept)                  0.181     0.099    0.011    0.374
    ## sd(sigma_logn2oeq_Intercept)                0.287     0.057    0.177    0.403
    ## cor(logn2o_Intercept,logn2o_mono3_cat)     -0.317     0.311   -0.810    0.405
    ## cor(logn2o_Intercept,logn2oeq_Intercept)    0.025     0.265   -0.497    0.557
    ## cor(logn2o_mono3_cat,logn2oeq_Intercept)    0.448     0.165    0.098    0.738
    ##                                           Rhat Bulk_ESS Tail_ESS
    ## sd(logn2o_Intercept)                     1.001      981     2017
    ## sd(logn2o_mono3_cat)                     1.001     3149     4612
    ## sd(logn2oeq_Intercept)                   1.000     3457     6360
    ## sd(sigma_logn2o_Intercept)               1.004      858     2373
    ## sd(sigma_logn2oeq_Intercept)             1.001     2470     3862
    ## cor(logn2o_Intercept,logn2o_mono3_cat)   1.003      945     1286
    ## cor(logn2o_Intercept,logn2oeq_Intercept) 1.005      556      823
    ## cor(logn2o_mono3_cat,logn2oeq_Intercept) 1.001     1518     3021
    ## 
    ## ~WSA9:state:size_cat (Number of levels: 352) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(logn2o_Intercept)                        0.065     0.011    0.042    0.087
    ## sd(logn2oeq_Intercept)                      0.004     0.002    0.000    0.008
    ## sd(sigma_logn2o_Intercept)                  0.539     0.055    0.432    0.647
    ## sd(sigma_logn2oeq_Intercept)                0.263     0.051    0.162    0.363
    ## cor(logn2o_Intercept,logn2oeq_Intercept)    0.395     0.348   -0.479    0.904
    ##                                           Rhat Bulk_ESS Tail_ESS
    ## sd(logn2o_Intercept)                     1.003     1478     2832
    ## sd(logn2oeq_Intercept)                   1.001     1506     2453
    ## sd(sigma_logn2o_Intercept)               1.002     1443     3611
    ## sd(sigma_logn2oeq_Intercept)             1.001     2372     4366
    ## cor(logn2o_Intercept,logn2oeq_Intercept) 1.001     4231     4956
    ## 
    ## Population-Level Effects: 
    ##                            Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS
    ## logn2o_Intercept              2.386     0.055    2.278    2.490 1.000     5971
    ## sigma_logn2o_Intercept       -1.855     0.281   -2.389   -1.294 1.000     5965
    ## logn2oeq_Intercept            3.115     0.051    3.016    3.217 1.001     7109
    ## sigma_logn2oeq_Intercept     -1.922     0.373   -2.634   -1.180 1.000     8445
    ## logn2o_surftemp              -0.021     0.002   -0.026   -0.017 1.000     5850
    ## sigma_logn2o_surftemp        -0.001     0.011   -0.023    0.021 1.001     6175
    ## logn2oeq_surftemp            -0.042     0.002   -0.046   -0.039 1.001     8285
    ## logn2oeq_log_elev            -0.080     0.008   -0.096   -0.065 1.001     7656
    ## logn2oeq_surftemp:log_elev    0.002     0.000    0.002    0.003 1.000     8407
    ## sigma_logn2oeq_surftemp      -0.065     0.010   -0.085   -0.046 1.000    10043
    ## sigma_logn2oeq_log_elev      -0.019     0.038   -0.095    0.054 1.000     6646
    ## logn2o_mono3_cat              0.225     0.058    0.108    0.340 1.000     4488
    ## sigma_logn2o_mono3_cat        0.256     0.037    0.187    0.331 1.000     6036
    ##                            Tail_ESS
    ## logn2o_Intercept               7057
    ## sigma_logn2o_Intercept         7704
    ## logn2oeq_Intercept             7366
    ## sigma_logn2oeq_Intercept       7581
    ## logn2o_surftemp                7677
    ## sigma_logn2o_surftemp          7768
    ## logn2oeq_surftemp              7516
    ## logn2oeq_log_elev              7564
    ## logn2oeq_surftemp:log_elev     7502
    ## sigma_logn2oeq_surftemp        8163
    ## sigma_logn2oeq_log_elev        7910
    ## logn2o_mono3_cat               5143
    ## sigma_logn2o_mono3_cat         7643
    ## 
    ## Simplex Parameters: 
    ##                            Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS
    ## logn2o_mono3_cat1[1]          0.018     0.012    0.001    0.045 1.000     5709
    ## logn2o_mono3_cat1[2]          0.093     0.027    0.041    0.149 1.000     5034
    ## logn2o_mono3_cat1[3]          0.231     0.060    0.126    0.366 1.000     4974
    ## logn2o_mono3_cat1[4]          0.659     0.062    0.520    0.763 1.000     4446
    ## sigma_logn2o_mono3_cat1[1]    0.107     0.066    0.007    0.255 1.000     6168
    ## sigma_logn2o_mono3_cat1[2]    0.129     0.087    0.006    0.328 1.000     8105
    ## sigma_logn2o_mono3_cat1[3]    0.459     0.151    0.168    0.758 1.000     7272
    ## sigma_logn2o_mono3_cat1[4]    0.304     0.139    0.038    0.569 1.000     6289
    ##                            Tail_ESS
    ## logn2o_mono3_cat1[1]           4429
    ## logn2o_mono3_cat1[2]           5645
    ## logn2o_mono3_cat1[3]           5478
    ## logn2o_mono3_cat1[4]           5210
    ## sigma_logn2o_mono3_cat1[1]     4471
    ## sigma_logn2o_mono3_cat1[2]     5010
    ## sigma_logn2o_mono3_cat1[3]     5802
    ## sigma_logn2o_mono3_cat1[4]     4264
    ## 
    ## Residual Correlations: 
    ##                         Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS
    ## rescor(logn2o,logn2oeq)    0.146     0.040    0.067    0.223 1.001     9839
    ##                         Tail_ESS
    ## rescor(logn2o,logn2oeq)     8508
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

### 3.4.2 Model checks

The same PPCs were employed for this model as above.

#### 3.4.2.1 Dissolved N2O

This model appeared to be an improvement on the previous model,
particularly with regard to the more constant variance indicated in the
observed *vs.* predicted plot (bottom, right panel).

<img src="NLA_N2O_models_files/figure-gfm/ppc_full_checks_mod_n2o4-1.png" style="display: block; margin: auto;" />

#### 3.4.2.2 Equilibrium N2O

This component of the model also seemed to be an improvement over model
3, with better representation in the tails as indicated in the skewness
*vs.* kurtosis PPC.

<img src="NLA_N2O_models_files/figure-gfm/ppc_full_checks_mod_n2oeq4-1.png" style="display: block; margin: auto;" />

#### 3.4.2.3 Bivariate

Again, an improvement over the previous model with a tighter fit of the
PPC to the observed bivariate density.

<img src="NLA_N2O_models_files/figure-gfm/ppc_bv_check_mod_n2o4-1.png" style="display: block; margin: auto;" />

#### 3.4.2.4 Saturation

This check also suggested an improvement over the previous models, with
better tail behavior and less bias in the proportion under-saturated
measure.

<img src="NLA_N2O_models_files/figure-gfm/ppc_sat_check_mod_n2o4-1.png" style="display: block; margin: auto;" />

#### 3.4.2.5 R-square

The Bayesian ![R^2](https://latex.codecogs.com/svg.image?R%5E2 "R^2")
estimates below indicated an improvement from the previous models.

    ##          Estimate Est.Error  Q2.5 Q97.5
    ## R2logn2o    0.606     0.025 0.551 0.651

    ##            Estimate Est.Error  Q2.5 Q97.5
    ## R2logn2oeq    0.875     0.005 0.865 0.883

### 3.4.3 Covariate effects

#### 3.4.3.1 Dissolved N2O

The conditional effects plots for the covariate effects on dissolved N2O
remained largely unchanged from the previous model.

<img src="NLA_N2O_models_files/figure-gfm/conditional_effects_mod_n2o4-1.png" style="display: block; margin: auto;" />

Below are estimates of the conditional effects of the covariates on
![\\sigma](https://latex.codecogs.com/svg.image?%5Csigma "\sigma") for
N2O. These plots suggested a large effect of NO3 on the variance of N2O,
but little to no effect of surface temperature.

<img src="NLA_N2O_models_files/figure-gfm/conditional_effects_sigma_n2o4-1.png" style="display: block; margin: auto;" />

#### 3.4.3.2 Equilibrium N2O

The covariate effects on equilibrium N2O remained largely the same as
for the previous model.

<img src="NLA_N2O_models_files/figure-gfm/conditional_effects_mod_n2oeq4-1.png" style="display: block; margin: auto;" />

The covariate effects on
![\\sigma](https://latex.codecogs.com/svg.image?%5Csigma "\sigma") for
N2O-eq suggested a negative effect of surface temperature and litte to
no effect of elevation.

<img src="NLA_N2O_models_files/figure-gfm/conditional_effects_sigma_n2oeq4-1.png" style="display: block; margin: auto;" />

## 3.5 Model 5

In the next model, additional complexity is added to the dissolved N2O
component by including a covariate for continuous lake surface area (log
scale) as well as interactions between NO3 and log(surface area) and
surface temperature.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/df_model.rda")

bf_n2o <- bf(log(n2o) ~ mo(no3_cat) +
               log_area +
               surftemp + 
               mo(no3_cat):log_area +
               mo(no3_cat):surftemp +
               (mo(no3_cat) | a | WSA9) + 
               (mo(no3_cat) | b | WSA9:state) + 
               (1 | c | WSA9:state:size_cat),
             sigma ~ log_area +
               mo(no3_cat) +
               (1 | WSA9) + 
               (1 | WSA9:state) + 
               (1 | WSA9:state:size_cat), 
             family = gaussian())

bf_n2oeq <- bf(log(n2o_eq) ~ surftemp +
                 log_elev +
                 surftemp:log_elev +
                 (1 | a | WSA9) + 
                 (1 | b | WSA9:state) +
                 (1 | c | WSA9:state:size_cat),
             sigma ~ surftemp +
               log_elev +
               (1 | WSA9) + 
               (1 | WSA9:state) + 
               (1 | WSA9:state:size_cat),
             family = gaussian())

priors <- c(
  prior(normal(2, 1), class = "Intercept", resp = "logn2o"),
  prior(normal(0, 1), class = "b", resp = "logn2o"),
  prior(exponential(2), class = "sd", resp = "logn2o"),
  prior(normal(-1, 2), class = "Intercept", dpar = "sigma", resp = "logn2o"),
  prior(normal(0, 1), class = "b", dpar = "sigma", resp = "logn2o"),
  prior(exponential(2), class = "sd", dpar = "sigma", resp = "logn2o"),
  
  prior(normal(2, 1), class = "Intercept", resp = "logn2oeq"), 
  prior(normal(0, 1), class = "b", resp = "logn2oeq"), 
  prior(exponential(2), class = "sd", resp = "logn2oeq"), 
  prior(normal(-1, 2), class = "Intercept", dpar = "sigma", resp = "logn2oeq"),
  prior(normal(0, 1), class = "b", dpar = "sigma", resp = "logn2oeq"),
  prior(exponential(2), class = "sd", dpar = "sigma", resp = "logn2oeq"),
  
  prior(lkj(2), class = "rescor"),
  prior(lkj(2), class = "cor")
  )

n2o_mod5 <- brm(bf_n2o + 
                  bf_n2oeq +
                  set_rescor(rescor = TRUE),
                data = df_model, 
                prior = priors,
  control = list(adapt_delta = 0.975, max_treedepth = 12),
  #sample_prior = "only",
  save_pars = save_pars(all = TRUE),
  seed = 54741,
  chains=4, 
  iter=5000, 
  cores=4)

save(n2o_mod5, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/n2o_mod5.rda")
```

### 3.5.1 Summarize fit

Below is a summary of the fitted parameters along with MCMC convergence
diagnostics.

    ##  Family: MV(gaussian, gaussian) 
    ##   Links: mu = identity; sigma = log
    ##          mu = identity; sigma = log 
    ## Formula: log(n2o) ~ mo(no3_cat) + log_area + surftemp + mo(no3_cat):log_area + mo(no3_cat):surftemp + (mo(no3_cat) | a | WSA9) + (mo(no3_cat) | b | WSA9:state) + (1 | c | WSA9:state:size_cat) 
    ##          sigma ~ log_area + mo(no3_cat) + (1 | WSA9) + (1 | WSA9:state) + (1 | WSA9:state:size_cat)
    ##          log(n2o_eq) ~ surftemp + log_elev + surftemp:log_elev + (1 | a | WSA9) + (1 | b | WSA9:state) + (1 | c | WSA9:state:size_cat) 
    ##          sigma ~ surftemp + log_elev + (1 | WSA9) + (1 | WSA9:state) + (1 | WSA9:state:size_cat)
    ##    Data: df_model (Number of observations: 984) 
    ##   Draws: 4 chains, each with iter = 5000; warmup = 2500; thin = 1;
    ##          total post-warmup draws = 10000
    ## 
    ## Priors: 
    ## b_logn2o ~ normal(0, 1)
    ## b_logn2o_sigma ~ normal(0, 1)
    ## b_logn2oeq ~ normal(0, 1)
    ## b_logn2oeq_sigma ~ normal(0, 1)
    ## Intercept_logn2o ~ normal(2, 1)
    ## Intercept_logn2o_sigma ~ normal(-1, 2)
    ## Intercept_logn2oeq ~ normal(2, 1)
    ## Intercept_logn2oeq_sigma ~ normal(-1, 2)
    ## L ~ lkj_corr_cholesky(2)
    ## Lrescor ~ lkj_corr_cholesky(2)
    ## <lower=0> sd_logn2o ~ exponential(2)
    ## <lower=0> sd_logn2o_sigma ~ exponential(2)
    ## <lower=0> sd_logn2oeq ~ exponential(2)
    ## <lower=0> sd_logn2oeq_sigma ~ exponential(2)
    ## simo_logn2o_mono3_cat:log_area1 ~ dirichlet(1)
    ## simo_logn2o_mono3_cat:surftemp1 ~ dirichlet(1)
    ## simo_logn2o_mono3_cat1 ~ dirichlet(1)
    ## simo_logn2o_sigma_mono3_cat1 ~ dirichlet(1)
    ## 
    ## Group-Level Effects: 
    ## ~WSA9 (Number of levels: 9) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(logn2o_Intercept)                        0.048     0.019    0.020    0.094
    ## sd(logn2o_mono3_cat)                        0.081     0.050    0.007    0.199
    ## sd(logn2oeq_Intercept)                      0.034     0.011    0.018    0.062
    ## sd(sigma_logn2o_Intercept)                  0.111     0.076    0.006    0.293
    ## sd(sigma_logn2oeq_Intercept)                0.209     0.102    0.036    0.445
    ## cor(logn2o_Intercept,logn2o_mono3_cat)     -0.056     0.360   -0.713    0.644
    ## cor(logn2o_Intercept,logn2oeq_Intercept)    0.464     0.284   -0.173    0.885
    ## cor(logn2o_mono3_cat,logn2oeq_Intercept)    0.259     0.341   -0.467    0.824
    ##                                           Rhat Bulk_ESS Tail_ESS
    ## sd(logn2o_Intercept)                     1.000     4413     4624
    ## sd(logn2o_mono3_cat)                     1.003     1412     2265
    ## sd(logn2oeq_Intercept)                   1.000     4061     6633
    ## sd(sigma_logn2o_Intercept)               1.000     2790     4064
    ## sd(sigma_logn2oeq_Intercept)             1.001     2176     1898
    ## cor(logn2o_Intercept,logn2o_mono3_cat)   1.000     5823     6525
    ## cor(logn2o_Intercept,logn2oeq_Intercept) 1.000     4935     6214
    ## cor(logn2o_mono3_cat,logn2oeq_Intercept) 1.000     4684     4770
    ## 
    ## ~WSA9:state (Number of levels: 96) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(logn2o_Intercept)                        0.046     0.014    0.016    0.072
    ## sd(logn2o_mono3_cat)                        0.097     0.024    0.053    0.146
    ## sd(logn2oeq_Intercept)                      0.033     0.003    0.027    0.040
    ## sd(sigma_logn2o_Intercept)                  0.207     0.088    0.024    0.369
    ## sd(sigma_logn2oeq_Intercept)                0.285     0.056    0.177    0.395
    ## cor(logn2o_Intercept,logn2o_mono3_cat)     -0.265     0.284   -0.762    0.336
    ## cor(logn2o_Intercept,logn2oeq_Intercept)    0.163     0.200   -0.232    0.559
    ## cor(logn2o_mono3_cat,logn2oeq_Intercept)    0.266     0.209   -0.155    0.649
    ##                                           Rhat Bulk_ESS Tail_ESS
    ## sd(logn2o_Intercept)                     1.001     1323     1765
    ## sd(logn2o_mono3_cat)                     1.001     2612     2569
    ## sd(logn2oeq_Intercept)                   1.001     3321     5286
    ## sd(sigma_logn2o_Intercept)               1.007      776     1082
    ## sd(sigma_logn2oeq_Intercept)             1.002     2464     3702
    ## cor(logn2o_Intercept,logn2o_mono3_cat)   1.001     1708     3060
    ## cor(logn2o_Intercept,logn2oeq_Intercept) 1.005      960     1715
    ## cor(logn2o_mono3_cat,logn2oeq_Intercept) 1.002      883     1843
    ## 
    ## ~WSA9:state:size_cat (Number of levels: 352) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(logn2o_Intercept)                        0.035     0.015    0.004    0.062
    ## sd(logn2oeq_Intercept)                      0.003     0.002    0.000    0.007
    ## sd(sigma_logn2o_Intercept)                  0.479     0.056    0.374    0.591
    ## sd(sigma_logn2oeq_Intercept)                0.260     0.052    0.160    0.360
    ## cor(logn2o_Intercept,logn2oeq_Intercept)    0.235     0.413   -0.658    0.872
    ##                                           Rhat Bulk_ESS Tail_ESS
    ## sd(logn2o_Intercept)                     1.004      819     1375
    ## sd(logn2oeq_Intercept)                   1.003     1569     3553
    ## sd(sigma_logn2o_Intercept)               1.006     1318     3346
    ## sd(sigma_logn2oeq_Intercept)             1.001     2326     3796
    ## cor(logn2o_Intercept,logn2oeq_Intercept) 1.002     2950     5425
    ## 
    ## Population-Level Effects: 
    ##                            Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS
    ## logn2o_Intercept              2.380     0.055    2.273    2.487 1.001     4020
    ## sigma_logn2o_Intercept       -1.596     0.097   -1.789   -1.408 1.002     4961
    ## logn2oeq_Intercept            3.116     0.051    3.019    3.218 1.000     6860
    ## sigma_logn2oeq_Intercept     -1.922     0.374   -2.639   -1.179 1.000     7865
    ## logn2o_log_area               0.029     0.003    0.023    0.034 1.000     9054
    ## logn2o_surftemp              -0.025     0.002   -0.029   -0.021 1.000     3996
    ## sigma_logn2o_log_area        -0.095     0.020   -0.135   -0.055 1.001     6068
    ## logn2oeq_surftemp            -0.042     0.002   -0.046   -0.039 1.000     7418
    ## logn2oeq_log_elev            -0.080     0.008   -0.097   -0.065 1.000     6834
    ## logn2oeq_surftemp:log_elev    0.003     0.000    0.002    0.003 1.000     7339
    ## sigma_logn2oeq_surftemp      -0.065     0.010   -0.085   -0.046 1.000     9671
    ## sigma_logn2oeq_log_elev      -0.018     0.038   -0.096    0.054 1.000     6265
    ## logn2o_mono3_cat              0.026     0.127   -0.226    0.275 1.004     1971
    ## logn2o_mono3_cat:log_area    -0.036     0.010   -0.054   -0.016 1.001     2602
    ## logn2o_mono3_cat:surftemp     0.014     0.006    0.003    0.026 1.004     1478
    ## sigma_logn2o_mono3_cat        0.246     0.036    0.179    0.321 1.001     4246
    ##                            Tail_ESS
    ## logn2o_Intercept               6632
    ## sigma_logn2o_Intercept         6145
    ## logn2oeq_Intercept             7546
    ## sigma_logn2oeq_Intercept       7774
    ## logn2o_log_area                8090
    ## logn2o_surftemp                7106
    ## sigma_logn2o_log_area          7491
    ## logn2oeq_surftemp              8290
    ## logn2oeq_log_elev              7865
    ## logn2oeq_surftemp:log_elev     7992
    ## sigma_logn2oeq_surftemp        8132
    ## sigma_logn2oeq_log_elev        7038
    ## logn2o_mono3_cat               3252
    ## logn2o_mono3_cat:log_area      4549
    ## logn2o_mono3_cat:surftemp      2416
    ## sigma_logn2o_mono3_cat         6738
    ## 
    ## Simplex Parameters: 
    ##                               Estimate Est.Error l-95% CI u-95% CI  Rhat
    ## logn2o_mono3_cat1[1]             0.026     0.024    0.001    0.089 1.001
    ## logn2o_mono3_cat1[2]             0.091     0.066    0.005    0.252 1.001
    ## logn2o_mono3_cat1[3]             0.250     0.146    0.034    0.615 1.002
    ## logn2o_mono3_cat1[4]             0.633     0.166    0.206    0.867 1.003
    ## logn2o_mono3_cat:log_area1[1]    0.060     0.041    0.005    0.161 1.000
    ## logn2o_mono3_cat:log_area1[2]    0.042     0.039    0.001    0.142 1.000
    ## logn2o_mono3_cat:log_area1[3]    0.297     0.161    0.047    0.681 1.000
    ## logn2o_mono3_cat:log_area1[4]    0.602     0.176    0.166    0.865 1.000
    ## logn2o_mono3_cat:surftemp1[1]    0.043     0.040    0.006    0.144 1.002
    ## logn2o_mono3_cat:surftemp1[2]    0.080     0.056    0.015    0.222 1.001
    ## logn2o_mono3_cat:surftemp1[3]    0.276     0.111    0.105    0.573 1.001
    ## logn2o_mono3_cat:surftemp1[4]    0.601     0.142    0.168    0.785 1.002
    ## sigma_logn2o_mono3_cat1[1]       0.131     0.074    0.010    0.288 1.000
    ## sigma_logn2o_mono3_cat1[2]       0.151     0.096    0.010    0.367 1.000
    ## sigma_logn2o_mono3_cat1[3]       0.444     0.149    0.149    0.734 1.001
    ## sigma_logn2o_mono3_cat1[4]       0.275     0.138    0.026    0.544 1.001
    ##                               Bulk_ESS Tail_ESS
    ## logn2o_mono3_cat1[1]              4261     5305
    ## logn2o_mono3_cat1[2]              1875     3916
    ## logn2o_mono3_cat1[3]              2027     2856
    ## logn2o_mono3_cat1[4]              1411     2439
    ## logn2o_mono3_cat:log_area1[1]     4395     3509
    ## logn2o_mono3_cat:log_area1[2]     7129     6046
    ## logn2o_mono3_cat:log_area1[3]     5220     5683
    ## logn2o_mono3_cat:log_area1[4]     4431     4505
    ## logn2o_mono3_cat:surftemp1[1]     3541     3002
    ## logn2o_mono3_cat:surftemp1[2]     3667     2580
    ## logn2o_mono3_cat:surftemp1[3]     3634     3352
    ## logn2o_mono3_cat:surftemp1[4]     2731     2289
    ## sigma_logn2o_mono3_cat1[1]        6840     4269
    ## sigma_logn2o_mono3_cat1[2]        7630     5795
    ## sigma_logn2o_mono3_cat1[3]        6331     6964
    ## sigma_logn2o_mono3_cat1[4]        4902     5797
    ## 
    ## Residual Correlations: 
    ##                         Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS
    ## rescor(logn2o,logn2oeq)    0.141     0.038    0.066    0.216 1.000    11829
    ##                         Tail_ESS
    ## rescor(logn2o,logn2oeq)     8204
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

### 3.5.2 Model checks

Again, the same PPCs as above were performed for this model.

#### 3.5.2.1 Dissolved N2O

This PPCs for dissolved N2O looked similar to the previous model.

<img src="NLA_N2O_models_files/figure-gfm/ppc_full_checks_mod_n2o5-1.png" style="display: block; margin: auto;" />

#### 3.5.2.2 Equilibrium N2O

Again, the PPCs for this model were similar to the previous model, which
was unsurprising given that it was the same model for N2O-eq.

<img src="NLA_N2O_models_files/figure-gfm/ppc_full_checks_mod_n2oeq5-1.png" style="display: block; margin: auto;" />

#### 3.5.2.3 Bivariate

This PPC was also similar to the previous model.

<img src="NLA_N2O_models_files/figure-gfm/ppc_bv_check_mod_n2o5-1.png" style="display: block; margin: auto;" />

#### 3.5.2.4 Saturation

This check was also similar to the prevoius model, with perhaps slightly
less bias in the proportion unsaturated estimates. There was also a
potentially concerning extreme prediction in the observed *vs* predicted
PPC.

<img src="NLA_N2O_models_files/figure-gfm/ppc_sat_check_mod_n2o5-1.png" style="display: block; margin: auto;" />

#### 3.5.2.5 R-square

    ##          Estimate Est.Error  Q2.5 Q97.5
    ## R2logn2o    0.629      0.03 0.563  0.68

    ##            Estimate Est.Error  Q2.5 Q97.5
    ## R2logn2oeq    0.874     0.005 0.864 0.882

### 3.5.3 Covariate effects

#### 3.5.3.1 Dissolved N2O

The conditional effects plots suggested a similar effect of NO3, but
interesting interactions between NO3 and lake area and NO3 and surface
temperature. For lake area, the effect was estimated to be larger and
more negative at the highest levels of NO3; and slightly negative at the
lowest level of NO3. For surface temperature, the effect was estimated
to be largest and positive at the highest level of NO3; and negative at
the lowest level of NO3.

<img src="NLA_N2O_models_files/figure-gfm/conditional_effects_mod_n2o5-1.png" style="display: block; margin: auto;" />

The estimated covariate effects on
![\\sigma](https://latex.codecogs.com/svg.image?%5Csigma "\sigma")
suggested a negative relationship with log(area) and a positive
relationship, again, with NO3.

<img src="NLA_N2O_models_files/figure-gfm/conditional_effects_sigma_n2o5-1.png" style="display: block; margin: auto;" />

#### 3.5.3.2 Equilibrium N2O

The estimated covariate effects on equilibrium N2O remained largely the
same as estimated in the previous model.

<img src="NLA_N2O_models_files/figure-gfm/conditional_effects_mod_n2oeq5-1.png" style="display: block; margin: auto;" />

<img src="NLA_N2O_models_files/figure-gfm/conditional_effects_sigma_n2oeq5-1.png" style="display: block; margin: auto;" />

## 3.6 A Final Model

As demonstrated above, models excluding the NO3 covariate consistently
resulted in poorer fits to to the observed dissolved N2O data and
potentially strongly biased estimates of the saturation ratio. Including
surface temperature and elevation in the equilibrium N2O part of the
model also resulted in substantially improved replication of key aspects
of the observed data. Likewise, added flexibility in the distributional
terms for both dissolved and equilibrium N2O led to improvements.

To make inferences to the population of interest from a model including
these covariates, however, however, the covariates needed to be (1)
fully observed across that population or (2) their missingness needed to
be modeled. For the lake area and elevation covariates, data *was*
available for all lakes from previously compiled geospatial databases.
However, neither surface temperature or NO3 were observed for lakes
outside of the sample. That is, they were partially observed with
respect to the target population. Therefore, a more complex model was
constructed below that included surface temperature and NO3 as
additional responses conditioned on the survey design variables and
fully observed covariates. This approach to inference for N2O was
similar to Bayesian structural equation modeling approaches ([Merkle et
al. 2021](#ref-Merkle_etal_2021); [Merkle and Rosseel
2018](#ref-Merkle_Rosseel_2018)). The logical dependence structure could
be characterized as:

![\\begin{aligned} 
{\\boldsymbol{N_2O\_{diss}}} &\\sim Survey + Area + {\\boldsymbol{NO_3}} + {\\boldsymbol{Temp}} \\\\ 
{\\boldsymbol{N_2O\_{equil}}} &\\sim Survey + Elev  + {\\boldsymbol{Temp}}\\\\
{\\boldsymbol{NO_3}} &\\sim Survey + Area + {\\boldsymbol{Temp}} \\\\ 
{\\boldsymbol{Temp}} &\\sim Survey + Lat + Elev + Day
\\end{aligned}](https://latex.codecogs.com/svg.image?%5Cbegin%7Baligned%7D%20%0A%7B%5Cboldsymbol%7BN_2O_%7Bdiss%7D%7D%7D%20%26%5Csim%20Survey%20%2B%20Area%20%2B%20%7B%5Cboldsymbol%7BNO_3%7D%7D%20%2B%20%7B%5Cboldsymbol%7BTemp%7D%7D%20%5C%5C%20%0A%7B%5Cboldsymbol%7BN_2O_%7Bequil%7D%7D%7D%20%26%5Csim%20Survey%20%2B%20Elev%20%20%2B%20%7B%5Cboldsymbol%7BTemp%7D%7D%5C%5C%0A%7B%5Cboldsymbol%7BNO_3%7D%7D%20%26%5Csim%20Survey%20%2B%20Area%20%2B%20%7B%5Cboldsymbol%7BTemp%7D%7D%20%5C%5C%20%0A%7B%5Cboldsymbol%7BTemp%7D%7D%20%26%5Csim%20Survey%20%2B%20Lat%20%2B%20Elev%20%2B%20Day%0A%5Cend%7Baligned%7D "\begin{aligned} 
{\boldsymbol{N_2O_{diss}}} &\sim Survey + Area + {\boldsymbol{NO_3}} + {\boldsymbol{Temp}} \\ 
{\boldsymbol{N_2O_{equil}}} &\sim Survey + Elev  + {\boldsymbol{Temp}}\\
{\boldsymbol{NO_3}} &\sim Survey + Area + {\boldsymbol{Temp}} \\ 
{\boldsymbol{Temp}} &\sim Survey + Lat + Elev + Day
\end{aligned}")

Variables in bold text above were treated as partially observed with
respect to the population of interest (i.e., observed only in the
sample), whereas variables not in bold were considered fully observed.
The partially observed variables, being dissolved and equilibrium N2O,
NO3, and surface temperature, were each modeled conditional on the
survey design variables and other partially and/or fully observed
covariates. This piece-wise approach required a more complex set of
post-processing steps compared to a typical MRP analysis. In order to
propagate estimates and uncertainty through the dependency structure and
make inferences, the fitted model was used to first predict surface
temperature in the target population, since it depended only on the
fully observed covariates. That predictive distribution was then used
alongside the relevant fully observed covariates to predict NO3 in the
target population. Finally, the predictive distributions for temperature
and NO3 were used to predict the N2O responses. These steps were carried
out in the “Predict to population” section to follow.

In the final model below, the sub-model for surface temperature assumed
a Gamma distributed error distribution and the linear predictor included
the survey design variables, latitude, elevation, and Julian date. The
shape parameter was also modeled as a function of latitude to address
increasing response variance along the latitudinal gradient. The NO3
sub-model was a cumulative logit model and the linear predictor included
all of the survey factors as well as surface temperature and lake area.

The dissolved and equilibrium N2O responses were each modeled with Gamma
distributed errors, but with the same covariate structure as in model 5.
The same structure was also employed for the shape terms in these
responses, corresponding to the
![\\sigma](https://latex.codecogs.com/svg.image?%5Csigma "\sigma") terms
in the previous models. Though not shown in this document, the Gamma
error structure appeared to result in slightly better performance in the
predictive checks compared to the Gaussian errors in previous models.
This was primarily apparent in the saturation ratio checks, which may
have been more sensitive to model performance in the tails of the N2O
responses. Others have also indicated that the Gamma error distribution
can work well for dissolved N2O data ([Webb et al.
2019](#ref-Webb_etal_2019)).

Note that there was no observation-level residual correlation term for
this model, since the residuals are undefined for the Gamma and
cumulative logit models. Dropping the observation-level residual
correlation term was deemed a reasonable compromise that enabled the
inclusion of NO3 as a covariate on dissolved N2O. The random intercepts,
however, still allowed for potential correlations between the four
responses at the group levels.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/df_model.rda")

bf_n2o <- bf(n2o ~ mo(no3_cat) +
               log_area +
               surftemp + 
               mo(no3_cat):log_area +
               mo(no3_cat):surftemp +
               (mo(no3_cat) | a | WSA9) + 
               (mo(no3_cat) | b | WSA9:state) + 
               (1 | c | WSA9:state:size_cat),
             shape ~ log_area +
               mo(no3_cat) +
               (1 | WSA9) + 
               (1 | WSA9:state) + 
               (1 | WSA9:state:size_cat),
             family = Gamma(link = "log"))

bf_n2oeq <- bf(n2o_eq ~ surftemp +
                 log_elev +
                 surftemp:log_elev +
                 (1 | a | WSA9) + 
                 (1 | b | WSA9:state) +
                 (1 | c | WSA9:state:size_cat),
             shape ~ surftemp +
               log_elev +
               (1 | WSA9) + 
               (1 | WSA9:state) + 
               (1 | WSA9:state:size_cat),
             family = Gamma(link = "log"))

bf_temp <- bf(surftemp ~ lat +
                s(log_elev) +
                s(jdate) +
                (1 | a | WSA9) + 
                (1 | b | WSA9:state) +
                (1 | c | WSA9:state:size_cat),
              shape ~ lat,
              family = Gamma(link = "log"))

bf_no3 <- bf(no3_cat ~ surftemp +
               log_area +
               (1 | a | WSA9) +
               (1 | b | WSA9:state) +
               (1 | c | WSA9:state:size_cat),
             family = cumulative(link = "logit", threshold="flexible"))

priors <- c(
  prior(normal(2, 1), class = "Intercept", resp = "n2o"),
  prior(normal(0, 1), class = "b", resp = "n2o"),
  prior(exponential(2), class = "sd", resp = "n2o"),
  prior(normal(5, 4), class = "Intercept", dpar = "shape", resp = "n2o"),
  prior(normal(0, 1), class = "b", dpar = "shape", resp = "n2o"),
  prior(exponential(2), class = "sd", dpar = "shape", resp = "n2o"),
  
  prior(normal(2, 1), class = "Intercept", resp = "n2oeq"), 
  prior(normal(0, 1), class = "b", resp = "n2oeq"),  
  prior(exponential(2), class = "sd", resp = "n2oeq"),
  prior(normal(5, 4), class = "Intercept", dpar = "shape", resp = "n2oeq"),
  prior(normal(0, 1), class = "b", dpar = "shape", resp = "n2oeq"),
  prior(exponential(2), class = "sd", dpar = "shape", resp = "n2oeq"),
  
  prior(normal(3, 1), class = "Intercept", resp = "surftemp"), 
  prior(normal(0, 1), class = "b", resp = "surftemp"), 
  prior(exponential(0.5), class = "sds", resp = "surftemp"),
  prior(exponential(2), class = "sd", resp = "surftemp"),
  prior(normal(5, 4), class = "Intercept", dpar = "shape", resp = "surftemp"),
  prior(normal(0, 1), class = "b", dpar = "shape", resp = "surftemp"),
  
  prior(normal(0, 3), class = "Intercept", resp = "no3cat"),
  prior(normal(0, 1), class = "b", resp = "no3cat"),
  prior(exponential(1), class = "sd", resp = "no3cat"),
  
  prior(lkj(2), class = "cor")
  )

n2o_mod6 <- brm(bf_n2o + 
                  bf_n2oeq + 
                  bf_temp + 
                  bf_no3 + 
                  set_rescor(rescor = FALSE),
                data = df_model, 
                prior = priors,
  control = list(adapt_delta = 0.975, max_treedepth = 14),
  #sample_prior = "only",
  save_pars = save_pars(all = TRUE),
  seed = 85132,#14548,
  #init = my_inits,
  init_r = 0.5,
  chains=4, 
  iter=5000, 
  cores=4)

save(n2o_mod6, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/n2o_mod6.rda")
```

### 3.6.1 Summarize fit

Below is a summary of the fitted parameters and MCMC diagnostics.

    ##  Family: MV(gamma, gamma, gamma, cumulative) 
    ##   Links: mu = log; shape = log
    ##          mu = log; shape = log
    ##          mu = log; shape = log
    ##          mu = logit; disc = identity 
    ## Formula: n2o ~ mo(no3_cat) + log_area + surftemp + mo(no3_cat):log_area + mo(no3_cat):surftemp + (mo(no3_cat) | a | WSA9) + (mo(no3_cat) | b | WSA9:state) + (1 | c | WSA9:state:size_cat) 
    ##          shape ~ log_area + mo(no3_cat) + (1 | WSA9) + (1 | WSA9:state) + (1 | WSA9:state:size_cat)
    ##          n2o_eq ~ surftemp + log_elev + surftemp:log_elev + (1 | a | WSA9) + (1 | b | WSA9:state) + (1 | c | WSA9:state:size_cat) 
    ##          shape ~ surftemp + log_elev + (1 | WSA9) + (1 | WSA9:state) + (1 | WSA9:state:size_cat)
    ##          surftemp ~ lat + s(log_elev) + s(jdate) + (1 | a | WSA9) + (1 | b | WSA9:state) + (1 | c | WSA9:state:size_cat) 
    ##          shape ~ lat
    ##          no3_cat ~ surftemp + log_area + (1 | a | WSA9) + (1 | b | WSA9:state) + (1 | c | WSA9:state:size_cat) 
    ##    Data: df_model (Number of observations: 984) 
    ##   Draws: 4 chains, each with iter = 5000; warmup = 2500; thin = 1;
    ##          total post-warmup draws = 10000
    ## 
    ## Priors: 
    ## b_n2o ~ normal(0, 1)
    ## b_n2o_shape ~ normal(0, 1)
    ## b_n2oeq ~ normal(0, 1)
    ## b_n2oeq_shape ~ normal(0, 1)
    ## b_no3cat ~ normal(0, 1)
    ## b_surftemp ~ normal(0, 1)
    ## b_surftemp_shape ~ normal(0, 1)
    ## Intercept_n2o ~ normal(2, 1)
    ## Intercept_n2o_shape ~ normal(5, 4)
    ## Intercept_n2oeq ~ normal(2, 1)
    ## Intercept_n2oeq_shape ~ normal(5, 4)
    ## Intercept_no3cat ~ normal(0, 3)
    ## Intercept_surftemp ~ normal(3, 1)
    ## Intercept_surftemp_shape ~ normal(5, 4)
    ## L ~ lkj_corr_cholesky(2)
    ## <lower=0> sd_n2o ~ exponential(2)
    ## <lower=0> sd_n2o_shape ~ exponential(2)
    ## <lower=0> sd_n2oeq ~ exponential(2)
    ## <lower=0> sd_n2oeq_shape ~ exponential(2)
    ## <lower=0> sd_no3cat ~ exponential(1)
    ## <lower=0> sd_surftemp ~ exponential(2)
    ## <lower=0> sds_surftemp ~ exponential(0.5)
    ## simo_n2o_mono3_cat:log_area1 ~ dirichlet(1)
    ## simo_n2o_mono3_cat:surftemp1 ~ dirichlet(1)
    ## simo_n2o_mono3_cat1 ~ dirichlet(1)
    ## simo_n2o_shape_mono3_cat1 ~ dirichlet(1)
    ## 
    ## Smooth Terms: 
    ##                           Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS
    ## sds(surftemp_slog_elev_1)    1.161     0.370    0.638    2.079 1.001     2264
    ## sds(surftemp_sjdate_1)       0.571     0.277    0.226    1.273 1.000     2679
    ##                           Tail_ESS
    ## sds(surftemp_slog_elev_1)     4174
    ## sds(surftemp_sjdate_1)        4853
    ## 
    ## Group-Level Effects: 
    ## ~WSA9 (Number of levels: 9) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(n2o_Intercept)                           0.048     0.018    0.020    0.091
    ## sd(n2o_mono3_cat)                           0.045     0.035    0.002    0.129
    ## sd(n2oeq_Intercept)                         0.033     0.011    0.018    0.061
    ## sd(surftemp_Intercept)                      0.031     0.014    0.011    0.064
    ## sd(no3cat_Intercept)                        0.690     0.256    0.296    1.308
    ## sd(shape_n2o_Intercept)                     0.211     0.141    0.011    0.535
    ## sd(shape_n2oeq_Intercept)                   0.402     0.183    0.075    0.814
    ## cor(n2o_Intercept,n2o_mono3_cat)           -0.049     0.336   -0.667    0.614
    ## cor(n2o_Intercept,n2oeq_Intercept)          0.393     0.268   -0.202    0.824
    ## cor(n2o_mono3_cat,n2oeq_Intercept)          0.102     0.328   -0.559    0.686
    ## cor(n2o_Intercept,surftemp_Intercept)      -0.350     0.296   -0.839    0.285
    ## cor(n2o_mono3_cat,surftemp_Intercept)       0.089     0.333   -0.564    0.703
    ## cor(n2oeq_Intercept,surftemp_Intercept)    -0.167     0.299   -0.705    0.437
    ## cor(n2o_Intercept,no3cat_Intercept)        -0.057     0.295   -0.605    0.522
    ## cor(n2o_mono3_cat,no3cat_Intercept)         0.141     0.333   -0.539    0.724
    ## cor(n2oeq_Intercept,no3cat_Intercept)       0.272     0.274   -0.305    0.740
    ## cor(surftemp_Intercept,no3cat_Intercept)    0.185     0.300   -0.436    0.718
    ##                                           Rhat Bulk_ESS Tail_ESS
    ## sd(n2o_Intercept)                        1.000     3011     2931
    ## sd(n2o_mono3_cat)                        1.003     1540     3502
    ## sd(n2oeq_Intercept)                      1.000     3147     4572
    ## sd(surftemp_Intercept)                   1.001     3705     4498
    ## sd(no3cat_Intercept)                     1.000     4184     5618
    ## sd(shape_n2o_Intercept)                  1.000     2379     3641
    ## sd(shape_n2oeq_Intercept)                1.001     2208     1617
    ## cor(n2o_Intercept,n2o_mono3_cat)         1.000     6558     6134
    ## cor(n2o_Intercept,n2oeq_Intercept)       1.000     4604     6010
    ## cor(n2o_mono3_cat,n2oeq_Intercept)       1.000     3241     5029
    ## cor(n2o_Intercept,surftemp_Intercept)    1.000     4744     5916
    ## cor(n2o_mono3_cat,surftemp_Intercept)    1.000     5246     6447
    ## cor(n2oeq_Intercept,surftemp_Intercept)  1.000     7342     7508
    ## cor(n2o_Intercept,no3cat_Intercept)      1.000     5152     6479
    ## cor(n2o_mono3_cat,no3cat_Intercept)      1.002     3435     5897
    ## cor(n2oeq_Intercept,no3cat_Intercept)    1.001     6449     7400
    ## cor(surftemp_Intercept,no3cat_Intercept) 1.000     6490     7928
    ## 
    ## ~WSA9:state (Number of levels: 96) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(n2o_Intercept)                           0.047     0.013    0.021    0.072
    ## sd(n2o_mono3_cat)                           0.101     0.019    0.068    0.142
    ## sd(n2oeq_Intercept)                         0.033     0.003    0.026    0.040
    ## sd(surftemp_Intercept)                      0.035     0.006    0.023    0.047
    ## sd(no3cat_Intercept)                        0.878     0.128    0.649    1.147
    ## sd(shape_n2o_Intercept)                     0.383     0.178    0.032    0.707
    ## sd(shape_n2oeq_Intercept)                   0.561     0.114    0.336    0.781
    ## cor(n2o_Intercept,n2o_mono3_cat)           -0.145     0.254   -0.602    0.376
    ## cor(n2o_Intercept,n2oeq_Intercept)          0.183     0.189   -0.180    0.549
    ## cor(n2o_mono3_cat,n2oeq_Intercept)          0.201     0.164   -0.122    0.519
    ## cor(n2o_Intercept,surftemp_Intercept)      -0.023     0.237   -0.476    0.438
    ## cor(n2o_mono3_cat,surftemp_Intercept)      -0.232     0.215   -0.644    0.194
    ## cor(n2oeq_Intercept,surftemp_Intercept)    -0.134     0.190   -0.494    0.246
    ## cor(n2o_Intercept,no3cat_Intercept)         0.462     0.211    0.018    0.822
    ## cor(n2o_mono3_cat,no3cat_Intercept)         0.145     0.200   -0.251    0.531
    ## cor(n2oeq_Intercept,no3cat_Intercept)       0.054     0.140   -0.220    0.329
    ## cor(surftemp_Intercept,no3cat_Intercept)   -0.231     0.189   -0.586    0.154
    ##                                           Rhat Bulk_ESS Tail_ESS
    ## sd(n2o_Intercept)                        1.004     1025     1337
    ## sd(n2o_mono3_cat)                        1.001     3366     4975
    ## sd(n2oeq_Intercept)                      1.001     2809     4611
    ## sd(surftemp_Intercept)                   1.000     4124     5116
    ## sd(no3cat_Intercept)                     1.000     4389     6247
    ## sd(shape_n2o_Intercept)                  1.004      649     1304
    ## sd(shape_n2oeq_Intercept)                1.002     1908     1810
    ## cor(n2o_Intercept,n2o_mono3_cat)         1.002     1069     1955
    ## cor(n2o_Intercept,n2oeq_Intercept)       1.003      794     1260
    ## cor(n2o_mono3_cat,n2oeq_Intercept)       1.002      833     1873
    ## cor(n2o_Intercept,surftemp_Intercept)    1.004     1849     3618
    ## cor(n2o_mono3_cat,surftemp_Intercept)    1.001     2626     4267
    ## cor(n2oeq_Intercept,surftemp_Intercept)  1.000     6340     6860
    ## cor(n2o_Intercept,no3cat_Intercept)      1.004      893     1944
    ## cor(n2o_mono3_cat,no3cat_Intercept)      1.000     1746     3085
    ## cor(n2oeq_Intercept,no3cat_Intercept)    1.000     6378     7983
    ## cor(surftemp_Intercept,no3cat_Intercept) 1.001     3573     5949
    ## 
    ## ~WSA9:state:size_cat (Number of levels: 352) 
    ##                                          Estimate Est.Error l-95% CI u-95% CI
    ## sd(n2o_Intercept)                           0.038     0.014    0.006    0.064
    ## sd(n2oeq_Intercept)                         0.004     0.002    0.000    0.008
    ## sd(surftemp_Intercept)                      0.010     0.007    0.000    0.024
    ## sd(no3cat_Intercept)                        0.308     0.182    0.016    0.673
    ## sd(shape_n2o_Intercept)                     0.895     0.110    0.684    1.114
    ## sd(shape_n2oeq_Intercept)                   0.511     0.104    0.309    0.715
    ## cor(n2o_Intercept,n2oeq_Intercept)          0.310     0.349   -0.486    0.841
    ## cor(n2o_Intercept,surftemp_Intercept)       0.015     0.365   -0.677    0.705
    ## cor(n2oeq_Intercept,surftemp_Intercept)    -0.102     0.381   -0.761    0.659
    ## cor(n2o_Intercept,no3cat_Intercept)        -0.090     0.342   -0.716    0.601
    ## cor(n2oeq_Intercept,no3cat_Intercept)      -0.143     0.359   -0.757    0.611
    ## cor(surftemp_Intercept,no3cat_Intercept)    0.173     0.378   -0.599    0.807
    ##                                           Rhat Bulk_ESS Tail_ESS
    ## sd(n2o_Intercept)                        1.011      580     1036
    ## sd(n2oeq_Intercept)                      1.009     1062     3160
    ## sd(surftemp_Intercept)                   1.003     1917     3889
    ## sd(no3cat_Intercept)                     1.002     1078     2386
    ## sd(shape_n2o_Intercept)                  1.002     1068     2786
    ## sd(shape_n2oeq_Intercept)                1.002     1397     2071
    ## cor(n2o_Intercept,n2oeq_Intercept)       1.004     2099     4705
    ## cor(n2o_Intercept,surftemp_Intercept)    1.000     4876     6414
    ## cor(n2oeq_Intercept,surftemp_Intercept)  1.001     4195     5741
    ## cor(n2o_Intercept,no3cat_Intercept)      1.001     3188     5334
    ## cor(n2oeq_Intercept,no3cat_Intercept)    1.001     2749     4946
    ## cor(surftemp_Intercept,no3cat_Intercept) 1.001     2607     5750
    ## 
    ## Population-Level Effects: 
    ##                          Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS
    ## n2o_Intercept               2.392     0.056    2.285    2.500 1.000     3063
    ## shape_n2o_Intercept         3.215     0.189    2.849    3.589 1.001     3997
    ## n2oeq_Intercept             3.111     0.053    3.011    3.218 1.000     5146
    ## shape_n2oeq_Intercept       3.873     0.751    2.385    5.326 1.001     5619
    ## surftemp_Intercept          3.791     0.060    3.672    3.906 1.000     7414
    ## shape_surftemp_Intercept    8.637     0.460    7.721    9.522 1.001    10817
    ## no3cat_Intercept[1]        -3.025     0.615   -4.258   -1.864 1.000     5932
    ## no3cat_Intercept[2]        -2.059     0.605   -3.274   -0.903 1.000     6278
    ## no3cat_Intercept[3]        -1.027     0.600   -2.235    0.115 1.000     6689
    ## no3cat_Intercept[4]        -0.028     0.602   -1.247    1.113 1.000     7002
    ## n2o_log_area                0.028     0.003    0.022    0.034 1.000     6958
    ## n2o_surftemp               -0.025     0.002   -0.029   -0.021 1.000     3145
    ## shape_n2o_log_area          0.190     0.041    0.110    0.271 1.000     4322
    ## n2oeq_surftemp             -0.042     0.002   -0.045   -0.038 1.000     5836
    ## n2oeq_log_elev             -0.080     0.008   -0.097   -0.064 1.000     5425
    ## n2oeq_surftemp:log_elev     0.002     0.000    0.002    0.003 1.000     5857
    ## shape_n2oeq_surftemp        0.131     0.020    0.092    0.171 1.001     6725
    ## shape_n2oeq_log_elev        0.030     0.077   -0.117    0.186 1.001     4914
    ## surftemp_lat               -0.016     0.001   -0.019   -0.013 1.000     7349
    ## shape_surftemp_lat         -0.105     0.011   -0.127   -0.083 1.000    11330
    ## no3cat_surftemp            -0.141     0.023   -0.187   -0.096 1.001     6713
    ## no3cat_log_area             0.068     0.035   -0.001    0.137 1.000    11692
    ## surftemp_slog_elev_1       -3.477     0.479   -4.414   -2.557 1.000     5891
    ## surftemp_sjdate_1          -0.008     0.528   -1.074    1.017 1.000     4010
    ## n2o_mono3_cat               0.007     0.087   -0.172    0.175 1.004     1650
    ## n2o_mono3_cat:log_area     -0.046     0.009   -0.063   -0.027 1.001     2079
    ## n2o_mono3_cat:surftemp      0.018     0.004    0.009    0.026 1.004     1417
    ## shape_n2o_mono3_cat        -0.510     0.070   -0.649   -0.375 1.001     4691
    ##                          Tail_ESS
    ## n2o_Intercept                5517
    ## shape_n2o_Intercept          6108
    ## n2oeq_Intercept              6075
    ## shape_n2oeq_Intercept        6535
    ## surftemp_Intercept           7632
    ## shape_surftemp_Intercept     7718
    ## no3cat_Intercept[1]          6622
    ## no3cat_Intercept[2]          6695
    ## no3cat_Intercept[3]          6957
    ## no3cat_Intercept[4]          6940
    ## n2o_log_area                 7668
    ## n2o_surftemp                 6260
    ## shape_n2o_log_area           5605
    ## n2oeq_surftemp               6494
    ## n2oeq_log_elev               6234
    ## n2oeq_surftemp:log_elev      6696
    ## shape_n2oeq_surftemp         7112
    ## shape_n2oeq_log_elev         6636
    ## surftemp_lat                 7598
    ## shape_surftemp_lat           7919
    ## no3cat_surftemp              7328
    ## no3cat_log_area              8177
    ## surftemp_slog_elev_1         6330
    ## surftemp_sjdate_1            5523
    ## n2o_mono3_cat                1855
    ## n2o_mono3_cat:log_area       2821
    ## n2o_mono3_cat:surftemp       1621
    ## shape_n2o_mono3_cat          6283
    ## 
    ## Simplex Parameters: 
    ##                            Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS
    ## n2o_mono3_cat1[1]             0.025     0.023    0.001    0.083 1.000     4899
    ## n2o_mono3_cat1[2]             0.183     0.111    0.014    0.425 1.002     1056
    ## n2o_mono3_cat1[3]             0.408     0.158    0.103    0.725 1.002     2167
    ## n2o_mono3_cat1[4]             0.384     0.169    0.058    0.713 1.003     1608
    ## n2o_mono3_cat:log_area1[1]    0.046     0.029    0.003    0.112 1.000     4018
    ## n2o_mono3_cat:log_area1[2]    0.033     0.028    0.001    0.104 1.001     6173
    ## n2o_mono3_cat:log_area1[3]    0.289     0.134    0.057    0.598 1.000     3286
    ## n2o_mono3_cat:log_area1[4]    0.631     0.140    0.301    0.864 1.000     3186
    ## n2o_mono3_cat:surftemp1[1]    0.028     0.019    0.003    0.060 1.001     3744
    ## n2o_mono3_cat:surftemp1[2]    0.066     0.033    0.011    0.128 1.001     3821
    ## n2o_mono3_cat:surftemp1[3]    0.281     0.079    0.122    0.431 1.003     2897
    ## n2o_mono3_cat:surftemp1[4]    0.625     0.088    0.461    0.798 1.003     2740
    ## shape_n2o_mono3_cat1[1]       0.116     0.067    0.009    0.263 1.001     5363
    ## shape_n2o_mono3_cat1[2]       0.149     0.096    0.008    0.365 1.001     3469
    ## shape_n2o_mono3_cat1[3]       0.401     0.151    0.116    0.703 1.001     4122
    ## shape_n2o_mono3_cat1[4]       0.334     0.143    0.051    0.604 1.001     3853
    ##                            Tail_ESS
    ## n2o_mono3_cat1[1]              5237
    ## n2o_mono3_cat1[2]              3343
    ## n2o_mono3_cat1[3]              4583
    ## n2o_mono3_cat1[4]              3021
    ## n2o_mono3_cat:log_area1[1]     3309
    ## n2o_mono3_cat:log_area1[2]     4948
    ## n2o_mono3_cat:log_area1[3]     4461
    ## n2o_mono3_cat:log_area1[4]     4184
    ## n2o_mono3_cat:surftemp1[1]     3410
    ## n2o_mono3_cat:surftemp1[2]     3213
    ## n2o_mono3_cat:surftemp1[3]     2971
    ## n2o_mono3_cat:surftemp1[4]     3106
    ## shape_n2o_mono3_cat1[1]        4297
    ## shape_n2o_mono3_cat1[2]        3637
    ## shape_n2o_mono3_cat1[3]        5848
    ## shape_n2o_mono3_cat1[4]        4561
    ## 
    ## Family Specific Parameters: 
    ##             Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## disc_no3cat    1.000     0.000    1.000    1.000   NA       NA       NA
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

### 3.6.2 Model checks

Below, the same PPCs for dissolved and equilibrium N2O were used as
before.

#### 3.6.2.1 Dissolved N2O

The PPCs for dissolved N2O were similar to those for models 4 and 5
above.

<img src="NLA_N2O_models_files/figure-gfm/ppc_full_checks_mod_n2o6-1.png" style="display: block; margin: auto;" />

#### 3.6.2.2 Equilibrium N2O

The PPCs for eqilibrium N2O were also similar to the same checks in
models 4 and 5.

<img src="NLA_N2O_models_files/figure-gfm/ppc_full_checks_mod_n2oeq6-1.png" style="display: block; margin: auto;" />

#### 3.6.2.3 Bivariate

This model provided a reasonable representation of the bivariate
relationship between the two N2O responses.

<img src="NLA_N2O_models_files/figure-gfm/ppc_bv_check_mod_n2o6-1.png" style="display: block; margin: auto;" />

#### 3.6.2.4 Saturation

The saturation ratio PPCs below suggested similar behavior as with
models 4 and 5 above, but with perhaps slightly less bias in the
predictions for the proportion of undersaturated waterbodies and fewer
extreme predictions for the means and standard deviations. The observed
*vs.* predicted PPC also appears to have a better behaved variance and
no extreme predictions, compared to models 4 and 5.

<img src="NLA_N2O_models_files/figure-gfm/ppc_sat_check_mod_n2o6-1.png" style="display: block; margin: auto;" />

The plot below shows the same PPC, but for the second-vist data.

<img src="NLA_N2O_models_files/figure-gfm/ppc_sat_check_testdata_mod_n2o6-1.png" style="display: block; margin: auto;" />

#### 3.6.2.5 R-square

Below are estimates for the Bayesian
![R^2](https://latex.codecogs.com/svg.image?R%5E2 "R^2"), which were
largely similar for the N2O responses as with models 4 and 5 above. The
![R^2](https://latex.codecogs.com/svg.image?R%5E2 "R^2") for the surface
temperature response also suggested a good fit.

    ##       Estimate Est.Error  Q2.5 Q97.5
    ## R2n2o    0.646     0.059 0.503 0.731

    ##         Estimate Est.Error  Q2.5 Q97.5
    ## R2n2oeq    0.863     0.006 0.851 0.874

    ##            Estimate Est.Error  Q2.5 Q97.5
    ## R2surftemp    0.744      0.01 0.723 0.763

Below are the ![R^2](https://latex.codecogs.com/svg.image?R%5E2 "R^2")
estimates for the second-visit data. That these estimates were similar
to those for the data used to fit the model was encouraging and
suggested that the model may perform reasonably well out-of-sample.

    ##       Estimate Est.Error  Q2.5 Q97.5
    ## R2n2o    0.607     0.137 0.325  0.85

    ##         Estimate Est.Error Q2.5 Q97.5
    ## R2n2oeq    0.857     0.008 0.84 0.872

    ##            Estimate Est.Error  Q2.5 Q97.5
    ## R2surftemp     0.75     0.018 0.715 0.783

### 3.6.3 Covariate effects

#### 3.6.3.1 Dissolved N2O

The conditional effects plot for the covariate effects on dissolved N2O
suggested a similar effect of NO3 as in previous models, but with
interesting potential interactions between NO3 and lake area and NO3 and
surface temperature. The lake area effect was estimated to be larger and
more negative at the highest levels of NO3 and slightly negative at the
lowest level of NO3. The surface temperature effect was estimated to be
largest and positive at the highest level of NO3 and negative at the
lowest level of NO3.

<img src="NLA_N2O_models_files/figure-gfm/conditional_effects_mod_n2oF-1.png" style="display: block; margin: auto;" />

The estimated covariate effects on
![\\sigma](https://latex.codecogs.com/svg.image?%5Csigma "\sigma") for
dissolved N2O suggested a negative relationship with log(area) and a
positive relationship with NO3.

<img src="NLA_N2O_models_files/figure-gfm/conditional_effects_sigma_n2oF-1.png" style="display: block; margin: auto;" />

#### 3.6.3.2 Equilibrium N2O

The estimated covariate effects on equilibrium N2O remained largely the
same as estimated in the previous model.

<img src="NLA_N2O_models_files/figure-gfm/conditional_effects_mod_n2oeqF-1.png" style="display: block; margin: auto;" />

<img src="NLA_N2O_models_files/figure-gfm/conditional_effects_sigma_n2oeqF-1.png" style="display: block; margin: auto;" />

# 4 Predict to population

As previously described, in order to make inferences to the population
of interest, the final model above was used to, first, predict surface
temperature in the target population, since it depended only on the
fully observed covariates. Next, the predictive distribution for surface
temperature was used, along with the relevant fully observed covariates,
to predict NO3 in the target population. Finally, the predictive
distributions for temperature and NO3 were used to predict the N2O
responses. The code for these steps is outlined in the following.

The first step used the final model to predict surface temperature to
the population:

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/sframe.rda")
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/n2o_mod6.rda")

predict_temp <- sframe %>%
  mutate(jdate = 205) %>%
  add_predicted_draws(n2o_mod6, resp=c("surftemp"), 
                      allow_new_levels = TRUE, 
                      cores =1, 
                      ndraws = 500) %>%
  mutate(surftemp = .prediction)

save(predict_temp, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/predict_temp.rda")
```

NO3 was next predicted using the posterior predictions for surface
temperature and the other fully observed covariates. Note that the
posterior predictive distribution for NO3 was subsampled in order to
minimize excess simulations.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/n2o_mod6.rda")
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/predict_temp.rda")

temp_X <- predict_temp %>% # select relevant columns as predictors
  ungroup() %>%
  select(WSA9,
         state,
         size_cat,
         log_area,
         .row,
         .draw,
         surftemp) %>%
  select(WSA9, state, size_cat, log_area, surftemp)


rm(predict_temp) # reduce memory
gc()

# set number of cores to use for parallel predictions
# and register the workers
cl <- parallel::makeCluster(5) 
doSNOW::registerDoSNOW(cl) 

# make a progress bar
pb <- txtProgressBar(max = 1500, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

system.time( # approx 26 hrs with 5 workers & 500 draws from PPD
predict_no3 <- foreach(sub_X = isplitRows(temp_X, chunkSize = 155299), 
                       .combine = 'c',
                       .packages = c("brms"),
                       .options.snow = opts
                       ) %dopar% {
                         apply(brms::posterior_predict(n2o_mod6,
                                                 newdata = sub_X,
                                                 resp = "no3cat",
                                                 allow_new_levels = T,
                                                 ndraws = 500,
                                                 cores = 1), 2, sample, 1)
                         }
)


close(pb)
parallel::stopCluster(cl)

save(predict_no3, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/predict_no3.rda")
```

Finally, dissolved and equilibrium N2O were predicted using the surface
temperature and NO3 predictions along with the survey variables and
known covariates. Again, the posterior was subsampled in order to reduce
excess simulations.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/predict_no3.rda")
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/predict_temp.rda")

# Assemble dataframe containing relevant covariates (known and predicted)
n2o_X <- predict_temp %>%
  ungroup() %>%
  mutate(no3_cat = predict_no3) %>%
  select(WSA9,
         state,
         size_cat,
         log_area,
         surftemp,
         log_elev,
         no3_cat)

# clear objects to reduce memory overhead
rm(predict_no3, predict_temp) 
gc()

# save the predictors for n2o and n2oeq
save(n2o_X, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/n2o_X.rda")
```

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/n2o_mod6.rda")
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/n2o_X.rda")

# set number of cores to use for parallel predictions
# and register the workers
cl <- parallel::makeCluster(6) 
doSNOW::registerDoSNOW(cl) 

# make a progress bar
pb <- txtProgressBar(max = 1500, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# make predictions in parallel
system.time(
predict_n2o <- foreach(sub_X = isplitRows(n2o_X, chunkSize = 155299),
                 .combine = rbind,
                 .options.snow = opts,
                 .packages = c("brms")) %dopar% {
  apply(posterior_predict(n2o_mod6,
                          newdata = sub_X,
                          resp = c("n2o", "n2oeq"),
                          allow_new_levels = T,
                          ndraws = 500,
                          cores = 1),
        2, sample, 1)
                   }
)

close(pb)
parallel::stopCluster(cl)

colnames(predict_n2o) <- c("n2o", "n2oeq")

save(predict_n2o, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/predict_n2o.rda")
```

Finally, the predictions for all four partially observed responses were
assembled into a new dataframe for use in inference:

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/predict_n2o.rda")
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/predict_no3.rda")
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/predict_temp.rda")

all_predictions <- predict_temp %>%
  ungroup() %>%
  mutate(no3cat = predict_no3) %>%
  bind_cols(predict_n2o) %>%
  mutate(n2osat = n2o / n2oeq, # calculate saturation ratio
         .row = rep(1:465897, each = 500),
         .draw = rep(seq(1,500, 1), 465897)) %>%
  mutate(area_ha = exp(log_area)) %>% # include area on ha scale
  select(WSA9,
         state,
         size_cat,
         area_ha,
         lat,
         lon,
         .row,
         .draw,
         surftemp,
         no3cat,
         n2o,
         n2oeq,
         n2osat)

rm(predict_n2o, predict_temp, predict_no3) # clean up workspace for RAM
gc()
 

save(all_predictions, file = "C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/all_predictions.rda")
```

# 5 Population estimates

A number of estimates for the target population were assembled and
presented below. First, the full posterior predictive distributions for
dissolved N2O, equilibrium N2O, and the saturation ratio were assessed.
These distributions summarized the predicted distribution of
concentrations or ratios for all lakes in the population of interest and
included parameter uncertainty propagated through the model. Next,
population means were assessed, followed by comparisons of some
model-based estimates to previously calculated design-based estimates.

## 5.1 Posterior predictive distributions

Below, a density plot summarized the posterior predictive distributions
across the target population of lakes. The PPDs consisted of based on
500 simulations for each variable. Note that the x-axis was truncated at
50 nmol/L for a clearer visualization of the bulk of the predictive
distribution. For reference, the max predicted value was 4403.2 nmol/L
for dissolved N2O, 20.4 nmol/L for dissolved N2O, and 793.5 for the
saturation ratio.

<img src="NLA_N2O_models_files/figure-gfm/plot_n2o_posterior_preds-1.png" style="display: block; margin: auto;" />

<img src="NLA_N2O_models_files/figure-gfm/plot_n2oeq_posterior_preds-1.png" style="display: block; margin: auto;" />

<img src="NLA_N2O_models_files/figure-gfm/plot_sat_posterior_preds-1.png" style="display: block; margin: auto;" />

## 5.2 Estimated means

### 5.2.1 National

Below are density plots summarizing the posterior distribution of
*means* for N2O concentrations and the saturation ratio for all US lakes
and reservoirs \> 4ha in the lower 48 states).

<img src="NLA_N2O_models_files/figure-gfm/plot_n2o_nat_posterior_mean-1.png" style="display: block; margin: auto;" />

<img src="NLA_N2O_models_files/figure-gfm/plot_n2oeq_nat_posterior_mean-1.png" style="display: block; margin: auto;" />

<img src="NLA_N2O_models_files/figure-gfm/plot_sat_nat_posterior_mean-1.png" style="display: block; margin: auto;" />

To illustrate the skewness in the predictive distribution for the
saturation ratio, an estimate for the median ratio is shown below. The
entire posterior distribution of the mean above was larger than 1, which
represents the boundary of under- *vs.* oversaturation. By comparison,
the posterior estimate of the median below only included values less
than one, suggesting that though the mean saturation ratio was likely
greater than 1, most lakes in the national population were
undersaturated (i.e., ratio less than 1). In distributions with strong
right-skew, the mean can often be considerably larger than the median.

<img src="NLA_N2O_models_files/figure-gfm/plot_sat_nat_posterior_median-1.png" style="display: block; margin: auto;" />

Below is a plot of the posterior mean estimate for the proportion of
unsaturated lakes at the national scale.

<img src="NLA_N2O_models_files/figure-gfm/plot_undersat_posterior_mean-1.png" style="display: block; margin: auto;" />

### 5.2.2 Ecoregion

Below are posterior estimates of the means for dissolved and equilibrium
N2O and the saturation ratio by WSA9 ecoregion.

<img src="NLA_N2O_models_files/figure-gfm/plot_n2o_wsa9_posterior_mean-1.png" style="display: block; margin: auto;" />

<img src="NLA_N2O_models_files/figure-gfm/plot_n2oeq_wsa9_posterior_mean-1.png" style="display: block; margin: auto;" />

<img src="NLA_N2O_models_files/figure-gfm/plot_sat_wsa9_posterior_mean-1.png" style="display: block; margin: auto;" />

A plot of the posterior estimates for the median saturation ratio below
indicated, again, that most lakes in each ecoregion were undersaturated
(i.e., median \<\< 1).

<img src="NLA_N2O_models_files/figure-gfm/plot_sat_wsa9_posterior_median-1.png" style="display: block; margin: auto;" />

A plot of the estimates of the proportion of under-saturated lakes by
ecoregion is below. These summaries again suggested that most lakes in
each ecoregion were likely undersaturated (i.e., median \<\< 1).

<img src="NLA_N2O_models_files/figure-gfm/plot_prop_sat_wsa9_posterior_median-1.png" style="display: block; margin: auto;" />

### 5.2.3 State

Comparisons of mean estimates by state are below. Density estimates were
not included to minimize the vertical plot space.

<img src="NLA_N2O_models_files/figure-gfm/plot_state_mean_n2o-1.png" style="display: block; margin: auto;" />

<img src="NLA_N2O_models_files/figure-gfm/plot_state_mean_n2oeq-1.png" style="display: block; margin: auto;" />

Below is a plot of estimates for the mean (black circles) and median
(grey circles) saturation ratio by state. A horizontal, dashed, black
line is shown at ratio = 1, indicating the boundary for under- *vs.*
oversaturation. Only a few states (e.g., NV, DE) had median estimates
that were 1 or greater, suggesting that, for most states, most lakes
were undersaturated.

<img src="NLA_N2O_models_files/figure-gfm/plot_state_mean_median_sat-1.png" style="display: block; margin: auto;" />

Finally, a plot of the estimated proportion of undersaturated lakes for
each state in the target population. Point estimates are the posterior
median of the proportion and bars are the upper and lower boundaries of
the central 95th percentile of the posterior distributions of
proportions.

<img src="NLA_N2O_models_files/figure-gfm/plot_state_prop_sat-1.png" style="display: block; margin: auto;" />

### 5.2.4 Size category

The estimated means and medians by size category are below for dissolved
and equilibrium N2O and the saturation ratio.

<img src="NLA_N2O_models_files/figure-gfm/plot_n2o_size_posterior_mean-1.png" style="display: block; margin: auto;" />

<img src="NLA_N2O_models_files/figure-gfm/plot_n2oeq_size_posterior_mean-1.png" style="display: block; margin: auto;" />

<img src="NLA_N2O_models_files/figure-gfm/plot_sat_size_posterior_mean-1.png" style="display: block; margin: auto;" />

<img src="NLA_N2O_models_files/figure-gfm/plot_sat_size_posterior_median-1.png" style="display: block; margin: auto;" />

Mean *vs.* median below.

<img src="NLA_N2O_models_files/figure-gfm/plot_size_cat_mean_median_sat-1.png" style="display: block; margin: auto;" />

And, finally, the estimated proportion of undersaturated lakes in the
target population by size category.

<img src="NLA_N2O_models_files/figure-gfm/plot_prop_sat_size_cat_posterior_median-1.png" style="display: block; margin: auto;" />

## 5.3 Model- *vs.* design-based

Below, estimates from the model-based approach are compared to
previously calculated, design-based estimates. In general, the
model-based estimates were similar to the design-based estimates. The
model-based estimates were typically within the confidence bounds of the
design-based estimates, but with much greater precision. Improved
precision was expected due to the “shrinkage” induced by the multilevel
parameterization, which affords some “borrowing” of information across
the various levels of the survey factors.

### 5.3.1 Dissolved N2O

Below, National mean estimates for dissolved N2O from the model and
design-based approaches were compared. The sample-based estimate was
also included as a reference. The black, vertical, dashed line
represents the mean of the sample.

``` r
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/df_model.rda")
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/n2o_survey_ests.rda")
#load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/all_predictions.rda")

all_predictions %>%
  group_by(.draw) %>%
  summarise(mean_n2o = mean(n2o)) %>%
  summarise(estimate = round(median(mean_n2o), 2), # posterior median
    LCL = round(quantile(mean_n2o, probs = 0.025), 2),
    UCL = round(quantile(mean_n2o, probs = 0.975), 2)) %>% 
  mutate(type = "model") %>%
  bind_rows(cbind(n2o_survey_ests[10, 2:4], type = rep("survey", 1))) %>%
  add_row(estimate = round(mean(df_model$n2o), 2),
          type = "sample") %>%
  print()
```

    ## # A tibble: 3 x 4
    ##   estimate   LCL   UCL type  
    ##      <dbl> <dbl> <dbl> <chr> 
    ## 1     7.51  7.36  7.63 model 
    ## 2     8.1   7     9.1  survey
    ## 3     8.72 NA    NA    sample

<img src="NLA_N2O_models_files/figure-gfm/plot_n2o_means_national-1.png" style="display: block; margin: auto;" />

Below, estimates were compared by ecoregion.

``` r
#load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/all_predictions.rda")
#load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/n2o_survey_ests.rda")

all_predictions %>%
  group_by(WSA9, .draw) %>%
  summarise(mean_n2o = mean(n2o)) %>%
  group_by(WSA9, .groups = "drop") %>%
  summarise(estimate = round(median(mean_n2o), 2),
    LCL = round(quantile(mean_n2o, probs = 0.025), 2),
    UCL = round(quantile(mean_n2o, probs = 0.975), 2),
    .groups = "drop") %>% 
  mutate(ecoregion = factor(WSA9)) %>%
  mutate(type = "model") %>%
  select(ecoregion, estimate, LCL, UCL, type) %>%
  mutate(ecoregion = forcats::fct_reorder(ecoregion, estimate)) %>%
  bind_rows(cbind(n2o_survey_ests[-10,], type = rep("survey", 9))) %>%
  arrange(ecoregion) %>%
  print()
```

    ## # A tibble: 18 x 5
    ##    ecoregion estimate   LCL   UCL type  
    ##    <fct>        <dbl> <dbl> <dbl> <chr> 
    ##  1 NPL           6.95  6.8   7.1  model 
    ##  2 NPL           6.9   6.4   7.4  survey
    ##  3 SPL           7.07  6.84  7.31 model 
    ##  4 SPL           6.5   4.9   8.1  survey
    ##  5 NAP           7.33  7.15  7.5  model 
    ##  6 NAP           7.7   7.2   8.1  survey
    ##  7 CPL           7.49  7.23  7.7  model 
    ##  8 CPL           8.4   5.1  11.7  survey
    ##  9 UMW           7.57  7.43  7.72 model 
    ## 10 UMW          10.8   6.4  15.2  survey
    ## 11 SAP           7.62  7.37  7.85 model 
    ## 12 SAP           7.1   5.9   8.2  survey
    ## 13 XER           7.65  7.41  7.85 model 
    ## 14 XER          10.6   7.5  13.7  survey
    ## 15 WMT           7.8   7.66  7.93 model 
    ## 16 WMT           7.8   7.1   8.4  survey
    ## 17 TPL           7.96  7.76  8.18 model 
    ## 18 TPL           7.8   5.9   9.6  survey

<img src="NLA_N2O_models_files/figure-gfm/plot_mean_n2o_wsa9-1.png" style="display: block; margin: auto;" />

Means were compared according to size categories below.

    ## # A tibble: 10 x 5
    ##    size   estimate   LCL   UCL type  
    ##    <ord>     <dbl> <dbl> <dbl> <chr> 
    ##  1 4_10        7.5   7.4   7.7 model 
    ##  2 4_10        7.6   6.5   8.8 survey
    ##  3 10_20       7.5   7.3   7.6 model 
    ##  4 10_20       7.6   7.1   8.1 survey
    ##  5 50_max      7.5   7.4   7.6 model 
    ##  6 50_max      8     7.4   8.5 survey
    ##  7 min_4       7.5   7.3   7.6 model 
    ##  8 min_4       8.2   6.4   9.9 survey
    ##  9 20_50       7.6   7.4   7.7 model 
    ## 10 20_50       8.6   7.7   9.5 survey

<img src="NLA_N2O_models_files/figure-gfm/plot_size_mean_n2o-1.png" style="display: block; margin: auto;" />

### 5.3.2 Saturation

Below, the same comparisons were made for the saturation estimates.

``` r
#load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/df_model.rda")
#load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/all_predictions.rda")
load("C:/Users/rmartin/OneDrive - Environmental Protection Agency (EPA)/Documents/AE_Reservoirs/DissolvedGasNla/modelFiles/sat_survey_ests.rda")

all_predictions %>%
  group_by(.draw) %>%
  summarise(mean_sat = mean(n2osat), .groups = "drop") %>%
  summarise(estimate = round(median(mean_sat), 3),
    LCL = round(quantile(mean_sat, probs = 0.025), 3),
    UCL = round(quantile(mean_sat, probs = 0.975), 3),
    .groups = "drop") %>% 
  mutate(type = "model") %>%
  bind_rows(cbind(sat_survey_ests[10, 2:4], type = rep("survey", 1))) %>%
  add_row(estimate = round(mean(df_model$n2o / df_model$n2o_eq), 3),
          type = "sample") %>%
  print()
```

    ## # A tibble: 3 x 4
    ##   estimate    LCL   UCL type  
    ##      <dbl>  <dbl> <dbl> <chr> 
    ## 1     1.11  1.10   1.11 model 
    ## 2     1.10  0.952  1.24 survey
    ## 3     1.17 NA     NA    sample

<img src="NLA_N2O_models_files/figure-gfm/plot_nat_sat_mean-1.png" style="display: block; margin: auto;" />

<img src="NLA_N2O_models_files/figure-gfm/plot_wsa9_sat_mean-1.png" style="display: block; margin: auto;" />

    ## # A tibble: 10 x 5
    ##    size   estimate   LCL   UCL type  
    ##    <ord>     <dbl> <dbl> <dbl> <chr> 
    ##  1 min_4      1.11 1.11   1.12 model 
    ##  2 min_4      1.12 0.874  1.37 survey
    ##  3 4_10       1.09 1.08   1.11 model 
    ##  4 4_10       1.02 0.889  1.15 survey
    ##  5 10_20      1.07 1.06   1.09 model 
    ##  6 10_20      1.02 0.956  1.08 survey
    ##  7 20_50      1.07 1.06   1.09 model 
    ##  8 20_50      1.14 1.02   1.27 survey
    ##  9 50_max     1.06 1.05   1.07 model 
    ## 10 50_max     1.06 0.987  1.12 survey

<img src="NLA_N2O_models_files/figure-gfm/plot_size_sat_mean-1.png" style="display: block; margin: auto;" />

# 6 Session Info

``` r
sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 22000)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] brms_2.18.0      Rcpp_1.0.9       tidybayes_3.0.2  bayesplot_1.9.0 
    ##  [5] itertools_0.1-3  iterators_1.0.14 foreach_1.5.2    future_1.28.0   
    ##  [9] forcats_0.5.2    stringr_1.5.0    purrr_0.3.5      readr_2.1.3     
    ## [13] tidyr_1.2.1      tibble_3.1.8     tidyverse_1.3.2  dplyr_1.0.10    
    ## [17] ggrepel_0.9.2    kableExtra_1.3.4 gridExtra_2.3    ggExtra_0.10.0  
    ## [21] moments_0.14.1   ggpubr_0.4.0     ggplot2_3.4.0   
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.4.1         backports_1.4.1      systemfonts_1.0.4   
    ##   [4] plyr_1.8.8           igraph_1.3.5         splines_4.1.2       
    ##   [7] svUnit_1.0.6         crosstalk_1.2.0      listenv_0.8.0       
    ##  [10] inline_0.3.19        rstantools_2.2.0     digest_0.6.31       
    ##  [13] htmltools_0.5.4      fansi_1.0.3          magrittr_2.0.3      
    ##  [16] checkmate_2.1.0      googlesheets4_1.0.1  tzdb_0.3.0          
    ##  [19] globals_0.16.1       modelr_0.1.9         RcppParallel_5.1.5  
    ##  [22] matrixStats_0.62.0   xts_0.12.2           svglite_2.1.0       
    ##  [25] timechange_0.1.1     prettyunits_1.1.1    colorspace_2.0-3    
    ##  [28] rvest_1.0.3          ggdist_3.2.0         haven_2.5.1         
    ##  [31] xfun_0.35            callr_3.7.3          crayon_1.5.2        
    ##  [34] jsonlite_1.8.4       zoo_1.8-11           glue_1.6.2          
    ##  [37] gtable_0.3.1         gargle_1.2.1         webshot_0.5.4       
    ##  [40] V8_4.2.1             distributional_0.3.1 pkgbuild_1.3.1      
    ##  [43] car_3.1-1            rstan_2.26.11        abind_1.4-5         
    ##  [46] scales_1.2.1         mvtnorm_1.1-3        DBI_1.1.3           
    ##  [49] rstatix_0.7.0        miniUI_0.1.1.1       viridisLite_0.4.1   
    ##  [52] xtable_1.8-4         diffobj_0.3.5        StanHeaders_2.26.11 
    ##  [55] stats4_4.1.2         DT_0.26              htmlwidgets_1.6.0   
    ##  [58] httr_1.4.4           threejs_0.3.3        arrayhelpers_1.1-0  
    ##  [61] posterior_1.3.1      ellipsis_0.3.2       pkgconfig_2.0.3     
    ##  [64] loo_2.5.1            farver_2.1.1         dbplyr_2.2.1        
    ##  [67] utf8_1.2.2           labeling_0.4.2       tidyselect_1.2.0    
    ##  [70] rlang_1.0.6          reshape2_1.4.4       later_1.3.0         
    ##  [73] munsell_0.5.0        cellranger_1.1.0     tools_4.1.2         
    ##  [76] cli_3.4.1            generics_0.1.3       broom_1.0.1         
    ##  [79] ggridges_0.5.4       evaluate_0.19        fastmap_1.1.0       
    ##  [82] yaml_2.3.6           processx_3.8.0       knitr_1.41          
    ##  [85] fs_1.5.2             nlme_3.1-161         mime_0.12           
    ##  [88] xml2_1.3.3           shinythemes_1.2.0    compiler_4.1.2      
    ##  [91] rstudioapi_0.14      curl_4.3.3           ggsignif_0.6.3      
    ##  [94] reprex_2.0.2         stringi_1.7.8        highr_0.9           
    ##  [97] ps_1.7.2             Brobdingnag_1.2-9    lattice_0.20-45     
    ## [100] Matrix_1.5-3         markdown_1.1         shinyjs_2.1.0       
    ## [103] tensorA_0.36.2       vctrs_0.5.1          pillar_1.8.1        
    ## [106] lifecycle_1.0.3      bridgesampling_1.1-2 httpuv_1.6.6        
    ## [109] R6_2.5.1             promises_1.2.0.1     parallelly_1.32.1   
    ## [112] codetools_0.2-18     colourpicker_1.1.1   gtools_3.9.4        
    ## [115] assertthat_0.2.1     withr_2.5.0          shinystan_2.6.0     
    ## [118] mgcv_1.8-41          parallel_4.1.2       hms_1.1.2           
    ## [121] grid_4.1.2           coda_0.19-4          rmarkdown_2.19      
    ## [124] carData_3.0-5        googledrive_2.0.0    shiny_1.7.1         
    ## [127] lubridate_1.9.0      base64enc_0.1-3      dygraphs_1.1.1.6

# 7 References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Burkner_2017" class="csl-entry">

Bürkner, Paul-Christian. 2017. “Brms: An r Package for Bayesian
Multilevel Models Using Stan.” Journal Article. *2017* 80 (1): 28.
<https://doi.org/10.18637/jss.v080.i01>.

</div>

<div id="ref-Gelman_etal_2014" class="csl-entry">

Gelman, Andrew, John B. Carlin, Hal S. Stern, David B. Dunson, Aki
Vehtari, and Donald B. Rubin. 2014. *Bayesian Data Analysis*. Book. 4th
ed. New York: CRC Press.

</div>

<div id="ref-Gelman_etal_2020" class="csl-entry">

Gelman, Andrew, Jennifer Hill, and Aki. Vehtari. 2020. *Regression and
Other Stories: Analytical Methods for Social Research*. Book. 1st ed.
Cambridge: CRC Press.

</div>

<div id="ref-Gelman_Little_1997" class="csl-entry">

Gelman, Andrew, and Thomas Little. 1997. “Postratification into Many
Categories Using Hierarchical Logistic Regression.” Journal Article.
*Survey Methodology* 23 (2): 127–35.

</div>

<div id="ref-Kennedy_Gelman_2021" class="csl-entry">

Kennedy, Lauren, and Andrew Gelman. 2021. “Know Your Population and Know
Your Model: Using Model-Based Regression and Poststratification to
Generalize Findings Beyond the Observed Sample.” Journal Article.
*Psychological Methods* 26 (5): 547–58.
<https://psycnet.apa.org/doi/10.1037/met0000362>.

</div>

<div id="ref-Link_Barker_2010" class="csl-entry">

Link, William A., and Richard J. Barker. 2010. *Bayesian Inference: With
Ecological Applications*. Book. 1st ed. Boston: Academic Press.

</div>

<div id="ref-McElreath_2020" class="csl-entry">

McElreath, Richard. 2020. *Statistical Rethinking: A Bayesian Course
with Examples in r and Stan*. Book. 2nd ed. Boca Raton: CRC Press.

</div>

<div id="ref-Merkle_etal_2021" class="csl-entry">

Merkle, Edgar C., Ellen Fitzsimmons, James Uanhoro, and Ben Goodrich.
2021. “Efficient Bayesian Structural Equation Modeling in Stan.” Journal
Article. *Journal of Statistical Software* 100 (6): 1–22.
https://doi.org/<https://doi.org/10.18637/jss.v100.i06>.

</div>

<div id="ref-Merkle_Rosseel_2018" class="csl-entry">

Merkle, Edgar C., and Yves Rosseel. 2018. “Blavaan: Bayesian Structural
Equation Models via Parameter Expansion.” Journal Article. *Journal of
Statistical Software* 85 (4): 1–30.
https://doi.org/<https://doi.org/10.18637/jss.v085.i04>.

</div>

<div id="ref-Park_etal_2004" class="csl-entry">

Park, David K., Andrew Gelman, and Joseph Bafumi. 2004. “Bayesian
Multilevel Estimation with Poststratification: State-Level Estimates
from National Polls.” Journal Article. *Political Analysis* 12 (4):
375–85. <https://doi.org/10.1093/pan/mph024>.

</div>

<div id="ref-Poggiato_etal_2021" class="csl-entry">

Poggiato, Giovanni, Tamara Munkemuller, Daria Bystrova, Julyan Arbel,
James S. Clark, and Wilfried Thuiller. 2021. “On the Interpretations of
Joint Modeling in Community Ecology.” Journal Article. *Trends in
Ecology & Evolution* 36 (5): 391–401.
https://doi.org/<https://doi.org/10.1016/j.tree.2021.01.002>.

</div>

<div id="ref-R_Core_Team_2021" class="csl-entry">

R Core Team. 2021. *R: A Language and Environment for Statistical
Computing*. Vienna, Austria: R Foundation for Statistical Computing.
<https://www.R-project.org/>.

</div>

<div id="ref-Stan_Development_Team_2018_c" class="csl-entry">

Stan Development Team. 2018a. “RStan: The r Interface to Stan.” Journal
Article. <http://mc-stan.org>.

</div>

<div id="ref-Stan_Development_Team_2018_a" class="csl-entry">

———. 2018b. “Stan Modeling Language Users Guide and Reference Manual,
Version 2.18.0.” Journal Article. <http://mc-stan.org>.

</div>

<div id="ref-Stan_Development_Team_2018_b" class="csl-entry">

———. 2018c. “The Stan Core Library, Version 2.18.0.” Journal Article.
<http://mc-stan.org>.

</div>

<div id="ref-Warton_etal_2015" class="csl-entry">

Warton, David I., Guillame F. Blanchet, Robert B. O’Hara, Otso
Ovaskainen, Sara Taskinen, Steven C. Walker, and Francis K. C. Hui.
2015. “So Many Variables: Joint Modeling in Community Ecology.” Journal
Article. *Trends in Ecology & Evolution* 30 (12): 766–79.
https://doi.org/<https://doi.org/10.1016/j.tree.2015.09.007.>

</div>

<div id="ref-Webb_etal_2019" class="csl-entry">

Webb, Jackie R., Nicole M. Hayes, Gavin L. Simpson, Peter R. Leavitt,
Helen M. Baulch, and Kerri Finlay. 2019. “Widespread Nitrous Oxide
Undersaturation in Farm Waterbodies Creates an Unexpected Greenhouse Gas
Sink.” Journal Article. *Proceedings of the National Academy of
Sciences* 116 (20): 9814–19. <https://doi.org/10.1073/pnas.1820389116>.

</div>

<div id="ref-Zachmann_etal_2022" class="csl-entry">

Zachmann, Luke J., Erin M. Borgman, Dana L. Witwicki, Megan C. Swan,
Cheryl McIntyre, and N. Thompson Hobbs. 2022. “Bayesian Models for
Analysis of Inventory and Monitoring Data with Non-Ignorable
Missingness.” Journal Article. *Journal of Agricultural, Biological and
Environmental Statistics* 27: 125–48.
https://doi.org/<https://doi.org/10.1007/s13253-021-00473-z>.

</div>

</div>
