Design-based estimates for N2O from the NLA 2017 Survey
================
Roy Martin, Jake Beaulieu, Michael McManus
2025-03-05

# 1 Purpose

Calculate estimates for the NLA 2017 survey with spsurvey and local mean
variance estimator for comparison to Bayesian MRP Estimates. Tom Kincaid
created a .csv file to be used for population estimates, which is loaded
and munged below.

# 2 Import and munge survey data

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
    ## 
    ## Attaching package: 'Hmisc'
    ## 
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     src, summarize
    ## 
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     format.pval, units
    ## 
    ## 
    ## Loading required package: sf
    ## 
    ## Linking to GEOS 3.10.2, GDAL 3.4.1, PROJ 8.2.1; sf_use_s2() is TRUE
    ## 
    ## Loading required package: survey
    ## 
    ## Loading required package: grid
    ## 
    ## Loading required package: Matrix
    ## 
    ## 
    ## Attaching package: 'Matrix'
    ## 
    ## 
    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack
    ## 
    ## 
    ## Loading required package: survival
    ## 
    ## 
    ## Attaching package: 'survey'
    ## 
    ## 
    ## The following object is masked from 'package:Hmisc':
    ## 
    ##     deff
    ## 
    ## 
    ## The following object is masked from 'package:graphics':
    ## 
    ##     dotchart
    ## 
    ## 
    ## spsurvey version 5.0.0 introduced significant changes to the inputs and outputs of many functions. Please review the updated materials, vignettes, and documentation by visiting 
    ##  https://cran.r-project.org/package=spsurvey

``` r
# modified file to include N2O emission rates.  See "Users\JBEAULIE\OneDrive - Environmental Protection Agency (EPA)\gitRepository\DissolvedGasNla\scripts\dataMunge.Rmd" for details
# updated file path
df_survey <- read.csv(file = "./../inputData/populationEstimates/NLA17_DissolvedGases_forPopEst.csv")
df_survey <- df_survey %>% 
  data.frame() %>%
  filter(!is.na(DISSOLVED_N2O))
nr <- nrow(df_survey) # 984

# Reorder levels for the lake size class (AREA_CAT6) variable
levels(df_survey$AREA_CAT6) <- list(
   "(1,4]" = "(1,4]",
   "(4,10]" = "(4,10]",
   "(10,20]" = "(10,20]",
   "(20,50]" = "(20,50]",
   ">50" = ">50")

head(df_survey)
```

    ## # A tibble: 6 × 81
    ##      UID SITE_ID       UNIQUE_ID AG_ECO3 AG_ECO3_NM AG_ECO9 AG_ECO9_NM AREA_CAT6
    ##    <int> <chr>         <chr>     <chr>   <chr>      <chr>   <chr>      <chr>    
    ## 1 225890 NLA17_NV-100… NLA_NV-1… WMTNS   West       XER     Xeric      (4,10]   
    ## 2 225930 NLA17_NV-100… NLA_NV-1… WMTNS   West       XER     Xeric      (1,4]    
    ## 3 225970 NLA17_ND-100… NLA_ND-1… PLNLOW  Plains an… NPL     Northern … (20,50]  
    ## 4 226000 NLA17_NV-102… NLA_NV-1… WMTNS   West       XER     Xeric      >50      
    ## 5 226021 NLA17_OH-100… NLA_OH-1… EHIGH   Eastern H… NAP     Northern … (1,4]    
    ## 6 226022 NLA17_TN-100… NLA_TN-1… EHIGH   Eastern H… SAP     Southern … (20,50]  
    ## # ℹ 73 more variables: AREA_HA <dbl>, BORD_LAKE <chr>, CNTYNAME <chr>,
    ## #   COMID <int>, DES_FTYPE <chr>, ELEVATION <int>, EPA_REG <chr>,
    ## #   EVAL_CAT <chr>, FCODE <int>, FEOW_ID <int>, FRAME07 <chr>, FRAME12 <chr>,
    ## #   FRAME17 <chr>, FRAME17_ID <int>, FS_EW <chr>, FTYPE <chr>, GNIS_ID <int>,
    ## #   GNIS_NAME <chr>, HUC2 <chr>, HUC2_NM <chr>, HUC8 <chr>, HUC8_NM <chr>,
    ## #   LAKE_ORGN <chr>, LAT_DD83 <dbl>, LON_DD83 <dbl>, MAJ_BAS_NM <chr>,
    ## #   MAJ_BASIN <chr>, NA_L1CODE <int>, NA_L1NAME <chr>, NA_L2CODE <dbl>, …

``` r
tail(df_survey)
```

    ## # A tibble: 6 × 81
    ##       UID SITE_ID      UNIQUE_ID AG_ECO3 AG_ECO3_NM AG_ECO9 AG_ECO9_NM AREA_CAT6
    ##     <int> <chr>        <chr>     <chr>   <chr>      <chr>   <chr>      <chr>    
    ## 1 1000619 NLA17_ND-10… NLA_ND-1… PLNLOW  Plains an… TPL     Temperate… (20,50]  
    ## 2 1000623 NLA17_CA-10… NLA_CA-1… WMTNS   West       WMT     Western M… (10,20]  
    ## 3 1000624 NLA17_MD-10… NLA_MD-1… PLNLOW  Plains an… CPL     Coastal P… (1,4]    
    ## 4 1000626 NLA17_UT-10… NLA_UT-1… WMTNS   West       XER     Xeric      (4,10]   
    ## 5 1000628 NLA17_AZ-10… NLA_AZ-1… WMTNS   West       WMT     Western M… (20,50]  
    ## 6 1000629 NLA17_MS-10… NLA_MS-1… PLNLOW  Plains an… CPL     Coastal P… >50      
    ## # ℹ 73 more variables: AREA_HA <dbl>, BORD_LAKE <chr>, CNTYNAME <chr>,
    ## #   COMID <int>, DES_FTYPE <chr>, ELEVATION <int>, EPA_REG <chr>,
    ## #   EVAL_CAT <chr>, FCODE <int>, FEOW_ID <int>, FRAME07 <chr>, FRAME12 <chr>,
    ## #   FRAME17 <chr>, FRAME17_ID <int>, FS_EW <chr>, FTYPE <chr>, GNIS_ID <int>,
    ## #   GNIS_NAME <chr>, HUC2 <chr>, HUC2_NM <chr>, HUC8 <chr>, HUC8_NM <chr>,
    ## #   LAKE_ORGN <chr>, LAT_DD83 <dbl>, LON_DD83 <dbl>, MAJ_BAS_NM <chr>,
    ## #   MAJ_BASIN <chr>, NA_L1CODE <int>, NA_L1NAME <chr>, NA_L2CODE <dbl>, …

``` r
# Create data frame for spsurvey estimates
dframe <- data.frame(
  siteID = df_survey$SITE_ID,
  N2O_SAT_RATIO=df_survey$N2O_SAT_RATIO,
  DISSOLVED_N2O=df_survey$DISSOLVED_N2O,
  National=rep("National", nr),
  LAKE_ORGN=df_survey$LAKE_ORGN,
  AREA_CAT6=df_survey$AREA_CAT6,
  AG_ECO9_NM=df_survey$AG_ECO9_NM,
  State = df_survey$STRATUM,
  wgt = df_survey$WGT_TP,
  xcoord = df_survey$XCOORD,
  ycoord = df_survey$YCOORD)
```

# 3 Survey estimates with spsurvey package

Estimates made using the local mean variance estimator

# 4 Summarize estimates and export

    ## # A tibble: 65 × 8
    ##    Type    Subpopulation nResp Estimate StdError MarginofError LCB95Pct UCB95Pct
    ##    <chr>   <chr>         <int>    <dbl>    <dbl>         <dbl>    <dbl>    <dbl>
    ##  1 Nation… National        984     8.05     0.53          1.03     7.02     9.09
    ##  2 LAKE_O… MAN_MADE        517     7.73     0.68          1.34     6.39     9.07
    ##  3 LAKE_O… NATURAL         467     8.62     0.8           1.56     7.06    10.2 
    ##  4 AREA_C… min_4           198     8.17     0.9           1.76     6.41     9.92
    ##  5 AREA_C… 10_20           174     7.6      0.24          0.48     7.13     8.08
    ##  6 AREA_C… 20_50           189     8.57     0.46          0.9      7.67     9.47
    ##  7 AREA_C… 4_10            158     7.64     0.58          1.13     6.51     8.78
    ##  8 AREA_C… 50_max          265     7.95     0.26          0.51     7.44     8.47
    ##  9 AG_ECO… Coastal Plai…   122     8.4      1.67          3.27     5.13    11.7 
    ## 10 AG_ECO… Northern App…    94     7.65     0.22          0.44     7.21     8.09
    ## # ℹ 55 more rows

``` r
write.csv(mean_dissolved, file = "./../inputData/populationEstimates/Survey_Ests_N20_dissolved.csv")
```

    ## # A tibble: 65 × 8
    ##    Type    Subpopulation nResp Estimate StdError MarginofError LCB95Pct UCB95Pct
    ##    <chr>   <chr>         <int>    <dbl>    <dbl>         <dbl>    <dbl>    <dbl>
    ##  1 Nation… National        984     1.1      0.07          0.14     0.95     1.24
    ##  2 LAKE_O… MAN_MADE        517     1.1      0.1           0.2      0.9      1.3 
    ##  3 LAKE_O… NATURAL         467     1.09     0.09          0.18     0.91     1.27
    ##  4 AREA_C… min_4           198     1.12     0.13          0.25     0.87     1.37
    ##  5 AREA_C… 10_20           174     1.02     0.03          0.06     0.96     1.08
    ##  6 AREA_C… 20_50           189     1.14     0.06          0.13     1.02     1.27
    ##  7 AREA_C… 4_10            158     1.02     0.07          0.13     0.89     1.15
    ##  8 AREA_C… 50_max          265     1.05     0.03          0.07     0.99     1.12
    ##  9 AG_ECO… Coastal Plai…   122     1.27     0.26          0.51     0.76     1.77
    ## 10 AG_ECO… Northern App…    94     0.97     0.03          0.05     0.91     1.02
    ## # ℹ 55 more rows

``` r
write.csv(mean_sat_ratio, file = "./../inputData/populationEstimates/Survey_Ests_N20_sat.csv")
```

# 5 Session Info

``` r
sessionInfo()
```

    ## R version 4.4.0 (2024-04-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 22.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Etc/UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] spsurvey_5.5.1  survey_4.4-2    survival_3.5-8  Matrix_1.7-0   
    ##  [5] sf_1.0-16       Hmisc_5.2-2     lubridate_1.9.3 forcats_1.0.0  
    ##  [9] stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5    
    ## [13] tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.5       xfun_0.43          htmlwidgets_1.6.4  lattice_0.22-6    
    ##  [5] tzdb_0.4.0         sampling_2.10      Rdpack_2.6.2       vctrs_0.6.5       
    ##  [9] tools_4.4.0        generics_0.1.3     proxy_0.4-27       fansi_1.0.6       
    ## [13] AlgDesign_1.2.1.1  cluster_2.1.6      pkgconfig_2.0.3    KernSmooth_2.23-22
    ## [17] data.table_1.15.4  checkmate_2.3.2    lifecycle_1.0.4    deldir_2.0-4      
    ## [21] compiler_4.4.0     munsell_0.5.1      mitools_2.4        htmltools_0.5.8.1 
    ## [25] class_7.3-22       yaml_2.3.8         htmlTable_2.4.3    Formula_1.2-5     
    ## [29] nloptr_2.1.1       pillar_1.9.0       MASS_7.3-60.2      classInt_0.4-10   
    ## [33] reformulas_0.4.0   boot_1.3-30        rpart_4.1.23       nlme_3.1-164      
    ## [37] gtools_3.9.5       tidyselect_1.2.1   digest_0.6.35      stringi_1.8.3     
    ## [41] splines_4.4.0      fastmap_1.1.1      colorspace_2.1-0   cli_3.6.2         
    ## [45] magrittr_2.0.3     base64enc_0.1-3    utf8_1.2.4         e1071_1.7-14      
    ## [49] foreign_0.8-86     withr_3.0.0        scales_1.3.0       backports_1.4.1   
    ## [53] timechange_0.3.0   rmarkdown_2.26     lme4_1.1-36        nnet_7.3-19       
    ## [57] gridExtra_2.3      hms_1.1.3          lpSolve_5.6.23     evaluate_0.23     
    ## [61] knitr_1.46         rbibutils_2.3      crossdes_1.1-2     rlang_1.1.3       
    ## [65] Rcpp_1.0.12        glue_1.7.0         DBI_1.2.2          minqa_1.2.8       
    ## [69] rstudioapi_0.16.0  R6_2.5.1           units_0.8-5
