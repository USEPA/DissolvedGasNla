# DissolvedGasNla

The purpose of this repository is to provide the code required to replicate the analysis presented in
Beaulieu, J.J., Martin, R. and McManus, M. Pervasive nitrous oxide undersaturation in U.S. lakes and reservoirs.

# Users Guide

The entire manuscript can be reproduced by runing the scripts in the following order: 
1. scripts/dataMunge.Rmd
2. scripts/DG_sensitivity_to_measurement_error.Rmd
3. scripts/Design_Based_Estimates_N2O.Rmd
4. manuscript/manuscript_file.Rmd

These scripts depend on data contained in this repository and a 45GB dataset of predicted N2O for the population of interest. This dataset can be created by running scripts/NLA17_N2O_models.Rmd, but the computations are resource intensive, requiring 48 hours of run time and ~ 256GB RAM. Alternatively, the dataset can be downloaded from Zenodo (https://zenodo.org/records/15159394) using `utils::download.file` in manuscript/manuscript_support.R. Download time using a laptop and home internet connection providing download speeds up to 3 Gig was 29 hours.

`renv` is used for package management. After forking the repo, run `renv::restore` to recreate the package depency environment used for this project. The code has been tested using R version 4.4.1. 

# EPA Disclaimer
The United States Environmental Protection Agency (EPA) GitHub project code is provided on an "as is" basis and the user assumes responsibility for its use. EPA has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality, or availability of the information. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by EPA. The EPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by EPA or the United States Government. 
