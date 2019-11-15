library(rgdal) # spatial data
library(sf) # spatial data
library(gstat) # lagged scatterplot
library(spdplyr) # apply dplyr verbs to sp objects
library(rmarkdown)


# ENFORCE EPA FORMAT ON NLA NAMES
toEPA <- function(X1){
  names(X1) = tolower(names(X1))
  names(X1) = gsub(pattern = c("\\(| |#|)|/|-|\\+|:|_"), replacement = ".", x = names(X1))
  X1
}
