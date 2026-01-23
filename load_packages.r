# load_packages.r

## Check required packages ------------------------
options(repos = c(CRAN = "https://cran.rstudio.com/"))
required_packages <- c("tidyverse", "rgbif", "maps", "rredlist", "taxize", "rfishbase", "readxl", "rotl", "phylobase", "data.table", "ape")
install_if_missing <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
        install.packages(package)
    }
}
invisible(lapply(required_packages, install_if_missing))

## Attach packages ------------------------
library(tidyverse)
library(rgbif)
library(maps)
library(rredlist)
library(taxize)
library(rfishbase)
library(readxl)
library(rotl)
library(phylobase)
library(data.table)
library(ape)
