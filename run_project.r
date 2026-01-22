# run_project.r

# This script should be used to run the project in the proper sequence to calculate wEDGE scores for Belize

## Source code ------------------------
source("load_packages.r") # attach required packages
source("load_redlist.r") # fetch redlist data from IUCN and create list of Belize species to analyze
source("calculate_w.r") # calculate weight for species (proportion of ocurrences or range in Belize)
source("calculate_GE.r") # define function to calculate GE for species that have no EDGE score
source("calculate_EDGE.r") # define function to calculate EDGE for species that have no EDGE score
source("load_EDGE.r") # load in existing EDGE scores for Belize species
source("calculate_wEDGE.r") # assemble final wEDGE scores
