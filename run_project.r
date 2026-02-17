# run_project.r

# This script should be used to run the project in the proper sequence to calculate wEDGE scores for Belize

## Source code ------------------------
source("load_packages.r") # attach required packages
source("config.r") # load configuration (root ages, synonyms, type labels, etc.)
source("load_redlist.r") # fetch redlist data from IUCN and create list of Belize species to analyze along with their redlist status
source("calculate_w.r") # calculate weight for species (proportion of occurrences or range in Belize)
source("calculate_GE.r") # define a function to calculate GE for species that have no EDGE score and run it
source("calculate_EDGE.r") # calculate EDGE values for species based on either Tree of Life phylogenies or published molecular dated phylogenies
source("load_EDGE.r") # load in published  EDGE scores for Belize species
source("calculate_wEDGE.r") # assemble final wEDGE scores, including those derived from published and calculated EDGE values
