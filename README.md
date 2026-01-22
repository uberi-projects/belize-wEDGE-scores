# Belize wEDGE Scores

This repository is intended to calculate wEDGE scores for target Belize species, to support development of a national priority species list.

To run this repository's project, follow these steps:
1. Ensure you have **access to IUCN data**. Either put an API key in your .Renviron folder (IUCN_REDLIST_KEY='yourkeyhere'), or drop existing taxonomy_batch.rds files into outputs to skip using the API.
2. Simply run the R script "run_project.r", which will result in an output file, "wEDGE.xlsx," for your use.

As you run the script, other outputs will be created, which are stored R objects (.rds files). These are temporary checkpoints throughout the script, so that if you must rerun later you don't need to wait for the IUCN, GBIF, or FishBase APIs again. If you have existing versions of these output files from an external source already, such as a UB-ERI collaborator, feel free to add those directly into the outputs folder yourself to greatly reduce time spent in the R scripts.


# Helper Scripts

"run_project.r" sources several helper scripts to calculate wEDGE scores for Belize:

1. "load_packages.r" attaches required packages
2. "load_redlist.r" fetches redlist data from IUCN and creates list of Belize species to analyze (this script requires IUCN API access or existing outputs to run)
3. "calculate_w.r" calculates weight for species (proportion of ocurrences or range in Belize)
4. "load_EDGE.r" loads in existing EDGE scores for Belize species
5. "calculate_wEDGE.r" assembles final wEDGE scores


## wEDGE Resources

[Birds EDGE List](https://conbio.onlinelibrary.wiley.com/doi/full/10.1111/cobi.14141)

[EDGE List](https://www.edgeofexistence.org/download-edge-lists/) on Edge of Existence

[EDGE2 Protocol](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001991)

[Raptor wEDGE Scores](https://www.researchgate.net/profile/Christopher-Mcclure-3/publication/387629803_Determining_global_and_regional_raptor_conservation_priorities/links/6820e3f960241d5140291f1d/Determining-global-and-regional-raptor-conservation-priorities.pdf)
