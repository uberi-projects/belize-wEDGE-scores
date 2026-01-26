# calculate_EDGE.r

## Source code from rgumbs/EDGE2 repository ------------------------
source("https://raw.githubusercontent.com/rgumbs/EDGE2/main/EDGE.2.calc")

## Define binomial synonyms ------------------------
synonyms_mammals <- c(
    "Physeter catodon" = "Physeter macrocephalus",
    "Tapirus bairdii" = "Tapirella bairdii",
    "Coendou mexicanus" = "Sphiggurus mexicanus",
    "Odocoileus pandora" = "Mazama pandora",
    "Antrozous dubiaquercus" = "Bauerus dubiaquercus",
    "Micronycteris nicefori" = "Trinycteris nicefori",
    "Micronycteris brachyotis" = "Lampronycteris brachyotis"
)
synonyms_amphibians <- c()

## Define function to create phylogenetic tree for a key taxa ------------------------
create_phylo <- function(redlist_species, ge_data, synonyms) {
    matches <- tnrs_match_names(redlist_species)
    ott_ids <- matches$ott_id[!is.na(matches$ott_id) & !matches$approximate_match]
    created_tree <- tol_induced_subtree(ott_ids = ott_ids)
    created_tree$tip.label <- created_tree$tip.label %>%
        sub("_ott[0-9]+$", "", .) %>%
        gsub("_", " ", .)
    created_tree$tip.label <- ifelse(
        created_tree$tip.label %in% names(synonyms),
        synonyms[created_tree$tip.label],
        created_tree$tip.label
    )
    intersecting_species <- intersect(created_tree$tip.label, ge_data$species)
    created_tree <- keep.tip(created_tree, intersecting_species)
    ge_data <- ge_data[ge_data$species %in% intersecting_species, ]
    ge_data <- ge_data[
        match(created_tree$tip.label, ge_data$species),
    ]
    created_tree <- created_tree %>%
        reorder.phylo(order = "cladewise")
    created_tree <- compute.brlen(created_tree, method = "Grafen")
    list(
        tree = created_tree,
        ge_data = ge_data
    )
}

## Create phylogenetic trees ------------------------
mammals_tree_list <- create_phylo(belize_redlist_mammals$species, GEs_mammals, synonyms_mammals)
amphibians_tree_list <- create_phylo(belize_redlist_amphibians$species, GEs_amphibians, synonyms_amphibians)

## Create directory for output trees ------------------------
directory_trees <- "outputs/trees"
if (!dir.exists(directory_trees)) {
    dir.create(directory_trees, recursive = TRUE)
}

## Calculate mammal EDGE scores ------------------------
EDGE_mammals_calculated_list <- EDGE2_mod(tree = mammals_tree_list$tree, pext = mammals_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_mammals.rda"))
EDGE_mammals_calculated <- EDGE_mammals_calculated_list[[1]]

## Calculate amphibians EDGE scores ------------------------
EDGE_amphibians_calculated_list <- EDGE2_mod(tree = amphibians_tree_list$tree, pext = amphibians_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_amphibians.rda"))
EDGE_amphibians_calculated <- EDGE_amphibians_calculated_list[[1]]
