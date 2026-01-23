# calculate_EDGE.r

## Source code from rgumbs/EDGE2 repository ------------------------
source("https://raw.githubusercontent.com/rgumbs/EDGE2/main/EDGE.2.calc")

## Create mammal phylogenetic tree ------------------------
matches <- tnrs_match_names(belize_redlist_mammals$species)
ott_ids <- matches$ott_id[!is.na(matches$ott_id) & !matches$approximate_match]
mammal_tree <- tol_induced_subtree(ott_ids = ott_ids)
mammal_tree$tip.label <- mammal_tree$tip.label %>%
    sub("_ott[0-9]+$", "", .) %>%
    gsub("_", " ", .)
synonyms_mammals <- c(
    "Physeter catodon" = "Physeter macrocephalus",
    "Tapirus bairdii" = "Tapirella bairdii",
    "Coendou mexicanus" = "Sphiggurus mexicanus",
    "Odocoileus pandora" = "Mazama pandora",
    "Antrozous dubiaquercus" = "Bauerus dubiaquercus",
    "Micronycteris nicefori" = "Trinycteris nicefori",
    "Micronycteris brachyotis" = "Lampronycteris brachyotis"
)
mammal_tree$tip.label <- ifelse(
    mammal_tree$tip.label %in% names(synonyms_mammals),
    synonyms_mammals[mammal_tree$tip.label],
    mammal_tree$tip.label
)
intersecting_species <- intersect(mammal_tree$tip.label, GEs_mammals$species)
mammal_tree <- keep.tip(mammal_tree, intersecting_species)
GEs_mammals <- GEs_mammals[GEs_mammals$species %in% intersecting_species, ]
GEs_mammals <- GEs_mammals[
    match(mammal_tree$tip.label, GEs_mammals$species),
]
mammal_tree <- mammal_tree %>%
    reorder.phylo(order = "cladewise")
mammal_tree <- compute.brlen(mammal_tree, method = "Grafen")

## Create mammal phylogenetic tree ------------------------
EDGE_mammals_calculated_list <- EDGE2_mod(
    tree = mammal_tree,
    pext = GEs_mammals
)
file.rename("tree.rda", "outputs/tree.rda")
EDGE_mammals_calculated <- EDGE_mammals_calculated_list[[1]]
