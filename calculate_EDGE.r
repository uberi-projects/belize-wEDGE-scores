# calculate_EDGE.r

## Source code from rgumbs/EDGE2 repository ------------------------
source("https://raw.githubusercontent.com/rgumbs/EDGE2/main/EDGE.2.calc")

## Define binomial synonyms ------------------------
synonyms_blank <- c()
synonyms_mammals <- c(
    "Physeter catodon" = "Physeter macrocephalus",
    "Tapirus bairdii" = "Tapirella bairdii",
    "Coendou mexicanus" = "Sphiggurus mexicanus",
    "Odocoileus pandora" = "Mazama pandora",
    "Antrozous dubiaquercus" = "Bauerus dubiaquercus",
    "Micronycteris nicefori" = "Trinycteris nicefori",
    "Micronycteris brachyotis" = "Lampronycteris brachyotis"
)

## Create directory for output trees ------------------------
directory_trees <- "outputs/trees"
if (!dir.exists(directory_trees)) {
    dir.create(directory_trees, recursive = TRUE)
}

## Define function to create phylogenetic tree for a key taxa ------------------------
create_phylo <- function(group, redlist_species, ge_data, synonyms = synonyms_blank) {
    file_name <- paste0("outputs/trees/tree_list_", group, ".rds")
    if (file.exists(file_name)) {
        created_tree_list <- readRDS(file_name)
        message(paste("Read existing", group, "tree list file (found in outputs)"))
    } else {
        redlist_species <- redlist_species[!is.na(redlist_species)]
        matches <- tnrs_match_names(redlist_species)
        ott_ids <- matches$ott_id[!is.na(matches$ott_id) & !matches$approximate_match]
        ott_ids <- unique(ott_ids)
        ott_ids <- ott_ids[is_in_tree(ott_ids)]
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
        created_tree_list <- list(
            tree = created_tree,
            ge_data = ge_data
        )
        saveRDS(created_tree_list, file_name)
    }
    created_tree_list
}

## Create phylogenetic trees ------------------------
birds_tree_list <- create_phylo("birds", belize_redlist_birds$species, GEs_birds)
mammals_tree_list <- create_phylo("mammals", belize_redlist_mammals$species, GEs_mammals, synonyms_mammals)
amphibians_tree_list <- create_phylo("amphibians", belize_redlist_amphibians$species, GEs_amphibians)
reptiles_tree_list <- create_phylo("reptiles", belize_redlist_reptiles$species, GEs_reptiles)
turtles_tree_list <- create_phylo("turtles", belize_redlist_turtles$species, GEs_turtles)
corals_tree_list <- create_phylo("corals", belize_redlist_corals$species, GEs_corals)
fish_freshwater_tree_list <- create_phylo("fish_freshwater", belize_redlist_fish_freshwater$species, GEs_fish_freshwater)
fish_marine_tree_list <- create_phylo("fish_marine", belize_redlist_fish_marine$species, GEs_fish_marine)
fish_mixed_tree_list <- create_phylo("fish_mixed", belize_redlist_fish_mixed$species, GEs_fish_mixed)

## Calculate birds EDGE scores ------------------------
EDGE_birds_calculated_list <- EDGE2_mod(tree = birds_tree_list$tree, pext = birds_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_birds.rda"))
EDGE_birds_calculated <- EDGE_birds_calculated_list[[1]]

## Calculate mammals EDGE scores ------------------------
EDGE_mammals_calculated_list <- EDGE2_mod(tree = mammals_tree_list$tree, pext = mammals_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_mammals.rda"))
EDGE_mammals_calculated <- EDGE_mammals_calculated_list[[1]]

## Calculate amphibians EDGE scores ------------------------
EDGE_amphibians_calculated_list <- EDGE2_mod(tree = amphibians_tree_list$tree, pext = amphibians_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_amphibians.rda"))
EDGE_amphibians_calculated <- EDGE_amphibians_calculated_list[[1]]

## Calculate reptiles EDGE scores ------------------------
EDGE_reptiles_calculated_list <- EDGE2_mod(tree = reptiles_tree_list$tree, pext = reptiles_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_reptiles.rda"))
EDGE_reptiles_calculated <- EDGE_reptiles_calculated_list[[1]]

## Calculate turtles EDGE scores ------------------------
EDGE_turtles_calculated_list <- EDGE2_mod(tree = turtles_tree_list$tree, pext = turtles_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_turtles.rda"))
EDGE_turtles_calculated <- EDGE_turtles_calculated_list[[1]]

## Calculate corals EDGE scores ------------------------
EDGE_corals_calculated_list <- EDGE2_mod(tree = corals_tree_list$tree, pext = corals_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_corals.rda"))
EDGE_corals_calculated <- EDGE_corals_calculated_list[[1]]

## Calculate fish_freshwater EDGE scores ------------------------
EDGE_fish_freshwater_calculated_list <- EDGE2_mod(tree = fish_freshwater_tree_list$tree, pext = fish_freshwater_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_fish_freshwater.rda"))
EDGE_fish_freshwater_calculated <- EDGE_fish_freshwater_calculated_list[[1]]

## Calculate fish_marine EDGE scores ------------------------
EDGE_fish_marine_calculated_list <- EDGE2_mod(tree = fish_marine_tree_list$tree, pext = fish_marine_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_fish_marine.rda"))
EDGE_fish_marine_calculated <- EDGE_fish_marine_calculated_list[[1]]

## Calculate fish_mixed EDGE scores ------------------------
EDGE_fish_mixed_calculated_list <- EDGE2_mod(tree = fish_mixed_tree_list$tree, pext = fish_mixed_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_fish_mixed.rda"))
EDGE_fish_mixed_calculated <- EDGE_fish_mixed_calculated_list[[1]]

## Compile calculated EDGE scores ------------------------
df_EDGE_all_calculated <- EDGE_birds_calculated %>%
    bind_rows(EDGE_mammals_calculated) %>%
    bind_rows(EDGE_amphibians_calculated) %>%
    bind_rows(EDGE_reptiles_calculated) %>%
    bind_rows(EDGE_turtles_calculated) %>%
    bind_rows(EDGE_corals_calculated) %>%
    bind_rows(EDGE_fish_freshwater_calculated) %>%
    bind_rows(EDGE_fish_marine_calculated) %>%
    bind_rows(EDGE_fish_mixed_calculated) %>%
    mutate(species = Species, Common.names = NA) %>%
    select(species, Common.names, EDGE)
