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
create_phylo <- function(group, redlist_species, ge_data, root_age_my, constraint_taxon, synonyms = synonyms_blank) {
    file_name <- paste0("outputs/trees/tree_list_", group, ".rds")
    if (file.exists(file_name)) {
        created_tree_list <- readRDS(file_name)
        message(paste("Read existing", group, "tree list file (found in outputs)"))
    } else {
        redlist_species <- redlist_species[!is.na(redlist_species)]
        species_for_tree <- unique(c(redlist_species, constraint_taxon))
        matches <- tnrs_match_names(species_for_tree)
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
        tips_to_keep <- unique(c(intersecting_species, constraint_taxon))
        created_tree <- keep.tip(created_tree, tips_to_keep)
        created_tree <- root(
            created_tree,
            outgroup = constraint_taxon,
            resolve.root = TRUE
        )
        ge_data <- ge_data[ge_data$species %in% intersecting_species, ]
        ge_data <- ge_data[
            match(created_tree$tip.label, ge_data$species),
        ]
        created_tree <- created_tree %>%
            reorder.phylo(order = "cladewise")
        created_tree <- compute.brlen(created_tree, method = "Grafen")
        current_height <- max(node.depth.edgelength(created_tree))
        scale_factor <- root_age_my / current_height
        created_tree$edge.length <- created_tree$edge.length * scale_factor
        created_tree_list <- list(tree = created_tree, ge_data = ge_data, root_age_my = root_age_my, constraint = constraint_taxon)
        saveRDS(created_tree_list, file_name)
    }
    created_tree_list
}

## Create phylogenetic trees ------------------------
birds_tree_list <- create_phylo("birds", belize_redlist_birds$species, GEs_birds, 110, "Struthio camelus")
mammals_tree_list <- create_phylo("mammals", belize_redlist_mammals$species, GEs_mammals, 160, "Ornithorhynchus anatinus", synonyms_mammals)
amphibians_tree_list <- create_phylo("amphibians", belize_redlist_amphibians$species, GEs_amphibians, 350, "Ambystoma mexicanum")
reptiles_tree_list <- create_phylo("reptiles", belize_redlist_reptiles$species, GEs_reptiles, 170, "Sphenodon punctatus")
turtles_tree_list <- create_phylo("turtles", belize_redlist_turtles$species, GEs_turtles, 220, "Chelydra serpentina")
corals_tree_list <- create_phylo("corals", belize_redlist_corals$species, GEs_corals, 450, "Nematostella vectensis")
fish_tree_list <- create_phylo(
    "fish",
    c(belize_redlist_fish_freshwater$species, belize_redlist_fish_marine$species, belize_redlist_fish_mixed$species),
    bind_rows(GEs_fish_freshwater, GEs_fish_marine, GEs_fish_mixed),
    420, "Carcharodon carcharias"
)

## Calculate EDGE scores ------------------------
EDGE_birds_calculated_list <- EDGE2_mod(tree = drop.tip(birds_tree_list$tree, birds_tree_list$constraint), pext = birds_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_birds.rda"))
EDGE_birds_calculated <- EDGE_birds_calculated_list[[1]]
EDGE_mammals_calculated_list <- EDGE2_mod(tree = drop.tip(mammals_tree_list$tree, mammals_tree_list$constraint), pext = mammals_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_mammals.rda"))
EDGE_mammals_calculated <- EDGE_mammals_calculated_list[[1]]
EDGE_amphibians_calculated_list <- EDGE2_mod(tree = drop.tip(amphibians_tree_list$tree, amphibians_tree_list$constraint), pext = amphibians_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_amphibians.rda"))
EDGE_amphibians_calculated <- EDGE_amphibians_calculated_list[[1]]
EDGE_reptiles_calculated_list <- EDGE2_mod(tree = drop.tip(reptiles_tree_list$tree, reptiles_tree_list$constraint), pext = reptiles_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_reptiles.rda"))
EDGE_reptiles_calculated <- EDGE_reptiles_calculated_list[[1]]
EDGE_turtles_calculated_list <- EDGE2_mod(tree = drop.tip(turtles_tree_list$tree, turtles_tree_list$constraint), pext = turtles_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_turtles.rda"))
EDGE_turtles_calculated <- EDGE_turtles_calculated_list[[1]]
EDGE_corals_calculated_list <- EDGE2_mod(tree = drop.tip(corals_tree_list$tree, corals_tree_list$constraint), pext = corals_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_corals.rda"))
EDGE_corals_calculated <- EDGE_corals_calculated_list[[1]]
EDGE_fish_calculated_list <- EDGE2_mod(tree = drop.tip(fish_tree_list$tree, fish_tree_list$constraint), pext = fish_tree_list$ge_data)
file.rename("tree.rda", file.path(directory_trees, "tree_fish.rda"))
EDGE_fish_calculated <- EDGE_fish_calculated_list[[1]]

## Compile calculated EDGE scores ------------------------
df_EDGE_all_calculated <- EDGE_birds_calculated %>%
    bind_rows(EDGE_mammals_calculated) %>%
    bind_rows(EDGE_amphibians_calculated) %>%
    bind_rows(EDGE_reptiles_calculated) %>%
    bind_rows(EDGE_turtles_calculated) %>%
    bind_rows(EDGE_corals_calculated) %>%
    bind_rows(EDGE_fish_calculated) %>%
    mutate(species = Species, Common.names = NA) %>%
    select(species, ED_Calc = ED, EDGE_Calc = EDGE)

## Load dated coral tree ------------------------
# https://doi.org/10.1038/s41586-025-09615-6
coral_tree_text <- readLines(
    "data_deposit/tree_mol_clock_penalized_likelihood_method_R.newick.rtf",
    warn = FALSE
)
coral_tree_text_newick <- coral_tree_text[grepl("\\(", coral_tree_text) & grepl(";", coral_tree_text)]
coral_tree_newick <- paste(coral_tree_text_newick, collapse = "")
writeLines(coral_tree_newick, "data_deposit/coral_tree.newick")
coral_tree_dated <- read.tree("data_deposit/coral_tree.newick")
coral_tree_dated <- coral_tree_dated
coral_tree_dated$tip.label <- str_replace(coral_tree_dated$tip.label, "^([A-Za-z]+_[a-z]+).*", "\\1")
coral_tree_dated$tip.label <- gsub("_", " ", coral_tree_dated$tip.label)
coral_tree_dated$tip.label[which(coral_tree_dated$tip.label == "Porites cf")[2]] <- "Porites cf 2"

## Add important missing coral species ------------------------
add_missing_tree_species <- function(tree, new_species, congener) {
    tip_idx <- which(tree$tip.label == congener)
    edge_row <- which(tree$edge[, 2] == tip_idx)
    orig_bl <- tree$edge.length[edge_row]
    tree <- bind.tip(tree, tip.label = new_species, where = tip_idx, position = orig_bl / 2)
    tree
}
coral_tree_dated <- add_missing_tree_species(coral_tree_dated, "Acropora cervicornis", "Acropora palmata")
coral_tree_dated <- add_missing_tree_species(coral_tree_dated, "Agaricia tenuifolia", "Agaricia lamarcki")
coral_tree_dated <- add_missing_tree_species(coral_tree_dated, "Agaricia agaricites", "Agaricia lamarcki")
coral_tree_dated <- add_missing_tree_species(coral_tree_dated, "Agaricia fragilis", "Agaricia lamarcki")

## Calculate derived coral EDGE scores ------------------------
all_coral_tree_species <- coral_tree_dated$tip.label
GEs_corals_full <- data.frame(
    species = all_coral_tree_species,
    GE = ifelse(all_coral_tree_species %in% belize_redlist_corals$species, GEs_corals$GE[match(all_coral_tree_species, GEs_corals$species)], 0)
)
GEs_corals_full <- GEs_corals_full[match(coral_tree_dated$tip.label, GEs_corals_full$species), ]
EDGE_corals_derived_list <- EDGE2_mod(tree = coral_tree_dated, pext = GEs_corals_full)
EDGE_corals_derived <- EDGE_corals_derived_list[[1]]
EDGE_corals_derived <- filter(EDGE_corals_derived, EDGE > 0)

## Compile derived EDGE scores ------------------------
df_EDGE_all_derived <- EDGE_corals_derived %>%
    mutate(species = Species) %>%
    select(species, ED_Deriv = ED, EDGE_Deriv = EDGE)
