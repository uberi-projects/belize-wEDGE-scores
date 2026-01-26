# calculate_wEDGE.r

## Finalize EDGE scores to be used ------------------------
df_EDGE_all <- mutate(df_EDGE_all, EDGE_Source = "Published")
df_EDGE_all_calculated <- mutate(df_EDGE_all_calculated, EDGE_Source = "Calculated")
df_EDGE_all_combined <- df_EDGE_all %>%
    full_join(df_EDGE_all_calculated, by = "species", suffix = c("_all", "_calc")) %>%
    mutate(
        species,
        EDGE_Source = coalesce(EDGE_Source_all, EDGE_Source_calc),
        Common.names = coalesce(Common.names_all, Common.names_calc),
        EDGE = coalesce(EDGE_all, EDGE_calc),
        .keep = "none"
    )

## Assemble wEDGE for birds ------------------------
wEDGE_birds <- weights_belize_birds %>%
    left_join(df_EDGE_all_combined, by = "species") %>%
    filter(!is.na(EDGE)) %>%
    mutate(wEDGE = EDGE * weight, Type = "Bird", gbif_id = as.numeric(gbif_id)) %>%
    left_join(select(belize_redlist_birds, red_list_category_code, original_species), by = join_by(species == original_species))

## Assemble wEDGE for amphibians ------------------------
wEDGE_amphibians <- weights_belize_amphibians %>%
    left_join(df_EDGE_all_combined, by = "species") %>%
    filter(!is.na(EDGE)) %>%
    mutate(wEDGE = EDGE * weight, Type = "Amphibian", gbif_id = as.numeric(gbif_id)) %>%
    left_join(select(belize_redlist_amphibians, red_list_category_code, original_species), by = join_by(species == original_species))

## Assemble wEDGE for mammals ------------------------
wEDGE_mammals <- weights_belize_mammals %>%
    left_join(df_EDGE_all_combined, by = "species") %>%
    filter(!is.na(EDGE)) %>%
    mutate(wEDGE = EDGE * weight, Type = "Mammal", gbif_id = as.numeric(gbif_id)) %>%
    left_join(select(belize_redlist_mammals, red_list_category_code, original_species), by = join_by(species == original_species))

## Assemble wEDGE for reptiles ------------------------
wEDGE_reptiles <- weights_belize_reptiles %>%
    left_join(df_EDGE_all_combined, by = "species") %>%
    filter(!is.na(EDGE)) %>%
    mutate(wEDGE = EDGE * weight, Type = "Reptile", gbif_id = as.numeric(gbif_id)) %>%
    left_join(select(belize_redlist_reptiles, red_list_category_code, original_species), by = join_by(species == original_species))

## Assemble wEDGE for turtles ------------------------
wEDGE_turtles <- weights_belize_turtles %>%
    left_join(df_EDGE_all_combined, by = "species") %>%
    filter(!is.na(EDGE)) %>%
    mutate(wEDGE = EDGE * weight, Type = "Turtle", gbif_id = as.numeric(gbif_id)) %>%
    left_join(select(belize_redlist_turtles, red_list_category_code, original_species), by = join_by(species == original_species))

## Assemble wEDGE for corals ------------------------
wEDGE_corals <- weights_belize_corals %>%
    left_join(df_EDGE_all_combined, by = "species") %>%
    filter(!is.na(EDGE)) %>%
    mutate(wEDGE = EDGE * weight, Type = "Coral", gbif_id = as.numeric(gbif_id)) %>%
    left_join(select(belize_redlist_corals, red_list_category_code, original_species), by = join_by(species == original_species))

## Assemble wEDGE for freshwater fish ------------------------
wEDGE_fish_freshwater <- weights_belize_fish_freshwater %>%
    left_join(df_EDGE_all_combined, by = "species") %>%
    filter(!is.na(EDGE)) %>%
    mutate(wEDGE = EDGE * weight, Type = "Freshwater Fish", gbif_id = as.numeric(gbif_id)) %>%
    left_join(select(belize_redlist_fish_freshwater, red_list_category_code, taxon_scientific_name), by = join_by(species == taxon_scientific_name))

## Assemble wEDGE for marine fish ------------------------
wEDGE_fish_marine <- distinct(weights_belize_fish_marine) %>%
    left_join(df_EDGE_all_combined, by = "species") %>%
    filter(!is.na(EDGE)) %>%
    mutate(wEDGE = EDGE * weight, Type = "Marine Fish", gbif_id = as.numeric(gbif_id)) %>%
    left_join(select(belize_redlist_fish_marine, red_list_category_code, taxon_scientific_name), by = join_by(species == taxon_scientific_name))

## Assemble wEDGE for mixed habitat fish ------------------------
wEDGE_fish_mixed <- weights_belize_fish_mixed %>%
    left_join(df_EDGE_all_combined, by = "species") %>%
    filter(!is.na(EDGE)) %>%
    mutate(wEDGE = EDGE * weight, Type = "Mixed Fish", gbif_id = as.numeric(gbif_id)) %>%
    left_join(select(belize_redlist_fish_mixed, red_list_category_code, taxon_scientific_name), by = join_by(species == taxon_scientific_name))

## Create directory for output wEDGE scores ------------------------
directory_wEDGE <- "outputs/wEDGE"
if (!dir.exists(directory_wEDGE)) {
    dir.create(directory_wEDGE, recursive = TRUE)
}

## Compile and export results ------------------------
wEDGE_all <- bind_rows(
    wEDGE_birds, wEDGE_amphibians, wEDGE_mammals, wEDGE_reptiles, wEDGE_turtles,
    wEDGE_corals, wEDGE_fish_freshwater, wEDGE_fish_marine, wEDGE_fish_mixed
) %>%
    select(Type, EDGE_Source, Species = species, Common_Names = Common.names, w = weight, EDGE, wEDGE, Redlist = red_list_category_code) %>%
    arrange(Type, -wEDGE)
write.csv(wEDGE_all, file.path(directory_wEDGE, "wEDGE_all.csv"))
