# calculate_wEDGE.r

## Finalize EDGE scores to be used ------------------------
df_EDGE_all_combined <- df_EDGE_all %>%
    full_join(df_EDGE_all_calculated, by = "species") %>%
    select(-Common.names)

## Assemble wEDGE for birds ------------------------
wEDGE_birds <- weights_belize_birds %>%
    left_join(df_EDGE_all_combined, by = "species") %>%
    mutate(wEDGE_Calc = EDGE_Calc * weight, Type = "Bird", gbif_id = as.numeric(gbif_id)) %>%
    mutate(wEDGE_Pub = EDGE_Pub * weight, Type = "Bird") %>%
    left_join(select(belize_redlist_birds, red_list_category_code, original_species), by = join_by(species == original_species))

## Assemble wEDGE for amphibians ------------------------
wEDGE_amphibians <- weights_belize_amphibians %>%
    left_join(df_EDGE_all_combined, by = "species") %>%
    mutate(wEDGE_Calc = EDGE_Calc * weight, Type = "Amphibian", gbif_id = as.numeric(gbif_id)) %>%
    mutate(wEDGE_Pub = EDGE_Pub * weight, Type = "Amphibian") %>%
    left_join(select(belize_redlist_amphibians, red_list_category_code, original_species), by = join_by(species == original_species))

## Assemble wEDGE for mammals ------------------------
wEDGE_mammals <- weights_belize_mammals %>%
    left_join(df_EDGE_all_combined, by = "species") %>%
    mutate(wEDGE_Calc = EDGE_Calc * weight, Type = "Mammal", gbif_id = as.numeric(gbif_id)) %>%
    mutate(wEDGE_Pub = EDGE_Pub * weight, Type = "Mammal") %>%
    left_join(select(belize_redlist_mammals, red_list_category_code, original_species), by = join_by(species == original_species))

## Assemble wEDGE for reptiles ------------------------
wEDGE_reptiles <- weights_belize_reptiles %>%
    left_join(df_EDGE_all_combined, by = "species") %>%
    mutate(wEDGE_Calc = EDGE_Calc * weight, Type = "Reptile", gbif_id = as.numeric(gbif_id)) %>%
    mutate(wEDGE_Pub = EDGE_Pub * weight, Type = "Reptile") %>%
    left_join(select(belize_redlist_reptiles, red_list_category_code, original_species), by = join_by(species == original_species))

## Assemble wEDGE for turtles ------------------------
wEDGE_turtles <- weights_belize_turtles %>%
    left_join(df_EDGE_all_combined, by = "species") %>%
    mutate(wEDGE_Calc = EDGE_Calc * weight, Type = "Turtle", gbif_id = as.numeric(gbif_id)) %>%
    mutate(wEDGE_Pub = EDGE_Pub * weight, Type = "Turtle") %>%
    left_join(select(belize_redlist_turtles, red_list_category_code, original_species), by = join_by(species == original_species))

## Assemble wEDGE for corals ------------------------
wEDGE_corals <- weights_belize_corals %>%
    left_join(df_EDGE_all_combined, by = "species") %>%
    mutate(wEDGE_Calc = EDGE_Calc * weight, Type = "Coral", gbif_id = as.numeric(gbif_id)) %>%
    mutate(wEDGE_Pub = EDGE_Pub * weight, Type = "Coral") %>%
    left_join(select(belize_redlist_corals, red_list_category_code, original_species), by = join_by(species == original_species))

## Assemble wEDGE for fish ------------------------
belize_redlist_fish <- bind_rows(belize_redlist_fish_freshwater, belize_redlist_fish_marine, belize_redlist_fish_mixed) %>% distinct()
wEDGE_fish <- bind_rows(weights_belize_fish_freshwater, weights_belize_fish_marine, weights_belize_fish_mixed) %>%
    distinct() %>%
    left_join(df_EDGE_all_combined, by = "species") %>%
    mutate(wEDGE_Calc = EDGE_Calc * weight, Type = "Fish", gbif_id = as.numeric(gbif_id)) %>%
    mutate(wEDGE_Pub = EDGE_Pub * weight, Type = "Fish") %>%
    left_join(select(belize_redlist_fish, red_list_category_code, taxon_scientific_name), by = join_by(species == taxon_scientific_name))

## Create directory for output wEDGE scores ------------------------
directory_wEDGE <- "outputs/wEDGE"
if (!dir.exists(directory_wEDGE)) {
    dir.create(directory_wEDGE, recursive = TRUE)
}

## Compile and export results ------------------------
wEDGE_all <- bind_rows(wEDGE_birds, wEDGE_amphibians, wEDGE_mammals, wEDGE_reptiles, wEDGE_turtles, wEDGE_corals, wEDGE_fish) %>%
    select(Type, Species = species, w = weight, ED, EDGE_Calc, EDGE_Pub, wEDGE_Calc, wEDGE_Pub, Redlist = red_list_category_code) %>%
    arrange(Type, -wEDGE_Calc, -wEDGE_Pub) %>%
    filter(!is.na(wEDGE_Calc) | !is.na(wEDGE_Pub))
write.csv(wEDGE_all, file.path(directory_wEDGE, "wEDGE_all.csv"))
