# calculate_wEDGE.r

## Finalize EDGE scores to be used ------------------------
df_EDGE_all_combined <- df_EDGE_all %>%
    full_join(df_EDGE_all_calculated, by = "species") %>%
    full_join(df_EDGE_all_derived, by = "species") %>%
    select(-Common.names)

## Helper to assemble wEDGE for a standard taxa group ------------------------
assemble_wEDGE <- function(weights_data, redlist_data, type_label, join_col = "original_species") {
    result <- weights_data %>%
        left_join(df_EDGE_all_combined, by = "species") %>%
        mutate(wEDGE_Calc = EDGE_Calc * weight, Type = type_label, gbif_id = as.numeric(gbif_id)) %>%
        mutate(wEDGE_Pub = EDGE_Pub * weight)
    if ("EDGE_Deriv" %in% names(result) && any(!is.na(result$EDGE_Deriv))) {
        result <- mutate(result, wEDGE_Deriv = EDGE_Deriv * weight)
    }
    result %>%
        left_join(
            select(redlist_data, red_list_category_code, !!sym(join_col)),
            by = setNames(join_col, "species")
        )
}

## Assemble wEDGE for each standard taxa group ------------------------
wEDGE_groups <- list(
    birds      = list(weights = weights_belize_birds,      redlist = belize_redlist_birds),
    amphibians = list(weights = weights_belize_amphibians,  redlist = belize_redlist_amphibians),
    mammals    = list(weights = weights_belize_mammals,     redlist = belize_redlist_mammals),
    reptiles   = list(weights = weights_belize_reptiles,    redlist = belize_redlist_reptiles),
    turtles    = list(weights = weights_belize_turtles,     redlist = belize_redlist_turtles),
    corals     = list(weights = weights_belize_corals,      redlist = belize_redlist_corals)
)
wEDGE_results <- list()
for (group_name in names(wEDGE_groups)) {
    grp <- wEDGE_groups[[group_name]]
    wEDGE_results[[group_name]] <- assemble_wEDGE(
        grp$weights, grp$redlist, type_labels[[group_name]]
    )
}

## Assemble wEDGE for fish (combines sub-groups) ------------------------
belize_redlist_fish <- bind_rows(belize_redlist_fish_freshwater, belize_redlist_fish_marine, belize_redlist_fish_mixed) %>% distinct()
wEDGE_results[["fish"]] <- bind_rows(weights_belize_fish_freshwater, weights_belize_fish_marine, weights_belize_fish_mixed) %>%
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
wEDGE_all <- bind_rows(wEDGE_results) %>%
    select(Type, Species = species, w = weight, ED_Calc, ED_Deriv, EDGE_Calc, EDGE_Pub, EDGE_Deriv, wEDGE_Calc, wEDGE_Pub, wEDGE_Deriv, Redlist = red_list_category_code) %>%
    arrange(Type, -wEDGE_Calc, -wEDGE_Pub, -wEDGE_Deriv) %>%
    filter(!is.na(wEDGE_Calc) | !is.na(wEDGE_Pub) | !is.na(wEDGE_Deriv))
write.csv(wEDGE_all, file.path(directory_wEDGE, "wEDGE_all.csv"))
