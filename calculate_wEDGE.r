# calculate_wEDGE.r

## Assemble wEDGE for amphibians ------------------------
wEDGE_amphibians <- weights_belize_amphibians %>%
    left_join(df_EDGE_all, by = "species") %>%
    filter(!is.na(EDGE)) %>%
    mutate(wEDGE = EDGE * weight, Type = "Amphibian", gbif_id = as.numeric(gbif_id))

## Assemble wEDGE for mammals ------------------------
wEDGE_mammals <- weights_belize_mammals %>%
    left_join(df_EDGE_all, by = "species") %>%
    filter(!is.na(EDGE)) %>%
    mutate(wEDGE = EDGE * weight, Type = "Mammal", gbif_id = as.numeric(gbif_id))

## Assemble wEDGE for reptiles ------------------------
wEDGE_reptiles <- weights_belize_reptiles %>%
    left_join(df_EDGE_all, by = "species") %>%
    filter(!is.na(EDGE)) %>%
    mutate(wEDGE = EDGE * weight, Type = "Reptile", gbif_id = as.numeric(gbif_id))

## Assemble wEDGE for turtles ------------------------
wEDGE_turtles <- weights_belize_turtles %>%
    left_join(df_EDGE_all, by = "species") %>%
    filter(!is.na(EDGE)) %>%
    mutate(wEDGE = EDGE * weight, Type = "Turtle", gbif_id = as.numeric(gbif_id))

## Assemble wEDGE for freshwater fish ------------------------
wEDGE_fish_freshwater <- weights_belize_fish_freshwater %>%
    left_join(df_EDGE_all, by = "species") %>%
    filter(!is.na(EDGE)) %>%
    mutate(wEDGE = EDGE * weight, Type = "Freshwater Fish", gbif_id = as.numeric(gbif_id))

## Assemble wEDGE for marine fish ------------------------
wEDGE_fish_marine <- weights_belize_fish_marine %>%
    left_join(df_EDGE_all, by = "species") %>%
    filter(!is.na(EDGE)) %>%
    mutate(wEDGE = EDGE * weight, Type = "Marine Fish", gbif_id = as.numeric(gbif_id))

## Assemble wEDGE for brackish fish ------------------------
wEDGE_fish_brackish <- weights_belize_fish_brackish %>%
    left_join(df_EDGE_all, by = "species") %>%
    filter(!is.na(EDGE)) %>%
    mutate(wEDGE = EDGE * weight, Type = "Brackish Fish", gbif_id = as.numeric(gbif_id))

## Assemble wEDGE for mixed habitat fish ------------------------
wEDGE_fish_mixed <- weights_belize_fish_mixed %>%
    left_join(df_EDGE_all, by = "species") %>%
    filter(!is.na(EDGE)) %>%
    mutate(wEDGE = EDGE * weight, Type = "Mixed Fish", gbif_id = as.numeric(gbif_id))

## Compile and export results ------------------------
wEDGE_all <- bind_rows(
    wEDGE_amphibians, wEDGE_mammals, wEDGE_reptiles, wEDGE_turtles,
    wEDGE_fish_freshwater, wEDGE_fish_marine, wEDGE_fish_brackish, wEDGE_fish_mixed
) %>%
    arrange(Type, -wEDGE)
write.csv(wEDGE_all, "outputs/wEDGE_all.csv")
