# load_redlist.r

## Use IUCN API to get Belize redlist (needs to be set in .Renviron) ------------------------
if (file.exists("outputs/belize_redlist_noDD.rds")) {
    belize_redlist_noDD <- readRDS("outputs/belize_redlist_noDD.rds")
    message("Read existing redlist file (found in outputs)")
} else {
    belize_redlist <- rl_countries("BZ", key = Sys.getenv("IUCN_REDLIST_KEY"), latest = TRUE)
    belize_redlist_noDD <- belize_redlist$assessments %>%
        filter(red_list_category_code != "DD") %>%
        select(taxon_scientific_name, red_list_category_code, year_published)
    saveRDS(belize_redlist_noDD, "outputs/belize_redlist_noDD.rds")
}

## Add in missing clades in batches ------------------------
batch_indices <- split(
    seq_along(belize_redlist_noDD$taxon_scientific_name),
    cut(seq_along(belize_redlist_noDD$taxon_scientific_name), 25, labels = FALSE)
)
for (i in seq_along(batch_indices)) {
    batch_file <- paste0("outputs/taxonomy_batch_", i, ".rds")
    if (file.exists(batch_file)) {
        message("Skipping batch ", i, " (found in outputs)")
    } else {
        taxa <- belize_redlist_noDD$taxon_scientific_name[batch_indices[[i]]]
        taxonomy <- classification(taxa, db = "gbif", ask = FALSE, rank = "species")
        taxonomy <- taxonomy[!sapply(taxonomy, is.logical)]
        tax_df <- imap_dfr(taxonomy, function(.x, .y) {
            species_id <- if ("species" %in% .x$rank) {
                .x$id[.x$rank == "species"]
            } else {
                NA
            }
            .x %>%
                mutate(
                    original_species = .y,
                    gbif_id = species_id
                )
        }) %>%
            select(name, rank, original_species, gbif_id) %>%
            pivot_wider(
                names_from = rank,
                values_from = name,
                values_fn = ~ .x[1]
            )
        saveRDS(tax_df, batch_file)
    }
}

## Combine outputs from batches ------------------------
files <- list.files("outputs", full.names = TRUE)
belize_redlist_taxa <- lapply(files, function(f) {
    readRDS(f) %>% mutate(across(everything(), as.character))
}) %>%
    bind_rows()

## Filter to only desired taxa (except fish) ------------------------
belize_redlist_mammals <- filter(belize_redlist_taxa, class == "Mammalia", !is.na(gbif_id))
belize_redlist_birds <- filter(belize_redlist_taxa, class == "Aves", !is.na(gbif_id))
belize_redlist_reptiles <- filter(belize_redlist_taxa, class == "Squamata", !is.na(gbif_id))
belize_redlist_turtles <- filter(belize_redlist_taxa, class == "Testudines", !is.na(gbif_id))
belize_redlist_amphibians <- filter(belize_redlist_taxa, class == "Amphibia", !is.na(gbif_id))

## Filter to only desired taxa (fish) through FishBase species list for Belize ------------------------
if (file.exists("outputs/fb_countries.rds")) {
    fb_countries <- readRDS("outputs/fb_countries.rds")
    message("Read existing fb countries file (found in outputs)")
} else {
    saveRDS(fb_countries, "outputs/fb_countries.rds")
}
if (file.exists("outputs/fb_species.rds")) {
    fb_species <- readRDS("outputs/fb_species.rds")
    message("Read existing fb species file (found in outputs)")
} else {
    saveRDS(fb_species, "outputs/fb_species.rds")
}
if (file.exists("outputs/fb_families.rds")) {
    fb_families <- readRDS("outputs/fb_families.rds")
    message("Read existing fb families file (found in outputs)")
} else {
    saveRDS(fb_families, "outputs/fb_families.rds")
}
fb_belize_species <- fb_countries %>%
    filter(C_Code == "084") %>%
    left_join(fb_species, by = "SpecCode") %>%
    mutate(taxon_scientific_name = paste(Genus, Species)) %>%
    left_join(fb_families, by = "FamCode") %>%
    select(SpecCode, Freshwater = Freshwater.x, Brackish = Brackish.x, Saltwater = Saltwater.x, Land, taxon_scientific_name, Family) %>%
    mutate(Habitats = Freshwater + Brackish + Saltwater) %>%
    rowwise() %>%
    mutate(
        gbif_match = list(name_backbone(name = taxon_scientific_name))
    ) %>%
    unnest_wider(gbif_match)

belize_redlist_taxa_fishbase <- fb_belize_species %>%
    mutate(gbif_id = speciesKey) %>%
    left_join(belize_redlist_taxa, by = "gbif_id") %>%
    filter(!is.na(Habitats))

belize_redlist_fish_freshwater <- filter(belize_redlist_taxa_fishbase, Freshwater.x == "1" & Habitats == 1, !is.na(gbif_id))
belize_redlist_fish_marine <- filter(belize_redlist_taxa_fishbase, Saltwater.x == "1" & Habitats == 1, !is.na(gbif_id))
belize_redlist_fish_brackish <- filter(belize_redlist_taxa_fishbase, Brackish.x == "1" & Habitats == 1, !is.na(gbif_id))
belize_redlist_fish_mixed <- filter(belize_redlist_taxa_fishbase, Habitats > 1, !is.na(gbif_id))
