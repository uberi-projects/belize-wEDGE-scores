# load_redlist.r

## Source code ------------------------
source("load_packages.r")

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
belize_redlist_taxa <- bind_rows(lapply(files, readRDS))

## Filter to only desired taxa ------------------------
belize_redlist_mammals <- filter(belize_redlist_taxa, class == "Mammalia", !is.na(gbif_id))
belize_redlist_birds <- filter(belize_redlist_taxa, class == "Aves", !is.na(gbif_id))
belize_redlist_reptiles <- filter(belize_redlist_taxa, class == "Squamata", !is.na(gbif_id))
belize_redlist_turtles <- filter(belize_redlist_taxa, class == "Testudines", !is.na(gbif_id))
belize_redlist_amphibians <- filter(belize_redlist_taxa, class == "Amphibia", !is.na(gbif_id))
belize_redlist_fish <- belize_redlist_taxa %>%
    filter(
        !is.na(gbif_id),
        order %in% c(
            "Acipenseriformes", "Albuliformes", "Alepocephaliformes", "Amiiformes", "Anabantiformes",
            "Ateleopodiformes", "Argentiniformes", "Batrachoidiformes", "Beloniformes", "Beryciformes", "Blenniiformes",
            "Caproiformes", "Carangiformes", "Carcharhiniformes", "Centrarchiformes", "Ceratodontiformes", "Chimaeriformes",
            "Chimaeriformes", "Clupeiformes", "Cichliformes", "Coelacanthiformes", "Cypriniformes", "Cyprinodontiformes", "Elopiformes",
            "Esociformes", "Gadiformes", "Galaxiiformes", "Gerreiformes", "Gasterosteiformes", "Gobiesociformes", "Gobiiformes",
            "Gonorynchiformes", "Gymnotiformes", "Heterodontiformes", "Hexanchiformes", "Holocentriformes", "Hiodontiformes",
            "Istiophoriformes", "Kurtiformes", "Lampriformes", "Lamniformes", "Lepisosteiformes", "Lophiiformes", "Mugiliformes",
            "Myliobatiformes", "Myxiniformes", "Notacanthiformes", "Ophidiiformes", "Orectolobiformes", "Osmeriformes", "Percopsiformes",
            "Perciformes", "Pleuronectiformes", "Polypteriformes", "Polymixiiformes", "Pristiophoriformes", "Rajiformes", "Rhinopristiformes",
            "Salmoniformes", "Scorpaeniformes", "Scombriformes", "Siluriformes", "Squatiniformes", "Stomiatiformes", "Stylephoriformes",
            "Syngnathiformes", "Synbranchiformes", "Tetraodontiformes", "Torpediniformes", "Zeiformes"
        )
    )
