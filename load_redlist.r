# load_redlist.r

## Source code ------------------------
source("load_packages.r")

## Use IUCN API to get Belize redlist (needs to be set in .Renviron) ------------------------
belize_redlist <- rl_countries("BZ", key = Sys.getenv("IUCN_REDLIST_KEY"), latest = TRUE)
belize_redlist_noDD <- belize_redlist$assessments %>%
    filter(red_list_category_code != "DD") %>%
    select(taxon_scientific_name, red_list_category_code, year_published)

## Add in missing clades in batches ------------------------
batch_indices <- split(
    seq_along(belize_redlist_noDD$taxon_scientific_name),
    cut(seq_along(belize_redlist_noDD$taxon_scientific_name), 25, labels = FALSE)
)
taxonomy_batches <- vector("list", length(batch_indices))
for (i in seq_along(batch_indices)) {
    taxa <- belize_redlist_noDD$taxon_scientific_name[batch_indices[[i]]]
    taxonomy <- classification(
        taxa,
        db = "gbif",
        ask = FALSE,
        rank = "species"
    )
    taxonomy <- taxonomy[!sapply(taxonomy, is.logical)]
    tax_df <- imap_dfr(
        taxonomy,
        ~ mutate(.x, original_species = .y)
    ) %>%
        select(name, rank, original_species) %>%
        pivot_wider(
            names_from = rank,
            values_from = name,
            values_fn = ~ .x[1]
        )
    saveRDS(
        tax_df,
        file = paste0("outputs/taxonomy_batch_", i, ".rds")
    )
    taxonomy_batches[[i]] <- tax_df
}

## Combine outputs from batches ------------------------
files <- list.files("outputs", full.names = TRUE)
belize_redlist_taxa <- bind_rows(lapply(files, readRDS))
