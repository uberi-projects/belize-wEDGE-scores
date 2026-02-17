# calculate_w.r

## Define function to calculate weights for Belize from GBIF occurrence data ------------------------
calculate_weights_belize <- function(id_vector, name_vector) {
    out <- list()
    for (i in seq_along(id_vector)) {
        id <- id_vector[i]
        name <- name_vector[i]
        global_n <- occ_search(taxonKey = id, hasCoordinate = TRUE, limit = 0)$meta$count
        if (global_n == 0) next
        belize_n <- occ_search(taxonKey = id, hasCoordinate = TRUE, country = "BZ", limit = 0)$meta$count
        weight <- belize_n / global_n
        out[[i]] <- data.frame(gbif_id = id, species = name, weight = weight)
        Sys.sleep(0.1)
    }
    bind_rows(out) %>%
        filter(weight > 0.0001) %>%
        arrange(-weight)
}

## Helper to load cached weights or calculate them ------------------------
load_or_calculate_weights <- function(group_name, redlist_data, species_col = "species") {
    file_path <- paste0("outputs/weights/weights_belize_", group_name, ".rds")
    if (file.exists(file_path)) {
        message(paste("Read existing", group_name, "weights file (found in outputs)"))
        readRDS(file_path)
    } else {
        result <- calculate_weights_belize(redlist_data$gbif_id, redlist_data[[species_col]])
        saveRDS(result, file_path)
        result
    }
}

## Create directory for output weights ------------------------
directory_weights <- "outputs/weights"
if (!dir.exists(directory_weights)) {
    dir.create(directory_weights, recursive = TRUE)
}

## Calculate weights for non-fish taxa using GBIF occurrence data ------------------------
non_fish_groups <- list(
    amphibians = belize_redlist_amphibians,
    mammals = belize_redlist_mammals,
    reptiles = belize_redlist_reptiles,
    turtles = belize_redlist_turtles,
    corals = belize_redlist_corals
)
for (group_name in names(non_fish_groups)) {
    assign(
        paste0("weights_belize_", group_name),
        load_or_calculate_weights(group_name, non_fish_groups[[group_name]]),
        envir = .GlobalEnv
    )
}

## Calculate weights for fish taxa using GBIF occurrence data ------------------------
belize_redlist_fish_freshwater <- distinct(belize_redlist_fish_freshwater)
belize_redlist_fish_marine <- distinct(belize_redlist_fish_marine)
belize_redlist_fish_mixed <- distinct(belize_redlist_fish_mixed)

fish_groups <- list(
    fish_freshwater = belize_redlist_fish_freshwater,
    fish_marine = belize_redlist_fish_marine,
    fish_mixed = belize_redlist_fish_mixed
)
for (group_name in names(fish_groups)) {
    assign(
        paste0("weights_belize_", group_name),
        load_or_calculate_weights(group_name, fish_groups[[group_name]], species_col = "species.x"),
        envir = .GlobalEnv
    )
}

## Load BirdLife International range maps ------------------------
files_dir <- "data_deposit/species"
full_path <- file.path(files_dir, "BOTW_2025.gpkg")

if (file.exists(full_path)) {
    global_birds_range <- st_read(dsn = full_path, layer = "all_species")
    message("Geopackage found. Read all_species layer (found in data_deposit)")
} else {
    warning(
        "Geopackage not found"
    )
}

## Read belize basemap ------------------------
belize_map <- st_read("basemap/Belize_Basemap.shp") %>%
    st_transform(4326) %>%
    st_make_valid()

## Calculate proportion of occurrences in Belize for birds using BirdLife International range maps ------------------------
