# calculate_w.r

## Source code ------------------------
source("load_packages.r")

## Read belize basemap ------------------------
belize_map <- st_read("shapefiles/Belize_Basemap.shp") %>%
    st_transform(4326) %>%
    st_make_valid()

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

## Calculate proportion of occurrences in Belize for amphibians using GBIF occurrence data ------------------------
if (file.exists("outputs/weights_belize_amphibians.rds")) {
    weights_belize_amphibians <- readRDS("outputs/weights_belize_amphibians.rds")
    message("Read existing amphibians weights file (found in outputs)")
} else {
    weights_belize_amphibians <- calculate_weights_belize(belize_redlist_amphibians$gbif_id, belize_redlist_amphibians$species)
    saveRDS(weights_belize_amphibians, "outputs/weights_belize_amphibians.rds")
}
weights_belize_amphibians

## Calculate proportion of occurrences in Belize for mammals using GBIF occurrence data ------------------------
if (file.exists("outputs/weights_belize_mammals.rds")) {
    weights_belize_mammals <- readRDS("outputs/weights_belize_mammals.rds")
    message("Read existing mammals weights file (found in outputs)")
} else {
    weights_belize_mammals <- calculate_weights_belize(belize_redlist_mammals$gbif_id, belize_redlist_mammals$species)
    saveRDS(weights_belize_mammals, "outputs/weights_belize_mammals.rds")
}
weights_belize_mammals

## Calculate proportion of occurrences in Belize for reptiles using GBIF occurrence data ------------------------
if (file.exists("outputs/weights_belize_reptiles.rds")) {
    weights_belize_reptiles <- readRDS("outputs/weights_belize_reptiles.rds")
    message("Read existing reptiles weights file (found in outputs)")
} else {
    weights_belize_reptiles <- calculate_weights_belize(belize_redlist_reptiles$gbif_id, belize_redlist_reptiles$species)
    saveRDS(weights_belize_reptiles, "outputs/weights_belize_reptiles.rds")
}
weights_belize_reptiles

## Calculate proportion of occurrences in Belize for turtles using GBIF occurrence data ------------------------
if (file.exists("outputs/weights_belize_turtles.rds")) {
    weights_belize_turtles <- readRDS("outputs/weights_belize_turtles.rds")
    message("Read existing turtles weights file (found in outputs)")
} else {
    weights_belize_turtles <- calculate_weights_belize(belize_redlist_turtles$gbif_id, belize_redlist_turtles$species)
    saveRDS(weights_belize_turtles, "outputs/weights_belize_turtles.rds")
}
weights_belize_turtles
