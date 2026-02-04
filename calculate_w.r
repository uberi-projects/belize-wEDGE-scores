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

## Create directory for output weights ------------------------
directory_weights <- "outputs/weights"
if (!dir.exists(directory_weights)) {
    dir.create(directory_weights, recursive = TRUE)
}

## Calculate proportion of occurrences in Belize for amphibians using GBIF occurrence data ------------------------
if (file.exists("outputs/weights/weights_belize_amphibians.rds")) {
    weights_belize_amphibians <- readRDS("outputs/weights/weights_belize_amphibians.rds")
    message("Read existing amphibians weights file (found in outputs)")
} else {
    weights_belize_amphibians <- calculate_weights_belize(belize_redlist_amphibians$gbif_id, belize_redlist_amphibians$species)
    saveRDS(weights_belize_amphibians, "outputs/weights/weights_belize_amphibians.rds")
}

## Calculate proportion of occurrences in Belize for mammals using GBIF occurrence data ------------------------
if (file.exists("outputs/weights/weights_belize_mammals.rds")) {
    weights_belize_mammals <- readRDS("outputs/weights/weights_belize_mammals.rds")
    message("Read existing mammals weights file (found in outputs)")
} else {
    weights_belize_mammals <- calculate_weights_belize(belize_redlist_mammals$gbif_id, belize_redlist_mammals$species)
    saveRDS(weights_belize_mammals, "outputs/weights/weights_belize_mammals.rds")
}

## Calculate proportion of occurrences in Belize for reptiles using GBIF occurrence data ------------------------
if (file.exists("outputs/weights/weights_belize_reptiles.rds")) {
    weights_belize_reptiles <- readRDS("outputs/weights/weights_belize_reptiles.rds")
    message("Read existing reptiles weights file (found in outputs)")
} else {
    weights_belize_reptiles <- calculate_weights_belize(belize_redlist_reptiles$gbif_id, belize_redlist_reptiles$species)
    saveRDS(weights_belize_reptiles, "outputs/weights/weights_belize_reptiles.rds")
}

## Calculate proportion of occurrences in Belize for turtles using GBIF occurrence data ------------------------
if (file.exists("outputs/weights/weights_belize_turtles.rds")) {
    weights_belize_turtles <- readRDS("outputs/weights/weights_belize_turtles.rds")
    message("Read existing turtles weights file (found in outputs)")
} else {
    weights_belize_turtles <- calculate_weights_belize(belize_redlist_turtles$gbif_id, belize_redlist_turtles$species)
    saveRDS(weights_belize_turtles, "outputs/weights/weights_belize_turtles.rds")
}

## Calculate proportion of occurrences in Belize for corals using GBIF occurrence data ------------------------
if (file.exists("outputs/weights/weights_belize_corals.rds")) {
    weights_belize_corals <- readRDS("outputs/weights/weights_belize_corals.rds")
    message("Read existing corals weights file (found in outputs)")
} else {
    weights_belize_corals <- calculate_weights_belize(belize_redlist_corals$gbif_id, belize_redlist_corals$species)
    saveRDS(weights_belize_corals, "outputs/weights/weights_belize_corals.rds")
}

## Calculate proportion of occurrences in Belize for freshwater fish using GBIF occurrence data ------------------------
belize_redlist_fish_freshwater <- distinct(belize_redlist_fish_freshwater)
if (file.exists("outputs/weights/weights_belize_fish_freshwater.rds")) {
    weights_belize_fish_freshwater <- readRDS("outputs/weights/weights_belize_fish_freshwater.rds")
    message("Read existing fish_freshwater weights file (found in outputs)")
} else {
    weights_belize_fish_freshwater <- calculate_weights_belize(belize_redlist_fish_freshwater$gbif_id, belize_redlist_fish_freshwater$species.x)
    saveRDS(weights_belize_fish_freshwater, "outputs/weights/weights_belize_fish_freshwater.rds")
}

## Calculate proportion of occurrences in Belize for marine fish using GBIF occurrence data ------------------------
belize_redlist_fish_marine <- distinct(belize_redlist_fish_marine)
if (file.exists("outputs/weights/weights_belize_fish_marine.rds")) {
    weights_belize_fish_marine <- readRDS("outputs/weights/weights_belize_fish_marine.rds")
    message("Read existing fish_marine weights file (found in outputs)")
} else {
    weights_belize_fish_marine <- calculate_weights_belize(belize_redlist_fish_marine$gbif_id, belize_redlist_fish_marine$species.x)
    saveRDS(weights_belize_fish_marine, "outputs/weights/weights_belize_fish_marine.rds")
}

## Calculate proportion of occurrences in Belize for brackish fish using GBIF occurrence data ------------------------
belize_redlist_fish_brackish <- distinct(belize_redlist_fish_brackish)
if (file.exists("outputs/weights/weights_belize_fish_brackish.rds")) {
    weights_belize_fish_brackish <- readRDS("outputs/weights/weights_belize_fish_brackish.rds")
    message("Read existing fish_brackish weights file (found in outputs)")
} else {
    weights_belize_fish_brackish <- calculate_weights_belize(belize_redlist_fish_brackish$gbif_id, belize_redlist_fish_brackish$species.x)
    saveRDS(weights_belize_fish_brackish, "outputs/weights/weights_belize_fish_brackish.rds")
}

## Calculate proportion of occurrences in Belize for mixed habitat fish using GBIF occurrence data ------------------------
belize_redlist_fish_mixed <- distinct(belize_redlist_fish_mixed)
if (file.exists("outputs/weights/weights_belize_fish_mixed.rds")) {
    weights_belize_fish_mixed <- readRDS("outputs/weights/weights_belize_fish_mixed.rds")
    message("Read existing fish_mixed weights file (found in outputs)")
} else {
    weights_belize_fish_mixed <- calculate_weights_belize(belize_redlist_fish_mixed$gbif_id, belize_redlist_fish_mixed$species.x)
    saveRDS(weights_belize_fish_mixed, "outputs/weights/weights_belize_fish_mixed.rds")
}

## Load BirdLife International range maps ------------------------
all_files <- c(
  "global_birds.cpg",
  "global_birds.dbf",
  "global_birds.prj", 
  "global_birds.sbn", 
  "global_birds.sbx", 
  "global_birds.shx", 
  "global_birds.shp"
)
files_dir <- "data_deposit/global_bird_range_maps"
full_path <- file.path (files_dir, all_files)
if (all(file.exists(full_path))) {
  global_birds_shp <- st_read(file.path(files_dir, "global_birds.shp"))
  message("All shapefile components found. Read global_birds shapefile (found in data_deposit)")
} else {
  missing_files <- all_files[!file.exists(full_path)]
  warning(
    "Missing shapefile components: ",
    paste(missing_files, collapse = ",")
    )
}

global_birds_shp <- st_make_valid(global_birds_shp)

#Read belize basemap ------------------------
belize_map <- st_read("basemap/Belize_Basemap.shp") %>%
  st_transform(4326) %>%
  st_make_valid()

## Calculate proportion of occurrences in Belize for birds using BirdLife International range maps ------------------------

