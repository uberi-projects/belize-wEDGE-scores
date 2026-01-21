# calculate_w.r

## Source code ------------------------
source("load_packages.r")

## Read belize basemap ------------------------
belize_map <- st_read("shapefiles/Belize_Basemap.shp") %>%
    st_transform(4326) %>%
    st_make_valid()

## Define function to fetch occurrence data for target species from GBIF------------------------
create_species_sf <- function(id_vector) {
    species_sf <- list()
    for (id in id_vector) {
        species_data <- occ_search(
            taxonKey = id,
            hasCoordinate = TRUE
        )
        sf_name <- paste0("species_sf_", id)
        species_sf[[sf_name]] <- st_as_sf(species_data$data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
    }
    return(species_sf)
}

## Define function to calculate weights for Belize from occurrence data ------------------------
calculate_weights_belize <- function(name_vector, species_sf_list) {
    weights <- c()
    names <- c()
    for (sf_name in names(species_sf_list)) {
        species_sf <- species_sf_list[[sf_name]]
        num_inside_belize <- st_within(species_sf, belize_map, sparse = FALSE)
        weight_belize <- sum(num_inside_belize) / nrow(species_sf)
        weights[sf_name] <- weight_belize
    }
    weights <- data.frame(Source = names(weights), Name = name_vector, Weight = unname(weights))
    weights <- weights %>% filter(Weight > 0.0001)
    return(weights)
}

## Calculate proportion of occurrences in Belize for mammals using GBIF occurence data ------------------------
species_sf <- create_species_sf(belize_redlist_fish$gbif_id)
weights_belize <- calculate_weights_belize(belize_redlist_fish$species, species_sf)
weights_belize
