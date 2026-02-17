# calculate_w.r
# Spatial weight calculation: proportion of species range area inside Belize
# Uses clustered convex hulls from GBIF occurrence points:
#   1. Download GBIF occurrences for each species
#   2. Cluster points geographically (to handle disjunct ranges)
#   3. Build a convex hull per cluster, buffer by 25km
#   4. Union all cluster hulls = estimated global range
#   5. weight = area of range in Belize / total range area
# For the legacy approach (proportion of GBIF occurrence counts), see legacy_calculate_w.r

## Read belize basemap (needed for spatial intersection below) ------------------------
belize_map <- st_read("basemap/Belize_Basemap.shp") %>%
    st_transform(4326) %>%
    st_make_valid()

## Merge Belize polygons into a single boundary for intersection ------------------------
belize_boundary <- st_union(belize_map)

## Equal-area projection for accurate area calculations ------------------------
mollweide_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

## Define function to build clustered convex hulls from occurrence points -----------------
build_clustered_hull <- function(occ_sf, cluster_distance_deg = 5, buffer_km = 25) {
    coords_matrix <- st_coordinates(occ_sf)

    if (nrow(coords_matrix) < 3) {
        return(NULL)
    }

    # Cluster points: groups within cluster_distance_deg of each other (~550km at equator)
    dist_matrix <- dist(coords_matrix)
    hc <- hclust(dist_matrix, method = "complete")
    clusters <- cutree(hc, h = cluster_distance_deg)

    # Build a convex hull per cluster, skip clusters with < 3 points
    hulls <- list()
    for (cl in unique(clusters)) {
        cluster_pts <- occ_sf[clusters == cl, ]
        if (nrow(cluster_pts) < 3) {
            # For tiny clusters (1-2 points), buffer the points instead
            buffered <- cluster_pts %>%
                st_transform(mollweide_crs) %>%
                st_buffer(buffer_km * 1000) %>%
                st_union() %>%
                st_transform(4326)
            hulls[[length(hulls) + 1]] <- buffered
        } else {
            hull <- st_convex_hull(st_combine(cluster_pts))
            # Buffer hull in equal-area projection
            hull_buffered <- hull %>%
                st_transform(mollweide_crs) %>%
                st_buffer(buffer_km * 1000) %>%
                st_transform(4326)
            hulls[[length(hulls) + 1]] <- hull_buffered
        }
    }

    # Union all cluster hulls into single range polygon
    range_polygon <- do.call(c, hulls) %>%
        st_union() %>%
        st_make_valid()

    return(range_polygon)
}

## Define function to calculate spatial weights for Belize ------------------------
calculate_spatial_weights_belize <- function(id_vector, name_vector) {
    out <- list()

    for (i in seq_along(id_vector)) {
        id <- id_vector[i]
        name <- name_vector[i]
        message(paste0("[", i, "/", length(id_vector), "] Processing ", name, "..."))

        # Download GBIF occurrences using taxon key
        global_occs <- try(
            occ_data(
                taxonKey = id,
                hasCoordinate = TRUE,
                hasGeospatialIssue = FALSE,
                limit = 300
            ),
            silent = TRUE
        )

        if (inherits(global_occs, "try-error") || is.null(global_occs$data)) {
            message(paste("  Skipping", name, "- GBIF query failed or returned no data"))
            next
        }

        coords <- global_occs$data %>%
            select(decimalLongitude, decimalLatitude) %>%
            filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
            distinct()

        if (nrow(coords) < 3) {
            message(paste("  Skipping", name, "- fewer than 3 unique coordinate points"))
            next
        }

        message(paste("  Retrieved", nrow(coords), "unique occurrence points"))

        # Convert to sf points
        occ_sf <- st_as_sf(coords, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

        # Build clustered convex hull range
        species_range <- try(build_clustered_hull(occ_sf), silent = TRUE)

        if (inherits(species_range, "try-error") || is.null(species_range)) {
            message(paste("  Skipping", name, "- hull construction failed"))
            next
        }

        # Calculate global range area (Mollweide equal-area projection)
        global_area_m2 <- species_range %>%
            st_transform(mollweide_crs) %>%
            st_area() %>%
            sum() %>%
            as.numeric()

        if (global_area_m2 == 0) {
            message(paste("  Skipping", name, "- global range area is zero"))
            next
        }

        # Intersect species range with Belize boundary
        belize_range <- try(
            st_intersection(species_range, belize_boundary),
            silent = TRUE
        )

        if (inherits(belize_range, "try-error") || length(belize_range) == 0) {
            belize_area_m2 <- 0
        } else {
            belize_area_m2 <- belize_range %>%
                st_transform(mollweide_crs) %>%
                st_area() %>%
                sum() %>%
                as.numeric()
        }

        weight <- belize_area_m2 / global_area_m2
        global_area_km2 <- global_area_m2 / 1e6
        belize_area_km2 <- belize_area_m2 / 1e6

        n_clusters <- length(unique(cutree(
            hclust(dist(st_coordinates(occ_sf)), method = "complete"),
            h = 5
        )))

        out[[i]] <- data.frame(
            gbif_id = id,
            species = name,
            weight = weight,
            n_points = nrow(coords),
            n_clusters = n_clusters,
            global_area_km2 = round(global_area_km2, 1),
            belize_area_km2 = round(belize_area_km2, 1)
        )

        message(paste0(
            "  Weight: ", round(weight, 6),
            " (", round(belize_area_km2, 1), " / ",
            round(global_area_km2, 1), " km2)",
            " [", n_clusters, " cluster(s)]"
        ))

        Sys.sleep(0.1)
    }

    result <- bind_rows(out)
    if (nrow(result) == 0) {
        warning("No species produced valid spatial weights in this group")
        return(data.frame(
            gbif_id = numeric(), species = character(), weight = numeric(),
            n_points = integer(), n_clusters = integer(),
            global_area_km2 = numeric(), belize_area_km2 = numeric()
        ))
    }
    result %>%
        filter(weight > 0.0001) %>%
        arrange(-weight)
}

## Helper to load cached spatial weights or calculate them ------------------------
load_or_calculate_spatial_weights <- function(group_name, redlist_data, species_col = "species") {
    file_path <- paste0("outputs/weights_spatial/weights_belize_", group_name, ".rds")
    if (file.exists(file_path)) {
        message(paste("Read existing", group_name, "spatial weights file (found in outputs)"))
        readRDS(file_path) %>% select(gbif_id, species, weight)
    } else {
        result <- calculate_spatial_weights_belize(
            id_vector = redlist_data$gbif_id,
            name_vector = redlist_data[[species_col]]
        )
        saveRDS(result, file_path)
        result %>% select(gbif_id, species, weight)
    }
}

## Create directory for output spatial weights ------------------------
directory_weights_spatial <- "outputs/weights_spatial"
if (!dir.exists(directory_weights_spatial)) {
    dir.create(directory_weights_spatial, recursive = TRUE)
}

## Calculate spatial weights for non-fish taxa ------------------------
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
        load_or_calculate_spatial_weights(group_name, non_fish_groups[[group_name]]),
        envir = .GlobalEnv
    )
}

## Calculate spatial weights for fish taxa ------------------------
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
        load_or_calculate_spatial_weights(group_name, fish_groups[[group_name]], species_col = "species.x"),
        envir = .GlobalEnv
    )
}
