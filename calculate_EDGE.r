# calculate_EDGE.r

## Perform Monte Carlo EDGE2 calculation for missing species ------------------------
test_data <- data.frame(
    species = c("Jaguar", "Hicatee", "Whatever else"),
    RL_cat  = c("NT", "CR", "CR"),
    ED      = c(50.2, 10.5, 5.3)
)
EDGE_matrix <- matrix(NA, nrow = nrow(test_data), ncol = 1000)
for (i in 1:1000) {
    GE_draw <- sapply(test_data$RL_cat, sample_GE)
    EDGE_matrix[, i] <- test_data$ED * GE_draw
}
species_EDGE <- data.frame(
    species = test_data$species,
    EDGE_median = apply(EDGE_matrix, 1, median),
    EDGE_CI_lower = apply(EDGE_matrix, 1, function(x) quantile(x, 0.025)),
    EDGE_CI_upper = apply(EDGE_matrix, 1, function(x) quantile(x, 0.975))
)
print(species_EDGE)

## Create mammal phylogenetic tree ------------------------
matches <- tnrs_match_names(belize_redlist_mammals$species)
ott_ids <- matches$ott_id[!is.na(matches$ott_id) & !matches$approximate_match]
mammal_tree <- tol_induced_subtree(ott_ids = ott_ids)
