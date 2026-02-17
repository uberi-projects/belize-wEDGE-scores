# calculate_GE.r

## Source code from rgumbs/EDGE2 repository ------------------------
source("https://raw.githubusercontent.com/rgumbs/EDGE2/main/GE.2.calc")

## Generate pext samples ------------------------
if (file.exists("outputs/pext_samples_EDGE2.rds")) {
    pext.samples <- readRDS("outputs/pext_samples_EDGE2.rds")
} else {
    pext.samples <- GE.2.calc(pext.vals)
    saveRDS(pext.samples, "outputs/pext_samples_EDGE2.rds")
}

## Get GE (pext mean) from IUCN categories ------------------------
GE_CR <- median(pext.samples$pext[pext.samples$RL.cat == "CR"])
GE_EN <- median(pext.samples$pext[pext.samples$RL.cat == "EN"])
GE_VU <- median(pext.samples$pext[pext.samples$RL.cat == "VU"])
GE_NT <- median(pext.samples$pext[pext.samples$RL.cat == "NT"])
GE_LC <- median(pext.samples$pext[pext.samples$RL.cat == "LC"])
iucn_categories <- data.frame(red_list_category_code = c("CR", "EN", "VU", "NT", "LC"), GE = c(GE_CR, GE_EN, GE_VU, GE_NT, GE_LC))

## Assign GEs for each Belize species ------------------------
assign_GE <- function(df) {
    df %>%
        select(species, red_list_category_code) %>%
        left_join(iucn_categories, by = "red_list_category_code") %>%
        na.omit() %>%
        select(-red_list_category_code)
}
GEs_birds <- assign_GE(belize_redlist_birds)
GEs_amphibians <- assign_GE(belize_redlist_amphibians)
GEs_mammals <- assign_GE(belize_redlist_mammals)
GEs_reptiles <- assign_GE(belize_redlist_reptiles)
GEs_turtles <- assign_GE(belize_redlist_turtles)
GEs_corals <- assign_GE(belize_redlist_corals)
GEs_fish_freshwater <- assign_GE(belize_redlist_fish_freshwater)
GEs_fish_marine <- assign_GE(belize_redlist_fish_marine)
GEs_fish_mixed <- assign_GE(belize_redlist_fish_mixed)
