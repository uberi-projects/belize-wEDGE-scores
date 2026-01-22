# calculate_GE.r

## Source code from rgumbs/EDGE2 repository ------------------------
source("https://raw.githubusercontent.com/rgumbs/EDGE2/main/GE.2.calc")
source("https://raw.githubusercontent.com/rgumbs/EDGE2/main/EDGE.2.calc")

## Generate pext samples ------------------------
if (file.exists("outputs/pext_samples_EDGE2.rds")) {
    pext.samples <- readRDS("outputs/pext_samples_EDGE2.rds")
} else {
    pext.samples <- GE.2.calc(pext.vals)
    saveRDS(pext.samples, "outputs/pext_samples_EDGE2.rds")
}

## Get GE sample from IUCN categories ------------------------
sample_GE <- function(RL_cat, pext_samples = pext.samples) {
    cat_pext <- pext.samples$pext[pext.samples$RL.cat == RL_cat]
    sample(cat_pext, 1)
}
