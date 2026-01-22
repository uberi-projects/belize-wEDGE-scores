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

## Get GE (pext mean) from IUCN categories ------------------------
GE_CR <- median(pext.samples$pext[pext.samples$RL.cat == "CR"])
GE_EN <- median(pext.samples$pext[pext.samples$RL.cat == "EN"])
GE_VU <- median(pext.samples$pext[pext.samples$RL.cat == "VU"])
GE_NT <- median(pext.samples$pext[pext.samples$RL.cat == "NT"])
GE_LC <- median(pext.samples$pext[pext.samples$RL.cat == "LC"])
