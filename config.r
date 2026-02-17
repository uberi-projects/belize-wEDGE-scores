# config.r

# Configuration for phylogenetic tree construction
tree_config <- list(
    birds = list(root_age = 110, constraint = "Struthio camelus"),
    mammals = list(root_age = 160, constraint = "Ornithorhynchus anatinus"),
    amphibians = list(root_age = 350, constraint = "Ambystoma mexicanum"),
    reptiles = list(root_age = 170, constraint = "Sphenodon punctatus"),
    turtles = list(root_age = 220, constraint = "Chelydra serpentina"),
    corals = list(root_age = 450, constraint = "Nematostella vectensis"),
    fish = list(root_age = 420, constraint = "Carcharodon carcharias")
)

# Binomial synonyms for matching species across datasets
synonyms_blank <- c()
synonyms_mammals <- c(
    "Physeter catodon" = "Physeter macrocephalus",
    "Tapirus bairdii" = "Tapirella bairdii",
    "Coendou mexicanus" = "Sphiggurus mexicanus",
    "Odocoileus pandora" = "Mazama pandora",
    "Antrozous dubiaquercus" = "Bauerus dubiaquercus",
    "Micronycteris nicefori" = "Trinycteris nicefori",
    "Micronycteris brachyotis" = "Lampronycteris brachyotis"
)

# Missing coral species to graft onto dated phylogeny (new_species -> congener)
missing_corals <- list(
    list(new = "Acropora cervicornis", congener = "Acropora palmata"),
    list(new = "Agaricia tenuifolia", congener = "Agaricia lamarcki"),
    list(new = "Agaricia agaricites", congener = "Agaricia lamarcki"),
    list(new = "Agaricia fragilis", congener = "Agaricia lamarcki")
)

# wEDGE type labels
type_labels <- list(
    birds = "Bird",
    amphibians = "Amphibian",
    mammals = "Mammal",
    reptiles = "Reptile",
    turtles = "Turtle",
    corals = "Coral",
    fish = "Fish"
)
