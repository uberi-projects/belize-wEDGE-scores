# load_EDGE.r

## Load data ------------------------
if (file.exists("data_deposit/EDGE_all.xlsx")) {
    df_EDGE_amphibians <- read_excel("data_deposit/EDGE_all.xlsx", sheet = 2)
    df_EDGE_mammals <- read_excel("data_deposit/EDGE_all.xlsx", sheet = 5)
    df_EDGE_fish1 <- read_excel("data_deposit/EDGE_all.xlsx", sheet = 6)
    df_EDGE_lepidosaurs <- read_excel("data_deposit/EDGE_all.xlsx", sheet = 7)
    df_EDGE_crocodylians <- read_excel("data_deposit/EDGE_all.xlsx", sheet = 8)
    df_EDGE_turtles <- read_excel("data_deposit/EDGE_all.xlsx", sheet = 9)
    df_EDGE_fish2 <- read_excel("data_deposit/EDGE_all.xlsx", sheet = 10)
} else {
    message("No data found for general EDGE scores; check README.md for instructions on including the data")
}
if (file.exists("data_deposit/EDGE_birds.xlsx")) {
    df_EDGE_birds <- read_excel("data_deposit/EDGE_birds.xlsx") %>%
        mutate(Common.names = NA) %>%
        select(Species, Common.names, EDGE)
} else {
    message("No data found for bird EDGE scores; check README.md for instructions on including the data")
}

## Generate list of EDGE scores by species ------------------------
df_EDGE_all_exceptbirds <- bind_rows(
    df_EDGE_amphibians, df_EDGE_mammals, df_EDGE_fish1, df_EDGE_lepidosaurs,
    df_EDGE_crocodylians, df_EDGE_turtles, df_EDGE_fish2
) %>%
    rename(EDGE = EDGE.median) %>%
    select(Species, Common.names, EDGE)
df_EDGE_all <- bind_rows(df_EDGE_all_exceptbirds, df_EDGE_birds) %>%
    rename(species = Species)
