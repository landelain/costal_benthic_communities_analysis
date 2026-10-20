# Script d'analyse des fichiers CSV
# Author : écrit par moi selon ce qu'on a fait dans les Lab 1 et 2, corrigé par LLM (Lumo AI par défaut, chatGPT si besoin) 
# Application : Benthic communities analysis

library(psych)
library(gplots)

library(readr)
library(dplyr)
library(sf)
library(ggplot2)


### PATH DEFINITION ###
WORKING_DIR <- "~/Documents/EPFL/2026_MA_3/MULTI_R/costal_benthic_communities_analysis/Cyrielle"
DATA_DIR <- "~/Documents/EPFL/2026_MA_3/MULTI_R/costal_benthic_communities_analysis/data_clean"
BEN_DATA_FILE <- "benthic_clean.csv"
ENV_DATA_FILE <- "env_clean.csv"
JOINED_BEN_ENV_DATA_FILE <- "benthic_env_clean.csv"

setwd(DATA_DIR)



### FONCTIONS ###

open_and_analyze_data <- function(csv_path) {
  
  # Import des données
  data <- read.csv(csv_path, header = TRUE, sep = ",")
  
  # Filtrer uniquement les colonnes numériques
  num_cols <- sapply(data, is.numeric) 
  data_num <- data[, num_cols, drop = FALSE]
  
  # Informations de base sur le dataframe

  str(data)
  cat("\nNom des colonnes:\n")
  print(colnames(data))
  cat("\nNom des lignes:\n")
  print(rownames(data))
  cat("\nDimensions (lignes, colonnes):\n")
  print(dim(data))
  
  
  # Verification des NA/NaN
  na_count <- sum(is.na(data_num))
  
  cat("\nValeurs NA:", na_count, "\n")

  
  ### Stats de base
  
  cat("\nSummary stats:\n")
  cat("\nboxplots\n")
  print(summary(data))
  cat("\ndistribution\n")
  print(describeBy(data[, num_cols]))
  cat("\nMIN-MAX\n")
  print(range(data_num, na.rm = TRUE))
  
  # Somme des lignes (toutes les espèces pour un site)
  line_sums <- rowSums(data[, num_cols], na.rm = TRUE)
  cat("\nSommes par ligne (toutes les espèces/parametres pour un site):\n")
  print(line_sums)
  cat("\nRange des sommes de lignes (min, max):\n")
  print(range(line_sums, na.rm = TRUE))
  
  # Somme des colonnes (une espèce sur tous les sites)
  col_sums <- colSums(data_num, na.rm = TRUE)
  cat("\nSommes par colonne (une espèce/parametre sur tous les sites):\n")
  print(col_sums)
  cat("\nRange des sommes de colonnes (min, max):\n")
  print(range(col_sums, na.rm = TRUE))
  
  
  # fin
  return(list(
    data = data,
    data_num = data_num,
    num_cols = num_cols,
    line_sums = line_sums,
    col_sums = col_sums
  ))
}


### MAIN ###
ben_data <- open_and_analyze_data(BEN_DATA_FILE)
#ben_data$data   # Le dataframe original
#ben_data_num$data_num   # Le dataframe ne contenant que les colonnes numeric
#ben_data$line_sums   # Sommes par ligne
#ben_data$col_sums    # Sommes par colonne




##############
### PART 0 ###
##############


# importe les données ;
# csélectionne les 27 catégories benthiques ;
# agrège les abondances par site en calculant leur moyenne ;
# conserve les coordonnées moyennes de chaque site ;
# calcule la couverture totale à partir des 27 catégories ;
# prépare une couche géographique pour la cartographie.

setwd(WORKING_DIR)
# 1. Importer les données
df <- read.csv(JOINED_BEN_ENV_DATA_FILE, header = TRUE, sep = ",")

# 2. Identifier les 27 catégories benthiques
cover_cols <- c(
  "ac_erect_single",
  "ag_articulated_calcareous",
  "ag_corticated_foliose",
  "ag_corticated_macrophyte",
  "cca_crustose",
  "cy_filamentous",
  "hc_arborescent",
  "hc_bushy",
  "hc_column",
  "hc_encrusting",
  "hc_foliose",
  "hc_massive",
  "hc_table",
  "oc_arborescent",
  "oc_bushy",
  "oc_cluster",
  "oc_digitate",
  "oc_encrusting",
  "oc_lobate",
  "oc_massive",
  "sg",
  "sp_encrusting",
  "sp_massive",
  "sp_repent",
  "turf_filamentous",
  "zo_encrusting",
  "zo_massive"
)

# Vérifier les colonnes
stopifnot(all(cover_cols %in% names(df)))

# 3. Créer un identifiant unique de site
df <- df %>%
  mutate(
    Site_ID = paste(Region, Site, sep = "_")
  )

# 4. Agréger les transects par site (somme par site, de toutes les moyennes de couverture par espèces entre les transects)
site_data <- df %>%
  group_by(Site_ID, Region, Site) %>%
  summarise(
    n_transects = n_distinct(Transect_2),
    Longitude = mean(Longitude, na.rm = TRUE),
    Latitude = mean(Latitude, na.rm = TRUE),
    
    across(
      all_of(cover_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    Total_cover = rowSums(
      across(all_of(cover_cols)),
      na.rm = FALSE
    )
  )

# Vérifier le nombre de sites
nrow(site_data)
length(unique(site_data$Site_ID))

# Résumé des abondances
summary(site_data$Total_cover)

# 5. Créer une couche géographique en WGS 84
sites_sf <- st_as_sf(
  site_data,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)

# Projection métrique pour Taïwan
sites_projected <- st_transform(sites_sf, 3826)

# 6. Cartographier la couverture totale
ggplot(sites_sf) +
  geom_sf(aes(color = Total_cover), size = 3) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(
    title = "Benthic cover across Taiwan",
    color = "Mean cover"
  )


st_write(
  sites_sf,
  "cover_per_sites.gpkg",
  layer = "total_cover",      # nom de la couche dans le gpkg
  append = FALSE,             # TRUE pour ajouter à un gpkg existant
  delete_dsn = FALSE,          # écrase le fichier s'il existe déjà
  delete_layer = TRUE         # écrase la couche si elle existe déjà
)
# Pour relire plus tard : site_sf <- st_read("cover_per_sites.gpkg")






##############
### PART 1 ###
##############


### Moran's I global ###






