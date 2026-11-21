# Script d'analyse des fichiers CSV
# Author : écrit par moi selon ce qu'on a fait dans les Lab 1 et 2, corrigé par LLM (Lumo AI par défaut, chatGPT si besoin) 
# Application : Benthic communities analysis

library(psych)
library(gplots)


### PATH DEFINITION ###
WORKING_DIR <- "~/Documents/EPFL/2026_MA_3/MULTI_R/costal_benthic_communities_analysis/Cyrielle"
DATA_DIR <- "~/Documents/EPFL/2026_MA_3/MULTI_R/costal_benthic_communities_analysis/data_clean"
BEN_DATA_FILE <- "benthic_clean.csv"
ENV_DATA_FILE <- "env_clean.csv"
JOINED_BEN_ENV_DATA_FILE <- "benthic_env_clean.csv"
JOINED_BEN_ENV_DATA_FILE_W_KMEANS_CLUSTERS <- "clusters_communities_kmeans.csv"

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


library(readr)
library(dplyr)
library(sf)
library(ggplot2)


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

# 5. Créer une couche géographique en WGS 84, coords : attention, c'est (X, Y) = (Longitude, Latitude), dans cet ordre.
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


### Moran's I global on total coverage per site ###


library(spdep)

# 1. Extraire les coordonnées projetées
coords <- st_coordinates(sites_sf)

# 2. Construire les voisins KNN
knn4 <- knearneigh(coords, k = 4)
nb4 <- knn2nb(knn4)

# 3. Construire la matrice de poids spatiaux
lw4 <- nb2listw(
  nb4,
  style = "W",
  zero.policy = TRUE
)

# 4. Moran's I global sur la couverture totale
x <- sites_sf$Total_cover

# Vérifier les valeurs manquantes
sum(is.na(x))

# Retirer les sites incomplets de manière cohérente
keep <- is.finite(x)

coords_ok <- coords[keep, , drop = FALSE]

nb_ok <- knn2nb(
  knearneigh(coords_ok, k = 4)
)

lw_ok <- nb2listw(
  nb_ok,
  style = "W",
  zero.policy = TRUE
)

x_ok <- x[keep]

# Moran's I avec approximation analytique
moran.test(
  x_ok,
  lw_ok,
  zero.policy = TRUE
)

# Test par permutations
set.seed(123)

moran.mc(
  x_ok,
  lw_ok,
  nsim = 9999,
  zero.policy = TRUE
)




##############
### PART 2 ###
##############

# L'objectif est de déterminer si les transects appartenant à une même communauté ont tendance à être géographiquement proches.

### La profondeur peut fortement influencer la composition benthique. De plus, certains sites possèdent plusieurs transects aux mêmes coordonnées géographiques, mais à différentes profondeurs.
#Je recommande donc de :
#  - regrouper les transects par localisation et profondeur ;
#- calculer les proportions des cinq communautés pour chaque localisation-profondeur ;
#- réaliser les analyses spatiales séparément pour chaque profondeur.
#Cela évite notamment de confondre un gradient bathymétrique avec une structure spatiale horizontale ###


# Packages
library(readr)
library(dplyr)
library(tidyr)
library(sf)
library(spdep)
library(ggplot2)
library(purrr)

setwd(WORKING_DIR)

df <- read.csv(
  JOINED_BEN_ENV_DATA_FILE_W_KMEANS_CLUSTERS,
  header = TRUE,
  sep = ","
)

# Check data
str(df)
table(df$community)
table(df$Depth)


#### Aggregation by site and depth
# On crée une observation par combinaison Region + Site + Depth.
# On calcule ensuite la proportion de transects appartenant à chacune des cinq communautés
# comm_1, comm_2, comm_3, comm_4, comm_5 représentent les proportions de chaque communauté au sein de chaque localisation-profondeur.

# Create unique location identifier
df <- df %>%
  mutate(
    Site_ID = paste(Region, Site, Depth, sep = "_"),
    community = as.integer(community)
  )


# 1. One geographical point per site and depth
site_coords <- df %>%
  group_by(Site_ID, Region, Site, Depth) %>%
  summarise(
    Latitude = mean(Latitude, na.rm = TRUE),
    Longitude = mean(Longitude, na.rm = TRUE),
    n_transects = n(),
    .groups = "drop"
  )


# 2. Count transects per community
community_data <- df %>%
  filter(!is.na(community)) %>%
  count(Site_ID, community) %>%
  
  # Include all five communities
  complete(
    Site_ID,
    community = 1:5,
    fill = list(n = 0)
  ) %>%
  
  group_by(Site_ID) %>%
  mutate(
    proportion = n / sum(n)
  ) %>%
  
  ungroup() %>%
  
  select(Site_ID, community, proportion) %>%
  
  pivot_wider(
    names_from = community,
    values_from = proportion,
    names_prefix = "comm_",
    values_fill = 0
  )


# 3. Add coordinates and metadata
community_data <- site_coords %>%
  left_join(community_data, by = "Site_ID")

# Check
nrow(community_data)

summary(
  rowSums(
    community_data %>%
      select(starts_with("comm_")),
    na.rm = TRUE
  )
)


#### Spatial data analysis
# Elle calcule :
# 1. Global Moran's I.
# 2. Local Moran's I.
# 3. P-values par permutations.
# 4. Correction des tests multiples.
# 5. Classification LISA.
# distances géodésiques pour construire les voisins KNN, plutôt que de choisir une projection métrique unique qui pourrait introduire des distorsions sur cette étendue géographique.

# Moran's I positif et significatif : la proportion de cette communauté tend à présenter une distribution spatiale regroupée.
# Moran's I négatif et significatif : les valeurs similaires tendent à être spatialement séparées.
# Non significatif : aucune autocorrélation globale n'est détectée avec la matrice de voisinage utilisée.


# Convert to spatial data
sites_sf <- st_as_sf(
  community_data,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)

# Check spatial distribution
plot(st_geometry(sites_sf))


# Function for spatial analysis
run_spatial_analysis <- function(data, k = 5) {
  
  # Coordinates
  coords <- st_coordinates(data)
  
  # Check number of observations
  if (nrow(data) < 3) {
    stop("Not enough sites for spatial analysis.")
  }
  
  # Number of neighbours
  k <- min(k, nrow(data) - 1)
  
  # Check duplicated coordinates
  if (any(duplicated(as.data.frame(coords)))) {
    stop(
      paste(
        "Duplicated coordinates detected at depth",
        unique(data$Depth),
        "- check site locations before calculating KNN."
      )
    )
  }
  
  # KNN neighbours (geodesic distances)
  knn <- knearneigh(
    coords,
    k = k,
    longlat = TRUE
  )
  
  nb <- knn2nb(knn)
  
  # Spatial weights
  lw <- nb2listw(
    nb,
    style = "W",
    zero.policy = TRUE
  )
  
  # Community variables
  communities <- paste0("comm_", 1:5)
  
  global_results <- list()
  local_results <- list()
  
  for (com in communities) {
    
    x <- data[[com]]
    
    # Skip variables without variation
    if (
      all(is.na(x)) ||
      is.na(sd(x, na.rm = TRUE)) ||
      sd(x, na.rm = TRUE) == 0
    ) {
      next
    }
    
    # -------------------------
    # GLOBAL MORAN'S I
    # -------------------------
    
    set.seed(123)
    
    global <- moran.mc(
      x,
      lw,
      nsim = 9999,
      zero.policy = TRUE
    )
    
    global_results[[com]] <- data.frame(
      Depth = data$Depth[1],
      Community = com,
      Moran_I = as.numeric(global$statistic),
      P_value = global$p.value
    )
    
    
    # -------------------------
    # LOCAL MORAN'S I
    # -------------------------
    
    set.seed(123)
    
    lisa <- localmoran_perm(
      x,
      lw,
      nsim = 9999,
      zero.policy = TRUE
    )
    
    # Extract permutation p-value
    p_col <- grep(
      "^Pr",
      colnames(lisa),
      value = TRUE
    )[1]
    
    if (is.na(p_col)) {
      stop("Could not identify the permutation p-value column.")
    }
    
    p_raw <- lisa[, p_col]
    
    # Spatial lag of standardized values
    z <- as.numeric(scale(x))
    
    lag_z <- lag.listw(
      lw,
      z,
      zero.policy = TRUE
    )
    
    # LISA quadrants
    local_results[[com]] <- data.frame(
      Site_ID = data$Site_ID,
      Depth = data$Depth,
      Community = com,
      Local_I = lisa[, "Ii"],
      P_value = p_raw,
      z = z,
      lag_z = lag_z
    )
    
  }
  
  # Combine global results
  global_df <- bind_rows(global_results)
  
  # Combine local results
  local_df <- bind_rows(local_results)
  
  # Correct multiple testing across all local tests
  # within the current depth
  local_df <- local_df %>%
    mutate(
      P_adjusted = p.adjust(
        P_value,
        method = "BH"
      ),
      
      LISA_cluster = case_when(
        
        P_adjusted >= 0.05 ~ "Not significant",
        
        z > 0 & lag_z > 0 ~ "High-High",
        
        z < 0 & lag_z < 0 ~ "Low-Low",
        
        z > 0 & lag_z < 0 ~ "High-Low",
        
        z < 0 & lag_z > 0 ~ "Low-High",
        
        TRUE ~ "Not classified"
      )
    )
  
  return(
    list(
      global = global_df,
      local = local_df
    )
  )
}


# Exécuter les analyses pour chaque profondeur
# Available depths
depths <- sort(unique(sites_sf$Depth))

# Run analysis separately for each depth
results <- list()

for (d in depths) {
  
  data_depth <- sites_sf %>%
    filter(Depth == d)
  
  # Skip groups with too few sites
  if (nrow(data_depth) < 10) {
    message("Skipping depth: ", d, " m (too few sites)")
    next
  }
  
  message("Analysing depth: ", d, " m")
  
  results[[as.character(d)]] <-
    run_spatial_analysis(
      data_depth,
      k = 5
    )
}


# Combine global Moran results
global_results <- map_dfr(
  results,
  "global"
)

# Combine LISA results
local_results <- map_dfr(
  results,
  "local"
)

# Display global Moran results
print(global_results)


# write.csv(
#   global_results,
#   "global_moran_results.csv",
#   row.names = FALSE
# )
# 
# write.csv(
#   local_results,
#   "local_moran_results.csv",
#   row.names = FALSE
# )



#### Visualisation of results

# Carte de distribution des cinq communautés (communauté dominante dans chaque localisation-profondeur)
# Attention : la communauté dominante ne représente qu'une simplification des proportions. Un site avec 40 % de communauté 1 et 35 % de communauté 2 sera représenté comme communauté 1.

# Identify dominant community
sites_map <- community_data %>%
  mutate(
    dominant_community = max.col(
      as.matrix(
        select(., starts_with("comm_"))
      )
    )
  )

# Convert to sf
sites_map_sf <- st_as_sf(
  sites_map,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)

# Plot
ggplot(sites_map_sf) +
  
  geom_sf(
    aes(color = factor(dominant_community)),
    size = 2.5
  ) +
  
  facet_wrap(~Depth) +
  
  theme_minimal() +
  
  labs(
    title = "Spatial distribution of benthic communities",
    color = "Community"
  )




# Carte des proportions de chaque communauté

# Convert to long format
community_long <- community_data %>%
  
  pivot_longer(
    cols = starts_with("comm_"),
    names_to = "Community",
    values_to = "Proportion"
  )

# Convert to sf
community_sf <- st_as_sf(
  community_long,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)

# Plot
ggplot(community_sf) +
  
  geom_sf(
    aes(color = Proportion),
    size = 2.5
  ) +
  
  facet_grid(
    Community ~ Depth
  ) +
  
  scale_color_viridis_c() +
  
  theme_minimal() +
  
  labs(
    title = "Relative abundance of benthic communities",
    color = "Proportion"
  )




# Carte LISA : localisation des clusters

# Merge LISA results with spatial data
lisa_map <- sites_sf %>%
  
  left_join(
    local_results,
    by = c("Site_ID", "Depth")
  )

# Plot
ggplot(lisa_map) +
  
  geom_sf(
    aes(color = LISA_cluster),
    size = 2.5
  ) +
  
  facet_grid(
    Community ~ Depth
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Local spatial autocorrelation of benthic communities",
    color = "LISA cluster"
  )




#####################################################
#### EXPORT ALL RESULTS FOR QGIS AND GEODA ##########
#####################################################

# Create output directory
output_dir <- "outputs_r_code"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


#### 1. EXPORT GLOBAL MORAN RESULTS

write.csv(
  global_results,
  file.path(output_dir, "global_moran_results.csv"),
  row.names = FALSE
)


#### 2. EXPORT LOCAL MORAN / LISA RESULTS

# Export as CSV
write.csv(
  local_results,
  file.path(output_dir, "lisa_results.csv"),
  row.names = FALSE
)


#### 3. EXPORT ALL SITES WITH COMMUNITY PROPORTIONS

# Export all sites and their five community proportions
st_write(
  sites_sf,
  file.path(output_dir, "sites_all_depths.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)


#### 4. EXPORT DOMINANT COMMUNITY MAP

# Identify dominant community
sites_map <- community_data %>%
  mutate(
    dominant_community = max.col(
      as.matrix(
        select(., starts_with("comm_"))
      )
    )
  )

# Convert to sf
sites_map_sf <- st_as_sf(
  sites_map,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)

# Export
st_write(
  sites_map_sf,
  file.path(output_dir, "dominant_community.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)


#### 5. EXPORT COMMUNITY PROPORTIONS

# Convert to long format
community_long <- community_data %>%
  pivot_longer(
    cols = starts_with("comm_"),
    names_to = "Community",
    values_to = "Proportion"
  )

# Convert to sf
community_sf <- st_as_sf(
  community_long,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)

# Export
st_write(
  community_sf,
  file.path(output_dir, "community_proportions.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)


#### 6. EXPORT LISA RESULTS WITH GEOMETRIES

# Merge LISA results with spatial data
lisa_map <- sites_sf %>%
  left_join(
    local_results,
    by = c("Site_ID", "Depth")
  )

# Export
st_write(
  lisa_map,
  file.path(output_dir, "lisa_results.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)


#### 7. EXPORT SPATIAL DATA SEPARATELY FOR EACH DEPTH

for (d in depths) {
  
  data_depth <- sites_sf %>%
    filter(Depth == d)
  
  st_write(
    data_depth,
    file.path(
      output_dir,
      paste0("community_spatial_", d, "m.gpkg")
    ),
    delete_dsn = TRUE,
    quiet = TRUE
  )
}


#### 8. EXPORT LISA RESULTS SEPARATELY FOR EACH DEPTH

for (d in depths) {
  
  lisa_depth <- lisa_map %>%
    filter(Depth == d)
  
  st_write(
    lisa_depth,
    file.path(
      output_dir,
      paste0("lisa_results_", d, "m.gpkg")
    ),
    delete_dsn = TRUE,
    quiet = TRUE
  )
}


#### 9. EXPORT COMMUNITY PROPORTIONS SEPARATELY FOR EACH DEPTH

for (d in depths) {
  
  community_depth <- community_sf %>%
    filter(Depth == d)
  
  st_write(
    community_depth,
    file.path(
      output_dir,
      paste0("community_proportions_", d, "m.gpkg")
    ),
    delete_dsn = TRUE,
    quiet = TRUE
  )
}


#### 10. EXPORT DOMINANT COMMUNITY SEPARATELY FOR EACH DEPTH

for (d in depths) {
  
  dominant_depth <- sites_map_sf %>%
    filter(Depth == d)
  
  st_write(
    dominant_depth,
    file.path(
      output_dir,
      paste0("dominant_community_", d, "m.gpkg")
    ),
    delete_dsn = TRUE,
    quiet = TRUE
  )
}


#### 11. EXPORT SUMMARY OF GLOBAL MORAN RESULTS

global_summary <- global_results %>%
  mutate(
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01 ~ "**",
      P_value < 0.05 ~ "*",
      TRUE ~ "Not significant"
    )
  )

write.csv(
  global_summary,
  file.path(output_dir, "global_moran_summary.csv"),
  row.names = FALSE
)


#### 12. CHECK EXPORTED FILES

cat("\nAll results exported successfully!\n\n")

print(
  list.files(
    output_dir,
    full.names = TRUE
  )
)

