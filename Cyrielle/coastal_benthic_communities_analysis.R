# Script d'analyse des fichiers CSV
# Author : écrit par moi selon ce qu'on a fait dans les Lab 1 et 2, corrigé par LLM (Lumo AI par défaut, chatGPT si besoin) 
# Application : Benthic communities analysis

library(psych)
library(gplots)



### PATH DEFINITION ###
WORKING_DIR <- "~/Documents/EPFL/2026_MA_3/MULTI_R/costal_benthic_communities_analysis/Cyrielle"
DATA_DIR <- "~/Documents/EPFL/2026_MA_3/MULTI_R/costal_benthic_communities_analysis/data"
BEN_DATA_FILE <- "Linetal_dataset_Benthic.csv"
ENV_DATA_FILE <- "Linetal_dataset_Env.csv"

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


env_data <- open_and_analyze_data(ENV_DATA_FILE)



### PLOTS ###

# heatmap2 des species
options(repr.plt.width=10, repr.plt.height=10)

heatmap.2(as.matrix(log10(ben_data$data_num + 1)),
          dendrogram = "none",
          trace = "none",
          main = "Benthic species coverage")

# OBSERVATION: filamentous est surdominant, appliquer une transformation log sur cette donnée ?
#              --> voir les differentes transformations de datasets dans le cours semaine 2




# radarchart des differents parametres enviro des differents sites
# cf Lab 1 part 7.5
options(repr.plt.width=10, repr.plt.height=10)




