# Packages
library(readr)
library(dplyr)
library(tidyr)
library(sf)
library(spdep)
library(ggplot2)
library(purrr)
library(psych)
library(gplots)
library(vegan)
library(ade4)


### PATH DEFINITION ###
WORKING_DIR <- "~/Documents/EPFL/2026_MA_3/MULTI_R/costal_benthic_communities_analysis/Cyrielle"
DATA_DIR <- "~/Documents/EPFL/2026_MA_3/MULTI_R/costal_benthic_communities_analysis/Cyrielle"
BEN_DATA_FILE <- "benthic_clean.csv"
ENV_DATA_FILE <- "env_clean.csv"
JOINED_BEN_ENV_DATA_FILE <- "benthic_env_clean.csv"
JOINED_BEN_ENV_DATA_FILE_W_KMEANS_CLUSTERS <- "clusters_communities_kmeans.csv"


setwd(WORKING_DIR)


### LECTURE CSV AND FOMATING INTO 6 DATAFRAMES
### df: le dataframe clean complet (sites, env variables, species abundances, cluster kmeans)  [transect v.s. all]
### env: le dataframe des variables environementales (mixed type)                              [transect v.s. var enviro]
### env_fac: le dataframe des variables environementales de type categorielles (str type)      [transect v.s. var enviro (cat)]
### env_num: le dataframe des variables environementales de type numeriques (double type)      [transect v.s. var enviro (num)]
### spe: le dataframe des abundances des species par transect                                  [transect v.s. species]
### xy: le dataframe des longitude-latitude des transects                                      [transect v.s. lon/lat]


df <- read.csv(
  file.path(DATA_DIR, JOINED_BEN_ENV_DATA_FILE_W_KMEANS_CLUSTERS),
  header = TRUE, sep = ",", stringsAsFactors = FALSE
)

str(df)
table(df$community)

### IDENTIFIANTS DES ÉCHANTILLONS ###
# vegan s'appuie sur les rownames : ils doivent être uniques
stopifnot(!anyDuplicated(df$Transect))
stopifnot(identical(df$Transect, df$Transect_2))   # Transect_2 est un doublon -> on l'ignore
rownames(df) <- df$Transect

### REPÉRAGE DES COLONNES ###
col_range <- function(data, from, to) {
  names(data)[which(names(data) == from):which(names(data) == to)]
}
sp_cols  <- col_range(df, "ac_erect_single", "zo_massive")   # 27 taxons / formes de croissance
env_cols <- col_range(df, "Depth", "Management_status")

### NETTOYAGE COMMUN (même lignes partout) ###
# Un site avec une abondance totale de 0 rend Bray-Curtis indéfini -> on le retire
row_tot <- rowSums(df[, sp_cols])
cat("Sites à abondance totale nulle :", sum(row_tot == 0), "\n")
cat("NA dans les abondances :", sum(is.na(df[, sp_cols])), "\n")
keep <- row_tot > 0 & complete.cases(df[, sp_cols])
df <- df[keep, ]

### 1) MATRICE SITES x ESPÈCES ###
spe <- df[, sp_cols]
stopifnot(all(sapply(spe, is.numeric)))

# Taxons jamais observés : inutiles dans la PCoA
zero_taxa <- names(spe)[colSums(spe) == 0]
cat("Taxons absents partout :", if (length(zero_taxa)) zero_taxa else "aucun", "\n")
spe <- spe[, colSums(spe) > 0]

summary(rowSums(spe))   # couverture totale par transect (ce n'est pas forcément 100)

### 2) VARIABLES ENVIRONNEMENTALES ###
env <- df[, env_cols]

# 2a. Coordonnées : à garder à part (spatial, pas un "facteur" environnemental)
xy <- env[, c("Longitude", "Latitude")]

# 2b. Variables quantitatives
env_num <- env[, setdiff(names(env), c("Longitude", "Latitude", "Management_status"))]
stopifnot(all(sapply(env_num, is.numeric)))

# Variables sans variance -> à retirer (envfit ne peut rien en faire)
sd_ok <- sapply(env_num, function(x) sd(x, na.rm = TRUE) > 0)
cat("Variables constantes retirées :", names(env_num)[!sd_ok], "\n")
env_num <- env_num[, sd_ok]
colSums(is.na(env_num))                     # NA restants par variable

# 2c. Variables qualitatives (projetées comme centroïdes par envfit)
env_fac <- data.frame(
  Management_status = factor(df$Management_status),
  Region            = factor(df$Region),
  Site              = factor(df$Site),
  row.names         = rownames(df)
)

### 3) CLUSTERS K-MEANS ###
clust <- factor(df$community)
names(clust) <- rownames(df)
table(clust)

### VÉRIFICATION FINALE : alignement des lignes ###
stopifnot(
  identical(rownames(spe), rownames(env_num)),
  identical(rownames(spe), rownames(env_fac)),
  identical(rownames(spe), rownames(xy)),
  identical(rownames(spe), names(clust))
)
dim(spe); dim(env_num); dim(env_fac)








### from ex5: PCoA sur dataframe spe [transect v.s. species], pour identifier les species qui déterminent les (dis)similitudes entre les transects (scores des espèces)
### et puis a partir de cette PCoA une projection des variables environnementales
### pour identifier les variables environnementales qui expliquent les (dis)similitudes dans la composition taxonomique entre les échantillons 


### Commentaires Claude:
# Transformation des abondances. Ce sont des couvertures avec des ordres de grandeur très différents : turf et hc_encrusting sont autour de 5, beaucoup de taxons sont à moins de 0,1. Demande-toi s’il faut transformer spe avant de calculer les distances. Regarde ?decostand et compare "hellinger", "sqrt" et aucune transformation.
# Choix de la distance. Il faut choisir la distance avant la PCoA. Regarde aussi ce que fait l’argument add si des valeurs propres négatives apparaissent.
# Variables au niveau du site. Presque toutes tes variables environnementales ont la même valeur pour tous les transects d’un même site. Lance sapply(env_num, function(x) length(unique(x))) pour le voir. Les p-values d’envfit seront donc optimistes, car il y a pseudo-réplication. Regarde l’argument strata d’envfit.
# Depth. Si elle n’a que 2 ou 3 niveaux, il vaut peut-être mieux la passer dans env_fac en facteur.
# Corrélation entre variables. Beaucoup de tes variables sont probablement corrélées (SST, DHW, typhons, nutriments). pairs.panels(env_num) de psych, déjà chargé, te le montrera avant l’interprétation.




# FROM EX5:
# Compute a Hellinger dissimilarity matrix from the species dataset, QUESTION TO MYSELF: would an other metric would be better ???
spe.bray <- vegdist(spe, dist="hellinger")

# Perform classical multidimensional scaling (PCoA) on the Bray-Curtis distance matrix
# Return coordinates in 2 dimensions and include eigenvalues
spe.b.pcoa <- cmdscale(spe.bray, k = 4, eig = TRUE)

#display the cmdscale output
spe.b.pcoa




# 1. Plot the ordination scores of ***SITES*** from the first two PCoA axes as text labels
ordiplot(scores(spe.b.pcoa, choices = c(3, 4)),
         type = "points",
         main = "PCoA with species weighted averages")



# 2. Calculate weighted average scores of ***SPECIES based on site coordinates*** in PCoA space
spe.wa <- wascores(spe.b.pcoa$points[, 1:2], spe)
# Add species names at their weighted average positions on the plot in red, smaller font
text(spe.wa, rownames(spe.wa), cex = 0.7, col = "red")



# 3.Perform envfit to ***FIT ENVIRONMENTAL VARIABLES*** onto the PCoA ordination

# Scale the data (standardization) and print the boxplot of each variable in the env dataset.
options(repr.plot.width = 18, repr.plot.height = 8)
env_num.scaled <- scale(env_num)

colMeans(env_num.scaled)  # check that we get mean of 0 and sd of 1 # faster version of apply(scaled.dat, 2, mean)
apply(env_num.scaled, 2, sd)
#boxplot(env_num.scaled, las=2)


# envfit to fit environmental variables onto the PCoA ordination
spe.b.pcoa.env_fac <- envfit(spe.b.pcoa, env_fac)
spe.b.pcoa.env_num.scaled <- envfit(spe.b.pcoa, env_num.scaled)
# Display the results of the envfit, including significance and vector directions
spe.b.pcoa.env_fac
spe.b.pcoa.env_num.scaled


# Plot only significant environmental variables (p-value < 0.05) onto the ordination plot
# Using color code 3 and 4 (typically green and blue)
plot(spe.b.pcoa.env_fac, p.max = 0.05, col = 3)
plot(spe.b.pcoa.env_num.scaled, p.max = 0.05, col = 4)



























