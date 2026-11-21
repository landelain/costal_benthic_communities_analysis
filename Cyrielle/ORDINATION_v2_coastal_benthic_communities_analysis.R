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
library(vegan)   # charge aussi 'permute'
library(ade4)


### PATH DEFINITION ###
WORKING_DIR <- "~/Documents/EPFL/2026_MA_3/MULTI_R/costal_benthic_communities_analysis/Cyrielle"
DATA_DIR <- "~/Documents/EPFL/2026_MA_3/MULTI_R/costal_benthic_communities_analysis/Cyrielle"
BEN_DATA_FILE <- "benthic_clean.csv"
ENV_DATA_FILE <- "env_clean.csv"
JOINED_BEN_ENV_DATA_FILE <- "benthic_env_clean.csv"
JOINED_BEN_ENV_DATA_FILE_W_KMEANS_CLUSTERS <- "clusters_communities_kmeans.csv"

setwd(WORKING_DIR)


### LECTURE CSV ET FORMATAGE EN DATAFRAMES
### df      : dataframe clean complet (sites, env, abondances, cluster kmeans)   [transect x all]
### env     : variables environnementales (types mixtes)                        [transect x env]
### env_fac : variables environnementales catégorielles (factor)                [transect x env (cat)]
### env_num : variables environnementales numériques                            [transect x env (num)]
### spe     : abondances des taxons par transect                                [transect x taxons]
### xy      : longitude-latitude des transects                                  [transect x lon/lat]
### station : unité d'échantillonnage environnementale (Region_Site_Depth)       [vecteur, 1 par transect]
### idx_lookup : table de correspondance index court <-> nom du transect

df <- read.csv(
  file.path(DATA_DIR, JOINED_BEN_ENV_DATA_FILE_W_KMEANS_CLUSTERS),
  header = TRUE, sep = ",", stringsAsFactors = FALSE
)

str(df)
table(df$community)

### IDENTIFIANTS DES ÉCHANTILLONS ###
stopifnot(!anyDuplicated(df$Transect))
stopifnot(identical(df$Transect, df$Transect_2))   # Transect_2 est un doublon -> ignoré
rownames(df) <- df$Transect

### REPÉRAGE DES COLONNES ###
col_range <- function(data, from, to) {
  names(data)[which(names(data) == from):which(names(data) == to)]
}
sp_cols  <- col_range(df, "ac_erect_single", "zo_massive")
env_cols <- col_range(df, "Depth", "Management_status")

### NETTOYAGE COMMUN (mêmes lignes partout) ###
row_tot <- rowSums(df[, sp_cols])
cat("Sites à abondance totale nulle :", sum(row_tot == 0), "\n")
cat("NA dans les abondances :", sum(is.na(df[, sp_cols])), "\n")
keep <- row_tot > 0 & complete.cases(df[, sp_cols])
df <- df[keep, ]

### NOUVEAU : INDEX COURT POUR L'AFFICHAGE ###
# Créé APRÈS le filtrage pour qu'il reste aligné. Les rownames restent le nom
# complet du transect (sécurité pour l'alignement) ; l'index sert uniquement
# d'étiquette sur les graphiques.
df$idx <- seq_len(nrow(df))
idx_lookup <- df[, c("idx", "Transect", "Region", "Site", "Depth", "community")]
# write.csv(idx_lookup, "transect_index_lookup.csv", row.names = FALSE)

### NOUVEAU : STATION = unité réelle de mesure de l'environnement ###
station <- factor(paste(df$Region, df$Site, df$Depth, sep = "_"))
cat("Nombre de stations :", nlevels(station), "\n")
summary(as.vector(table(station)))   # nb de transects par station

### 1) MATRICE SITES x ESPÈCES ###
spe <- df[, sp_cols]
stopifnot(all(sapply(spe, is.numeric)))
zero_taxa <- names(spe)[colSums(spe) == 0]
cat("Taxons absents partout :", if (length(zero_taxa)) zero_taxa else "aucun", "\n")
spe <- spe[, colSums(spe) > 0]
summary(rowSums(spe))

### 2) VARIABLES ENVIRONNEMENTALES ###
env <- df[, env_cols]
xy  <- env[, c("Longitude", "Latitude")]

env_num <- env[, setdiff(names(env), c("Longitude", "Latitude", "Management_status"))]
stopifnot(all(sapply(env_num, is.numeric)))
sd_ok <- sapply(env_num, function(x) sd(x, na.rm = TRUE) > 0)
cat("Variables constantes retirées :", names(env_num)[!sd_ok], "\n")
env_num <- env_num[, sd_ok]
colSums(is.na(env_num))

env_fac <- data.frame(
  Management_status = factor(df$Management_status),
  Region            = factor(df$Region),
  Site              = factor(df$Site),
  row.names         = rownames(df)
)

### NOUVEAU : quelles variables varient À L'INTÉRIEUR d'une station ? ###
# TRUE  = variable mesurée au niveau station (même valeur pour tous ses transects)
# FALSE = variable mesurée au niveau transect
is_station_level <- sapply(env_num, function(x)
  all(tapply(x, station, function(v) length(unique(v)) <= 1)))
print(is_station_level)
vars_station  <- names(env_num)[is_station_level]
vars_transect <- names(env_num)[!is_station_level]

# Management_status doit lui aussi être constant par station
stopifnot(all(tapply(as.character(env_fac$Management_status), station,
                     function(v) length(unique(v)) == 1)))

### 3) CLUSTERS K-MEANS ###
clust <- factor(df$community)
names(clust) <- rownames(df)
table(clust)

### VÉRIFICATION FINALE : alignement des lignes ###
stopifnot(
  identical(rownames(spe), rownames(env_num)),
  identical(rownames(spe), rownames(env_fac)),
  identical(rownames(spe), rownames(xy)),
  identical(rownames(spe), names(clust)),
  length(station) == nrow(spe)
)
dim(spe); dim(env_num); dim(env_fac)


################################################################################
### EXPLORATION DES VARIABLES ENVIRONNEMENTALES (figures séparées de la PCoA)
################################################################################

# Boxplot des variables standardisées : uniquement pour l'exploration visuelle.
# (envfit n'a PAS besoin de données standardisées : son r² et ses p-values
#  sont identiques avec ou sans scale().)
env_num.scaled <- scale(env_num)
par(mar = c(10, 4, 2, 1))
boxplot(env_num.scaled, las = 2, cex.axis = 0.7, main = "Variables env. standardisées")
par(mar = c(5, 4, 4, 2) + 0.1)

# Colinéarité, au niveau station (une ligne par station = la vraie réplication)
env_st <- aggregate(env_num, by = list(station = station), FUN = mean)
rownames(env_st) <- env_st$station
env_st$station <- NULL
cor_env <- cor(env_st, use = "pairwise.complete.obs")
high_cor <- which(abs(cor_env) > 0.8 & upper.tri(cor_env), arr.ind = TRUE)
data.frame(var1 = rownames(cor_env)[high_cor[, 1]],
           var2 = colnames(cor_env)[high_cor[, 2]],
           r    = round(cor_env[high_cor], 2))
# pairs.panels(env_st)   # lisible seulement par groupes de ~8 variables


################################################################################
### PCoA sur spe [transect x taxons]
################################################################################

# --- Distance --------------------------------------------------------------
# Distance de Hellinger = distance euclidienne sur les données transformées
# Hellinger. Elle est euclidienne : pas de valeurs propres négatives, et la
# PCoA donne le même résultat qu'une PCA (tb-PCA) sur spe.hel.
spe.hel  <- decostand(spe, method = "hellinger")
spe.dist <- dist(spe.hel)
# Alternative à comparer : Bray-Curtis sur couvertures racine-carrées
# spe.dist <- vegdist(sqrt(spe), method = "bray")   # -> ajouter add = TRUE dans cmdscale

# --- Ordination ------------------------------------------------------------
spe.pcoa <- cmdscale(spe.dist, k = 2, eig = TRUE)

eig <- spe.pcoa$eig
cat("Valeurs propres négatives :", sum(eig < -1e-8), "\n")   # 0 attendu avec Hellinger
var_exp <- 100 * eig / sum(eig[eig > 0])
round(var_exp[1:6], 1)
barplot(var_exp[1:10], names.arg = 1:10, ylab = "% variation",
        main = "Variation expliquée par axe PCoA")

site_sc <- spe.pcoa$points[, 1:2]
colnames(site_sc) <- c("PCoA1", "PCoA2")
ax_lab <- sprintf("PCoA%d (%.1f %%)", 1:2, var_exp[1:2])

# --- Fond de carte réutilisable -------------------------------------------
# Chaque figure repart de zéro : plus de superposition avec le boxplot.
cols <- hcl.colors(nlevels(clust), "Dark 3")

base_plot <- function(main = "", labels = FALSE) {
  plot(site_sc, type = "n", asp = 1,           # asp = 1 : obligatoire en PCoA
       xlab = ax_lab[1], ylab = ax_lab[2], main = main)
  abline(h = 0, v = 0, lty = 3, col = "grey60")
  points(site_sc, pch = 19, cex = 0.6,
         col = adjustcolor(cols[clust], alpha.f = 0.6))
  if (labels) text(site_sc, labels = df$idx, cex = 0.35, pos = 3, col = "grey30")
  legend("topright", legend = levels(clust), col = cols, pch = 19,
         title = "k-means", bty = "n", cex = 0.8)
}


### FIGURE 1 : transects -------------------------------------------------------
base_plot("PCoA (Hellinger) – transects")
# Avec les index courts (utiliser idx_lookup pour retrouver le transect) :
# base_plot("PCoA (Hellinger) – transects", labels = TRUE)


### FIGURE 2 : scores des taxons ------------------------------------------------
# Moyennes pondérées des taxons, calculées avec les MÊMES données que la distance
spe.wa <- wascores(site_sc, spe.hel)

# Taxons les mieux corrélés aux axes (r² de la régression taxon ~ PCoA1 + PCoA2)
spe.fit <- envfit(site_sc, spe.hel, permutations = 0)
r2_sp   <- sort(spe.fit$vectors$r, decreasing = TRUE)
round(r2_sp, 2)
top_sp  <- names(r2_sp)[1:10]                # à ajuster

base_plot("PCoA (Hellinger) – taxons principaux")
text(spe.wa[top_sp, ], labels = top_sp, cex = 0.7, col = "red", font = 2)


### FIGURE 3 : projection des variables environnementales ----------------------
# Problème de pseudo-réplication : la plupart des variables sont mesurées une
# fois par station mais répétées sur chaque transect. Un envfit sur 433 lignes
# compte ces répétitions comme des observations indépendantes -> p-values
# beaucoup trop optimistes.
# Solution : tester les variables station au niveau station
# (moyenne des scores des transects de chaque station ; avec Hellinger/PCA, c'est
#  exactement le score de la composition moyenne de la station).

st_sc <- aggregate(as.data.frame(site_sc), by = list(station = station), FUN = mean)
rownames(st_sc) <- st_sc$station
st_sc$station <- NULL
stopifnot(identical(rownames(st_sc), rownames(env_st)))

set.seed(1)
fit_env_st <- envfit(st_sc, env_st[, vars_station, drop = FALSE],
                     permutations = 999, na.rm = TRUE)
fit_env_st

# Variable catégorielle (une valeur par station)
fac_st <- data.frame(
  Management_status = factor(tapply(as.character(env_fac$Management_status),
                                    station, function(v) v[1])[rownames(st_sc)])
)
set.seed(1)
fit_fac_st <- envfit(st_sc, fac_st, permutations = 999)
fit_fac_st
# Region et Site ne sont pas testés : ce sont la structure d'échantillonnage
# elle-même (ils "expliquent" toujours la composition, par construction).

# Variables qui varient entre transects d'une même station : test au niveau
# transect, permutations restreintes à l'intérieur des stations
if (length(vars_transect) > 0) {
  set.seed(1)
  fit_env_tr <- envfit(site_sc, env_num[, vars_transect, drop = FALSE],
                       permutations = how(blocks = station, nperm = 999),
                       na.rm = TRUE)
  print(fit_env_tr)
}

# --- Choix des vecteurs à afficher --------------------------------------------
# p.max = 0.1 afficherait PLUS de vecteurs que 0.05, pas moins.
# Critère combiné : significatif ET r² suffisant.
r2_env <- fit_env_st$vectors$r
p_env  <- fit_env_st$vectors$pvals
show_env <- names(r2_env)[p_env <= 0.05 & r2_env >= 0.3]   # seuil à ajuster
show_env <- head(show_env[order(r2_env[show_env], decreasing = TRUE)], 8)
show_env

fit_env_plot <- envfit(st_sc, env_st[, show_env, drop = FALSE], permutations = 0)

base_plot("PCoA (Hellinger) – variables environnementales")
points(st_sc, pch = 3, cex = 0.6, col = "grey20")    # centroïdes des stations
plot(fit_env_plot, add = TRUE, col = "blue", cex = 0.7)
plot(fit_fac_st, add = TRUE, col = "darkgreen", cex = 0.7, p.max = 0.05)