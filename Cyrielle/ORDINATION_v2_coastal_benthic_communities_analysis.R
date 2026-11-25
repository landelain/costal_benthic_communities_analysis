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
library(reshape2)
library(pheatmap)


### PATH DEFINITION ###
WORKING_DIR <- "~/Documents/EPFL/2026_MA_3/MULTI_R/costal_benthic_communities_analysis/Cyrielle"
DATA_DIR <- "~/Documents/EPFL/2026_MA_3/MULTI_R/costal_benthic_communities_analysis/Cyrielle"
BEN_DATA_FILE <- "benthic_clean.csv"
ENV_DATA_FILE <- "env_clean.csv"
JOINED_BEN_ENV_DATA_FILE <- "benthic_env_clean.csv"
JOINED_BEN_ENV_DATA_FILE_W_KMEANS_CLUSTERS <- "clusters_communities_kmeans.csv"

setwd(WORKING_DIR)


### LECTURE CSV ET FORMATAGE EN DATAFRAMES
### df      : dataframe clean complet (sites, env, abondances, cluster kmeans)  [transect x all]
### env     : variables environnementales (types mixtes)                        [transect x env]
### env_fac : variables environnementales catégorielles (factor)                [transect x env (cat)]
### env_num : variables environnementales numériques                            [transect x env (num)]
### spe     : abondances des taxons par transect                                [transect x taxons]
### xy      : longitude-latitude des transects                                  [transect x lon/lat]
### station : unité d'échantillonnage environnementale (Region_Site_Depth)      [vecteur, 1 par transect]
### idx_lookup : table de correspondance index court <-> nom du transect

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
sp_cols  <- col_range(df, "ac_erect_single", "zo_massive")
env_cols <- col_range(df, "Depth", "Management_status")

### NETTOYAGE COMMUN (mêmes lignes partout) ###
# Un site avec une abondance totale de 0 rend Bray-Curtis indéfini -> on le retire
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
zero_taxa <- names(spe)[colSums(spe) == 0] # Taxons jamais observés : inutiles dans la PCoA
cat("Taxons absents partout :", if (length(zero_taxa)) zero_taxa else "aucun", "\n")
spe <- spe[, colSums(spe) > 0]
summary(rowSums(spe))

### 2) VARIABLES ENVIRONNEMENTALES ###
env <- df[, env_cols]

# 2a. Coordonnées : à garder à part (spatial, pas un "facteur" environnemental)
xy  <- env[, c("Longitude", "Latitude")]

# 2b. Variables quantitatives
env_num <- env[, setdiff(names(env), c("Longitude", "Latitude", "Management_status"))]
stopifnot(all(sapply(env_num, is.numeric)))
sd_ok <- sapply(env_num, function(x) sd(x, na.rm = TRUE) > 0) # Variables sans variance -> à retirer (envfit ne peut rien en faire)
cat("Variables constantes retirées :", names(env_num)[!sd_ok], "\n")
env_num <- env_num[, sd_ok]
colSums(is.na(env_num))

# 2c. Variables qualitatives (projetées comme centroïdes par envfit)
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
par(mar = c(10, 4, 2, 1)) #règle les marges extérieures d'un graphique en nombre de lignes de texte, dans l'ordre bas, gauche, haut, droite
boxplot(env_num.scaled, las = 2, cex.axis = 0.7, main = "Variables env. standardisées")
par(mar = c(5, 4, 4, 2) + 0.1)


# Colinéarité, au niveau station (une ligne par station = la vraie réplication)
env_st <- aggregate(env_num, by = list(station = station), FUN = mean) #on prend la moyenne de tous les transects par station, on garde séparé par profondeur
rownames(env_st) <- env_st$station
env_st$station <- NULL
cor_env <- cor(env_st, use = "pairwise.complete.obs")
high_cor <- which(abs(cor_env) > 0.6 & upper.tri(cor_env), arr.ind = TRUE)
data.frame(var1 = rownames(cor_env)[high_cor[, 1]],
           var2 = colnames(cor_env)[high_cor[, 2]],
           r    = round(cor_env[high_cor], 2))
# pairs.panels(env_st)   # lisible seulement par groupes de ~8 variables

# ordre par clustering hiérarchique
pheatmap(cor_env,
         main = "Correlation entre variables env. et clustering hierarchique",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101),
         cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = TRUE, number_round = 2)






################################################################################
### PCoA sur spe [transect x taxons]
################################################################################

### from ex5: PCoA sur dataframe spe [transect v.s. species], pour identifier les species qui déterminent les (dis)similitudes entre les transects (scores des espèces)
### et puis a partir de cette PCoA une projection des variables environnementales
### pour identifier les variables environnementales qui expliquent les (dis)similitudes dans la composition taxonomique entre les échantillons 


### Commentaires Claude:
# Transformation des abondances. Ce sont des couvertures avec des ordres de grandeur très différents : turf et hc_encrusting sont autour de 5, beaucoup de taxons sont à moins de 0,1. Demande-toi s’il faut transformer spe avant de calculer les distances. Regarde ?decostand et compare "hellinger", "sqrt" et aucune transformation.
# Choix de la distance. Il faut choisir la distance avant la PCoA. Regarde aussi ce que fait l’argument add si des valeurs propres négatives apparaissent.
# Variables au niveau du site. Presque toutes tes variables environnementales ont la même valeur pour tous les transects d’un même site. Lance sapply(env_num, function(x) length(unique(x))) pour le voir. Les p-values d’envfit seront donc optimistes, car il y a pseudo-réplication. Regarde l’argument strata d’envfit.
# Depth. Si elle n’a que 2 ou 3 niveaux, il vaut peut-être mieux la passer dans env_fac en facteur.
# Corrélation entre variables. Beaucoup de tes variables sont probablement corrélées (SST, DHW, typhons, nutriments). pairs.panels(env_num) de psych, déjà chargé, te le montrera avant l’interprétation.




# STEPS:
# 1. Compute a dissimilarity matrix from the **SPECIES** dataset
# 2. Perform classical multidimensional scaling (PCoA) on the distance matrix
      # Return coordinates in 2 dimensions and include eigenvalues

# 3. Plot the ordination scores of ***SITES*** from the first two PCoA axes as text labels (projection des transects expliqués par les deux premieres PCoA axes)

# 4. Calculate weighted average scores of ***SPECIES based on site coordinates*** in PCoA space (projection des species qui expliquent ces deux premieres PCoA axes)

# 5.Perform envfit to ***FIT ENVIRONMENTAL VARIABLES*** onto the PCoA ordination
      # Scale the data (standardization) and print the boxplot of each variable in the env dataset.
      # envfit to fit environmental variables onto the PCoA ordination
      # Plot only significant environmental variables (p-value < 0.05) onto the ordination plot


# --- Distance --------------------------------------------------------------
# Distance de Hellinger = distance euclidienne sur les données transformées
# Hellinger. Elle est euclidienne : pas de valeurs propres négatives, et la
# PCoA donne le même résultat qu'une PCA (tb-PCA) sur spe.hel.
# Alternative à comparer : Bray-Curtis sur couvertures racine-carrées

DIST <- "hellinger"          # "hellinger" ou "bray"

if (DIST == "hellinger") {
  spe.tr    <- decostand(spe, method = "hellinger")    # données transformées
  spe.dist  <- dist(spe.tr)                            # distance de Hellinger
  add_const <- FALSE                                   # déjà euclidienne
  DIST_NAME <- "Hellinger"
} else if (DIST == "bray") {
  spe.tr    <- sqrt(spe)                               # données transformées
  spe.dist  <- vegdist(spe.tr, method = "bray")        # Bray-Curtis sur sqrt (turf, hc_encrusting et CCA dominent les abondances.
                                                       # La distance entre deux transects reflète alors surtout leurs différences de turf et de CCA.
                                                       # Une transformation qui compresse les grandes valeurs rééquilibre les contributions.)
  add_const <- TRUE                                    # corrige les valeurs propres négatives -> ajouter add = TRUE dans cmdscale
  DIST_NAME <- "Bray-Curtis (sqrt)"
}



# --- Ordination ------------------------------------------------------------

# Diagnostic (juste print, non utilisé par la suite): y a-t-il des valeurs propres négatives SANS correction ?
spe.pcoa0 <- cmdscale(spe.dist, k = 4, eig = TRUE)     # sans correction
eig_raw   <- spe.pcoa0$eig
tol       <- 1e-8 * max(eig_raw)
eig_pos   <- eig_raw[eig_raw >  tol]
eig_neg   <- eig_raw[eig_raw < -tol] # 0 attendu avec Hellinger

cat(DIST_NAME, "- valeurs propres négatives :", length(eig_neg),
    "| poids relatif :", round(100 * sum(abs(eig_neg)) / sum(eig_pos), 2), "%\n")
# Si ce poids est faible (quelques %), les négatives sont négligeables :
# pas besoin de correction. S'il est élevé, la représentation euclidienne
# de Bray-Curtis est mauvaise -> préférer Hellinger.


# La PCoA reellement effectuée:
USE_CORRECTION <- FALSE    # TRUE seulement si le poids des négatives est important

spe.pcoa <- if (USE_CORRECTION) {
  cmdscale(spe.dist, k = 4, eig = TRUE, add = TRUE)
} else spe.pcoa0

eig <- spe.pcoa$eig


# Nombre d'axes "réels" : avec Hellinger, il y a au plus (nb de species - 1) donc (27 species - 1) = 26 axes.
# La correction de Cailliez add = TRUE rend toute la matrice euclidienne, donc p ne vaut plus ~26 mais peut monter jusqu’à ~430.
# cmdscale renvoie quand même n valeurs propres, et les suivantes ne sont que du
# bruit numérique (~1e-15). Il faut les exclure, sinon le broken stick est faussé.


# Broken stick sur les valeurs propres positives NON corrigées
# (Hellinger : identique à avant ; Bray-Curtis : n'inclut pas la variation
#  artificielle ajoutée par la correction)
p       <- length(eig_pos)
var_exp <- 100 * eig_pos / sum(eig_pos)
round(var_exp[1:6], 1)
round(cumsum(var_exp)[1:6], 1) # % cumulé


# Broken stick : % de variation attendu pour chaque axe si la variation totale
# était répartie au hasard entre les p axes
bs <- bstick(p, tot.var = 100)
cat("Nombre d'axes positifs :", p, "\n")

# On garde les axes consécutifs dont la variation observée dépasse le broken stick
n_keep <- if (all(var_exp > bs)) p else which(var_exp <= bs)[1] - 1
cat(paste(DIST_NAME, "- Axes retenus (broken stick) :"), n_keep, "\n")

n_show <- min(10, p)
mids <- barplot(var_exp[1:n_show], names.arg = 1:n_show,
                xlab = "Axe PCoA", ylab = "% variation",
                ylim = c(0, max(var_exp[1], bs[1]) * 1.1),
                col = ifelse(seq_len(n_show) <= n_keep, "steelblue", "grey80"),
                border = NA,
                main = paste("Variation expliquée par axe PCoA", DIST_NAME,"vs broken stick"))
lines(mids, bs[1:n_show], type = "b", pch = 19, col = "red", lwd = 2)
legend("topright", legend = c("Observé (retenu)", "Observé", "Broken stick"),
       col = c("steelblue", "grey80", "red"), pch = c(15, 15, 19),
       lty = c(NA, NA, 1), bty = "n")

components2plot <- c(3,4)
site_sc <- spe.pcoa$points[, components2plot]
colnames(site_sc) <- c("PCoA1", "PCoA2")
ax_lab <- sprintf("PCoA%d (%.1f %%)", components2plot, var_exp[components2plot])




# --- Fond de carte réutilisable -------------------------------------------
# Chaque figure repart de zéro : plus de superposition avec le boxplot.
cols <- hcl.colors(nlevels(clust), "Dark 3") # couleurs des 5 comunities cluster

# sites     : afficher les points des transects
# spider    : relier chaque transect au centroïde de son cluster
# hull      : enveloppe convexe (polygone transparent) de chaque cluster
# centroids : afficher les centroïdes des clusters (par défaut : si spider ou hull)
# lims      : matrice de points (2 colonnes) qui définit le cadre du graphique
base_plot <- function(main = "", labels = FALSE, spider = FALSE, hull = FALSE,
                      sites = TRUE, centroids = spider | hull, lims = NULL) {
  # lims : matrice de points (2 colonnes) qui définit le cadre du graphique.
  # Par défaut, tous les transects (même cadrage pour toutes les figures) ;
  # sinon, par ex. les taxons + centroïdes pour zoomer.
  # type = "n" : on ne dessine rien, on fixe seulement le cadre.
  if (is.null(lims)) lims <- site_sc
  plot(lims, type = "n", asp = 1,              # asp = 1 : obligatoire en PCoA
       xlab = ax_lab[1], ylab = ax_lab[2], main = main)
  abline(h = 0, v = 0, lty = 3, col = "grey60")
  
  # Enveloppes convexes : polygone transparent par cluster.
  # alpha = transparence du remplissage (0 = invisible, 255 = opaque).
  # Tracées EN PREMIER pour rester en fond, sous l'araignée et les points.
  if (hull) {
    for (i in seq_len(nlevels(clust))) {
      ordihull(site_sc, groups = clust, show.groups = levels(clust)[i],
               draw = "polygon", col = cols[i], alpha = 50,
               border = cols[i], lwd = 1, label = FALSE)
    }
  }
  
  # Araignée : relie chaque transect au centroïde de son cluster k-means.
  # Tracée AVANT les points pour que les points restent visibles par-dessus.
  if (spider) {
    for (i in seq_len(nlevels(clust))) {
      ordispider(site_sc, groups = clust, show.groups = levels(clust)[i],
                 col = adjustcolor(cols[i], alpha.f = 0.4), lwd = 0.5)
    }
  }
  
  if (sites) {
    points(site_sc, pch = 19, cex = 0.6,
           col = adjustcolor(cols[clust], alpha.f = 0.6))
    if (labels) text(site_sc, labels = df$idx, cex = 0.5, pos = 3, col = "grey30")
  }
  
  # Centroïdes des clusters (après les points pour rester lisibles)
  if (centroids) {
    cent <- apply(site_sc, 2, function(x) tapply(x, clust, mean))
    points(cent, pch = 21, bg = cols, col = "black", cex = 1.8)
    text(cent, labels = levels(clust), cex = 0.7, font = 2, col = "white")
  }
  
  legend("topright", legend = levels(clust), col = cols, pch = 19,
         title = "k-means", bty = "n", cex = 0.8)
}

### FIGURE 1 : transects -------------------------------------------------------
#base_plot(paste("PCoA", DIST_NAME, "– transects"), spider = TRUE)               # araignée seule
base_plot(paste("PCoA", DIST_NAME, "– transects"), hull = TRUE, spider = TRUE)  # les deux

# Avec les index courts (utiliser idx_lookup pour retrouver le transect) :
#base_plot(paste("PCoA", DIST_NAME, "– transects"), labels = TRUE, spider = TRUE)



### FIGURE 2 : scores des taxons ------------------------------------------------
# Moyennes pondérées des taxons, calculées avec les MÊMES données que la distance
spe.wa  <- wascores(site_sc, spe.tr)

# Taxons les mieux corrélés aux axes (r² de la régression taxon ~ PCoA1 + PCoA2)
spe.fit <- envfit(site_sc, spe.tr, permutations = 0)
r2_sp   <- sort(spe.fit$vectors$r, decreasing = TRUE)
round(r2_sp, 2)
top_sp  <- names(r2_sp)[1:min(27, length(r2_sp))]   # à ajuster

# Seulement les centroïdes des clusters, sans transects ni araignée,
# pour que les noms des taxons restent lisibles
# Cadre ajusté aux taxons affichés + centroïdes des clusters, avec une marge de
# 10 % pour que les noms ne soient pas coupés au bord
cent  <- apply(site_sc, 2, function(x) tapply(x, clust, mean))
pts   <- rbind(spe.wa[top_sp, ], cent)
marge <- 0.10 * apply(pts, 2, function(x) diff(range(x)))
lims2 <- rbind(apply(pts, 2, min) - marge,
               apply(pts, 2, max) + marge)



base_plot(paste("PCoA", DIST_NAME, "– taxons principaux"),
          sites = FALSE, spider = FALSE, hull=TRUE, centroids = FALSE, lims = lims2)
#text(spe.wa[top_sp, ], labels = top_sp, cex = 0.8, col = "red", font = 1)
ordipointlabel(spe.wa[top_sp, ], display = "species", add = TRUE, col = "red", pch = c("+"), cex = 0.8)



### FIGURE 3 : projection des variables environnementales ----------------------
# Problème de pseudo-réplication : la plupart des variables sont mesurées une
# fois par station mais répétées sur chaque transect. Un envfit sur 433 lignes
# compte ces répétitions comme des observations indépendantes -> p-values
# beaucoup trop optimistes.
# Solution : tester les variables station au niveau station
# (moyenne des scores des transects de chaque station ; avec Hellinger/PCA, c'est
#  exactement le score de la composition moyenne de la station). Avec Bray-Curtis,
# ce n’est plus qu’une approximation, qui reste raisonnable pour placer les stations
# et tester les variables.

st_sc <- aggregate(as.data.frame(site_sc), by = list(station = station), FUN = mean)
rownames(st_sc) <- st_sc$station
st_sc$station <- NULL
stopifnot(identical(rownames(st_sc), rownames(env_st)))


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


# --- Fonction commune aux 3 versions ------------------------------------------
# 1. transforme en log(x+1) les variables de 'vars_log' (aucune par défaut)
# 2. envfit au niveau station sur toutes les variables station (999 permutations)
# 3. choix des vecteurs à afficher
#    p.max = 0.1 afficherait PLUS de vecteurs que 0.05, pas moins.
#    Critère combiné : significatif ET r² suffisant.
# 4. figure (étiquettes "log(...+1)" pour les variables transformées)
run_envfit_version <- function(env_data, vars_log = character(0), titre,
                               r2_min = 0.01, p_max = 0.05, n_max = 21) {
  if (length(vars_log) > 0)
    env_data[vars_log] <- lapply(env_data[vars_log], log1p)
  
  set.seed(1)
  fit <- envfit(st_sc, env_data[, vars_station, drop = FALSE],
                permutations = 999, na.rm = TRUE)
  
  r2   <- fit$vectors$r
  p    <- fit$vectors$pvals
  show <- names(r2)[p <= p_max & r2 >= r2_min]                 # seuils à ajuster
  show <- head(show[order(r2[show], decreasing = TRUE)], n_max)
  
  base_plot(paste("PCoA", DIST_NAME, "–", titre,
                  sprintf("(r² ≥ %.2f, p ≤ %.2f)", r2_min, p_max)), sites = FALSE, centroids=FALSE, spider=FALSE, hull=TRUE)
  if (length(show) > 0) {
    fit_plot <- envfit(st_sc, env_data[, show, drop = FALSE],
                       permutations = 0, na.rm = TRUE)
    lab <- ifelse(show %in% vars_log, paste0("log(", show, "+1)"), show)
    plot(fit_plot, add = TRUE, col = "blue", cex = 0.7,
         labels = list(vectors = lab))
  } else {
    message("Aucune variable ne passe les seuils pour : ", titre)
  }
  #plot(fit_fac_st, add = TRUE, col = "darkgreen", cex = 0.7, p.max = 0.05)
  #points(st_sc, pch = 3, cex = 0.6, col = "grey20")    # centroïdes des transects au sein de chaque site
  
  invisible(list(fit = fit, show = show, data = env_data, vars_log = vars_log))
}


################################################################################
### VERSION 1 : variables brutes (classique, sans transformation)

res_raw <- run_envfit_version(env_st, titre = "variables brutes")
res_raw$fit
res_raw$show
fit_env_st <- res_raw$fit   # nom conservé pour la suite du script

# Preuve que standardiser ne change rien
fit_chk <- envfit(st_sc, scale(env_st[, vars_station]), permutations = 0)
all.equal(unname(fit_chk$vectors$r), unname(fit_env_st$vectors$r))   # TRUE attendu




# ################################################################################
# ### VERSION 2 : TOUTES les variables transformées log(x+1) (test de robustesse)
# 
# # log1p n'est défini que pour x > -1 : on exclut les variables avec des valeurs négatives
# min_val      <- sapply(env_st[, vars_station], min, na.rm = TRUE)
# VARS_LOG_ALL <- vars_station[min_val >= 0]
# cat("Non transformées (valeurs négatives) :",
#     if (any(min_val < 0)) vars_station[min_val < 0] else "aucune", "\n")
# # Remarque : Depth est aussi transformée ici ; c'est voulu pour un test de
# # robustesse « tout transformé », mais ce n'est pas un choix à interpréter.
# 
# res_all <- run_envfit_version(env_st, vars_log = VARS_LOG_ALL,
#                               titre = "toutes variables log(x+1)")
# res_all$fit
# res_all$show
# 
# 
# ################################################################################
# ### VERSION 3 : seules les variables choisies sont transformées
# 
# # c) Asymétrie des variables (certains sites beaucoup plus exposés que d'autres)
# #    -> diagnostic qui a guidé le choix ci-dessous
# sk <- sapply(env_st[, vars_station], psych::skew, na.rm = TRUE)
# round(sort(sk, decreasing = TRUE), 2)
# 
# # Variables transformées log(x+1), choisies selon leur nature:
# # variables de pression humaine, qui sont des comptages et des densités, plus DHW et Nitrite, qui sont nettement asymétriques
# VARS_LOG <- c("Population_density", "Tourist_visitors",
#               "Anthropogenic_land_use", "Forest_land_use",
#               "DHW", "Nitrite")
# stopifnot(all(VARS_LOG %in% vars_station))
# 
# res_sel <- run_envfit_version(env_st, vars_log = VARS_LOG,
#                               titre = "variables choisies log(x+1)")
# res_sel$fit
# res_sel$show
# env_st_tr <- res_sel$data   # tableau transformé (version 3), ex. pour la heatmap
# 
# 
# ################################################################################
# ### ROBUSTESSE : comparaison des 3 versions
# ################################################################################
# cmp <- data.frame(
#   r2_brut = res_raw$fit$vectors$r, p_brut = res_raw$fit$vectors$pvals,
#   r2_all  = res_all$fit$vectors$r, p_all  = res_all$fit$vectors$pvals,
#   r2_sel  = res_sel$fit$vectors$r, p_sel  = res_sel$fit$vectors$pvals
# )
# cmp$affichee_dans <- sapply(rownames(cmp), function(v)
#   paste(c("brut", "all", "sel")[c(v %in% res_raw$show,
#                                   v %in% res_all$show,
#                                   v %in% res_sel$show)], collapse = "+"))
# cmp[order(-cmp$r2_sel), ] |> transform(r2_brut = round(r2_brut, 3),
#                                        r2_all  = round(r2_all, 3),
#                                        r2_sel  = round(r2_sel, 3))
# 
# 



# Pour le rapport:
# Présente la version transformée comme analyse principale, avec le tableau cmp en annexe.
# Distingue les variables robustes, significatives dans les deux versions (Light_intensity, Wave_exposure, Typhoon_frequency, Nitrate, Nitrite, Forest_land_use, et dans une moindre mesure Depth et DHW), des variables sensibles à la transformation (Tourist_visitors, Population_density, Anthropogenic_land_use). Pour ces dernières, dis clairement que le signal dépend de quelques stations extrêmes.
# Les r² restent modestes : 0,26 au maximum, sur seulement 2 axes. Chaque variable explique au mieux un quart de la position des stations dans ce plan, et beaucoup sont corrélées entre elles (ta heatmap). Ces flèches montrent des associations, pas quelle variable « détermine » la composition. Pour aller jusque-là, il faudrait la db-RDA avec sélection de variables dont on a parlé.


#Notes:
# Asymétrie des variables (certains sites beaucoup plus exposés que d'autres) --> Variables transformées log(x+1), choisies selon leur nature:
# variables de pression humaine ("Population_density", "Tourist_visitors", "Anthropogenic_land_use", "Forest_land_use"), qui sont des comptages et des densités
# plus "DHW" et "Nitrite", qui sont nettement asymétriques



### FIGURE 4 : heatmap des abondances – clusters choisis ------------------------
# Transects des clusters choisis, ordonnés le long d'un axe (gradient),
# taxons ordonnés par leur position moyenne le long de ce même ordre
# -> fait apparaître le remplacement des taxons le long du gradient (structure diagonale)
#
# Arguments :
# GROUPS   : clusters à comparer (dans l'ordre d'affichage si ORDER_BY = "Cluster")
# ORDER_BY : ordre des transects : "Cluster" (cluster puis profondeur),
#            "Depth", "PCoA3" ou "PCoA4"
# SCALE    : "max"  : chaque taxon ramené à son maximum (taxons rares visibles)
#            "none" : valeurs transformées utilisées pour la distance (spe.tr)
#
# Objets globaux utilisés : spe.dist, USE_CORRECTION, spe.tr, clust, df,
#                           DIST, DIST_NAME, cols

heatmap_clusters <- function(GROUPS = c("4", "5"), ORDER_BY = "Cluster", SCALE = "max") {
  
  # Coordonnées des transects sur les 4 premiers axes (mêmes options que la PCoA principale)
  pcoa4 <- cmdscale(spe.dist, k = 4, eig = TRUE, add = USE_CORRECTION)
  sc4   <- pcoa4$points
  colnames(sc4) <- paste0("PCoA", 1:4)
  
  # --- Sélection des transects (TOUS les taxons sont gardés) --------------------
  sel <- clust %in% GROUPS
  stopifnot(any(sel))
  M0  <- as.matrix(spe.tr[sel, , drop = FALSE])   # plus de filtre sur les taxons
  present <- colSums(M0) > 0                       # taxons présents dans ces clusters
  cat("Taxons absents des clusters", paste(GROUPS, collapse = " & "), ":",
      if (any(!present)) colnames(M0)[!present] else "aucun", "\n")
  dep <- df$Depth[sel]
  cl  <- clust[sel]
  
  # --- Ordre des transects -------------------------------------------------------
  # PCoA3 sert à départager les transects de même profondeur
  ord_s <- switch(ORDER_BY,
                  Cluster = order(match(as.character(cl), GROUPS), dep, sc4[sel, "PCoA3"]),  # dans l'ordre de GROUPS
                  Depth   = order(dep, sc4[sel, "PCoA3"]),
                  order(sc4[sel, ORDER_BY])                                                  # "PCoA3" ou "PCoA4"
  )
  M0  <- M0[ord_s, , drop = FALSE]
  dep <- dep[ord_s]
  cl  <- cl[ord_s]
  
  # --- Échelle des abondances ----------------------------------------------------
  if (SCALE == "max") {
    # division par le max de chaque taxon ; les taxons absents (max = 0)
    # restent à 0 au lieu de donner NaN
    mx     <- apply(M0, 2, max)
    M0     <- sweep(M0, 2, ifelse(mx > 0, mx, 1), "/")
    lab_ab <- "Abondance\n(/ max du taxon)"
  } else {
    lab_ab <- if (DIST == "hellinger") "Hellinger" else "sqrt(couverture)"
  }
  
  # --- Ordre des taxons : position moyenne pondérée le long de l'ordre des transects
  pos    <- seq_len(nrow(M0))
  sp_pos <- colSums(M0 * pos) / colSums(M0)
  sp_pos[!present] <- NA                           # pas de position pour les absents
  # na.last = FALSE -> taxons absents en premier = en BAS de la heatmap
  ord_sp      <- order(sp_pos, decreasing = TRUE, na.last = FALSE)
  M0          <- M0[, ord_sp, drop = FALSE]
  present_ord <- present[ord_sp]
  n_absent    <- sum(!present)
  
  M <- M0
  M[M == 0] <- NA                                 # absences en blanc
  ns  <- nrow(M)
  nsp <- ncol(M)
  
  # --- Couleurs ------------------------------------------------------------------
  pal_ab  <- hcl.colors(100, "YlOrRd", rev = TRUE)
  zlim_ab <- range(M, na.rm = TRUE)
  
  pal_d  <- hcl.colors(100, "Blues 3", rev = TRUE)   # clair = peu profond, foncé = profond
  zlim_d <- range(dep)
  if (diff(zlim_d) == 0) zlim_d <- zlim_d + c(-0.5, 0.5)
  depth_col <- function(v) pal_d[pmax(1, pmin(100,
                                              1 + floor(99 * (v - zlim_d[1]) / diff(zlim_d))))]
  d_un <- sort(unique(dep))
  dleg <- if (length(d_un) <= 6) d_un else {
    pr <- pretty(zlim_d, n = 4); pr[pr >= zlim_d[1] & pr <= zlim_d[2]]
  }
  
  # --- Séparateurs ---------------------------------------------------------------
  # Blocs = suites consécutives de transects de même cluster ET même profondeur
  cl_chr     <- as.character(cl)
  chg_clu    <- cl_chr[-1] != cl_chr[-ns]
  chg_dep    <- dep[-1] != dep[-ns]
  clu_breaks <- which(chg_clu) + 0.5                 # changements de cluster
  dep_breaks <- which(chg_dep | chg_clu) + 0.5       # changements de profondeur ou de cluster
  blk        <- cumsum(c(TRUE, chg_dep | chg_clu))   # identifiant de bloc
  show_blocks <- ORDER_BY %in% c("Cluster", "Depth") # avec PCoA3/4, trop de blocs pour les afficher
  
  # --- Mise en page --------------------------------------------------------------
  # 1 heatmap | 2 échelle d'abondance
  # 3 bande cluster | 4 légendes cluster + profondeur
  # 5 bande profondeur | 4
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); layout(1) })                  # réinitialise la mise en page, même en cas d'erreur
  layout(matrix(c(1, 2,
                  3, 4,
                  5, 4), ncol = 2, byrow = TRUE),
         widths = c(6, 1), heights = c(6, 0.45, 0.9))
  ML <- 9                                          # marge gauche (noms des taxons) : identique partout pour l'alignement
  
  # 1. Heatmap
  par(mar = c(0.3, ML, 3, 1))
  image(1:ns, 1:nsp, M, col = pal_ab, zlim = zlim_ab,
        axes = FALSE, xlab = "", ylab = "",
        main = sprintf("PCoA %s – clusters %s : abondances (transects ordonnés selon %s)",
                       DIST_NAME, paste(GROUPS, collapse = " & "), ORDER_BY))
  abline(h = 0.5 + 0:nsp, col = "grey90", lwd = 0.5)
  if (show_blocks) abline(v = dep_breaks, col = "grey40", lwd = 0.8)
  if (ORDER_BY == "Cluster") abline(v = clu_breaks, col = "black", lwd = 2.5)  # séparation des clusters
  if (n_absent > 0) abline(h = n_absent + 0.5, col = "grey40", lty = 2)
  # noms des taxons absents des clusters en gris
  mtext(colnames(M), side = 2, at = 1:nsp, las = 2, line = 0.3, adj = 1,
        cex = 0.45, col = ifelse(present_ord, "black", "grey60"))
  box()
  
  # 2. Échelle de couleur des abondances
  par(mar = c(0.3, 1, 3, 4))
  zk <- seq(zlim_ab[1], zlim_ab[2], length.out = 100)
  image(1, zk, matrix(zk, nrow = 1), col = pal_ab, axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1, cex.axis = 0.7)
  box()
  mtext(lab_ab, side = 3, line = 0.3, cex = 0.6)
  
  # 3. Bande des clusters
  par(mar = c(0.2, ML, 0.2, 1))
  image(1:ns, 1, matrix(match(cl_chr, levels(clust)), ncol = 1),
        col = cols, breaks = seq(0.5, nlevels(clust) + 0.5, 1),
        axes = FALSE, xlab = "", ylab = "")
  if (ORDER_BY == "Cluster") {
    abline(v = clu_breaks, col = "white", lwd = 2.5)
    clu_mid <- tapply(seq_len(ns), factor(cl_chr, levels = GROUPS), mean)
    text(clu_mid, 1, labels = names(clu_mid), cex = 0.7, font = 2, col = "white")
  }
  axis(2, at = 1, labels = "Cluster", las = 2, cex.axis = 0.7, tick = FALSE)
  box()
  
  # 4. Légendes (cluster + profondeur)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("topleft", legend = GROUPS, fill = cols[match(GROUPS, levels(clust))],
         title = "Cluster", bty = "n", cex = 0.7, xpd = NA)
  legend("bottomleft", legend = paste(dleg, "m"), fill = depth_col(dleg),
         title = "Depth", bty = "n", cex = 0.7, xpd = NA)
  
  # 5. Bande de profondeur
  par(mar = c(2, ML, 0.2, 1))
  image(1:ns, 1, matrix(dep, ncol = 1), col = pal_d, zlim = zlim_d,
        axes = FALSE, xlab = "", ylab = "")
  if (show_blocks) {
    abline(v = dep_breaks, col = "white", lwd = 0.8)
    # valeur de profondeur au centre de chaque bloc
    blk_mid <- tapply(seq_len(ns), blk, mean)
    blk_dep <- tapply(dep, blk, function(v) v[1])
    text(blk_mid, 1, labels = paste(blk_dep, "m"), cex = 0.6,
         col = ifelse(blk_dep > mean(zlim_d), "white", "black"))
  }
  if (ORDER_BY == "Cluster") abline(v = clu_breaks, col = "black", lwd = 2.5)
  axis(2, at = 1, labels = "Depth", las = 2, cex.axis = 0.7, tick = FALSE)
  box()
  mtext(sprintf("%d transects, ordonnés selon %s", ns, ORDER_BY),
        side = 1, line = 0.7, cex = 0.7)
  
  # renvoie (sans l'afficher) la matrice affichée, pour vérification
  invisible(list(M = M0, depth = dep, cluster = cl))
}


### Utilisation ---------------------------------------------------------------
#heatmap_clusters(GROUPS = c("4", "5"), ORDER_BY = "Cluster")   # cluster, puis profondeur
heatmap_clusters(GROUPS = c("4", "5"), ORDER_BY = "Depth")     # profondeur seule
#heatmap_clusters(GROUPS = c("4", "5"), ORDER_BY = "PCoA3")     # le long de PCoA3
#heatmap_clusters(GROUPS = c("5", "4"), ORDER_BY = "Cluster", SCALE = "none")

# Pour enregistrer (évite l'erreur « figure margins too large » dans RStudio) :
#png("heatmap_clusters_4_5.png", width = 12, height = 8, units = "in", res = 200)
#heatmap_clusters(GROUPS = c("4", "5"), ORDER_BY = "Cluster")
#dev.off()