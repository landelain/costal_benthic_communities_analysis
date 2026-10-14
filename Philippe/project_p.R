# =========================================================================
# 1. SETUP AND DATA LOADING
# =========================================================================
library(vegan)
library(gplots)
library(gclus)

# Set your working directory
setwd("C:/Users/phili/OneDrive/Desktop/EPFL/ma3/multivar_r/project/costal_benthic_communities_analysis/data")

# Load datasets
benthic_data <- read.csv("Linetal_dataset_Benthic.csv", head=TRUE, sep=",", row.names = 1)
cor_data <- read.csv("Linetal_dataset_Cor.csv", head=TRUE, sep=",")
env_data <- read.csv("Linetal_dataset_Env.csv", head=TRUE, sep=",")

# Extract spatial coordinates (assuming they are in the Env dataset)
spa <- env_data[, c("Longitude", "Latitude")]

# =========================================================================
# 2. BASIC EXPLORATION & SUMMARY
# =========================================================================
head(benthic_data)
head(env_data)

# Benthic data is zero-inflated and right-skewed. 
# Calculate exact proportion of zeros to justify transformation later:
prop_zeros <- sum(benthic_data == 0) / (nrow(benthic_data) * ncol(benthic_data))
print(paste("Proportion of zeros in Benthic dataset:", round(prop_zeros, 2)))

# env_data 
# Calculate exact proportion of zeros to justify transformation later:
prop_zeros <- sum(env_data == 0) / (nrow(env_data) * ncol(env_data))
print(paste("Proportion of zeros in Benthic dataset:", round(prop_zeros, 2)))

summary(benthic_data)  
summary(env_data)
# =========================================================================
# 3. DISTRIBUTION CHECKS & CRITICAL TRANSFORMATIONS
# =========================================================================
# Plot raw histograms to identify skewness

options(repr.plot.width = 15, repr.plot.height = 12) # Slightly taller for 5 rows
par(mfrow = c(5, 5), mar = c(2, 2, 2, 1)) # 5 rows by 5 columns

env_numeric <- env_data[sapply(env_data, is.numeric)]

for (col_name in colnames(env_numeric)) {
  hist(env_numeric[[col_name]], main = col_name, xlab = "")
}

# distribution for benthic communities
# Find the top 20 most abundant benthic categories
top_species <- names(sort(apply(benthic_data, 2, sum), decreasing = TRUE)[1:20])

# Set up a 4x5 grid to fit exactly 20 plots
options(repr.plot.width = 16, repr.plot.height = 12)
par(mfrow = c(4, 5), mar = c(3, 3, 2, 1)) 

for (sp in top_species) {
  hist(benthic_data[[sp]], 
       main = sp, 
       xlab = "", 
       ylab = "",
       breaks = 20)
}

# Critical Step: Log-transformation of highly right-skewed variables.
# Variables like Nitrite and Anthropogenic_land_use show severe skew, 
# which violates the linearity assumptions of many multivariate models.
env_numeric$Nitrite_log <- log1p(env_numeric$Nitrite)
env_numeric$Anth_land_use_log <- log1p(env_numeric$Anthropogenic_land_use)

# Standardize all numerical environmental variables (z-scores) for distance calculations
env_scaled <- as.data.frame(scale(env_numeric))

# =========================================================================
# 4. BIVARIATE RELATIONSHIPS (PAIRS PLOT)
# =========================================================================
# Custom functions for pairs plot
panel.cor2 <- function(x, y){
  usr <- par("usr")
  on.exit(par(usr = usr))
  par(usr = c(0, 1, 0, 1))
  
  # Using Kendall tau because environmental data often has outliers / non-normality
  r <- round(cor(x, y, method = "kendall", use = "pairwise.complete.obs"), digits=2)
  txt <- paste0("tau = ", r)
  text(0.5, 0.5, txt, cex = 2)
}

panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr = usr))
  par(usr = c(usr[1:2], 0, 1.5))
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y / max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

# Select a subset of standardized, critical variables to visualize
# (e.g., Depth, Light, Wave Exposure, and our newly transformed variables)
vars_to_plot <- c("Depth", "Light_intensity", "Wave_exposure", "Nitrite_log", "Anth_land_use_log")
options(repr.plot.width = 20, repr.plot.height = 14)
pairs(env_scaled[, vars_to_plot], 
      panel = panel.smooth, 
      diag.panel = panel.hist, 
      upper.panel = panel.cor2, 
      cex.labels = 3, 
      main = "Bivariate Plots with Histograms and Smooth Curves (Standardized)"
)

# =========================================================================
# 5. SPATIAL VISUALIZATION (RICHNESS)
# =========================================================================
# Compute species richness per site
sit.pres <- apply(benthic_data > 0, 1, sum)

options(repr.plot.width = 10, repr.plot.height = 10)
par(mfrow = c(1, 1))
plot(spa$Longitude, spa$Latitude, 
     asp = 1, 
     main = "Geographic Distribution of Benthic Richness", 
     pch = 21, 
     col = "darkgrey", 
     bg = "orange", 
     cex = 3 * sit.pres / max(sit.pres), 
     xlab = "Longitude", 
     ylab = "Latitude"
)

# =========================================================================
# 6. DISSIMILARITY MATRICES & VISUALIZATION
# =========================================================================
image.real <- function(mat) { 
  mat <- t(mat)[,nrow(mat):1]
  image(mat, axes = FALSE, col = hcl.colors(15, palette="viridis"))
  # Suppressing axis labels if the matrix is huge (e.g., 400+ rows) to avoid black overlapping text
  if(nrow(mat) < 100) {
    axis(1, at = seq(0, 1, length = nrow(mat)), labels = rownames(mat), las=2)
    axis(2, at = seq(0, 1, length = ncol(mat)), labels = colnames(mat))
  }
  box() 
}

# Critical Step: Because of the high proportion of zeros (double-zero problem), 
# directly applying Bray-Curtis to raw abundance can falsely equate sites that share absences.
# We apply a Hellinger transformation first, which gives low weights to rare species.
spe.hel <- decostand(benthic_data, "hel")

# -------------------------------------------------------------------------
# 6A. Q-MODE: SITE-BY-SITE DISSIMILARITY
# -------------------------------------------------------------------------
# Compute Bray-Curtis on Hellinger-transformed data
benthic_data.db <- vegdist(spe.hel, method = "bray")

options(repr.plot.width = 8, repr.plot.height = 8)
image.real(as.matrix(benthic_data.db))
title("Bray-Curtis Dissimilarity Across Sites")

# -------------------------------------------------------------------------
# 6B. R-MODE: SPECIES CO-OCCURRENCE (CLUSTERED HEATMAP)
# -------------------------------------------------------------------------
# Transpose the Hellinger-transformed data so species are rows
spe_t <- t(spe.hel)

# Calculate the Euclidean distance between species profiles
species_dist <- dist(spe_t)

# Create a Clustered Heatmap for the 27 species
options(repr.plot.width = 12, repr.plot.height = 12)
heatmap.2(as.matrix(species_dist), 
          Rowv = TRUE, Colv = TRUE,      
          dendrogram = "both",           
          trace = "none",                
          col = hcl.colors(50, palette="viridis"), 
          margins = c(14, 14),           
          cexRow = 1.2, cexCol = 1.2,    
          main = "Clustered Species Dissimilarity Matrix\n(Co-occurrence Profiles)")

# -------------------------------------------------------------------------
# Calculate Environmental Euclidean Distance (Prep for Mantel Test)
env.de <- dist(env_scaled)

# =========================================================================
# 7. MANTEL TEST
# =========================================================================
# Does environmental similarity dictate biological similarity?
mantel_result <- mantel(benthic_data.db, env.de, method = "spearman")
print(mantel_result)

# =========================================================================
# 8. GEO-ECOLOGY: SPATIAL AUTOCORRELATION & DISTANCE DECAY
# =========================================================================
library(vegan)

# 1. Create a spatial distance matrix (using Longitude and Latitude)
# Note: For strict geographic accuracy over large scales, you'd use 
# Haversine distance, but Euclidean is standard for basic numerical ecology classes.
spa <- env_data[, c("Longitude", "Latitude")]
spa_dist <- dist(spa)

# 2. Ensure you have your biological distance matrix ready 
# (Bray-Curtis on Hellinger-transformed data)
spe_hel <- decostand(benthic_data, "hel")
spe_dist <- vegdist(spe_hel, method = "bray")

# -------------------------------------------------------------------------
# Plot 1: The Classic Distance-Decay Scatterplot
# -------------------------------------------------------------------------
options(repr.plot.width = 8, repr.plot.height = 6)
par(mfrow = c(1, 1))

# Plot biological dissimilarity against geographic distance
plot(spa_dist, spe_dist, 
     pch = 16, col = rgb(0.1, 0.5, 0.8, alpha = 0.3), # semi-transparent dots
     xlab = "Geographic Distance (Euclidean)", 
     ylab = "Bray-Curtis Dissimilarity",
     main = "Distance-Decay Relationship of Benthic Communities")

# Add a smoothing line (lowess) to highlight the trend
lines(lowess(spa_dist, spe_dist), col = "red", lwd = 3)

# -------------------------------------------------------------------------
# Plot 2: Mantel Correlogram (The Numerical Ecology Standard)
# -------------------------------------------------------------------------
# This groups the distances into classes and tests if the correlation 
# is significant at each specific distance scale.
correlog <- mantel.correlog(spe_dist, spa_dist)

options(repr.plot.width = 10, repr.plot.height = 6)
plot(correlog, main = "Mantel Correlogram: Spatial Autocorrelation")

# =========================================================================
# 9. SPECIES-ENVIRONMENT CROSS-CORRELATION
# =========================================================================
# Note: We use the Hellinger-transformed species (spe.hel) and the 
# scaled environmental data (env_scaled) to ensure fair comparisons.

# 1. Calculate the cross-correlation matrix (Kendall tau)
# This correlates every species (columns of spe.hel) against every env variable (columns of env_scaled)
spe_env_cor <- cor(spe.hel, env_scaled, method = "kendall", use = "pairwise.complete.obs")

# (Optional) Print the top correlations for a specific species, e.g., Turf Algae
print("Top environmental correlations for Filamentous Turf:")
sort(spe_env_cor["turf_filamentous", ], decreasing = TRUE)[1:5]

# 1. Augmentez considérablement la hauteur du graphique (passez à 16 ou 18)
options(repr.plot.width = 14, repr.plot.height = 18)

# 2. Tracez la carte de chaleur avec une police plus petite pour les lignes
heatmap.2(spe_env_cor, 
          Rowv = TRUE, Colv = TRUE,       
          dendrogram = "both",            
          trace = "none",                 
          col = hcl.colors(50, palette = "Blue-Red 2"), 
          margins = c(12, 18),           # Gardez une grande marge à droite
          cexRow = 0.6,                  # <-- RÉDUISEZ CETTE VALEUR (ex: 0.5 ou 0.6)
          cexCol = 1.0,                  # Taille de la police pour les variables environnementales
          main = "Kendall Correlation: Species vs. Environment")
