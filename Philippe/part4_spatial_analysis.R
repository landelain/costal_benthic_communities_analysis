library(vegan)     # For Hellinger transform and nMDS
library(NbClust)   # For optimal k-means partitioning
library(ape)       # For Moran's I spatial autocorrelation
library(party)     # For Conditional Random Forests (CRF)
library(dplyr)     # For column selection
library(spdep)     # For Local Moran's I (LISA)
library(ggplot2)   # For mapping
library(spgwr)
library(sp)

setwd("C:/Users/phili/OneDrive/Desktop/EPFL/ma3/multivar_r/project/costal_benthic_communities_analysis/Philippe")

dat <- read.csv("benthic_env_clean.csv")

# Extract species columns 
species_data <- dat %>% select(ac_erect_single:zo_massive)

# Extract environmental variables
env_data <- dat %>% select(Depth:Management_status)


# 1. Hellinger + K-Means
# Hellinger transformation
species_hel <- decostand(species_data, method = "hellinger")

# Run NbClust
cluster_eval <- NbClust(species_hel, distance = "euclidean", 
                        min.nc = 2, max.nc = 10, method = "kmeans", index = "all")

# Extract the optimal cluster assignments and add to the main dataframe
dat$Community <- as.factor(cluster_eval$Best.partition)

# 2. Moran's I
# Inverse distance matrix
coords <- cbind(dat$Longitude, dat$Latitude)
dist_matrix <- as.matrix(dist(coords))
inv_dist_matrix <- 1 / (dist_matrix+0.000001)
diag(inv_dist_matrix) <- 0

moran_result <- Moran.I(as.numeric(dat$Community), weight = inv_dist_matrix)
print(moran_result)

# K-means: Group sites with similar species
# Moran's I: Do the sites that share similar species also happen to be located physically next to each other?
# RES: since p~=0 -> Yes it does -> Heavy spatial structure on the data!!

# 3. nMDS (Non-metric Multidimensional Scaling) -> Apparently it is a good way to visualise a clustering in very high dimension for ecological data
# Run nMDS on the Hellinger-transformed Euclidean distance matrix
nmds_result <- metaMDS(species_hel, distance = "euclidean", k = 2, trymax = 100)

# Plot nMDS colored by the k-means communities
plot(nmds_result$points, col = dat$Community, pch = 16, 
     main = "nMDS of Benthic Communities", xlab = "NMDS1", ylab = "NMDS2")
legend("topright", legend = levels(dat$Community), 
       col = 1:length(levels(dat$Community)), pch = 16, title = "Community")

# 4.LISA for only 1 species, maybe would be ok to do it with NMDS1 since it holds the principal info for the entire community?
# Need to convert the matrix to a list of weight for the localmoran
listw <- mat2listw(inv_dist_matrix, style = "W")

# Run Local Moran's I on the most abundant species: turf_filamentous
lisa_result <- localmoran(dat$turf_filamentous, listw)

dat$Local_I <- lisa_result[, "Ii"]
dat$Local_P <- lisa_result[, "Pr(z != E(Ii))"]

# Categorize sites for Turf Algae
dat$Turf_Hotspot <- ifelse(dat$Local_P < 0.05, "Significant Cluster", "Not Significant")

# 5. LISA Map

map_turf_lisa <- ggplot(dat, aes(x = Longitude, y = Latitude, 
                                 color = Turf_Hotspot, 
                                 size = turf_filamentous)) + 
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Not Significant" = "grey70", 
                                "Significant Cluster" = "firebrick")) +
  theme_minimal() +
  labs(
    title = "Spatial Hotspots of Turf Algae (turf_filamentous)",
    subtitle = "Red points indicate significant spatial clusters based on Local Moran's I",
    x = "Longitude", 
    y = "Latitude", 
    size = "Turf Abundance",
    color = "LISA Significance"
  )

print(map_turf_lisa)

# RES: Spatial distribution is heterogeneous, turf algae forms localized, significant patches (clusters) rather than a uniform spread

# 6. GWR with anthropogenic variables -> mangement_status, anthropogenic_land_use, Population_density, Tourist_visitors

# Problem with Management_status since it is a categorical var -> change it into a binary 0/1
# If it is "no-take", assign 1, otherwise 0
dat$is_notake <- ifelse(dat$Management_status == "no-take", 1, 0)

# If it is "restricted", assign 1, otherwise 0
dat$is_restricted <- ifelse(dat$Management_status == "restricted", 1, 0)
# (Note: We don't need an 'is_open' column, because if both of the above are 0, the model knows it is "open" by default), it would introduce redundancy

# copy of data
dat_sp <- dat
coordinates(dat_sp) <- ~Longitude+Latitude

# Param GWR model -> need to know which bandwith we can chose to have better results
gwr_bw <- gwr.sel(turf_filamentous ~ Anthropogenic_land_use + Population_density + 
                    Tourist_visitors + is_notake + is_restricted, 
                  data = dat_sp, adapt = TRUE)

# Run GWR model
gwr_model <- gwr(turf_filamentous ~ Anthropogenic_land_use + Population_density + 
                   Tourist_visitors + is_notake + is_restricted, 
                 data = dat_sp, adapt = gwr_bw, se.fit = TRUE, hatmatrix = TRUE)

# Res model
gwr_results <- as.data.frame(gwr_model$SDF)
gwr_results$Longitude <- dat$Longitude
gwr_results$Latitude <- dat$Latitude

# Map
map_gwr_notake <- ggplot(gwr_results, aes(x = Longitude, y = Latitude, color = Population_density)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_viridis_c(option = "magma") +
  theme_minimal() +
  labs(
    title = "GWR: Local Effect of 'Population_density' on Turf Algae",
    subtitle = "Ligther values mean the Population_density strongly increases turf algae in that area",
    color = "Population_density\nCoefficient",
    x = "Longitude", 
    y = "Latitude"
  )

print(map_gwr_notake)

# Statistical significance check
# Local t-statistic 
gwr_results$t_pop <- gwr_results$Population_density / gwr_results$Population_density_se

gwr_results$Significant_pop <- ifelse(abs(gwr_results$t_pop) > 1.96, 
                                      "Significant", "Not Significant")

# MAP: Highlight only the statistically significant areas
map_gwr_sig_pop <- ggplot(gwr_results, aes(x = Longitude, y = Latitude, 
                                           color = Significant_pop, 
                                           size = abs(Population_density))) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Not Significant" = "grey80", 
                                "Significant" = "darkorange")) +
  theme_minimal() +
  labs(
    title = "GWR Significance: Where does population density actually matter?",
    subtitle = "Orange points show where the effect of population on turf algae is statistically significant (p < 0.05)",
    color = "Statistical\nSignificance",
    size = "Strength of\nEffect",
    x = "Longitude", 
    y = "Latitude"
  )
print(map_gwr_sig_pop)

# 7. Distance-Decay Scatterplot

# Calculate real geographic distances (in km)
# spDists with longlat = TRUE accounts for the Earth's curvature
coords_matrix <- as.matrix(cbind(dat$Longitude, dat$Latitude))
geo_dist_matrix <- spDists(coords_matrix, longlat = TRUE) 

# Calculate ecological similarity
# Using Bray-Curtis similarity (1 - Bray-Curtis distance)
bray_dist <- vegdist(species_data, method = "bray")
eco_sim_matrix <- 1 - as.matrix(bray_dist)

# Only lower triangle of matrices to not take into account pairs and self-comparisons
geo_dist_vec <- geo_dist_matrix[lower.tri(geo_dist_matrix)]
eco_sim_vec <- eco_sim_matrix[lower.tri(eco_sim_matrix)]

# Dataframe for plotting
decay_df <- data.frame(
  Geographic_Distance_km = geo_dist_vec,
  Similarity = eco_sim_vec
)

# Distance decay plot
plot_decay <- ggplot(decay_df, aes(x = Geographic_Distance_km, y = Similarity)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  # Linear trend line (you can change method to "loess" for non-linear)
  geom_smooth(method = "lm", color = "firebrick", se = TRUE) + 
  theme_minimal() +
  labs(
    title = "Distance Decay of Benthic Communities",
    subtitle = "Community similarity as a function of geographic distance",
    x = "Geographic Distance (km)",
    y = "Community Similarity (1 - Bray-Curtis)"
  )

print(plot_decay)

# 8. Overall Variability of Environmental Parameters

# Only the numeric columns
env_num <- env_data %>% select(where(is.numeric))

# Min-Max scaling
min_max_scale <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

env_minmax <- as.data.frame(lapply(env_num, min_max_scale))

# Boxplot
dev.off() # Clear any overlapping ggplot canvas

par(mar = c(5, 12, 4, 2)) 

boxplot(env_minmax, 
        horizontal = TRUE, 
        las = 1,                 
        col = "darkseagreen", 
        main = "Variability of Environmental Variables (Min-Max Scaled)",
        xlab = "Relative Spread (0 = Min, 1 = Max)",
        cex.axis = 0.7) # Shrinks the axis text by 30% so they all fit
# CLUSTERING (Hellinger + K-Means): Groups sites into distinct ecological 
# communities based on the similarity of their species composition.
#
# GLOBAL SPATIAL STRUCTURE (Moran's I & nMDS): Tests and visualizes whether 
# these biological communities physically cluster together on the map.
#
# LOCAL SPATIAL STRUCTURE (LISA): Identifies the exact, localized geographic patches 
# where specific species (like turf algae) form statistically significant clusters.
#
# LOCAL SPATIAL STRUCTURE DRIVERS (GWR): Local distance-weighted regressions. This maps exactly where specific 
# anthropogenic factors (like population density) drive those algae clusters.
#
# DISTANCE DECAY SCATTERPLOT: Shows how distance reduce the similarity.
# Not super strong but still it is visible. One thing that is important to see
# is the variability at each distance, it is quite high. 
#
# VARIABILTY OF ENV PARAM: see which variable has high variability to see which
# of the them could explain the clustering 
