# ENV 513, Multivariate Statistics

# summary: benthic species at different transects, plus environmental and spatial data

# 0) preparation 

# packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(vegan)
library(ade4)
library(pastecs) 
library(psych)
library(fmsb)
library(NbClust)


# reading in data
setwd("C:/Users/rona.schenk/OneDrive - HESSO/7_Courses/2026_F_Statistics/Group Project/costal_benthic_communities_analysis/Rona")
benthic <- read.csv("benthic_clean.csv")
com <- read.csv("communities_site.csv")
env <- read.csv("env_clean.csv")


# 1) exploration

# objects = benthic species, descriptors = environmental data

# descriptive stats
summary(benthic)

# 2) clustering (compare with paper)

# take transect as row names
rownames(benthic) <- benthic$Transect   # use transects as rownames 
benthic$Transect <- NULL                # delete transect column

dim(benthic)                            # control: 433 x 27 (# species)?

# control before Hellinger
all(sapply(benthic, is.numeric))   # Numerical in all columns? Should be true. 
anyNA(benthic)                     # Any NA in the columns? Should be false. 
sum(rowSums(benthic) == 0)         # Is any row sum = 0? Problematic bc of division / 0 in Hellinger transformation. Should be 0. 

# Hellinger-transformation
benthic_hel <- decostand(benthic, "hel")

# control after Hellinger
dim(benthic_hel)                   # Still 433 x 27? Nothing lost?
summary(rowSums(benthic_hel^2))    # Should be 1's. Hellinger = Square root of relative values -> square should equal 1.

# Euclidean distance (prep for NbClust)
benthic_dist <- dist(benthic_hel, method = "euclidean")

# Clustering with NbClust, for all indices
nb <- NbClust(benthic_hel,
              diss     = benthic_dist,
              distance = NULL,
              min.nc   = 2,
              max.nc   = 10,
              method   = "kmeans",
              index    = "all") # all results in plot, with 5 as best number. Also tested for "ch" (), "gap" (2) and "silhouette" (5)

nb$Best.nc                 # best k per index
table(nb$Best.nc[1, ])     # votes per k

# k-means (benthic_hel or benthic_dist?)
km <- kmeans(benthic_hel, centers = 5, nstart = 100, iter.max = 100)
community <- factor(km$cluster)
table(community)

km.2 <- kmeans(benthic_dist, centers = 5, nstart = 100, iter.max = 100)
community.2 <- factor(km.2$cluster)
table(community.2)


# save clusters in df (for both versions of k-means clustering)
community <- data.frame(
  transect  = rownames(benthic_hel),
  community = factor(km$cluster)
)

community.2 <- data.frame(
  transect  = rownames(benthic_dist),
  community.2 = factor(km.2$cluster.2)
)

# plot k-means clustering
env_comm <- merge(env, community, by.x = "Transect", by.y = "transect")

write.csv(env_comm, "env_comm_k5.csv", row.names = FALSE)

summary(env_comm[, c("Longitude", "Latitude")])

cols <- c("red", "blue", "darkgreen", "orange", "purple")

par(mar = c(5, 4, 4, 7), xpd = TRUE)   # rechter Rand breiter
plot(env_comm$Longitude, env_comm$Latitude,
     col = cols[env_comm$community], pch = 19,
     xlab = "Longitude", ylab = "Latitude", main = "Communities (k = 5)")

set.seed(1)
plot(jitter(env_comm$Longitude, amount = 0.03),
     jitter(env_comm$Latitude,  amount = 0.03),
     col  = adjustcolor(cols[env_comm$community], alpha.f = 0.6),
     pch  = 19, cex = 0.8,
     xlab = "Longitude", ylab = "Latitude",
     main = "Communities (k = 5), mit Jitter")

legend("topright", inset = c(-0.2, 0), legend = 1:5, col = cols,
       pch = 19, title = "Community", bty = "n")
par(mar = c(5, 4, 4, 2) + 0.1, xpd = FALSE)  # zurücksetzen

table(env_comm$Region, env_comm$community) # schauen, welche Community in welche Region gehört


# Fazit Clustering: If Hellinger transformation is done as in the paper (justified by many 0's), 
# then I get the same optimal number of clusters (5 = best number for most 
# indices, including ch and sihouette that we've seen in the course, 
# but not for gap). Via k-means I get the clusters. Which transect belongs
# to which cluster can be seen in the community dataframe.
# Question: k-means on distance or or Hellinger-matrix?

write.csv(community, "community.csv", row.names = FALSE)



# X) old



# subset bc only a few species are present
benthic_pr <- benthic[,c(1, 6,11,26)] # 6,11,26 are columns I want

# separate the transect information
benthic_pr <- benthic_pr %>%
  separate(
    transect,
    into = c("transect_1", "transect_2", "transect_3", "transect_4"),
    sep = "_",
    remove = FALSE
  )

# datatype
str(benthic_pr)

# crustose and encrusting are left-skewed -> power- transformation 

# -> sum up per column, is it one? 100 %?
sum(benthic$turf_filamentous)


# 2) visualization 

# histogram of one species

  # Set plot size in Jupyter Notebook (10 x 10 inches)
  options(repr.plot.width = 10, repr.plot.height = 10)
  
  # Plot a basic histogram of the variable 'Pato' from the 'spe' data frame
  # This displays the frequency (number of observations) in each bin
  hist(benthic$sp_encrusting)
  
  
  # Plot the same histogram but with probability = TRUE
  # This scales the y-axis to show density (so that the total area under the histogram = 1)
  hist(benthic$sp_encrusting,probability = TRUE)           
  
  
# heatmap
  # Set plot size in Jupyter Notebook (10 x 10 inches)
  options(repr.plot.width = 10, repr.plot.height = 10)
  
  # Plot heatmap for the spe dataset
  heatmap(as.matrix(benthic[,-1]),         # Convert spe to matrix
            dendrogram = "none",    # No dendrogram
            trace = "none",         # No trace lines
            col = heat.colors(20),  # Color gradient
            main = "Species Abundance Heatmap")  # Title
  
  
# radarplot
  # Set plot size in Jupyter Notebook (10 x 10 inches)
  options(repr.plot.width = 10, repr.plot.height = 10)
  
  # Radar chart with customized colors and styles
  radarchart(
    benthic,     # data frame where the first two rows contain the max and min values
    axistype = 1, # changes the appearance of the concentric circles
    pcol = c("red", "blue"),               # Line colors for each observation
    pfcol = c(rgb(1,0,0,0.3), rgb(0,0,1,0.3)), # Transparent fill colors
    plwd = 1,                             # Line width
    cglcol = "grey",                      # Color of concentric grid lines
    cglty = 1,                           # Grid line type (1 = solid)
    axislabcol = "darkgrey",             # Axis label color
    vlcex = 1.2                          # Variable label size
  )
  
  