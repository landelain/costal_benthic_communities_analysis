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
              index    = "all") # all results in plot, with 5 as best number. Also tested for "ch" (2), "gap" (2) and "silhouette" (5)

nb$Best.nc                 # best k per index
table(nb$Best.nc[1, ])     # votes per k

# k-means (?)
km <- kmeans(benthic_hel, centers = 5, nstart = 100, iter.max = 100)
community <- factor(km$cluster)
table(community)

# save clusters in df
community <- data.frame(
  transect  = rownames(benthic_hel),
  community = factor(km$cluster)
)


# Fazit Clustering: If Hellinger transformation is done as in the paper, 
# then I get the same optimal number of clusters (5 = best number for most 
# indices, including ch and sihouette that we've seen in the course, 
# but not for gap). Via k-means I get the clusters. Which transect belongs
# to which cluster can be seen in the community dataframe.



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
  
  