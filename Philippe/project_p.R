library(vegan)
library(gplots)
library(gclus)

getwd()
setwd("C:/Users/phili/OneDrive/Desktop/EPFL/ma3/multivar_r/project/costal_benthic_communities_analysis/data")
benthic_data <- read.csv("Linetal_dataset_Benthic.csv", head=TRUE, sep=",", row.names = 1)
cor_data <- read.csv("Linetal_dataset_Cor.csv", head=TRUE, sep=",")
env_data <- read.csv("Linetal_dataset_Env.csv", head=TRUE, sep=",")

head(benthic_data)
head(cor_data)
head(env_data)

summary(benthic_data)  # Not normally distributed since score of presence (right skewed since lot of 0 and some values for some species)
summary(cor_data)  
summary(env_data)  # More or less normally distributed

#Histograms
for (col_name in colnames(env_data)) {
  hist(env_data[[col_name]], main = paste("Histogram of", col_name), xlab = col_name)
}

#Histograms, correlation and linear plots
# Set plot size in Jupyter Notebook (20 x 14 inches)

panel.cor2 <- function(x, y){
  usr <- par("usr")            # Save current plotting parameters (user coordinates)
  on.exit(par(usr = usr))            # Restore plotting parameters when function exits
  par(usr = c(0, 1, 0, 1))    # Set user coordinates to [0,1] for x and y axes (normalized plot region)
  
  r <- round(cor(x, y), digits=2)        # Calculate correlation between x and y, rounded to 2 decimals
  txt <- paste0("R = ", r)                # Create text string to display correlation coefficient
  
  text(0.5, 0.5, txt, cex = 2)           # Add the correlation text at center (0.5, 0.5) with font size 2
}

panel.hist <- function(x, ...) {
  usr <- par("usr")             # Save current plotting parameters (user coordinates)
  on.exit(par(usr = usr))             # Ensure to restore these parameters when function exits
  
  par(usr = c(usr[1:2], 0, 1.5)) # Set new y-axis limits from 0 to 1.5, keep x-axis limits unchanged
  
  h <- hist(x, plot = FALSE)    # Compute histogram data for x without plotting
  
  breaks <- h$breaks            # Extract histogram bin boundaries
  nB <- length(breaks)          # Get number of bin boundaries
  
  y <- h$counts                 # Extract counts for each bin
  y <- y / max(y)               # Normalize counts to max = 1 for scaling
  
  rect(breaks[-nB], 0,         # Draw rectangles (bars) for histogram:
       breaks[-1],             # from each bin start (except last) ...
       y,                     # ... up to the normalized count height ...
       col = "cyan", ...)      # ... fill color cyan, plus other graphical args
}

# Select only numerical variables
env_numeric <- env_data[sapply(env_data, is.numeric)]
options(repr.plot.width = 20, repr.plot.height = 14)
pairs(env_numeric[13:17],                     # Create a pairs plot for the first 5 columns of the data frame 'env'
      panel = panel.smooth,        # Use 'panel.smooth' function to draw scatter plots with smooth curves in the lower panels
      diag.panel = panel.hist,     # Use 'panel.hist' function to draw histograms on the diagonal panels
      upper.panel = panel.cor2,    # Use 'panel.cor2' function to display correlation coefficients in the upper panels
      cex.labels = 3,              # Set the size of the variable labels (axis labels) to 3 (larger text)
      main = "Bivariate Plots with Histograms and Smooth Curves"  # Set the main title of the plot
)

#Dissimilarity
#imaeg.real.fct
image.real <- function(mat) { 
  mat <- t(mat)[,nrow(mat):1]
  image(mat, axes = FALSE, col = hcl.colors(15, palette="viridis"))
  axis(1, at = seq(0, 1, length = nrow(mat)), labels = rownames(mat), las=2)
  axis(2, at = seq(0, 1, length = ncol(mat)), labels = colnames(mat))
  box() 
}

# Compute Bray-Curtis dissimilarity matrix from the 'spe' dataset
benthic_data.db <- vegdist(benthic_data)

# Set plot size in Jupyter Notebook (8 x 8 inches)
options(repr.plot.width = 8, repr.plot.height = 8)

# Display the dissimilarity matrix as a color image
image.real(as.matrix(benthic_data.db))
title("Bray-Curtis Dissimilarity Matrix Across Samples")
