# ENV 513, Multivariate Statistics

# summary: benthic species at different transects, plus environmental and spatial data

# 0) preparation 

# packages
library(tidyr)
library(dplyr)
library(ggplot2)

# reading in data
setwd("C:/Users/rona.schenk/OneDrive - HESSO/7_Courses/2026_F_Statistics/Group Project/costal_benthic_communities_analysis/Rona")
benthic <- read.csv("Linetal_dataset_Benthic.csv")
cor <- read.csv("Linetal_dataset_Cor.csv")
env <- read.csv("Linetal_dataset_Env.csv")


# 1) exploration 

# objects = benthic species, descriptors = environmental data

# descriptive stats
summary(benthic)

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

benthic_pr_long <- benthic_pr %>%
  pivot_longer(
    cols = c(cca_crustose, hc_encrusting, turf_filamentous),
    names_to = "species",
    values_to = "abundance"
  ) %>%
  pivot_longer(
    cols = starts_with("transect_"),
    names_to = "transect_level",
    values_to = "transect_value"
  )
                     