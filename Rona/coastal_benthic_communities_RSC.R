# ENV 513, Multivariate Statistics

# summary: benthic species at different 

# reading in data
setwd("C:/Users/rona.schenk/OneDrive - HESSO/7_Courses/2026_F_Statistics/Group Project/costal_benthic_communities_analysis/Rona")
benthic <- read.csv("Linetal_dataset_Benthic.csv")
cor <- read.csv("Linetal_dataset_Cor.csv")
env <- read.csv("Linetal_dataset_Env.csv")

# exploring -> sum up per column, is it one? 100 %?
sum(benthic$turf_filamentous)

# descriptive stats
summary(benthic)


# subset bc only a few species are present
benthic_pr <- benthic[,c(6,11,26)] # 6,11,26 are columns I want

# crustose and encrusting are left-skewed -> power- transformation 

                     