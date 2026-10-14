
library(ggplot2)

# species

df_spec <- read.csv("data/Linetal_dataset_Benthic.csv")
colnames(df_spec)

df_spec$sp_repent

h <- hist(log1p(df_spec$sp_repent), freq=FALSE)
print(h)


# environment

df_env <- read.csv("data/Linetal_dataset_Env.csv")
colnames(df_env)

tabledata <- table(df_env$Region)
barplot(tabledata, las = 2, cex.names = 0.7)

ggplot(df_env, aes(x=Longitude, y=Latitude)) + geom_point()
