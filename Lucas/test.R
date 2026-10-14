
library(ggplot2)
library(dplyr)
library(tidyr)

# species

df_spec <- read.csv("data/Linetal_dataset_Benthic.csv")
colnames(df_spec)

df_spec$sp_repent

# h <- hist(log1p(df_spec$sp_repent), freq=FALSE)
# print(h)


df_spec$depth <- gsub("m", "", sapply(strsplit(as.character(df_spec$transect), "_"), `[`, 3))
depth_list <- split(df_spec, df_spec$depth)
species_cols <- setdiff(names(df_spec), c("transect", "depth"))

plots <- list()
for (i in 1:5) {
  totals <- data.frame(
    species = species_cols,
    total   = log1p(colSums(depth_list[[i]][, species_cols], na.rm = TRUE)/length(depth_list[[i]]$depth))
  )

  plots[[i]] <- ggplot(totals, aes(x = reorder(species, -total), y = total)) +
    geom_col(fill = "steelblue") +
    labs(x = NULL, y = "LOG Sum", title = paste("Depth", names(depth_list)[i])) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  ggsave(paste0("Lucas/bar_plot_depth_", names(depth_list)[i], ".png"), plots[[i]], width = 8, height = 5)

}

# for (i in 1:5) {
#   totals <- data.frame(
#     species = factor(species_cols, levels = species_cols),
#     total   = log1p(colSums(depth_list[[i]][, species_cols], na.rm = TRUE))
#   )

#   plots[[i]] <- ggplot(totals, aes(x = species, y = total)) +
#   geom_col(fill = "steelblue") +
#   labs(x = NULL, y = "log(1 + sum)", title = paste("Depth", names(depth_list)[i])) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
#   ggsave(paste0("Lucas/bar_plot_depth_", names(depth_list)[i], ".png"), plots[[i]], width = 8, height = 5)
# }


species_cols <- setdiff(names(df_spec), c("transect", "depth"))

plot_df <- df_spec %>%
  pivot_longer(cols = all_of(species_cols),
               names_to = "species", values_to = "value") %>%
  group_by(species, depth) %>%
  summarise(mean_value = sum(value, na.rm = TRUE) / n(), .groups = "drop") %>%
  mutate(
    species = factor(species, levels = species_cols),
    depth   = factor(depth, levels = sort(unique(as.numeric(as.character(depth)))))
  )

p <- ggplot(plot_df, aes(x = depth, y = mean_value)) +
  geom_col(fill = "steelblue") +
  facet_wrap(~ species, ncol = 6, scales = "free_y", axes = "all_x") +
  labs(x = "Depth", y = "Sum / n") +
  theme_minimal() +
  theme(
    strip.text  = element_text(size = 7),
    axis.text.x = element_text(size = 6, angle = 0, hjust = 1, vjust = 0.5),
    axis.title.x = element_text(size = 9)
  )

print(p)
ggsave("Lucas/facet_species_by_depth.png", p, width = 14, height = 10)




# environment

df_env <- read.csv("data/Linetal_dataset_Env.csv")
colnames(df_env)

tabledata <- table(df_env$Region)
barplot(tabledata, las = 2, cex.names = 0.7)

ggplot(df_env, aes(x=Longitude, y=Latitude)) + geom_point()



