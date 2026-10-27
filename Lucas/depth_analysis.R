library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

getwd()

df_spec <- read.csv("data_clean/benthic_clean.csv")
colnames(df_spec)

df_spec$depth <- gsub("m", "", sapply(strsplit(as.character(df_spec$Transect), "_"), `[`, 3))
depth_list <- split(df_spec, df_spec$depth)
species_cols <- setdiff(names(df_spec), c("Transect", "depth"))


# species histogram for each depth

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


# depth histograms for each species

species_cols <- setdiff(names(df_spec), c("Transect", "depth"))

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



# site names parsing 

parts <- strsplit(df_spec$Transect, "_")
df_spec$site     <- paste(sapply(parts, `[`, 1), sapply(parts, `[`, 2), sep = "_")
df_spec$depth    <- as.numeric(gsub("m", "", sapply(parts, `[`, 3)))
df_spec$transect <- sapply(parts, `[`, 4)

species_cols <- setdiff(names(df_spec), c("Transect", "site", "depth", "transect"))

long <- df_spec %>%
  pivot_longer(all_of(species_cols), names_to = "species", values_to = "cover") %>%
  mutate(present = as.integer(cover > 0))

df_spec %>% count(site, depth) %>% pivot_wider(names_from = depth, values_from = n)


# model fittings 

fit_one <- function(d) {
  # need some variation to fit anything
  if (sum(d$present) < 3 || sum(d$present) == nrow(d)) return(NULL)

  pres <- glm(present ~ depth, family = binomial, data = d)
  abun <- lm(log1p(cover) ~ depth, data = d)   # abundance (log-transformed)

  tibble(
    pres_slope = coef(pres)["depth"],
    pres_p     = summary(pres)$coefficients["depth", 4],
    abun_slope = coef(abun)["depth"],
    abun_p     = summary(abun)$coefficients["depth", 4]
  )
}

site_models <- long %>%
  group_by(site, species) %>%
  group_modify(~ { r <- tryCatch(fit_one(.x), error = function(e) NULL)
                   if (is.null(r)) tibble() else r }) %>%
  ungroup()


site_models

# model fitting plots

specie_of_interest <- "hc_massive"

p <- ggplot(filter(long, species == specie_of_interest),
       aes(depth, log1p(cover))) +
  geom_jitter(width = 0.3, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ site)

ggsave(paste0("Lucas/site_models_", specie_of_interest, ".png"), p, width = 14, height = 10)
