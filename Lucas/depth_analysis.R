library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(magrittr)
library(vegan)
library(gclus)



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



# --------------- build transect depth profile ---------------- 


# find the clusters depth profiles

# data formatting

getwd()
df_spec <- read.csv("costal_benthic_communities_analysis/data_clean/benthic_clean.csv")
df_env <- read.csv("costal_benthic_communities_analysis/data_clean/env_clean.csv")

parts <- strsplit(df_spec$Transect, "_")
df_spec$site     <- paste(sapply(parts, `[`, 1), sapply(parts, `[`, 2), sep = "_")
df_spec$depth    <- as.numeric(gsub("m", "", sapply(parts, `[`, 3)))
df_spec$transect <- sapply(parts, `[`, 4)

colnames(df_spec)

avg_transect <- df_spec %>%
  group_by(site, depth) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

colnames(avg_transect)


avg_hell <- avg_transect %>%
  mutate(across(-c(site, depth), as.numeric)) %>% 
  { bind_cols(select(., site, depth),
              as_tibble(decostand(select(., -c(site, depth)), method = "hellinger"))) }


head(avg_hell)
avg_hell$depth

# creation of the site depth profiles (hellinger distance)

transect5 <- filter(avg_hell, depth == 5)
transect10 <- filter(avg_hell, depth == 10)
transect15 <- filter(avg_hell, depth == 15)
transect20 <- filter(avg_hell, depth == 20)
transect40 <- filter(avg_hell, depth == 40)

length(transect5)
length(transect10)
length(transect15)
length(transect20)
length(transect40)


# dissimilarity of depths


image.real <- function(mat) { 
  mat <- t(mat)[,nrow(mat):1]
  image(mat, axes = FALSE, col = hcl.colors(15, palette="viridis"))
  axis(1, at = seq(0, 1, length = nrow(mat)), labels = rownames(mat), las=2)
  axis(2, at = seq(0, 1, length = ncol(mat)), labels = colnames(mat))
  box() 
}


# coldiss()
# Color plots of a dissimilarity matrix, without and with ordering
#
# License: GPL-2
# Author:  Francois Gillet
#          23 August 2012 - rev. 07 June 2016

"coldiss" <- function(D,
                      nc = 4,
                      byrank = TRUE,
                      diag = FALSE) {
  require(gclus)
  
  D <- as.dist(as.matrix(D))
  
  if (max(D) > 1)
    D <- D / max(D)
  
  if (byrank) {
    spe.color <- dmat.color(1 - D, cm.colors(nc))
  }
  else {
    spe.color <- dmat.color(1 - D, byrank = FALSE, cm.colors(nc))
  }
  
  spe.o <- order.single(1 - D)
  speo.color <- spe.color[spe.o, spe.o]
  
  op <- par(mfrow = c(1, 2), pty = "s")
  
  if (diag) {
    plotcolors(
      spe.color,
      rlabels = attributes(D)$Labels,
      main = "Dissimilarity Matrix",
      dlabels = attributes(D)$Labels
    )
    plotcolors(
      speo.color,
      rlabels = attributes(D)$Labels[spe.o],
      main = "Ordered Dissimilarity Matrix",
      dlabels = attributes(D)$Labels[spe.o]
    )
  }
  else {
    plotcolors(spe.color, rlabels = attributes(D)$Labels,
               main = "Dissimilarity Matrix")
    plotcolors(speo.color,
               rlabels = attributes(D)$Labels[spe.o],
               main = "Ordered Dissimilarity Matrix")
  }
  
  par(op)
}

# Usage:
# coldiss(D = dissimilarity.matrix, nc = 4, byrank = TRUE, diag = FALSE)
# If D is not a dissimilarity matrix (max(D) > 1), then D is divided by max(D)
# nc 							number of colours (classes)
# byrank = TRUE		equal-sized classes
# byrank = FALSE	equal-length intervals
# diag = TRUE			print object labels also on the diagonal

# Example:
# coldiss(spe.dj, nc = 9, byrank = FALSE, diag = TRUE)


test <- dist(transect5)
image.real(as.matrix(test))
coldiss(test, nc = 15, diag = TRUE)

colnames(df_spec)
rowSums(df_spec[, -c(1, 29, 30, 31)])
df_spec

colnames(df_env)


# ------------------------

# set up and overview of depths 
df_spec <- read.csv("costal_benthic_communities_analysis/data_clean/benthic_clean.csv")
df_env <- read.csv("costal_benthic_communities_analysis/data_clean/env_clean.csv")

parts <- strsplit(df_spec$Transect, "_")
df_spec$site     <- paste(sapply(parts, `[`, 1), sapply(parts, `[`, 2), sep = "_")
df_spec$depth    <- as.numeric(gsub("m", "", sapply(parts, `[`, 3)))
df_spec$transect <- sapply(parts, `[`, 4)


avg_transect <- df_spec %>%
  group_by(site, depth) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")


table(avg_transect$site, avg_transect$depth)

# how many depths per site
avg_transect %>% count(site, name = "n_depths") %>% arrange(n_depths) %>% print(n=65)





# total abundance vs depth
species_cols <- setdiff(colnames(avg_transect), c("site", "depth", "region"))
avg_transect$abund <- rowSums(avg_transect[, species_cols], na.rm = TRUE)

ggplot(avg_transect, aes(factor(depth), abund)) +
  geom_boxplot() + geom_jitter(width = 0.1, alpha = 0.5)

kruskal.test(abund ~ factor(depth), data = avg_transect)


#set up
avg_transect$region <- sub("_.*", "", avg_transect$site)
avg_transect$total_cover <- rowSums(avg_transect[, species_cols], na.rm = TRUE)
species_cols <- setdiff(colnames(avg_transect), c("site", "depth", "region"))

# effect of depth alone 
ggplot(avg_transect, aes(factor(depth), total_cover)) +
  geom_boxplot() + geom_jitter(width = 0.1, alpha = 0.5)
summary(lm(log1p(total_cover) ~ depth + region, data = avg_transect))

comm <- decostand(avg_transect[, species_cols], "hellinger")
d <- dist(comm)
adonis2(d ~ depth, data = avg_transect,
        permutations = how(nperm = 999, blocks = avg_transect$region))


# human activity

df_env$Anthropogenic_land_use

df_env$site <- paste(gsub(" ", "", df_env$Region),
                     gsub(" ", "", df_env$Site), sep = "_")

df_env$site <- sub("^EastTaiwan_",  "East_",  df_env$site)
df_env$site <- sub("^NorthTaiwan_", "North_", df_env$site)
df_env$site <- gsub("\\.", "", df_env$site)        # 82.5 -> 825
df_env$site <- sub("Guishandao_", "Guishandao", df_env$site)
df_env$site <- sub("Hujingyu_",   "Hujingyu",   df_env$site)
df_env$site <- sub("Tiezhenyu_",  "Tiezhenyu",  df_env$site)
df_env$site <- sub("Zhongao_",    "Zhongao",    df_env$site)

# case-insensitive match handles GuianHarbor vs Guianharbor
key <- function(x) tolower(x)
df_env$key       <- key(df_env$site)
avg_transect$key <- key(avg_transect$site)

# check if it worked
setdiff(avg_transect$key, df_env$key)  
setdiff(df_env$key, avg_transect$key)  



hum_vars <- c("Tourist_visitors", "Anthropogenic_land_use", "Population_density") 


env_site <- df_env %>%
  group_by(key) %>%
  summarise(across(all_of(hum_vars), \(x) mean(x, na.rm = TRUE)),
            mgmt = first(Management_status))

pca <- prcomp(env_site[, hum_vars], scale. = TRUE)
env_site$human <- pca$x[, 1]
if (mean(pca$rotation[, 1]) < 0) env_site$human <- -env_site$human 

meta <- left_join(avg_transect, env_site, by = "key")
sum(complete.cases(meta$human)) # check that we get 87
meta$region <- factor(meta$region)

# Does depth effect depend on humana activity ?
species_cols <- setdiff(colnames(avg_transect),
                        c("site", "depth", "region", "key", "total_cover"))

# --- diagnostics (replaces sum(ok), table(...[ok]), setdiff(...) ) ---
table(meta$region)                                  # sites x depths per region
setdiff(avg_transect$key, env_site$key)             # want character(0)
summary(meta$human)

# --- Does the depth effect depend on human activity? ---

# total cover
m1 <- lm(log1p(total_cover) ~ depth * human + region, data = meta)
summary(m1)

# composition
comm <- decostand(as.matrix(meta[, species_cols]), "hellinger")
adonis2(dist(comm) ~ depth * human, data = meta,
        by = "terms",
        permutations = how(nperm = 999, blocks = meta$region))
