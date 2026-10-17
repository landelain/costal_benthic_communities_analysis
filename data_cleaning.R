# ============================================================
# 00_prep.R : nettoyage des données du groupe
# Lin et al. 2024, communautés benthiques de Taïwan
#
# Entrée  : les 3 CSV bruts du paper (Benthic, Env, Cor)
# Sortie  : 3 CSV propres dans le dossier data_clean/
#   benthic_clean.csv      433 transects x 27 groupes benthiques
#   env_clean.csv          433 transects x variables environnementales
#   communities_site.csv   87 sites x part de chaque communauté du paper
# Les lignes de benthic_clean et env_clean sont dans le MEME ordre.
# ============================================================

# Seule ligne à adapter : le dossier qui contient les CSV bruts
data_dir <- "data"
out_dir  <- file.path(dirname(data_dir), "data_clean")
dir.create(out_dir, showWarnings = FALSE)


# 1. Charger ####################################################
benthic <- read.csv(file.path(data_dir, "Linetal_dataset_Benthic.csv"), row.names = 1)
env     <- read.csv(file.path(data_dir, "Linetal_dataset_Env.csv"))
cor     <- read.csv(file.path(data_dir, "Linetal_dataset_Cor.csv"), check.names = FALSE)


# 2. Vérifier ###################################################
stopifnot(nrow(benthic) == 433, nrow(env) == 433)            # 433 transects des deux côtés
stopifnot(sum(is.na(benthic)) == 0, sum(is.na(env)) == 0)    # pas de valeurs manquantes
stopifnot(min(benthic) >= 0)                                  # pas de recouvrement négatif
stopifnot(!anyDuplicated(rownames(benthic)))                  # pas de transect en double


# 3. Nettoyer les textes ########################################
env$Region            <- trimws(env$Region)
env$Site              <- trimws(env$Site)
env$Management_status <- trimws(env$Management_status)

# Un site est écrit "Gongguan" dans Benthic mais "Gonggang" dans Env et Cor :
# on garde l'orthographe de Env et Cor.
rownames(benthic) <- sub("Gongguan", "Gonggang", rownames(benthic))

# Erreur dans Env : à Dongsha, Southwest Outer Atoll, Env a 4 transects à 10 m et 5 à 20 m,
# alors que Benthic (et le fichier Cor) en ont 5 à 10 m et 4 à 20 m.
# Le transect en trop à 20 m dans Env est le 5e. Son benthos (Benthic, 10m_T5) ressemble
# nettement aux transects à 10 m et pas à ceux à 20 m : c'est un transect à 10 m mal étiqueté.
# On corrige sa profondeur et sa lumière (seules variables qui changent avec la profondeur).
sw20 <- which(env$Site == "Southwest Outer Atoll" & env$Depth == 20)
sw10 <- which(env$Site == "Southwest Outer Atoll" & env$Depth == 10)
if (length(sw20) == 5) {
  env$Depth[sw20[5]]           <- 10
  env$Light_intensity[sw20[5]] <- env$Light_intensity[sw10[1]]
  cat("Correction : 1 transect de Southwest Outer Atoll repassé de 20 m à 10 m\n")
}


# 4. Remettre benthic dans l'ordre de env #######################
# Benthic est trié par ordre alphabétique, Env ne l'est pas : la ligne 6
# n'est pas le même transect dans les deux fichiers. On les relie grâce
# au nom du transect, de la forme Region_Site_10m_T1.
parts   <- strsplit(rownames(benthic), "_")
site_b  <- sapply(parts, function(p) p[2])                                  # "EastOuterAtoll"
depth_b <- as.numeric(sub("m$", "", sapply(parts, function(p) p[3])))       # 10
num_b   <- as.numeric(sub("^T", "", sapply(parts, function(p) p[4])))       # 1 à 5

# Même clé "site profondeur" des deux côtés : sans espaces, tirets ni majuscules
# ("Guian harbor" et "GuianHarbor" deviennent "guianharbor", "82.5" et "825" deviennent "825")
cle <- function(site, depth) paste(tolower(gsub("[^A-Za-z0-9]", "", site)), depth)
key_b <- cle(site_b, depth_b)
key_e <- cle(env$Site, env$Depth)

# Chaque site doit exister des deux côtés, avec le même nombre de transects
manquants <- c(setdiff(key_b, key_e), setdiff(key_e, key_b))
if (length(manquants) > 0) stop("Sites présents dans un seul fichier : ", paste(manquants, collapse = ", "))
stopifnot(all(table(key_b)[names(table(key_e))] == table(key_e)))

# Dans un site, les lignes de env sont dans l'ordre T1, T2, ... (leur latitude augmente à chaque transect)
rang_b <- ave(num_b, key_b, FUN = rank)
rang_e <- ave(seq_along(key_e), key_e, FUN = seq_along)
ordre  <- match(paste(key_e, rang_e), paste(key_b, rang_b))
stopifnot(!anyNA(ordre), !anyDuplicated(ordre))
benthic <- benthic[ordre, ]

# Exception : Unstable_substrate_cover (seule variable mesurée transect par transect)
# est déjà rangée dans l'ordre de Benthic dans le fichier brut. Preuve : avec l'ordre brut,
# corrélation de Spearman de -0.77 avec le recouvrement vivant (logique : plus de sable et
# de débris = moins de vivant), contre -0.29 si on la déplace avec le reste de env.
# On la déplace donc comme benthic, pour qu'elle reste avec son transect.
env$Unstable_substrate_cover <- env$Unstable_substrate_cover[ordre]

# Contrôle : la profondeur écrite dans le nom correspond maintenant à Depth partout
stopifnot(all(depth_b[ordre] == env$Depth))
cat(sum(ordre != 1:433), "lignes de benthic ont été remises dans l'ordre de env\n")


# 5. Identifiant commun ##########################################
env$Transect <- rownames(benthic)
env     <- env[, c("Transect", setdiff(names(env), "Transect"))]
benthic <- cbind(Transect = rownames(benthic), benthic)


# 6. Fichier Cor : communautés du paper par site ################
# Colonnes 1 à 5 = part des transects du site dans chaque communauté du paper
# (vérifié avec les effectifs du paper : 100, 128, 155, 27 et 23 transects).
# y et x = position du site et de la profondeur sur la Figure 1b ; Label1 et Label2 = étiquettes de la figure.
names(cor)[names(cor) %in% as.character(1:5)] <-
  c("CCA", "Turf_StonyCoral", "Turf", "Digitate_Octocoral", "Bushy_Octocoral")
cor <- cor[, c("Region", "Site", "Depth", "Latitude", "Longitude",
               "CCA", "Turf_StonyCoral", "Turf", "Digitate_Octocoral", "Bushy_Octocoral")]
cor$Region <- trimws(cor$Region)
cor$Site   <- trimws(cor$Site)
# Communauté majoritaire du site
cor$Main_community <- names(cor)[6:10][apply(cor[, 6:10], 1, which.max)]
stopifnot(nrow(cor) == 87, all(abs(rowSums(cor[, 6:10]) - 1) < 1e-6))


# 7. Exporter ###################################################
write.csv(benthic, file.path(out_dir, "benthic_clean.csv"),    row.names = FALSE)
write.csv(env,     file.path(out_dir, "env_clean.csv"),        row.names = FALSE)
write.csv(cor,     file.path(out_dir, "communities_site.csv"), row.names = FALSE)
cat("OK : fichiers écrits dans", out_dir, "\n")


