
# ============================================================
# NEET-BREAST-CANCER
# 04 - DATA VISUALIZATION
# ============================================================

# 1. Charger les packages
packages <- c("dplyr", "ggplot2")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(dplyr)
library(ggplot2)

# 2. Charger les données
neet_long <- readRDS("data/neet_long.rds")
clinical_tumor <- readRDS("data/clinical_tumor.rds")

# 3. Préparer les données tumorales avec les sous-types PAM50
neet_subtype_clean <- neet_long %>%
  filter(group == "Tumor") %>%
  left_join(
    clinical_tumor %>%
      select(sampleID, PAM50Call_RNAseq),
    by = c("patient" = "sampleID")
  ) %>%
  filter(!is.na(PAM50Call_RNAseq)) %>%
  mutate(
    PAM50Call_RNAseq = factor(
      PAM50Call_RNAseq,
      levels = c("LumA", "LumB", "Her2", "Basal", "Normal")
    )
  ) %>%
  filter(!is.na(PAM50Call_RNAseq))

# 4. Créer le dossier des figures
dir.create("figures", showWarnings = FALSE)

# 5. Figure 1 : comparaison Tumor vs Normal
figure1 <- ggplot(
  neet_long,
  aes(x = group, y = expression, fill = group)
) +
  geom_boxplot() +
  facet_wrap(~sample, scales = "free_y") +
  labs(
    title = "Expression of NEET genes in breast cancer",
    x = "Sample group",
    y = "Expression (log2 normalized counts + 1)"
  ) +
  theme_classic() +
  theme(legend.position = "none")

print(figure1)

ggsave(
  "figures/Figure1_Tumor_vs_Normal.png",
  plot = figure1,
  width = 9,
  height = 5,
  dpi = 300
)

# 6. Figure 2 : expression selon les sous-types PAM50
figure2 <- ggplot(
  neet_subtype_clean,
  aes(x = PAM50Call_RNAseq, y = expression,
      fill = PAM50Call_RNAseq)
) +
  geom_boxplot() +
  facet_wrap(~sample, scales = "free_y") +
  labs(
    title = "NEET gene expression across PAM50 subtypes",
    x = "PAM50 molecular subtype",
    y = "Expression (log2 normalized counts + 1)"
  ) +
  theme_classic() +
  theme(legend.position = "none")

print(figure2)

ggsave(
  "figures/Figure2_PAM50_subtypes.png",
  plot = figure2,
  width = 10,
  height = 5,
  dpi = 300
)

cat("Figures générées avec succès.\n")
