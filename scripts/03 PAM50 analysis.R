
# ============================================================
# NEET-BREAST-CANCER
# 03 - PAM50 SUBTYPE ANALYSIS
# ============================================================

# 1. Charger les packages
packages <- c("dplyr", "tidyr", "FSA")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(dplyr)
library(tidyr)
library(FSA)

# 2. Charger les données préparées
neet_long <- readRDS("data/neet_long.rds")
clinical_tumor <- readRDS("data/clinical_tumor.rds")

# 3. Garder uniquement les échantillons tumoraux
neet_tumor <- neet_long %>%
  filter(group == "Tumor")

# 4. Associer les données d'expression aux sous-types PAM50
neet_subtype <- neet_tumor %>%
  left_join(
    clinical_tumor %>%
      select(sampleID, PAM50Call_RNAseq),
    by = c("patient" = "sampleID")
  )

# 5. Retirer les annotations manquantes
neet_subtype_clean <- neet_subtype %>%
  filter(!is.na(PAM50Call_RNAseq)) %>%
  mutate(
    PAM50Call_RNAseq = factor(
      PAM50Call_RNAseq,
      levels = c("LumA", "LumB", "Her2", "Basal", "Normal")
    )
  ) %>%
  filter(!is.na(PAM50Call_RNAseq))

print(table(neet_subtype_clean$PAM50Call_RNAseq))

# 6. Résumé de l'expression par sous-type
subtype_summary <- neet_subtype_clean %>%
  group_by(sample, PAM50Call_RNAseq) %>%
  summarise(
    median_expression = median(expression, na.rm = TRUE),
    mean_expression = mean(expression, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

print(subtype_summary, n = Inf)

# 7. Test de Kruskal-Wallis pour chaque gène
kruskal_results <- neet_subtype_clean %>%
  group_by(sample) %>%
  summarise(
    p_value = kruskal.test(
      expression ~ PAM50Call_RNAseq
    )$p.value,
    .groups = "drop"
  )

kruskal_results$p_adjusted <- p.adjust(
  kruskal_results$p_value,
  method = "BH"
)

print(kruskal_results)

# 8. Comparaisons deux à deux avec le test de Dunn
dunn_results <- neet_subtype_clean %>%
  group_by(sample) %>%
  group_modify(~ {
    dunnTest(
      expression ~ PAM50Call_RNAseq,
      data = .x,
      method = "bh"
    )$res
  }) %>%
  ungroup()

# Correction globale des p-values des comparaisons de Dunn
dunn_results$P.adj_global <- p.adjust(
  dunn_results$P.unadj,
  method = "BH"
)

print(dunn_results, n = Inf)

# 9. Calcul de la taille d'effet epsilon carré
effect_size <- neet_subtype_clean %>%
  group_by(sample) %>%
  summarise(
    H = as.numeric(
      kruskal.test(
        expression ~ PAM50Call_RNAseq
      )$statistic
    ),
    k = n_distinct(PAM50Call_RNAseq),
    n = n(),
    epsilon_squared = (H - k + 1) / (n - k),
    .groups = "drop"
  )

print(effect_size)

# 10. Sauvegarder les résultats
write.csv(
  subtype_summary,
  "results/NEET_PAM50_summary.csv",
  row.names = FALSE
)

write.csv(
  kruskal_results,
  "results/NEET_Kruskal_Wallis.csv",
  row.names = FALSE
)

write.csv(
  dunn_results,
  "results/NEET_Dunn_tests.csv",
  row.names = FALSE
)

write.csv(
  effect_size,
  "results/NEET_effect_sizes.csv",
  row.names = FALSE
)

cat("Analyse PAM50 terminée.\n")
