
# ============================================================
# NEET-BREAST-CANCER
# 02 - TUMOR VS NORMAL COMPARISON
# ============================================================

# 1. Charger les packages
packages <- c("dplyr", "tidyr", "ggplot2")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(dplyr)
library(tidyr)
library(ggplot2)

# 2. Charger les données préparées
neet_long <- readRDS("data/neet_long.rds")

# 3. Résumé descriptif de l'expression
expression_summary <- neet_long %>%
  group_by(sample, group) %>%
  summarise(
    median_expression = median(expression, na.rm = TRUE),
    mean_expression = mean(expression, na.rm = TRUE),
    .groups = "drop"
  )

print(expression_summary, n = Inf)

# 4. Comparaison statistique : test de Wilcoxon
wilcox_results <- neet_long %>%
  group_by(sample) %>%
  summarise(
    p_value = wilcox.test(
      expression[group == "Tumor"],
      expression[group == "Normal"]
    )$p.value,
    .groups = "drop"
  )

# Correction pour les tests multiples
wilcox_results$p_adjusted <- p.adjust(
  wilcox_results$p_value,
  method = "BH"
)

print(wilcox_results)

# 5. Calcul des différences d'expression
effect_results <- expression_summary %>%
  select(sample, group, median_expression) %>%
  pivot_wider(
    names_from = group,
    values_from = median_expression
  ) %>%
  mutate(
    log2_difference = Tumor - Normal,
    fold_change = 2^log2_difference
  )

print(effect_results)

# 6. Sauvegarder les résultats
write.csv(
  expression_summary,
  "results/NEET_expression_summary.csv",
  row.names = FALSE
)

write.csv(
  wilcox_results,
  "results/NEET_wilcoxon_results.csv",
  row.names = FALSE
)

write.csv(
  effect_results,
  "results/NEET_expression_comparison.csv",
  row.names = FALSE
)

cat("Comparaison tumorale et normale terminée.\n")
