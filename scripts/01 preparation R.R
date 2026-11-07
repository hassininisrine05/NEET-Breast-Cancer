
# ============================================================
# NEET-BREAST-CANCER
# 01 - DATA PREPARATION
# ============================================================

# 1. Informations sur le projet
project_name <- "Analyse intégrative des protéines NEET"
genes <- c("CISD1", "CISD2", "CISD3")

print(project_name)
print(genes)

# 2. Préparer les dossiers de résultats
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# 3. Installer et charger les packages nécessaires
packages <- c("UCSCXenaTools", "tidyr")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(UCSCXenaTools)
library(tidyr)

# 4. Explorer les données UCSC Xena
data("XenaData")

brca_data_info <- XenaData[
  XenaData$XenaHostNames == "tcgaHub" &
    grepl("TCGA Breast Cancer",
          XenaData$XenaCohorts,
          ignore.case = TRUE) &
    grepl("gene expression",
          XenaData$DataSubtype,
          ignore.case = TRUE),
]

print(
  brca_data_info[, c(
    "XenaDatasets",
    "DataSubtype",
    "SampleCount",
    "Unit"
  )],
  n = Inf
)

# 5. Sélectionner les données d'expression TCGA-BRCA
brca_expression <- XenaGenerate(
  subset = XenaHostNames == "tcgaHub" &
    XenaDatasets == "TCGA.BRCA.sampleMap/HiSeqV2"
)

brca_query <- XenaQuery(brca_expression)
brca_download <- XenaDownload(brca_query)
brca_data <- XenaPrepare(brca_download)

# Vérifier les dimensions
dim(brca_data)

# 6. Extraire les trois gènes NEET
neet_data <- brca_data[
  brca_data$sample %in% genes,
]

print(neet_data)

# 7. Identifier les échantillons tumoraux et normaux
sample_ids <- colnames(brca_data)[-1]
sample_types <- substr(sample_ids, 14, 15)

print(table(sample_types))

tumor_samples <- sample_ids[sample_types == "01"]
normal_samples <- sample_ids[sample_types == "11"]

# 8. Construire le tableau d'expression
neet_comparison <- neet_data[
  ,
  c("sample", tumor_samples, normal_samples)
]

# 9. Transformer les données au format long
neet_long <- pivot_longer(
  neet_comparison,
  cols = -sample,
  names_to = "patient",
  values_to = "expression"
)

# Identifier le groupe de chaque échantillon
neet_long$group <- ifelse(
  substr(neet_long$patient, 14, 15) == "01",
  "Tumor",
  "Normal"
)

print(table(neet_long$group))

# 10. Télécharger les annotations cliniques PAM50
clinical_query <- XenaGenerate(
  subset = XenaHostNames == "tcgaHub" &
    XenaDatasets ==
    "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"
)

clinical_query <- XenaQuery(clinical_query)
clinical_download <- XenaDownload(clinical_query)
brca_clinical <- XenaPrepare(clinical_download)

# Vérifier les annotations
print(table(
  brca_clinical$PAM50Call_RNAseq,
  useNA = "ifany"
))

# 11. Garder les annotations des tumeurs primaires
clinical_tumor <- brca_clinical[
  brca_clinical$sampleID %in% tumor_samples,
]

# Vérifier les sous-types
print(table(
  clinical_tumor$PAM50Call_RNAseq,
  useNA = "ifany"
))

# 12. Sauvegarder les données préparées
saveRDS(neet_long, "data/neet_long.rds")
saveRDS(clinical_tumor, "data/clinical_tumor.rds")

write.csv(
  neet_long,
  "data/NEET_expression_long.csv",
  row.names = FALSE
)

# 13. Message de fin
cat("Préparation des données terminée.\n")

