# ============================================================
# Alzheimer snRNA-seq Reanalysis
# Sex-linked gene expression validation
#
# This script should be run after:
# R/01_preprocessing_harmony.R
# ============================================================


# ------------------------------------------------------------
# 1. Required libraries
# ------------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)


# ------------------------------------------------------------
# 2. Check that the combined Seurat object exists
# ------------------------------------------------------------

if (!exists("combined")) {
  stop(
    "The 'combined' Seurat object was not found. ",
    "Run R/01_preprocessing_harmony.R first."
  )
}


# ------------------------------------------------------------
# 3. Sex-linked marker genes
# ------------------------------------------------------------

DefaultAssay(combined) <- "RNA"

sex_markers <- c(
  "XIST",
  "RPS4Y1",
  "KDM5D",
  "UTY",
  "EIF1AY",
  "DDX3Y",
  "ZFY"
)


# Keep only markers present in the dataset
available_sex_markers <- intersect(
  sex_markers,
  rownames(combined)
)

print(available_sex_markers)


# ------------------------------------------------------------
# 4. Visualize sex-linked gene expression by sample
# ------------------------------------------------------------

sex_dotplot <- DotPlot(
  combined,
  features = available_sex_markers,
  group.by = "sample_id",
  assay = "RNA"
) +
  RotatedAxis() +
  ggtitle("Sex-linked Gene Expression by Sample")

print(sex_dotplot)


VlnPlot(
  combined,
  features = available_sex_markers,
  group.by = "sample_id",
  pt.size = 0,
  ncol = 2
)


# ------------------------------------------------------------
# 5. Average expression by sample
# ------------------------------------------------------------

sex_average <- AverageExpression(
  combined,
  assays = "RNA",
  features = available_sex_markers,
  group.by = "sample_id"
)

print(sex_average$RNA)


# ------------------------------------------------------------
# 6. Add expression-inferred sex to metadata
# ------------------------------------------------------------

# Based on the observed expression pattern:
# A1 and A3: high XIST and negligible Y-linked expression
# A2 and A4: Y-linked marker expression and negligible XIST

combined$inferred_sex <- dplyr::case_when(
  combined$sample_id %in% c("A1", "A3") ~ "Female-like",
  combined$sample_id %in% c("A2", "A4") ~ "Male-like",
  TRUE ~ NA_character_
)


# ------------------------------------------------------------
# 7. Inspect disease-sex relationship
# ------------------------------------------------------------

table(
  combined$sample_id,
  combined$condition,
  combined$inferred_sex
)


sample_metadata_check <- combined@meta.data %>%
  select(
    sample_id,
    condition,
    inferred_sex
  ) %>%
  distinct() %>%
  arrange(sample_id)

print(sample_metadata_check)


# ------------------------------------------------------------
# Important interpretation
# ------------------------------------------------------------

message(
  "A1 and A3 (AD) show a female-like expression pattern, ",
  "whereas A2 and A4 (Control) show a male-like expression pattern. ",
  "Disease status and expression-inferred sex are therefore fully confounded ",
  "in this four-donor dataset."
)
