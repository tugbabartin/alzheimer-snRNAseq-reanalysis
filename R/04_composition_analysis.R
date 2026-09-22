# ============================================================
# Alzheimer snRNA-seq Reanalysis
# Cell-type composition analysis
#
# This script should be run after:
# R/01_preprocessing_harmony.R
# R/03_celltype_annotation.R
#
# IMPORTANT:
# Disease status is fully confounded with expression-inferred
# sex in this dataset. Therefore, AD vs Control differences
# shown here must be interpreted as descriptive/exploratory
# rather than disease-specific effects.
# ============================================================


# ------------------------------------------------------------
# 1. Required libraries
# ------------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)


# ------------------------------------------------------------
# 2. Check required metadata
# ------------------------------------------------------------

if (!exists("combined")) {
  stop(
    "The 'combined' Seurat object was not found. ",
    "Run the previous analysis scripts first."
  )
}

if (!"cell_type" %in% colnames(combined@meta.data)) {
  stop(
    "The 'cell_type' metadata column was not found. ",
    "Run R/03_celltype_annotation.R first."
  )
}


# ------------------------------------------------------------
# 3. Cell counts by condition
# ------------------------------------------------------------

composition_table <- table(
  combined$cell_type,
  combined$condition
)

print(composition_table)


# ------------------------------------------------------------
# 4. Cell-type percentages within each condition
# ------------------------------------------------------------

composition_percent <- prop.table(
  composition_table,
  margin = 2
) * 100

print(
  round(
    composition_percent,
    2
  )
)


# ------------------------------------------------------------
# 5. Cell counts by individual sample
# ------------------------------------------------------------

sample_composition <- table(
  combined$cell_type,
  combined$sample_id
)

print(sample_composition)


# ------------------------------------------------------------
# 6. Cell-type percentages within each sample
# ------------------------------------------------------------

sample_composition_percent <- prop.table(
  sample_composition,
  margin = 2
) * 100

print(
  round(
    sample_composition_percent,
    2
  )
)


# ------------------------------------------------------------
# 7. Total number of cells per sample
# ------------------------------------------------------------

sample_cell_counts <- table(
  combined$sample_id
)

print(sample_cell_counts)


# ------------------------------------------------------------
# 8. Convert condition-level percentages to data frame
# ------------------------------------------------------------

condition_composition_df <- as.data.frame(
  composition_percent
)

colnames(condition_composition_df) <- c(
  "cell_type",
  "condition",
  "percent"
)

condition_composition_df$percent <- round(
  condition_composition_df$percent,
  2
)

print(condition_composition_df)


# ------------------------------------------------------------
# 9. Convert sample-level percentages to data frame
# ------------------------------------------------------------

sample_composition_df <- as.data.frame(
  sample_composition_percent
)

colnames(sample_composition_df) <- c(
  "cell_type",
  "sample_id",
  "percent"
)

sample_composition_df$percent <- round(
  sample_composition_df$percent,
  2
)

print(sample_composition_df)


# ------------------------------------------------------------
# 10. Add condition information to sample-level table
# ------------------------------------------------------------

sample_condition_map <- combined@meta.data %>%
  dplyr::select(
    sample_id,
    condition
  ) %>%
  dplyr::distinct()


sample_composition_df <- sample_composition_df %>%
  dplyr::left_join(
    sample_condition_map,
    by = "sample_id"
  )


print(sample_composition_df)


# ------------------------------------------------------------
# 11. Sample-level composition plot
# ------------------------------------------------------------

ggplot(
  sample_composition_df,
  aes(
    x = sample_id,
    y = percent,
    fill = cell_type
  )
) +
  geom_col() +
  labs(
    title = "Cell-type Composition by Sample",
    x = "Sample",
    y = "Cell-type Percentage",
    fill = "Cell Type"
  ) +
  theme_classic()


# ------------------------------------------------------------
# 12. Condition-level composition plot
# ------------------------------------------------------------

ggplot(
  condition_composition_df,
  aes(
    x = condition,
    y = percent,
    fill = cell_type
  )
) +
  geom_col() +
  labs(
    title = "Cell-type Composition by Condition",
    x = "Condition",
    y = "Cell-type Percentage",
    fill = "Cell Type"
  ) +
  theme_classic()


# ------------------------------------------------------------
# 13. Interpretation reminder
# ------------------------------------------------------------

message(
  "Condition-level differences should not be interpreted ",
  "as Alzheimer-specific effects because disease status and ",
  "expression-inferred sex are fully confounded. ",
  "Sample-level percentages should be inspected to assess ",
  "donor-to-donor variability."
)
