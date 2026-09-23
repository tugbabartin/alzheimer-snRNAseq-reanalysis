# Alzheimer’s Disease snRNA-seq Reanalysis

Reanalysis of a public Alzheimer’s disease single-nucleus RNA sequencing (snRNA-seq) dataset using R, Seurat, Harmony, and complementary single-cell analysis tools.

## Overview

This project reanalyzes a publicly available Alzheimer’s disease snRNA-seq dataset containing four human brain samples:

| Sample | Reported Group |
|---|---|
| A1 | Alzheimer’s disease |
| A2 | Control |
| A3 | Alzheimer’s disease |
| A4 | Control |

The main objectives were to:

- process and quality-control raw gene-expression matrices,
- identify and remove potential doublets,
- integrate samples using Harmony,
- perform dimensionality reduction and clustering,
- identify major brain cell populations using marker genes,
- examine cell-type composition across samples,
- independently investigate sample sex using sex-linked gene expression,
- evaluate whether disease-associated expression analysis could be performed reliably.

## Analysis Workflow

The analysis was performed in R using a Seurat-based workflow.

### 1. Preprocessing and Quality Control

Raw count matrices were imported separately for each sample and converted into Seurat objects.

Quality-control filtering included:

- minimum detected features per nucleus,
- minimum RNA counts,
- mitochondrial RNA percentage filtering,
- sample-specific metadata assignment,
- doublet detection using `scDblFinder`.

Samples were subsequently merged for downstream analysis.

### 2. Normalization and Integration

The merged dataset was processed using:

- `NormalizeData`
- `FindVariableFeatures`
- `ScaleData`
- Principal Component Analysis (PCA)
- Harmony batch integration

Harmony was performed using sample identity as the integration variable.

The Harmony representation was then used for:

- nearest-neighbor graph construction,
- UMAP visualization,
- unsupervised clustering.

### 3. Cell-Type Annotation

Cluster-specific marker genes were identified using Seurat.

Cell identities were evaluated using canonical markers for major brain cell populations, including:

- excitatory neurons,
- inhibitory neurons,
- astrocytes,
- oligodendrocytes,
- oligodendrocyte precursor cells (OPCs),
- microglia,
- endothelial cells.

Ambiguous clusters were further evaluated using additional lineage-specific marker panels before final annotation.

## Sex-Linked Expression Validation

During the analysis, sex-linked genes were independently examined to investigate inconsistencies in the sample metadata.

Genes evaluated included:

- `XIST`
- `RPS4Y1`
- `KDM5D`
- `UTY`
- `EIF1AY`
- `DDX3Y`
- `ZFY`

The observed expression pattern indicated:

| Sample | Reported Group | Expression-Inferred Sex Pattern |
|---|---|---|
| A1 | AD | Female-like |
| A2 | Control | Male-like |
| A3 | AD | Female-like |
| A4 | Control | Male-like |

This revealed an important limitation of the experimental design:

**Disease status and expression-inferred sex are completely confounded in this dataset.**

Both Alzheimer’s disease-labelled samples show a female-like expression pattern, while both control-labelled samples show a male-like expression pattern.

As a result, an observed AD-vs-Control expression difference cannot be reliably separated into disease-associated and sex-associated effects.

## Differential Expression Considerations

Exploratory differential-expression analyses were investigated during development of the project.

Sex-linked genes such as `XIST`, `UTY`, `NLGN4Y`, `TTTY14`, and `USP9Y` appeared among the strongest group-associated signals.

Because the dataset contains only four biological donors and disease condition is completely confounded with expression-inferred sex, disease-specific differential-expression inference is not considered reliable.

For this reason, the repository focuses primarily on:

- preprocessing,
- data integration,
- cell-type identification,
- sample-level validation,
- cell-type composition,
- identification of experimental-design limitations.

This limitation is treated as an analytical result rather than attempting to force a disease-specific interpretation from the dataset.

## Repository Structure

```text
alzheimer-snRNAseq-reanalysis/
│
├── R/
│   ├── 01_preprocessing_harmony.R
│   ├── 02_sex_validation.R
│   ├── 03_celltype_annotation.R
│   └── 04_composition_analysis.R
│
├── .gitignore
└── README.md
```

### Scripts

**`01_preprocessing_harmony.R`**  
Data import, quality control, doublet detection, normalization, PCA, Harmony integration, UMAP, and clustering.

**`02_sex_validation.R`**  
Evaluation of XIST and Y-linked gene expression across samples.

**`03_celltype_annotation.R`**  
Cluster marker analysis and annotation of major brain cell populations.

**`04_composition_analysis.R`**  
Descriptive comparison of cell-type abundance across individual samples and reported disease groups.

## Main R Packages

The workflow uses packages including:

- Seurat
- Harmony
- SingleR
- celldex
- scDblFinder
- SingleCellExperiment
- Matrix
- dplyr
- ggplot2
- clusterProfiler
- org.Hs.eg.db

## Data

Raw sequencing matrices are not stored in this repository because of their size.

The analysis uses publicly available samples corresponding to:

- `GSM5348374` – A1
- `GSM5348375` – A2
- `GSM5348376` – A3
- `GSM5348377` – A4

Local file paths in the preprocessing script may need to be changed before running the analysis on another system.

## Important Limitations

This dataset contains only four biological donors:

- 2 Alzheimer’s disease-labelled samples
- 2 control-labelled samples

In addition, expression-inferred sex is perfectly aligned with the reported disease groups.

Therefore:

> Group-associated differences observed in this dataset should not be interpreted as Alzheimer’s disease-specific effects.

The small number of biological replicates also substantially limits statistical inference at the donor level.

## Reference

**Single cell RNA-Seq Alzheimer’s brain**  
PLOS ONE, 2023  
DOI: `10.1371/journal.pone.0277630`

## Notes

This repository was created as a reanalysis and learning project focused on single-nucleus RNA-seq analysis, quality control, cell-type annotation, and critical evaluation of experimental design.
