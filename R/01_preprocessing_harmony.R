# ============================================================
# Alzheimer snRNA-seq Reanalysis
# Initial processing, QC, Harmony integration and clustering
# Samples:
# A1 = AD
# A2 = Control
# A3 = AD
# A4 = Control
# ============================================================


# ------------------------------------------------------------
# 1. Libraries
# ------------------------------------------------------------

library(SingleR)
library(celldex)
library(Matrix)
library(Seurat)
library(dplyr)
library(harmony)
library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGGREST)
library(ggplot2)


# ------------------------------------------------------------
# 2. Data paths
# ------------------------------------------------------------

# A1 - AD
mtx_file_A1 <- "C:/R/GSM5348374_A1matrix.mtx"
barcodes_file_A1 <- "C:/R/GSM5348374_A1barcodes.tsv"
features_file_A1 <- "C:/R/GSM5348374_A1features.tsv"

# A2 - Control
mtx_file_A2 <- "C:/R/GSM5348375_A2matrix.mtx"
barcodes_file_A2 <- "C:/R/GSM5348375_A2barcodes.tsv"
features_file_A2 <- "C:/R/GSM5348375_A2features.tsv"

# A3 - AD
mtx_file_A3 <- "C:/R/GSM5348376_A3matrix.mtx"
barcodes_file_A3 <- "C:/R/GSM5348376_A3barcodes.tsv"
features_file_A3 <- "C:/R/GSM5348376_A3features.tsv"

# A4 - Control
mtx_file_A4 <- "C:/R/GSM5348377_A4matrix.mtx"
barcodes_file_A4 <- "C:/R/GSM5348377_A4barcodes.tsv"
features_file_A4 <- "C:/R/GSM5348377_A4features.tsv"


# ------------------------------------------------------------
# 3. Read raw matrices
# ------------------------------------------------------------

# A1
matrix_A1 <- readMM(mtx_file_A1)

barcodes_A1 <- read.table(
  barcodes_file_A1,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

features_A1 <- read.table(
  features_file_A1,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)


# A2
matrix_A2 <- readMM(mtx_file_A2)

barcodes_A2 <- read.table(
  barcodes_file_A2,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

features_A2 <- read.table(
  features_file_A2,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)


# A3
matrix_A3 <- readMM(mtx_file_A3)

barcodes_A3 <- read.table(
  barcodes_file_A3,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

features_A3 <- read.table(
  features_file_A3,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)


# A4
matrix_A4 <- readMM(mtx_file_A4)

barcodes_A4 <- read.table(
  barcodes_file_A4,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

features_A4 <- read.table(
  features_file_A4,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 4. Add gene and cell names
# ------------------------------------------------------------

# A1
rownames(matrix_A1) <- make.unique(features_A1[, 2])
colnames(matrix_A1) <- barcodes_A1[, 1]

# A2
rownames(matrix_A2) <- make.unique(features_A2[, 2])
colnames(matrix_A2) <- barcodes_A2[, 1]

# A3
rownames(matrix_A3) <- make.unique(features_A3[, 2])
colnames(matrix_A3) <- barcodes_A3[, 1]

# A4
rownames(matrix_A4) <- make.unique(features_A4[, 2])
colnames(matrix_A4) <- barcodes_A4[, 1]


# ------------------------------------------------------------
# 5. Create Seurat objects
# ------------------------------------------------------------

seurat_A1 <- CreateSeuratObject(
  counts = matrix_A1,
  project = "A1"
)

seurat_A2 <- CreateSeuratObject(
  counts = matrix_A2,
  project = "A2"
)

seurat_A3 <- CreateSeuratObject(
  counts = matrix_A3,
  project = "A3"
)

seurat_A4 <- CreateSeuratObject(
  counts = matrix_A4,
  project = "A4"
)


# ------------------------------------------------------------
# 6. Calculate mitochondrial percentage
# ------------------------------------------------------------

seurat_A1[["percent.mt"]] <- PercentageFeatureSet(
  seurat_A1,
  pattern = "^MT-"
)

seurat_A2[["percent.mt"]] <- PercentageFeatureSet(
  seurat_A2,
  pattern = "^MT-"
)

seurat_A3[["percent.mt"]] <- PercentageFeatureSet(
  seurat_A3,
  pattern = "^MT-"
)

seurat_A4[["percent.mt"]] <- PercentageFeatureSet(
  seurat_A4,
  pattern = "^MT-"
)


# ------------------------------------------------------------
# 7. Quality control
# ------------------------------------------------------------

seurat_A1 <- subset(
  seurat_A1,
  subset =
    nFeature_RNA > 200 &
    nCount_RNA > 10 &
    percent.mt < 15
)

seurat_A2 <- subset(
  seurat_A2,
  subset =
    nFeature_RNA > 200 &
    nCount_RNA > 10 &
    percent.mt < 15
)

seurat_A3 <- subset(
  seurat_A3,
  subset =
    nFeature_RNA > 200 &
    nCount_RNA > 10 &
    percent.mt < 15
)

seurat_A4 <- subset(
  seurat_A4,
  subset =
    nFeature_RNA > 200 &
    nCount_RNA > 10 &
    percent.mt < 15
)


# ------------------------------------------------------------
# 8. Sample metadata
# ------------------------------------------------------------

# Sample ID
seurat_A1$sample_id <- "A1"
seurat_A2$sample_id <- "A2"
seurat_A3$sample_id <- "A3"
seurat_A4$sample_id <- "A4"

# Disease condition
seurat_A1$condition <- "AD"
seurat_A2$condition <- "Control"
seurat_A3$condition <- "AD"
seurat_A4$condition <- "Control"


# ------------------------------------------------------------
# 9. Make cell names unique across samples
# ------------------------------------------------------------

seurat_A1 <- RenameCells(
  seurat_A1,
  add.cell.id = "A1"
)

seurat_A2 <- RenameCells(
  seurat_A2,
  add.cell.id = "A2"
)

seurat_A3 <- RenameCells(
  seurat_A3,
  add.cell.id = "A3"
)

seurat_A4 <- RenameCells(
  seurat_A4,
  add.cell.id = "A4"
)


# ------------------------------------------------------------
# 10. Merge samples
# ------------------------------------------------------------

combined <- merge(
  seurat_A1,
  y = list(
    seurat_A2,
    seurat_A3,
    seurat_A4
  )
)


# Check sample numbers
table(combined$sample_id)
table(combined$condition)


# ------------------------------------------------------------
# 11. Normalization and feature selection
# ------------------------------------------------------------

DefaultAssay(combined) <- "RNA"

combined <- NormalizeData(combined)

combined <- FindVariableFeatures(combined)

combined <- ScaleData(combined)


# ------------------------------------------------------------
# 12. PCA
# ------------------------------------------------------------

combined <- RunPCA(
  combined,
  npcs = 10
)


# ------------------------------------------------------------
# 13. Harmony integration
# ------------------------------------------------------------

combined <- RunHarmony(
  combined,
  group.by.vars = "sample_id"
)


# ------------------------------------------------------------
# 14. UMAP and clustering
# ------------------------------------------------------------

combined <- RunUMAP(
  combined,
  reduction = "harmony",
  dims = 1:10
)

combined <- FindNeighbors(
  combined,
  reduction = "harmony",
  dims = 1:10
)

combined <- FindClusters(
  combined,
  resolution = 0.5
)

Idents(combined) <- "seurat_clusters"


# ------------------------------------------------------------
# 15. UMAP plots
# ------------------------------------------------------------

DimPlot(
  combined,
  group.by = "sample_id",
  label = TRUE
) +
  ggtitle("Sample-based UMAP - Harmony")


DimPlot(
  combined,
  group.by = "condition",
  label = TRUE
) +
  ggtitle("Condition-based UMAP - Harmony")


DimPlot(
  combined,
  label = TRUE
) +
  ggtitle("Cluster-based UMAP - Harmony")


# ------------------------------------------------------------
# 16. Cluster composition
# ------------------------------------------------------------

table(
  combined$seurat_clusters,
  combined$sample_id
)

table(
  combined$seurat_clusters,
  combined$condition
)


# ------------------------------------------------------------
# 17. Prepare reference for cell-type annotation
# ------------------------------------------------------------

ref <- HumanPrimaryCellAtlasData()


# ------------------------------------------------------------
# 18. Inspect Seurat v5 layers
# ------------------------------------------------------------

Layers(combined)


# Inspect A1 normalized expression layer
layer_A1 <- LayerData(
  combined,
  layer = "data.A1"
)

dim(layer_A1)

layer_A1[1:5, 1:5]


# ------------------------------------------------------------
# 19. Combine normalized expression layers
# ------------------------------------------------------------

samples <- unique(combined$sample_id)

mat_list <- list()

for (sample in samples) {
  
  layer_name <- paste0(
    "data.",
    sample
  )
  
  mat_part <- LayerData(
    combined,
    layer = layer_name
  )
  
  mat_list[[sample]] <- mat_part
}


expr_matrix <- do.call(
  cbind,
  mat_list
)


dim(expr_matrix)

