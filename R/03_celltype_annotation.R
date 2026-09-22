# ============================================================
# Alzheimer snRNA-seq Reanalysis
# Cell-type annotation
#
# This script should be run after:
# R/01_preprocessing_harmony.R
#
# Cell types are assigned using cluster-specific marker genes
# and canonical brain cell-type markers.
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
# 3. Create a separate object for marker analysis
# ------------------------------------------------------------

# Keep the original combined object unchanged during
# Seurat v5 layer joining and marker discovery.

marker_obj <- combined

DefaultAssay(marker_obj) <- "RNA"

marker_obj <- JoinLayers(marker_obj)

print(Layers(marker_obj))


# ------------------------------------------------------------
# 4. Identify cluster-specific marker genes
# ------------------------------------------------------------

Idents(marker_obj) <- "seurat_clusters"

all_markers <- FindAllMarkers(
  marker_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


top10_markers <- all_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::slice_head(n = 10)


print(
  top10_markers %>%
    dplyr::select(
      cluster,
      gene,
      avg_log2FC
    ),
  n = Inf
)


# ------------------------------------------------------------
# 5. Canonical brain cell-type marker validation
# ------------------------------------------------------------

canonical_markers <- c(

  # Excitatory neurons
  "SLC17A7", "CAMK2A", "SATB2",

  # Inhibitory neurons
  "GAD1", "GAD2", "SLC6A1",

  # Astrocytes
  "AQP4", "SLC1A2", "SLC1A3", "GFAP",

  # Oligodendrocytes
  "MBP", "PLP1", "MOG", "MOBP", "MYRF",

  # OPC
  "PDGFRA", "CSPG4", "VCAN", "GPR17",

  # Microglia
  "P2RY12", "CX3CR1", "C1QA", "C1QB", "TYROBP",

  # Activated microglia
  "ITGAX", "CD86", "APOE", "LPL",

  # Endothelial cells
  "CLDN5", "VWF", "PECAM1", "FLT1"
)


canonical_markers <- intersect(
  canonical_markers,
  rownames(marker_obj)
)


DotPlot(
  marker_obj,
  features = canonical_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  ggtitle("Canonical Cell-type Markers")


# ------------------------------------------------------------
# 6. Detailed validation of ambiguous clusters
# ------------------------------------------------------------

ambiguous_clusters <- c(
  "2", "11", "14", "15", "16", "17", "18"
)


validation_markers <- c(

  # Pan-neuronal
  "RBFOX3", "SNAP25", "SYT1",

  # Excitatory neuron
  "SLC17A7", "CAMK2A", "SATB2",

  # Inhibitory neuron
  "GAD1", "GAD2", "SLC6A1",

  # Astrocyte
  "AQP4", "ALDH1L1", "SLC1A2", "SLC1A3",
  "GFAP", "SOX9",

  # OPC
  "PDGFRA", "CSPG4", "VCAN", "OLIG1", "OLIG2",

  # Oligodendrocyte lineage
  "GPR17", "SOX10", "MBP", "PLP1",
  "MOG", "MOBP", "MYRF",

  # Microglia / myeloid
  "P2RY12", "CX3CR1", "C1QA", "C1QB",
  "AIF1", "TYROBP", "PTPRC",

  # Activated microglia
  "ITGAX", "CD86", "FCGR1A", "HLA-DRA",

  # Endothelial
  "CLDN5", "PECAM1", "VWF"
)


validation_markers <- intersect(
  validation_markers,
  rownames(marker_obj)
)


DotPlot(
  marker_obj,
  features = validation_markers,
  group.by = "seurat_clusters",
  idents = ambiguous_clusters
) +
  RotatedAxis() +
  ggtitle("Detailed Validation of Ambiguous Clusters")


# ------------------------------------------------------------
# 7. Final cluster annotations
# ------------------------------------------------------------

final_cluster_names <- c(

  `0`  = "Excitatory neuron (CBLN2+)",

  `1`  = "Mature oligodendrocyte (ST18+/MOBP+)",

  `2`  = "Excitatory neuron",

  `3`  = "Astrocyte (SLCO1C1+/ETNPPL+)",

  `4`  = "Excitatory neuron (CPNE4+)",

  `5`  = "Mature oligodendrocyte (MYRF+/FA2H+)",

  `6`  = "Reactive astrocyte (GFAP+)",

  `7`  = "OPC / oligodendrocyte precursor (GPR17+)",

  `8`  = "Activated microglia (CD86+/TLR2+)",

  `9`  = "Inhibitory neuron (LHX6+/TACR1+)",

  `10` = "Inhibitory neuron (DLX6-AS1+/GPR149+)",

  `11` = "Excitatory neuron (RXFP1+)",

  `12` = "Endothelial cell (CLDN5+/VWF+)",

  `13` = "ITGAX+ activated microglia",

  `14` = "Astrocyte (PAX6+/TNC+)",

  `15` = "OPC / oligodendrocyte precursor",

  `16` = "Excitatory neuron",

  `17` = "OPC / oligodendrocyte precursor (VCAN+)",

  `18` = "Activated microglia (IRF8+/FCGR1A+)"
)


# ------------------------------------------------------------
# 8. Add annotations to Seurat metadata
# ------------------------------------------------------------

cell_type_vector <- unname(
  final_cluster_names[
    as.character(marker_obj$seurat_clusters)
  ]
)

names(cell_type_vector) <- colnames(marker_obj)

marker_obj$cell_type <- cell_type_vector


# Add the same annotations to the main combined object
combined_cell_types <- unname(
  final_cluster_names[
    as.character(combined$seurat_clusters)
  ]
)

names(combined_cell_types) <- colnames(combined)

combined$cell_type <- combined_cell_types


# ------------------------------------------------------------
# 9. Check cluster-to-cell-type mapping
# ------------------------------------------------------------

print(
  table(
    marker_obj$seurat_clusters,
    marker_obj$cell_type
  )
)


# ------------------------------------------------------------
# 10. Annotated UMAP
# ------------------------------------------------------------

Idents(marker_obj) <- "cell_type"

DimPlot(
  marker_obj,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) +
  NoLegend() +
  ggtitle("Annotated Cell Types")
