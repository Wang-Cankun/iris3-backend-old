library(Seurat)
library(patchwork)
library(loomR)
library(bayNorm)
library(tidyverse)

setwd('C:/Users/flyku/Documents/GitHub/iris3-backend/plumber')


mca.matrix <- readRDS(file = "./data/MCA_merged_mat.rds")
mca.metadata <- read.csv(file = "./data/MCA_All-batch-removed-assignments.csv", row.names = 1)


mca <- CreateSeuratObject(counts = mca.matrix, meta.data = mca.metadata, project = "MouseCellAtlas")
# Only keep annotated cells
mca <- subset(mca, cells = names(which(!is.na(mca$ClusterID))))
# Leaves us with 242k cells
mca

mca <- NormalizeData(mca, normalization.method = "LogNormalize", scale.factor = 10000)
mca <- FindVariableFeatures(mca)


levels(as.factor(mca$ClusterID))


lfile <- connect(filename = "./data/mouse_brain_zeisel_2018.loom", mode = "r+")
lfile

lfile[["col_attrs"]]

lfile[['matrix']]

lfile$row.attrs$Gene




dim(full.matrix)

lfile[["row_attrs/Gene"]]
lfile[["row_attrs/Gene"]]

lfile$row.attrs$Gene[]
lfile$col.attrs$CellID[]

lfile$col.attrs$Class[]

meta <- data.frame(cell = lfile$col.attrs$CellID[], class = lfile$col.attrs$Class[])
meta <- meta[which(!duplicated(meta$cell)),]
rownames(meta) <- NULL
meta <- column_to_rownames(meta,"cell")

full.matrix <- lfile$matrix[, ]

full.matrix <- as.sparse(full.matrix)
full.matrix <- t_sp(full.matrix)
rownames(full.matrix) <- lfile$row.attrs$Gene[]
colnames(full.matrix) <- lfile$col.attrs$CellID[]
full.matrix <- which(!duplicated(rownames(full.matrix)))

atlas <- CreateSeuratObject(counts = full.matrix, meta.data = meta, project = "MouseBrainAtlas")
atlas <- AddMetaData(atlas, meta, col.name = "cell_type")
atlas.obj <- FindVariableFeatures(atlas.obj, selection.method = "vst", nfeatures = 2000)

saveRDS(atlas,"mouse_brain_atlas.rds")
