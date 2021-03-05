

options(future.globals.maxSize = 8000 * 1024^2)
suppressPackageStartupMessages(library(fst))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(Polychrome))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(future))


## Do not use it, not working in OSC clusters
## Set multi-thread for Seurat
#plan("multiprocess", workers = 16)
#plan()


args <- commandArgs(TRUE)
wd <- args[1] # working directory
atlas_filename <- args[2] # rds seurat object
query_filename <- args[3] # raw filename
query_data_id <- args[4] # query data ID

load_test_data <- function(){
  # This function is used for testing, set wd to your working directory
  rm(list = ls(all = TRUE))
  wd <- 'C:/Users/flyku/Documents/GitHub/iris3-backend/plumber/data'
  atlas_filename <- "mouse_brain_atlas.rds"
  query_filename <- "Zeisel_expression.fst"
  query_data_id <- "query_example"
}


setwd(wd)
source("functions.R")

####### Load raw files
#atlas.obj <- atlas
atlas.obj <- read_rds(atlas_filename)
query_matrix <- read.fst(query_filename)

rownames(query_matrix) <- NULL
query_matrix <- column_to_rownames(query_matrix, var = "X1")
query.obj <- CreateSeuratObject(query_matrix, project = "all", min.cells = 5)

# Preview atlas object cell types
#Idents(atlas.obj) <- atlas.obj$predicted.id
#Plot.cluster2D(atlas.obj,reduction.method = "umap",pt_size = 0.1, txt = "Predicted.id")

query.obj <- FindVariableFeatures(query.obj, selection.method = "vst", nfeatures = 2000)
query.obj.gene <- rownames(query.obj)
query.obj <- ScaleData(query.obj, features = query.obj.gene)
query.obj <- RunPCA(query.obj, features = VariableFeatures(object = query.obj))
query.obj <- RunUMAP(query.obj, reduction = "pca", dims = 1:25)


## FindTransferAnchors: We recommend using PCA when reference and query datasets are from scRNA-seq
transfer.anchors <- FindTransferAnchors(reference = atlas.obj, query = query.obj, features = VariableFeatures(object = atlas.obj), reduction = "pcaproject",verbose = TRUE)

if(nrow(transfer.anchors@anchors) > 30) {
  celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = atlas.obj$cell_type, weight.reduction = query.obj[["pca"]],l2.norm = FALSE,dims = 1:25, k.weight = 30)
} else{
  celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = atlas.obj$cell_type, weight.reduction = query.obj[["pca"]],l2.norm = FALSE,dims = 1:25, k.weight = (nrow(transfer.anchors@anchors)-1))
}

query.obj <- AddMetaData(query.obj, metadata = celltype.predictions)

#Idents(atlas.obj) <- atlas.obj$cell_type
#p1 <- Plot.cluster2D(atlas.obj, reduction.method = "umap",pt_size = 0.4,txt = "atlas cell type")

Idents(query.obj) <- query.obj$predicted.id
p2 <- Plot.cluster2D(query.obj, reduction.method = "umap",pt_size = 0.4,txt = "query cell type")

png(paste(query_data_id,"_transfer_umap.png",sep = ""),width=2000, height=2000,res=300)
plot_grid(p2)
dev.off()

# Save Seurat object
Idents(query.obj) <- query.obj$predicted.id


#saveRDS(query.obj, paste0(query_data_id,".rds"))

# Save cell type labels
cell_info <- query.obj$predicted.id
cell_label <- cbind(colnames(query.obj),as.character(cell_info))
colnames(cell_label) <- c("cell_name","label")
cell_label <- cell_label[order(cell_label[,1]),]
write.table(cell_label,paste(query_data_id,"_cell_label.txt",sep = ""),quote = F,row.names = F,sep = "\t")

meta <- read.csv('Zeisel_index_label.csv')
query.obj <<- AddMetaData(query.obj, meta$Label, col.name = "cell_type")
query.obj <<- AddMetaData(query.obj, meta$Sex, col.name = "sex")
query.obj <<- AddMetaData(query.obj, meta$Sample, col.name = "sample")
print(ident_idx)

Idents(query.obj) <- query.obj$cell_type
p1 <- Plot.cluster2D(query.obj, reduction.method = "umap",pt_size = 0.4,txt = "provided cell type")

png(paste(query_data_id,"_provided_umap.png",sep = ""),width=2000, height=2000,res=300)
plot_grid(p1)
dev.off()


# Session Infomation
#sessionInfo()

