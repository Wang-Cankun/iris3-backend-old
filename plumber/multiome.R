# app.R
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(hdf5r))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggplot2))
library(patchwork)
library(fst)
library(dplyr)
library(tibble)
library(scales)
library(tinytex)
library(hdf5r)
library(patchwork)
library(dplyr)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
#install.packages("rmarkdown")
library(rmarkdown)
#BiocManager::install("rtracklayer")
library(rtracklayer)
library(ggplot2)
#install.packages('tinytex')

library(BSgenome.Hsapiens.UCSC.hg38)


# Global variable for each docker session
raw_expr_data <- NULL
obj <- readRDS("pbmc2.rds")
da_peaks <- readRDS("da_peaks.rds")
variable_genes <- NULL
filename <- 'pbmc.rds'
Idents(obj) <- obj$orig.ident
meta <- read.csv('Zeisel_index_label.csv')
ident_idx <- 1

#' Echo the parameter that was sent in
#' @param msg The message to echo back.
#' @get /echo
function(msg=""){
  return( 
    list(
      wd = getwd()
    ))
}

#' Plot out data from the iris dataset
#' @param spec If provided, filter the data to only this species (e.g. 'setosa')
#' @get /plot
#' @serializer png
function(spec){
  myData <- iris
  title <- "All Species"
  
  # Filter if the species was specified
  if (!missing(spec)){
    title <- paste0("Only the '", spec, "' Species")
    myData <- subset(iris, Species == spec)
  }
  print("Run plot")
  plot <- plot(myData$Sepal.Length, myData$Petal.Length,
               main=title, xlab="Sepal Length", ylab="Petal Length")
  #png("plot.png", plot)
  #return(base64enc::dataURI(file = "plot.png", mime = "image/png"))
  return(plot)
}

#' Read data into Seurat object
#' @param jobid
#' @param type Upload expression file type, CellGene, 10X h5, 10X folder
#' @post /job
function(jobid){
  dir.create(jobid,showWarnings = F)
}

#' Read data into Seurat object
# @param filename
# @param type Upload expression file type, CellGene, 10X h5, 10X folder
#' @post /load
function(req,filename,min_cells=1,min_genes=200, nVariableFeatures=3000, percentMt=5, removeRibosome=FALSE){
  print(obj)
  DefaultAssay(obj) <- "RNA"
  filter_n_genes <- dim(obj)[1]
  DefaultAssay(obj) <- "ATAC"
  return(list(
    raw_n_genes = dim(obj)[1], #peaks
    raw_n_cells = dim(obj)[2],
    raw_percent_zero = 0,
    raw_mean_expr = 0,
    filter_n_genes = filter_n_genes,
    filter_n_cells = dim(obj)[2],
    filter_percent_zero = 0,
    filter_mean_expr = 0
  ))
}

#' Get all gene names
#' @get /genes
function(){
  return(rownames(obj))
}

#' Get all Seurat Idents names
#' @get /idents
function(){
  return(colnames(obj@meta.data))
}

#' Set Seurat Idents by name
#' @post /idents
#' 
function(req,name='orig.ident'){
  ident_idx <<- which(colnames(obj@meta.data) == name)
  print(ident_idx)
  return(ident_idx)
}



#' Read data into Seurat object
#' @param filename
#' @param type Upload expression file type, CellGene, 10X h5, 10X folder
#' @post /cluster
function(req, nPCs=20, resolution=0.5){
  print('run cluster calculation')
  print(nPCs)
  nPCs <- as.numeric(nPCs)
  resolution <- as.numeric(resolution)
  obj <- NormalizeData(obj,verbose = F)
  obj <- ScaleData(obj, features = rownames(obj),verbose = F)
  variable_genes <<- VariableFeatures(obj)
  #Idents(obj) <- obj$cell_type
  #print(obj)
  if(nPCs > ncol(obj)){
    nPC <- ncol(obj)
  }
  obj <- RunPCA(obj, features = variable_genes, npcs = nPCs,verbose = F)
  obj <- FindNeighbors(obj, dims = 1:nPCs,verbose = F)
  obj <- FindClusters(obj, resolution = resolution,verbose = F)
  obj <<- RunUMAP(obj, dims = 1:nPCs,n.neighbors = 15,verbose = F)
  return(list(
    n_seurat_clusters = length(levels(obj$seurat_clusters))
  ))
}


#' Plot QC
#' @get /qcplot
#' @serializer png list(width = 1200, height =500)
function(){
  plot <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
  return(print(plot))
}

#' Plot nFeature_RNA
#' @get /qcplot1
#' @serializer png list(width = 500, height =500)
function(){
  plot <- VlnPlot(obj, features = c("nFeature_RNA"), ncol = 1)
  return(print(plot))
}

#' Plot nCount_RNA
#' @get /qcplot2
#' @serializer png list(width = 500, height =500)
function(){
  plot <- VlnPlot(obj, features = c("nCount_RNA"), ncol = 1)
  return(print(plot))
}

#' Plot percent.mt
#' @get /qcplot3
#' @serializer png list(width = 500, height =500)
function(){
  plot <- VlnPlot(obj, features = c("percent.mt"), ncol = 1)
  return(print(plot))
}

#' Plot percent.ribo
#' @get /qcplot4
#' @serializer png list(width = 500, height =500)
function(){
  plot <- VlnPlot(object = obj,features ='nucleosome_signal')
  return(print(plot))
}

#' Plot Variable genes
#' @get /var-genes-plot
#' @serializer png list(width = 1200, height =500)
function(){
  plot <- TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()
  return(print(plot))
}

#' Plot frag hist
#' @get /qcplot5
#' @serializer png list(width = 500, height =500)
function(){
  plot <- FragmentHistogram(object = obj, group.by = "nucleosome_group")
  return(print(plot))
}

#' Plot gene count
#' @get /qcplot6
#' @serializer png list(width = 500, height =500)
function(){
  plot <- VlnPlot(object = obj,features ='gex_genes_count')
  return(print(plot))
}



#' Get Variable genes list
#' @get /var-genes-list
function(){
  vargenes <- VariableFeatures(obj)
  Idents(obj) <- obj@meta.data$empty_ident
  this_obj <- as.matrix(GetAssayData(subset(obj, features = vargenes), slot = "data"))
  row_min <- apply(this_obj, 1, min)
  row_sd <- apply(this_obj, 1, sd)
  row_max <- apply(this_obj, 1, max)
  result <- data.frame(gene=rownames(this_obj), mean=rowMeans(this_obj), std=row_sd, min=row_min, max=row_max)
  return(list(result))
}


#' Plot umap
#' @get /umap-cluster
#' @serializer png list(width = 700, height = 600)
function(){
  plot <- VlnPlot(pbmc, features ="nCount_ATAC",log = TRUE)
  return(print(plot))
}

#'  UMAP RNA
#' @get /umap-rna
#' @serializer png list(width = 700, height = 600)
function(){
  plot <- DimPlot(obj, reduction = "umap.rna",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
  return(print(plot))
}

#'  UMAP atac
#' @get /umap-atac
#' @serializer png list(width = 700, height = 600)
function(){
  plot <- DimPlot(pbmc, reduction = "umap.atac",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
  return(print(plot))
}

#'  UMAP integrated
#' @get /umap-integrated
#' @serializer png list(width = 700, height = 600)
function(){
  plot <- DimPlot(pbmc, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  return(print(plot))
}


#' Plot cell type metadata pie chart
#' @get /pie-meta
#' @serializer png list(width = 600, height = 600)
function(){
  #plot <- VlnPlot(object = pbmc,features ='blacklist_ratio')
  
  meta <- read.csv('Zeisel_index_label.csv')
  #obj <<- AddMetaData(obj, meta$Label, col.name = "cell_type")
  #print(ident_idx)
  #Idents(obj) <- obj@meta.data[,7]
  
  #Idents(obj) <- obj@meta.data[,ident_idx]
  mytable <- table(meta$Label)
  
  x <- as.vector(mytable)
  mypct <- percent(x/sum(x))
  labels <- paste(names(mytable), "\n", mytable,"(",mypct,")", sep="")
  plot <- pie(x, labels, main = "Cell type metadata", col = rainbow(length(x)))
  
  return(print(plot))
}

#' Calculate DEG for two selections
#' @post /deg

function(req, ident1=1, ident2=2, min_pct=0.2, min_lfc=0.5){
  #ident_idx=8
  print('run deg')
  Idents(obj) <- obj@meta.data[,ident_idx]
  this_markers <- FindMarkers(obj, ident.1 = ident1, ident.2 = ident2, min.pct =min_pct, logfc.threshold = min_lfc)
  this_markers <- rownames_to_column(this_markers,"gene")
  return(list(this_markers))
}


#' @post /cluster
function(req, nPCs=20, resolution=0.5){
  print('run cluster calculation')
  print(nPCs)
  nPCs <- as.numeric(nPCs)
  resolution <- as.numeric(resolution)
  obj <- NormalizeData(obj,verbose = F)
  obj <- ScaleData(obj, features = rownames(obj),verbose = F)
  variable_genes <<- VariableFeatures(obj)
  #Idents(obj) <- obj$cell_type
  #print(obj)
  if(nPCs > ncol(obj)){
    nPC <- ncol(obj)
  }
  obj <- RunPCA(obj, features = variable_genes, npcs = nPCs,verbose = F)
  obj <- FindNeighbors(obj, dims = 1:nPCs,verbose = F)
  obj <- FindClusters(obj, resolution = resolution,verbose = F)
  obj <<- RunUMAP(obj, dims = 1:nPCs,n.neighbors = 15,verbose = F)
  return(list(
    n_seurat_clusters = length(levels(obj$seurat_clusters))
  ))
}


#' Plot features in umap
#' @post /umap-gene
#' @serializer png list(width = 600, height =600)
function(req, gene="Gad1"){
  Idents(obj) <- obj@meta.data[,ident_idx]
  plot <- FeaturePlot(obj,gene)
  return(print(plot))
}

#' Plot features in violin
#' @post /violin-gene
#' @serializer png list(width = 600, height =600)
function(req, gene="Gad1"){
  Idents(obj) <- obj@meta.data[,ident_idx]
  plot <- VlnPlot(obj, gene)
  return(print(plot))
}
