# app.R
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(hdf5r))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggplot2))
library(patchwork)
library(fst)
library(tidyverse)

# Global variable for each docker session
raw_expr_data <- NULL
obj <- NULL
filename <- 'Zeisel_expression.fst'
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
function(req,filename,min_cells=1,min_genes=200, nVariableFeatures=3000){
  raw_expr_data <<- read.fst(filename)
  rownames(raw_expr_data) <- NULL
  raw_expr_data <- column_to_rownames(raw_expr_data, "X1")
  raw_obj <- CreateSeuratObject(raw_expr_data)
  obj <<- CreateSeuratObject(raw_expr_data,min.cells = min_cells, min.features = min_genes)
  empty_ident <- as.factor(obj$orig.ident)
  levels(empty_ident) <- rep("empty.ident",length(levels(empty_ident)))
  obj <<- AddMetaData(obj, metadata = empty_ident,col.name = "empty.ident")
  Idents(obj) <- obj$orig.ident
  obj[["percent.mt"]] <<- PercentageFeatureSet(obj, pattern = "^MT-")
  rb.genes <- rownames(obj)[grep("^RP[SL]",rownames(obj))]
  percent.ribo <- Matrix::colSums(obj[rb.genes,])/Matrix::colSums(obj)*100
  obj <<- AddMetaData(obj, percent.ribo, col.name = "percent.ribo")
  raw_percent_zero <- length(which((as.matrix(GetAssayData(raw_obj)) > 0)))/length(GetAssayData(raw_obj))
  filter_percent_zero <- length(which((as.matrix(GetAssayData(obj)) > 0)))/length(GetAssayData(obj))
  raw_mean_expr <- mean(as.matrix(GetAssayData(raw_obj)))
  filter_mean_expr <- mean(as.matrix(GetAssayData(obj)))
  obj <<- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nVariableFeatures,verbose = F)
  return(list(
    raw_n_genes = dim(raw_obj)[1],
    raw_n_cells = dim(raw_obj)[2],
    raw_percent_zero = raw_percent_zero,
    raw_mean_expr = raw_mean_expr,
    filter_n_genes = dim(obj)[1],
    filter_n_cells = dim(obj)[2],
    filter_percent_zero = filter_percent_zero,
    filter_mean_expr = filter_mean_expr
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
# @param filename
# @param type Upload expression file type, CellGene, 10X h5, 10X folder
#' @post /cluster
function(req, nPCs=20, resolution=0.2){
  nPCs <- as.numeric(nPCs)
  resolution <- as.numeric(resolution)
  obj <<- NormalizeData(obj,verbose = F)
  obj <<- ScaleData(obj, features = rownames(obj),verbose = F)
  variable_genes <- VariableFeatures(obj)
  obj <<- AddMetaData(obj, meta$Label, col.name = "cell_type")
  #Idents(obj) <- obj$cell_type
  #print(obj)
  if(nPCs > ncol(obj)){
    nPC <- ncol(obj)
  }
  obj <<- RunPCA(obj, features = variable_genes, npcs = nPCs,verbose = F)
  obj <<- FindNeighbors(obj, dims = 1:nPCs,verbose = F)
  obj <<- FindClusters(obj, resolution = resolution,verbose = F)
  obj <<- RunUMAP(obj, dims = 1:nPCs,n.neighbors = 15,verbose = F)
  return(list(
    n_seurat_clusters = length(levels(obj$seurat_clusters))
  ))
}

#' Plot QC
#' @get /qcplot
#' @serializer png list(width = 800, height =500)
function(){
  Idents(obj) <- obj@meta.data$empty.ident
  plot <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  return(print(plot))
}

#' Plot Variable genes
#' @get /var-genes-plot
#' @serializer png list(width = 800, height =500)
function(){
  top10 <- head(VariableFeatures(obj), 10)
  plot1 <- VariableFeaturePlot(obj)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot <- CombinePlots(plots = list(plot1, plot2))
  plot <- plot1+plot2
  return(print(plot))
}

#' Plot umap
#' @get /umap-cluster
#' @serializer png list(width = 500, height =500)
function(){
  Idents(obj) <- obj@meta.data[,ident_idx]
  plot <- DimPlot(obj,reduction = "umap")
  return(print(plot))
}

#' Plot features in umap
#' @post /umap-gene
#' @serializer png list(width = 500, height =500)
function(req, gene="BRCA1"){
  Idents(obj) <- obj@meta.data[,ident_idx]
  plot <- FeaturePlot(obj,gene)
  return(print(plot))
}

#' Plot features in violin
#' @post /violin-gene
#' @serializer png list(width = 500, height =500)
function(req, gene="BRCA1"){
  Idents(obj) <- obj@meta.data[,ident_idx]
  plot <- VlnPlot(obj, gene)
  return(print(plot))
}


#' Plot features in violin
#' @post /violin-gene
#' @serializer png list(width = 500, height =500)
function(req, gene="BRCA1"){
  Idents(obj) <- obj@meta.data[,ident_idx]
  plot <- VlnPlot(obj, gene)
  return(print(plot))
}