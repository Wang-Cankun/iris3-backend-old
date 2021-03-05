# app.R
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(hdf5r))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggplot2))
library(patchwork)
library(fst)
library(tidyverse)
library(tibble)
library(scales)
library(Polychrome)
library(RColorBrewer)
#library(ComplexHeatmap)
# Global variable for each docker session
raw_expr_data <- NULL
obj <- NULL
variable_genes <- NULL
filename <- 'Zeisel_expression.fst'
meta <- read.csv('Zeisel_index_label.csv')
combine_regulon <- read.table("2020041684528_combine_regulon.txt", sep = "\t", header = T)
ident_idx <- 1
workdir <- "/var/www/nodejs/iris3-backend/plumber/"


## test
min_cells=1
min_genes=200
nVariableFeatures=3000
percentMt=5
removeRibosome=FALSE
nPCs=20
resolution=0.5
k_arg=20
f_arg=0.7
o_arg=500
promoter_arg=1000

source("functions.R")




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
  raw_expr_data <<- read.fst(filename)
  rownames(raw_expr_data) <- NULL
  raw_expr_data <- column_to_rownames(raw_expr_data, "X1")
  raw_obj <- CreateSeuratObject(raw_expr_data)
  obj <<- CreateSeuratObject(raw_expr_data,min.cells = min_cells, min.features = min_genes)
  empty_ident <- as.factor(obj$orig.ident)
  levels(empty_ident) <- rep("empty_ident",length(levels(empty_ident)))
  obj <<- AddMetaData(obj, metadata = empty_ident,col.name = "empty_ident")
  Idents(obj) <- obj$orig.ident
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  #obj[["percent.ribo"]] <<- PercentageFeatureSet(obj, pattern = "^MT-")
  #obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  rb.genes <- rownames(obj)[grep("^Rp[sl][[:digit:]]",rownames(obj))]
  percent.ribo <- Matrix::colSums(obj[rb.genes,])/Matrix::colSums(obj)*100
  obj <- AddMetaData(obj, percent.ribo, col.name = "percent.ribo")
  #obj <- subset(obj, subset = percent.mt < as.numeric(percentMt))
  raw_percent_zero <- length(which((as.matrix(GetAssayData(raw_obj)) > 0)))/length(GetAssayData(raw_obj))
  filter_percent_zero <- length(which((as.matrix(GetAssayData(obj)) > 0)))/length(GetAssayData(obj))
  raw_mean_expr <- mean(as.matrix(GetAssayData(raw_obj)))
  filter_mean_expr <- mean(as.matrix(GetAssayData(obj)))
  obj <<- FindVariableFeatures(obj, selection.method = "vst", nfeatures = as.numeric(nVariableFeatures),verbose = F)

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
  obj <- RunUMAP(obj, dims = 1:nPCs,n.neighbors = 15,verbose = F)
  return(list(
    n_seurat_clusters = length(levels(obj$seurat_clusters)),
    umap_pts = data.frame(umap1=as.vector(Embeddings(obj, reduction = "umap")[,1]), umap2=as.vector(Embeddings(obj, reduction = "umap")[,2]))
  ))
}


#' Plot QC
#' @get /qcplot
#' @serializer png list(width = 1200, height =500)
function(){
  Idents(obj) <- obj@meta.data$empty_ident
  plot <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
  return(print(plot))
}

#' Plot nFeature_RNA
#' @get /qcplot1
#' @serializer png list(width = 500, height =500)
function(){
  Idents(obj) <- obj@meta.data$empty_ident
  plot <- VlnPlot(obj, features = c("nFeature_RNA"), ncol = 1)
  return(print(plot))
}

#' Plot nCount_RNA
#' @get /qcplot2
#' @serializer png list(width = 500, height =500)
function(){
  Idents(obj) <- obj@meta.data$empty_ident
  plot <- VlnPlot(obj, features = c("nCount_RNA"), ncol = 1)
  return(print(plot))
}

#' Plot percent.mt
#' @get /qcplot3
#' @serializer png list(width = 500, height =500)
function(){
  Idents(obj) <- obj@meta.data$empty_ident
  plot <- VlnPlot(obj, features = c("percent.mt"), ncol = 1)
  return(print(plot))
}

#' Plot percent.ribo
#' @get /qcplot4
#' @serializer png list(width = 500, height =500)
function(){
  Idents(obj) <- obj@meta.data$empty_ident
  plot <- VlnPlot(obj, features = c("percent.ribo"), ncol = 1)
  return(print(plot))
}

#' Plot Variable genes
#' @get /var-genes-plot
#' @serializer png list(width = 1200, height =500)
function(){
  top10 <- head(VariableFeatures(obj), 10)
  plot1 <- VariableFeaturePlot(obj)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot <- CombinePlots(plots = list(plot1, plot2))
  plot <- plot1+plot2
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
  n_genes_per_cell <- obj$nFeature_RNA
  n_reads_per_cell <- obj$nCount_RNA 
  n_cells_per_gene <- apply(this_obj, 1, function(x) sum(x>0))
  pct_ribo_per_gene <- obj$percent.ribo
  pct_mito_per_gene <- obj$percent.mt
  gene_result <-
    data.frame(
      gene = rownames(this_obj),
      mean = rowMeans(this_obj),
      std = row_sd,
      min = row_min,
      max = row_max,
      n_cells_per_gene = n_cells_per_gene
    )
  cell_result <- 
    data.frame(
      n_reads_per_cell = n_reads_per_cell,
      n_genes_per_cell = n_genes_per_cell,
      pct_ribo_per_gene = pct_ribo_per_gene,
      pct_mito_per_gene = pct_mito_per_gene
    )
  hist_genes_per_cell <- data.frame(
    breaks=hist(n_genes_per_cell)$breaks[-1],
    counts=hist(n_genes_per_cell)$counts
  )
  hist_reads_per_cell <- data.frame(
    breaks=hist(n_reads_per_gene)$breaks[-1],
    counts=hist(n_reads_per_gene)$counts
  )
  hist_cells_per_gene <- data.frame(
    breaks=hist(n_cells_per_gene)$breaks[-1],
    counts=hist(n_cells_per_gene)$counts
  )
  return(
    list(
      gene_result,
      cell_result,
      hist_genes_per_cell,
      hist_reads_per_cell,
      hist_cells_per_gene
    )
  )
}

#' Get all metadata in a list
#' @get /meta-data
function(){
  
  n_count_rna <- obj@meta.data$nCount_RNA
  n_feature_rna <- obj@meta.data$nFeature_RNA
  pct_mito <- obj@meta.data$percent.mt
  pct_ribo <- obj@meta.data$percent.ribo
  
  meta <- read.csv('Zeisel_index_label.csv')
  obj <<- AddMetaData(obj, meta$Label, col.name = "cell_type")
  obj <<- AddMetaData(obj, meta$Sex, col.name = "sex")
  obj <<- AddMetaData(obj, meta$Sample, col.name = "sample")
  
  meta1_title <- "Cell type"
  meta1_name <- names(table(meta$Label))
  meta1_val <- as.vector(table(meta$Label))
  
  meta2_title <- "Sex"
  meta2_name <- names(table(meta$Sex))
  meta2_val <- as.vector(table(meta$Sex))
  
  meta3_title <- "Sample"
  meta3_name <- names(table(meta$Sample))
  meta3_val <- as.vector(table(meta$Sample))
  
  result <- list(meta1_title = meta1_title,
                 meta1_name=meta1_name, 
                 meta1_val=data.frame(name=meta1_name, value=meta1_val), 
                 meta2_title = meta2_title,
                 meta2_name=meta2_name, 
                 meta2_val=data.frame(name=meta2_name, value=meta2_val), 
                 meta3_title = meta3_title,
                 meta3_name=meta3_name, 
                 meta3_val=data.frame(name=meta3_name, value=meta3_val), 
                 n_count_rna = n_count_rna,
                 n_feature_rna = n_feature_rna,
                 pct_mito = pct_mito,
                 pct_ribo = pct_ribo
                 )
  message("meta-data send")
  return(result)
}



#' Plot umap
#' @get /umap-cluster
#' @serializer png list(width = 700, height = 600)
function(){
  print(ident_idx)
  #Idents(obj) <- obj$seurat_clusters
  Idents(obj) <- obj@meta.data[,ident_idx]
  this_text <- colnames(obj@meta.data)[ident_idx]
  #plot <- DimPlot(obj,reduction = "umap")
  plot <- Plot.cluster2D(obj, txt = this_text, pt_size = 0.5)
  return(print(plot))
}

#' Plot cell type metadata pie chart
#' @get /pie-meta
#' @serializer png list(width = 600, height = 600)
function(){
  meta <- read.csv('Zeisel_index_label.csv')
  obj <<- AddMetaData(obj, meta$Label, col.name = "cell_type")
  obj <<- AddMetaData(obj, meta$Sex, col.name = "sex")
  obj <<- AddMetaData(obj, meta$Sample, col.name = "sample")
  print(ident_idx)
  Idents(obj) <- obj@meta.data[,7]

  #Idents(obj) <- obj@meta.data[,ident_idx]
  mytable <- table(Idents(obj))
  
  x <- as.vector(mytable)
  mypct <- percent(x/sum(x))
  labels <- paste(names(mytable), "\n", mytable,"(",mypct,")", sep="")
  plot <- pie(x, labels, main = "Cell type metadata", col = rainbow(length(x)))
  return(print(plot))
}

#' Calculate DEG for two selections
#' @post /deg
function(req, ident1=4, ident2=5, min_pct=0.2, min_lfc=0.5){
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
    n_seurat_clusters = length(levels(obj$seurat_clusters)),
    umap_pts = data.frame(umap1=as.vector(Embeddings(obj, reduction = "umap")[,1]), umap2=as.vector(Embeddings(obj, reduction = "umap")[,]))
  ))
}

#' Get regulon results
#' @post /combine-regulon
function(req, jobid){
  print('return combine regulons')
  #jobid=2020041684528
  combine_regulon=read.table(url(paste0("https://bmbl.bmi.osumc.edu/iris3/data/",jobid,"/",jobid,"_combine_regulon.txt")), header = T)
  print(head(combine_regulon))
  return(list(
    list(combine_regulon)
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
function(req, gene="Gad1", split="sex"){
  #ident_idx=9
  if (split == "NULL") {
    Idents(obj) <- obj@meta.data[,ident_idx]
    plot <- VlnPlot(obj, gene, group.by = colnames(obj@meta.data)[ident_idx])
  } else{
    idx <- which(colnames(obj@meta.data) == split)
    Idents(obj) <- obj@meta.data[,ident_idx]
    plot <- VlnPlot(obj, gene, split.by = colnames(obj@meta.data)[idx], group.by = colnames(obj@meta.data)[ident_idx])
  }
  
  return(print(plot))
}


#' Plot gene expression in feature plot
#' @post /feature-gene
#' @serializer png list(width = 600, height =600)
function(req, gene="Gad1"){
  Idents(obj) <- obj@meta.data[,ident_idx]
  plot <- FeaturePlot(obj, gene)
  return(print(plot))
}


#' Plot dot plot
#' @post /dot-plot
#' @serializer png list(width = 1000, height =750)
function(req, top="3"){
  #ident_idx=9
  #cluster_markers <- FindAllMarkers(obj,logfc.threshold = 0.7)
  #write.csv(cluster_markers,"cluster_markers.csv")
  #top = 3
  top = as.numeric(top)
  m1 <- read_csv("cluster_markers.csv")%>%
    group_by(cluster) %>%
    top_n(top) %>%
    pull(gene) %>%
    unique()
    
  #idx <- which(colnames(obj@meta.data) == split)
  #Idents(obj) <- obj@meta.data[,ident_idx]
  plot <- DotPlot(obj, features = m1) + RotatedAxis()

  return(print(plot))
}


#' Plot dot plot
#' @post /annotate-cell-type
function(req){
  #ident_idx=9
  #cluster_markers <- FindAllMarkers(obj,logfc.threshold = 0.7)
  #write.csv(cluster_markers,"cluster_markers.csv")
  print('run cell type annotation')
  meta <- read.csv('Zeisel_index_label.csv')
  obj <<- AddMetaData(obj, meta$Label, col.name = "annotate_cell_type")
  Sys.sleep(5)
  return(list(
    n_annotate_cell_type = length(levels(as.factor(obj$annotate_cell_type))),
    annotate_cell_type=data.frame(annotate_cell_type = levels(as.factor(obj$annotate_cell_type)))
    
  ))
  
}

#' @post /transfer-cell-type
function(req){
  #ident_idx=9
  #cluster_markers <- FindAllMarkers(obj,logfc.threshold = 0.7)
  #write.csv(cluster_markers,"cluster_markers.csv")
  print('run cell type transfer')
  Sys.sleep(5)
  meta <- read.table('query_example_cell_label.txt',sep = "\t",header = T)
  meta <- meta[match(colnames(obj), meta$cell_name),]
  obj <<- AddMetaData(obj, meta$label, col.name = "transfer_cell_type")
  
  return(list(
    n_transfer_cell_type = length(levels(as.factor(obj$transfer_cell_type))),
    transfer_cell_type=data.frame(transfer_cell_type = levels(as.factor(obj$transfer_cell_type)))
    
  ))
  
}


#' @post /run-v1
function(req, jobid="1612737955509", k_arg=20, f_arg=0.7, o_arg=500,promoter_arg=1000){
  #ident_idx=9
  #cluster_markers <- FindAllMarkers(obj,logfc.threshold = 0.7)
  #write.csv(cluster_markers,"cluster_markers.csv")
  
  url_data <- httr::GET(paste0("https://bmbl.bmi.osumc.edu/iris3/data/",jobid))
  url_data <- httr::content(url_data)
  is_job_exist <- str_detect(url_data,"location.href")
  
  if(!is_job_exist) {
    current_wd <- getwd()
    wd=paste0(workdir,jobid)
    dir.create(wd, showWarnings = F)
    setwd(wd)
    
    #wd=paste0("/var/www/html/iris3/data/",jobid)
    
    cat(paste0("#!/bin/bash\necho 'test run ",jobid,"'\nsleep 10"), file="qsub2.sh",sep = "\n") 
    label_use_predict <- "2"
    
    
    cat(paste0("#!/bin/bash"), file="qsub.sh",sep = "\n") 
    cat(paste0("wd=",wd), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("exp_file=Zeisel_expression.csv"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("label_file=Zeisel_index_label.csv"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("gene_module_file="), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("jobid=",jobid), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("motif_min_length=12"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("motif_max_length=12"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("perl /var/www/html/iris3/program/prepare_email1.pl $jobid"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("Rscript /var/www/html/iris3/program/genefilter.R $jobid $wd$exp_file , $label_file , No ",resolution," ",nPCs," 5000 ",label_use_predict," No"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("echo gene_module_detection > running_status.txt"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("/var/www/html/iris3/program/qubic2/qubic -i $wd$jobid\\_filtered_expression.txt -q 0.06 -c 1.0 -k ",k_arg,"-o ",o_arg," -f ",f_arg), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("for file in *blocks"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("do"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("grep Conds $file |cut -d ':' -f2 >'$(basename $jobid_blocks.conds.txt)'"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("done"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("for file in *blocks"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("do"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("grep Genes $file |cut -d ':' -f2 >'$(basename $jobid_blocks.gene.txt)'"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("done"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("Rscript /var/www/html/iris3/program/ari_score.R $label_file $jobid , ",label_use_predict," "), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("echo gene_module_assignment > running_status.txt"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("Rscript /var/www/html/iris3/program/cts_gene_list.R $wd $jobid 1000   "), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("echo motif_finding_and_comparison > running_status.txt"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("/var/www/html/iris3/program/get_motif.sh $wd $motif_min_length $motif_max_length 1"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("Rscript /var/www/html/iris3/program/convert_meme.R $wd $motif_min_length"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("/var/www/html/iris3/program/get_motif.sh $wd $motif_min_length $motif_max_length 0"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("wait"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("cd $wd"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("Rscript /var/www/html/iris3/program/prepare_bbc.R $jobid $motif_min_length"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("mkdir tomtom"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("mkdir logo_tmp"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("mkdir logo"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("mkdir regulon_id"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("/var/www/html/iris3/program/get_logo.sh $wd"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("/var/www/html/iris3/program/get_tomtom.sh $wd"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("echo active_regulon_determination > running_status.txt"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("Rscript /var/www/html/iris3/program/merge_tomtom.R $wd $jobid $motif_min_length"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("echo regulon_inference > running_status.txt"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("Rscript /var/www/html/iris3/program/sort_regulon.R $wd $jobid"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("/var/www/html/iris3/program/get_atac_overlap.sh $wd"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("Rscript /var/www/html/iris3/program/prepare_heatmap.R $wd $jobid ",label_use_predict," "), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("Rscript /var/www/html/iris3/program/get_alternative_regulon.R $jobid"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("Rscript /var/www/html/iris3/program/generate_rss_scatter.R $jobid"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("Rscript /var/www/html/iris3/program/process_tomtom_result.R $jobid"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("mkdir json"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("/var/www/html/iris3/program/build_clustergrammar.sh $wd $jobid ",label_use_predict," "), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("zip -R $wd$jobid '*.regulon_gene_id.txt' '*.regulon_gene_symbol.txt' '*.regulon_rank.txt' '*_silh.txt' '*umap_embeddings.txt' '*.regulon_activity_score.txt' '*_cell_label.txt' '*.blocks' '*_blocks.conds.txt' '*_blocks.gene.txt' '*_filtered_expression.txt' '*_gene_id_name.txt' '*_marker_genes.txt' 'cell_cluster_unique_diffrenetially_expressed_genes.txt' '*_combine_regulon.txt'"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("perl /var/www/html/iris3/program/prepare_email.pl $jobid"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("echo 'finish'> done"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("chmod -R 777 ."), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("rm $wd$jobid\\_filtered_expression.txt"), file="qsub.sh",append = T,sep = "\n") 
    cat(paste0("rm $wd$jobid\\_filtered_expression.txt.chars"), file="qsub.sh",append = T,sep = "\n") 
    
    
    
    ### info.txt
    
    cat(paste0("is_load_exp,0"), file="info.txt",append = F,sep = "\n") 
    cat(paste0("k_arg,",k_arg), file="info.txt",append = T,sep = "\n") 
    cat(paste0("f_arg,",f_arg), file="info.txt",append = T,sep = "\n") 
    cat(paste0("o_arg,",o_arg), file="info.txt",append = T,sep = "\n") 
    cat(paste0("resolution_seurat,",resolution), file="info.txt",append = T,sep = "\n") 
    cat(paste0("n_variable_features,",nVariableFeatures), file="info.txt",append = T,sep = "\n") 
    cat(paste0("n_pca,",nPCs), file="info.txt",append = T,sep = "\n") 
    cat(paste0("label_use_predict,",label_use_predict), file="info.txt",append = T,sep = "\n") 
    cat(paste0("expfile,","Zeisel_expression.csv"), file="info.txt",append = T,sep = "\n") 
    cat(paste0("labelfile,","Zeisel_index_label.csv"), file="info.txt",append = T,sep = "\n") 
    cat(paste0("is_c,"), file="info.txt",append = T,sep = "\n") 
    cat(paste0("promoter_arg,",promoter_arg), file="info.txt",append = T,sep = "\n") 
    cat(paste0("bic_inference,2"), file="info.txt",append = T,sep = "\n") 
    cat(paste0("gene_module_file,"), file="info.txt",append = T,sep = "\n") 
    cat(paste0("is_imputation,No"), file="info.txt",append = T,sep = "\n") 
    
    
    message('Start IRIS3 v1')
    
    
    
    system(paste0("cp -r ",wd," /var/www/html/iris3/data"))
    system(paste0("cp -R /var/www/nodejs/iris3-backend/plumber/template/mouse/*"," /var/www/html/iris3/data/",jobid))
    system("nohup sh qsub2.sh > output.txt &")
    system(paste0("chmod -R 777 /var/www/html/iris3/data/",jobid))
    setwd(current_wd)
    return(jobid)
  } else {
    return (F)
  }
 
  
}

