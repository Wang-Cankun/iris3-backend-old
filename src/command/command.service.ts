import { Injectable } from '@nestjs/common'
import { DockerService } from 'src/docker/docker.service'

@Injectable()
export class CommandService {
  constructor(private dockerservice: DockerService) {}

  async runPreprocessing(option) {
    const cmd = `library(Seurat)\nlibrary(ggplot2)\n
    obj <- CreateSeuratObject(read.csv('Yan_2013_expression.csv',row.names = 1), min.cells = ${option.cellFilter},min.features = ${option.geneFilter})\n
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")\n
    ggsave("qc_plot.png",height = 8, width = 12,VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))\n
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = ${option.nVariableFeatures})\n
    top10 <- head(VariableFeatures(obj), 10)\n
    plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")\n
    plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")\n
    ggsave("feature_scatter.png",height = 8, width = 12, CombinePlots(plots = list(plot1, plot2)))+theme(aspect.ratio=2)\n
    plot1 <- VariableFeaturePlot(obj)\n
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)\n
    ggsave("top_variable_genes.png",height = 8, width = 12, CombinePlots(plots = list(plot1, plot2)))+theme(aspect.ratio=2)\n
    `

    /*
    const cmd2 = `
    
    obj <- ScaleData(obj, features = rownames(obj))\n
    obj <- RunPCA(obj, features = VariableFeatures(object = obj))\n
    obj <- FindNeighbors(obj, dims = 1:${option.nPCs})\n
    obj <- FindClusters(obj, resolution = ${option.nResolution})\n
    obj <- RunUMAP(obj, dims = 1:${option.nPCs})\n 
    ggsave("umap_cluster_plot.png",DimPlot(obj, reduction = "umap"))\n
    `*/
    return this.dockerservice.runCmd(cmd)
  }
}
