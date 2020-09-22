import { Injectable } from '@nestjs/common'
import { DockerService } from 'src/docker/docker.service'

@Injectable()
export class CommandService {
  constructor(private dockerservice: DockerService) {}

  async runPreprocessing(option) {
    const cmd = `library(Seurat)\nlibrary(ggplot2)\n
    obj <- CreateSeuratObject(read.csv('Yan_2013_expression.csv',row.names = 1))\n
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")\n
    ggsave("qc_plot.png",VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))\n
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)\n
  obj <- ScaleData(obj, features = rownames(obj))\n
  obj <- RunPCA(obj, features = VariableFeatures(object = obj))\n
  obj <- FindNeighbors(obj, dims = 1:${option.nPCs})\n
  obj <- FindClusters(obj, resolution = ${option.nResolution})\n
  obj <- RunUMAP(obj, dims = 1:${option.nPCs})\n
  ggsave("umap_cluster_plot.png",DimPlot(obj, reduction = "umap"))\n
    `
    return this.dockerservice.runCmd(cmd)
  }
}
