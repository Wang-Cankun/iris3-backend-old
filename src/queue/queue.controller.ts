import { InjectQueue, OnGlobalQueueCompleted } from '@nestjs/bull'
import {
  Body,
  Controller,
  Get,
  HttpService,
  Param,
  Post,
  UseInterceptors,
  UploadedFiles,
  Req
} from '@nestjs/common'
import {
  AnyFilesInterceptor,
  FileFieldsInterceptor,
  FileInterceptor
} from '@nestjs/platform-express'
import { Queue } from 'bull'

@Controller('queue')
export class QueueController {
  constructor(
    @InjectQueue('task')
    private readonly jobQueue: Queue,
    private httpService: HttpService
  ) {}

  @Post('transcode')
  async transcode() {
    console.log('Add transcode event')
    await this.jobQueue.add('transcode', {
      file: 'audio.mp3'
    })
  }

  @Get(':id')
  async getJob(@Param('id') id: string) {
    return await this.jobQueue.getJob(id)
  }

  @Post('test')
  async index() {
    console.log('Add test event')
    this.jobQueue.add('test', { data: 1, somedata: 2 })
    return 1
  }

  @Post('upload')
  @UseInterceptors(
    FileFieldsInterceptor(
      [
        { name: 'expFile', maxCount: 3 },
        { name: 'labelFile', maxCount: 1 }
      ],
      { dest: './tmp' }
    )
  )
  async uploadFile(@UploadedFiles() files, @Body() body) {
    const jobInfo = await this.jobQueue.add('upload', {
      file: files,
      body: body
    })
    //return { file, body, jobid }
    return jobInfo
  }

  @Post('load')
  async loadExpression(@Body() body) {
    const jobInfo = this.jobQueue.add('load', body)
    return jobInfo
  }

  @Post('load-multi-rna')
  async loadMultiRna(@Body() body) {
    const jobInfo = this.jobQueue.add('load-multi-rna', body)
    return jobInfo
  }

  @Post('load-multiome')
  async loadMultiome(@Body() body) {
    const jobInfo = this.jobQueue.add('load-multiome', body)
    return jobInfo
  }

  @Post('cluster')
  async cluster(@Body() body) {
    const jobInfo = this.jobQueue.add('cluster', body)
    return jobInfo
  }

  @Post('cluster-multiome')
  async clusterMultiome(@Body() body) {
    const jobInfo = this.jobQueue.add('cluster-multiome', body)
    return jobInfo
  }

  @Post('annotate-cell-type')
  async annotateCellType(@Body() body) {
    const jobInfo = this.jobQueue.add('annotate-cell-type', body)
    return jobInfo
  }
  @Post('transfer-cell-type')
  async transferCellType(@Body() body) {
    const jobInfo = this.jobQueue.add('transfer-cell-type', body)
    return jobInfo
  }

  @Post('run-v1')
  async runV1(@Body() body) {
    console.log(body)
    const jobInfo = this.jobQueue.add('run-v1', body)
    return jobInfo
  }

  @Post('genes')
  async genes() {
    const result = await this.httpService
      .get('http://localhost:8000/genes')
      .toPromise()
      .then((response) => response.data)
    return result
  }

  @Post('idents')
  async getIdents() {
    const result = await this.httpService
      .get('http://localhost:8000/idents')
      .toPromise()
      .then((response) => response.data)

    return result
  }

  @Post('set-idents')
  async setIdents(@Body() body) {
    const jobInfo = this.jobQueue.add('set-ident', body)
    return jobInfo
  }

  @Post('assays')
  async getAssay(@Body() body) {
    const result = this.jobQueue.add('get-assay', body)
    return result
  }

  @Post('set-assay')
  async setAssay(@Body() body) {
    const jobInfo = this.jobQueue.add('set-assay', body)
    return jobInfo
  }

  @Post('embeddings')
  async getEmbedding(@Body() body) {
    const result = this.jobQueue.add('get-embedding', body)
    return result
  }

  @Post('set-embedding')
  async setEmbedding(@Body() body) {
    const jobInfo = this.jobQueue.add('set-embedding', body)
    return jobInfo
  }

  @Post('merge-idents')
  async mergeIdents(@Body() body) {
    const jobInfo = this.jobQueue.add('merge-idents', body)
    return jobInfo
  }

  @Post('rename-idents')
  async renameIdents(@Body() body) {
    const jobInfo = this.jobQueue.add('rename-idents', body)
    return jobInfo
  }

  @Post('select-category')
  async selectCategory(@Body() body) {
    const jobInfo = this.jobQueue.add('select-category', body)
    return jobInfo
  }

  @Post('select-cells')
  async selectCells(@Body() body) {
    const jobInfo = this.jobQueue.add('select-cells', body)
    return jobInfo
  }

  @Post('subset-cells')
  async subsetCells(@Body() body) {
    const jobInfo = this.jobQueue.add('subset-cells', body)
    return jobInfo
  }

  @Post('set-obj')
  async setObj(@Body() body) {
    const jobInfo = this.jobQueue.add('set-obj', body)
    return jobInfo
  }

  @Post('qcplot')
  async qcPlot() {
    const jobInfo = await this.jobQueue.add('qcplot', { jobid: 1 })
    return jobInfo
  }

  @Post('qcplot1')
  async qcPlot1() {
    const jobInfo = await this.jobQueue.add('qcplot1', { jobid: 1 })
    return jobInfo
  }
  @Post('qcplot2')
  async qcPlot2() {
    const jobInfo = await this.jobQueue.add('qcplot2', { jobid: 1 })
    return jobInfo
  }
  @Post('qcplot3')
  async qcPlot3() {
    const jobInfo = await this.jobQueue.add('qcplot3', { jobid: 1 })
    return jobInfo
  }
  @Post('qcplot4')
  async qcPlot4() {
    const jobInfo = await this.jobQueue.add('qcplot4', { jobid: 1 })
    return jobInfo
  }
  @Post('qcplot5')
  async qcPlot5() {
    const jobInfo = await this.jobQueue.add('qcplot5', { jobid: 1 })
    return jobInfo
  }
  @Post('qcplot6')
  async qcPlot6() {
    const jobInfo = await this.jobQueue.add('qcplot6', { jobid: 1 })
    return jobInfo
  }
  @Post('var-genes-plot')
  async varGenesPlot() {
    const jobInfo = await this.jobQueue.add('varGenesPlot', { jobid: 1 })
    return jobInfo
  }

  @Post('var-genes-list')
  async varGenesList() {
    const jobInfo = await this.jobQueue.add('varGenesList', { jobid: 1 })
    return jobInfo
  }

  @Post('atac-qc-list')
  async atacQcList() {
    const jobInfo = await this.jobQueue.add('atacQcList', { jobid: 1 })
    return jobInfo
  }

  @Post('meta-data')
  async metaData() {
    const jobInfo = await this.jobQueue.add('metaData', { jobid: 1 })
    return jobInfo
  }

  @Post('combine-regulon')
  async combineRegulon(@Body() body) {
    console.log(body)
    const jobInfo = await this.jobQueue.add('combineRegulon', body)
    return jobInfo
  }

  @Post('deg')
  async deg(@Body() body) {
    const jobInfo = this.jobQueue.add('deg', body)
    return jobInfo
  }

  @Post('gsea-table')
  async gseaTable(@Body() body) {
    const jobInfo = this.jobQueue.add('gsea-table', body)
    return jobInfo
  }

  @Post('gsea-plot')
  async gseaPlot(@Body() body) {
    const jobInfo = this.jobQueue.add('gsea-plot', body)
    return jobInfo
  }

  @Post('gene-correlation-plot')
  async geneCorPlot(@Body() body) {
    const jobInfo = this.jobQueue.add('gene-correlation-plot', body)
    return jobInfo
  }

  @Post('coverage-plot')
  async coveragePlot(@Body() body) {
    const jobInfo = this.jobQueue.add('coverage-plot', body)
    return jobInfo
  }

  @Post('pie-meta')
  async pieMetaPlot() {
    console.log('pie meta plot')
    const jobInfo = await this.jobQueue.add('pieMetaPlot', { jobid: 1 })
    return jobInfo
  }
  @Post('umap-cluster')
  async umapClusterPlot() {
    console.log('umap plot')
    const jobInfo = await this.jobQueue.add('umapClusterPlot', { jobid: 1 })
    return jobInfo
  }

  @Post('umap-static')
  async umapStatic(@Body() body) {
    console.log('umap static')
    const jobInfo = await this.jobQueue.add('umapStatic', body)
    return jobInfo
  }

  @Post('umap-rna')
  async umapRNAPlot() {
    console.log('umap plot')
    const jobInfo = await this.jobQueue.add('umapRNAPlot', { jobid: 1 })
    return jobInfo
  }

  @Post('umap-atac')
  async umapATACPlot() {
    console.log('umap plot')
    const jobInfo = await this.jobQueue.add('umapATACPlot', { jobid: 1 })
    return jobInfo
  }

  @Post('umap-atac')
  async umapIntegratedPlot() {
    console.log('umap plot')
    const jobInfo = await this.jobQueue.add('umapIntegratedPlot', { jobid: 1 })
    return jobInfo
  }

  @Post('umap-gene')
  async umapGenePlot(@Body() body) {
    console.log(body)
    const jobInfo = await this.jobQueue.add('umapGenePlot', body)
    return jobInfo
  }

  @Post('violin-gene')
  async violinGenePlot(@Body() body) {
    const jobInfo = await this.jobQueue.add('violinGenePlot', body)
    return jobInfo
  }

  @Post('feature-gene')
  async featureGenePlot(@Body() body) {
    const jobInfo = await this.jobQueue.add('featureGenePlot', body)
    return jobInfo
  }

  @Post('dot-plot')
  async dotPlot(@Body() body) {
    const jobInfo = await this.jobQueue.add('dotPlot', body)
    return jobInfo
  }

  @Post('bicluster')
  async runBicluster(@Body() body) {
    const jobInfo = await this.jobQueue.add('bicluster', body)
    return jobInfo
  }

  @Post('motif')
  async runMotif(@Body() body) {
    const jobInfo = await this.jobQueue.add('motif', body)
    return jobInfo
  }

  @Post('regulon')
  async runRegulon(@Body() body) {
    const jobInfo = await this.jobQueue.add('regulon', body)
    return jobInfo
  }

  @Post('save')
  async saveSession(@Body() body) {
    const jobInfo = await this.jobQueue.add('save', body)
    return jobInfo
  }
}
