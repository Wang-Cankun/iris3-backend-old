import {
  InjectQueue,
  OnGlobalQueueCompleted,
  OnQueueActive,
  OnQueueCompleted,
  Process,
  Processor
} from '@nestjs/bull'
import { HttpService, Logger } from '@nestjs/common'
import { Job, JobCounts, Queue } from 'bull'
import { Job as JobEntity } from 'src/job/entities/job.entity'
import { JobService } from 'src/job/job.service'
import { FileService } from '../file/file.service'
import { PlumberService } from '../plumber/plumber.service'

@Processor('task')
export class QueueProcessor {
  private readonly logger = new Logger(QueueProcessor.name)
  constructor(
    @InjectQueue('task')
    private readonly jobQueue: Queue,
    private readonly httpService: HttpService,
    private readonly jobService: JobService,
    private readonly plumberService: PlumberService,
    private readonly fileService: FileService
  ) {}

  @OnQueueActive()
  onActive(job: Job) {
    this.logger.log(
      `Processing job ${job.id} of type ${job.name} with data ${job.data}...`
    )
  }

  @OnQueueCompleted()
  onComplete(job: Job) {
    this.logger.log(
      `Completed job ${job.id} of type ${job.name} with data ${job.data}...`
    )
  }

  @OnGlobalQueueCompleted()
  async onGlobalCompleted(jobId: number, result: any) {
    const job = await this.jobQueue.getJob(jobId)
    // console.log('(Global) on completed: job ', job.id, ' -> result: ', result)
  }

  @Process('transcode')
  async transcode(job: Job<unknown>) {
    let progress = 0
    for (let i = 0; i < 100; i++) {
      await console.log(job.id)
      progress += 1
      job.progress(progress)
    }
    return {}
  }

  @Process('test')
  handleTest(job: Job) {
    this.logger.log('Start transcoding...')
    this.logger.log(job.data)
    this.logger.log('Transcoding completed')
  }

  @Process('load')
  async loadExpression(job: Job) {
    const uploadFiles = await this.fileService.findAll(job.data.jobid)

    const rnaFile = uploadFiles.filter((file) => file.fieldname === 'singleRna')
    const labelFile = uploadFiles.filter(
      (file) => file.fieldname === 'labelFile-singleRna'
    )
    const loadPayload = { ...job.data, expr: rnaFile, label: labelFile }
    console.log(loadPayload)
    const result = await this.plumberService.runCommand('load', loadPayload)
    return result
  }

  @Process('load-multi-rna')
  async loadMultiRna(job: Job) {
    const result = await this.plumberService.runCommand(
      'load-multi-rna',
      job.data
    )
    return result
  }

  @Process('load-multiome')
  async loadMultiome(job: Job) {
    const result = await this.plumberService.runCommand(
      'load-multiome',
      job.data
    )
    return result
  }

  @Process('get-ident')
  async getIdent(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/idents', job.data)
      .toPromise()
      .then((response) => response.data)
    return result
  }

  @Process('set-ident')
  async setIdent(job: Job) {
    const result = await this.httpService
      .post('http://localhost:8000/idents', job.data)
      .toPromise()
      .then((response) => response.data)
    return result
  }

  @Process('get-assay')
  async getAssay(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/assays', job.data)
      .toPromise()
      .then((response) => response.data)
    return result
  }

  @Process('set-assay')
  async setAssay(job: Job) {
    return await this.plumberService.runCommand('set-assay', job.data)
  }

  @Process('get-embedding')
  async getEmbedding(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/embeddings', job.data)
      .toPromise()
      .then((response) => response.data)
    return result
  }

  @Process('set-embedding')
  async setEmbedding(job: Job) {
    return await this.plumberService.runCommand('set-embedding', job.data)
  }

  @Process('cluster')
  async cluster(job: Job) {
    const result = await this.plumberService.runCommand('cluster', job.data)
    return result
  }

  @Process('cluster-multiome')
  async clusterMultiome(job: Job) {
    const result = await this.plumberService.runCommand(
      'cluster-multiome',
      job.data
    )
    return result
  }

  @Process('merge-idents')
  async mergeIdents(job: Job) {
    const result = await this.plumberService.runCommand(
      'merge-idents',
      job.data
    )
    return result
  }

  @Process('rename-idents')
  async renameIdents(job: Job) {
    const result = await this.plumberService.runCommand(
      'rename-idents',
      job.data
    )
    return result
  }

  @Process('select-category')
  async selectCategory(job: Job) {
    const result = await this.plumberService.runCommand(
      'select-category',
      job.data
    )
    return result
  }

  @Process('select-cells')
  async selectCells(job: Job) {
    const result = await this.httpService
      .post('http://localhost:8000/select-cells', job.data)
      .toPromise()
      .then((res) => res.data)
    return result
  }
  @Process('subset-cells')
  async subsetCells(job: Job) {
    const result = await this.httpService
      .post('http://localhost:8000/subset-cells', job.data)
      .toPromise()
      .then((res) => res.data)
    return result
  }
  @Process('set-obj')
  async setObj(job: Job) {
    const result = await this.httpService
      .post('http://localhost:8000/set-obj', job.data)
      .toPromise()
      .then((res) => res.data)
    return result
  }
  @Process('transfer-cell-type')
  async transferCellType(job: Job) {
    const result = await this.httpService
      .post('http://localhost:8000/transfer-cell-type', job.data)
      .toPromise()
      .then((response) => response.data)
    return result
  }

  @Process('run-v1')
  async runV1(job: Job) {
    console.log(job.data)
    const result = await this.httpService
      .post('http://localhost:8000/run-v1', job.data)
      .toPromise()
      .then((response) => response.data)
    return result
  }

  @Process('annotate-cell-type')
  async annotateCellType(job: Job) {
    const result = await this.httpService
      .post('http://localhost:8000/annotate-cell-type', job.data)
      .toPromise()
      .then((response) => response.data)
    return result
  }

  @Process('varGenesPlot')
  async varsGenesPlot(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/var-genes-plot', {
        responseType: 'arraybuffer'
      })
      .toPromise()
      .then(
        (response) =>
          'data:image/png;base64,' +
          Buffer.from(response.data, 'binary').toString('base64')
      )
    return result
  }

  @Process('pieMetaPlot')
  async pieMetaPlot(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/pie-meta', {
        responseType: 'arraybuffer'
      })
      .toPromise()
      .then(
        (response) =>
          'data:image/png;base64,' +
          Buffer.from(response.data, 'binary').toString('base64')
      )
    return result
  }

  @Process('varGenesList')
  async varGenesList(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/var-genes-list')
      .toPromise()
      .then((res) => res.data)
    return result
  }

  @Process('atacQcList')
  async atacQcList(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/atac-qc-list')
      .toPromise()
      .then((res) => res.data)
    return result
  }

  @Process('metaData')
  async metaData(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/meta-data')
      .toPromise()
      .then((res) => res.data)
    return result
  }
  @Process('combineRegulon')
  async combineRegulon(job: Job) {
    console.log('running combine regulon')
    const result = await this.httpService
      .post('http://localhost:8000/combine-regulon', job.data)
      .toPromise()
      .then((res) => res.data)
    return result
  }

  @Process('deg')
  async deg(job: Job) {
    const result = await this.plumberService.runCommand('deg', job.data)
    return result
  }

  @Process('gsea-table')
  async gseaTable(job: Job) {
    const result = await this.plumberService.runCommand('gsea-table', job.data)
    return result
  }

  @Process('gsea-plot')
  async gseaPlot(job: Job) {
    const result = await this.plumberService.runStaticImage(
      'gsea-plot',
      job.data
    )
    return result
  }
  @Process('coverage-plot')
  async coveragePlot(job: Job) {
    const result = await this.plumberService.runStaticImage(
      'coverage-plot',
      job.data
    )
    return result
  }

  @Process('static-heatmap')
  async staticHeatmap(job: Job) {
    const result = await this.plumberService.runStaticImage(
      'static-heatmap',
      job.data
    )
    return result
  }

  @Process('gene-correlation-plot')
  async geneCorPlot(job: Job) {
    const result = await this.plumberService.runStaticImage(
      'gene-correlation-plot',
      job.data
    )
    return result
  }

  @Process('umapClusterPlot')
  async umapClusterPlot(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/umap-cluster', {
        responseType: 'arraybuffer'
      })
      .toPromise()
      .then((res) => res.data)
    return result
  }

  @Process('umapStatic')
  async umapStatic(job: Job) {
    const result = await this.httpService
      .post('http://localhost:8000/umap-static', job.data, {
        responseType: 'arraybuffer'
      })
      .toPromise()
      .then(
        (response) =>
          'data:image/png;base64,' +
          Buffer.from(response.data, 'binary').toString('base64')
      )
    return result
  }

  @Process('umapRNAPlot')
  async umapRNAPlot(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/umap-rna', {
        responseType: 'arraybuffer'
      })
      .toPromise()
      .then(
        (response) =>
          'data:image/png;base64,' +
          Buffer.from(response.data, 'binary').toString('base64')
      )
    return result
  }

  @Process('umapATACPlot')
  async umapATACPlot(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/umap-atac', {
        responseType: 'arraybuffer'
      })
      .toPromise()
      .then(
        (response) =>
          'data:image/png;base64,' +
          Buffer.from(response.data, 'binary').toString('base64')
      )
    return result
  }

  @Process('umapIntegratedPlot')
  async umapIntegratedPlot(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/umap-integrated', {
        responseType: 'arraybuffer'
      })
      .toPromise()
      .then(
        (response) =>
          'data:image/png;base64,' +
          Buffer.from(response.data, 'binary').toString('base64')
      )
    return result
  }

  @Process('umapGenePlot')
  async umapGenePlot(job: Job) {
    const result = await this.httpService
      .post('http://localhost:8000/gene-umap-static', job.data, {
        responseType: 'arraybuffer'
      })
      .toPromise()
      .then(
        (response) =>
          'data:image/png;base64,' +
          Buffer.from(response.data, 'binary').toString('base64')
      )
    return result
  }

  @Process('violinGenePlot')
  async violinGenePlot(job: Job) {
    const result = await this.httpService
      .post('http://localhost:8000/violin-gene', job.data, {
        responseType: 'arraybuffer'
      })
      .toPromise()
      .then(
        (response) =>
          'data:image/png;base64,' +
          Buffer.from(response.data, 'binary').toString('base64')
      )
    return result
  }

  @Process('featureGenePlot')
  async featureGenePlot(job: Job) {
    const result = await this.httpService
      .post('http://localhost:8000/gene-umap-static', job.data, {
        responseType: 'arraybuffer'
      })
      .toPromise()
      .then(
        (response) =>
          'data:image/png;base64,' +
          Buffer.from(response.data, 'binary').toString('base64')
      )
    return result
  }
  @Process('dotPlot')
  async dotPlot(job: Job) {
    const result = await this.httpService
      .post('http://localhost:8000/dot-plot', job.data, {
        responseType: 'arraybuffer'
      })
      .toPromise()
      .then(
        (response) =>
          'data:image/png;base64,' +
          Buffer.from(response.data, 'binary').toString('base64')
      )
    return result
  }
}
