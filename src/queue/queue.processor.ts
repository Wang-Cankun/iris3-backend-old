import {
  InjectQueue,
  OnGlobalQueueCompleted,
  OnQueueActive,
  OnQueueCompleted,
  Process,
  Processor
} from '@nestjs/bull'
import { HttpService, Logger } from '@nestjs/common'
import { InjectRepository } from '@nestjs/typeorm'
import { Job, JobCounts, Queue } from 'bull'
import { response } from 'express'
import moveFile from 'move-file'
import { Job as JobEntity } from 'src/job/entities/job.entity'
import { JobService } from 'src/job/job.service'
import { Repository } from 'typeorm'

@Processor('task')
export class QueueProcessor {
  private readonly logger = new Logger(QueueProcessor.name)
  constructor(
    @InjectQueue('task')
    private readonly jobQueue: Queue,
    private httpService: HttpService,
    private jobService: JobService
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

  @Process('upload')
  async submit(job: Job) {
    const dto = {
      ...job.data.body,
      expFile: job.data.file.expFile[0].filename,
      labelFile: job.data.file.labelFile[0].filename
    }
    console.log(dto)
    await this.jobService.create(dto)
    return 1
  }

  @Process('get-ident')
  async getIdent(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/ident', job.data)
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

  @Process('load')
  async loadExpression(job: Job) {
    const result = await this.httpService
      .post('http://localhost:8000/load', job.data)
      .toPromise()
      .then((response) => response.data)
    return result
  }

  @Process('cluster')
  async cluster(job: Job) {
    const result = await this.httpService
      .post('http://localhost:8000/cluster', job.data)
      .toPromise()
      .then((response) => response.data)
    return result
  }
  @Process('qcplot')
  async qcplot(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/qcplot', {
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
  @Process('qcplot1')
  async qcplot1(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/qcplot1', {
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

  @Process('qcplot2')
  async qcplot2(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/qcplot2', {
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
  @Process('qcplot3')
  async qcplot3(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/qcplot3', {
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
  @Process('qcplot4')
  async qcplot4(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/qcplot4', {
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

  @Process('varGenesList')
  async varGenesList(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/var-genes-list')
      .toPromise()
      .then((res) => res.data)
    return result
  }
  @Process('umapClusterPlot')
  async umapClusterPlot(job: Job) {
    const result = await this.httpService
      .get('http://localhost:8000/umap-cluster', {
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
      .post('http://localhost:8000/umap-gene', job.data, {
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
}
