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

  @Post('cluster')
  async cluster(@Body() body) {
    console.log(body)
    const jobInfo = this.jobQueue.add('cluster', body)
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
  @Post('umap-cluster')
  async umapClusterPlot() {
    const jobInfo = await this.jobQueue.add('umapClusterPlot', { jobid: 1 })
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
