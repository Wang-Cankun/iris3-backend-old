import {
  InjectQueue,
  OnGlobalQueueCompleted,
  OnQueueActive,
  OnQueueCompleted,
  Process,
  Processor
} from '@nestjs/bull'
import { Logger } from '@nestjs/common'
import { Job, JobCounts, Queue } from 'bull'

@Processor('file')
export class FileProcessor {
  private readonly logger = new Logger(FileProcessor.name)
  constructor(
    @InjectQueue('file')
    private readonly fileQueue: Queue
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
    const job = await this.fileQueue.getJob(jobId)
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
    await 1
    return 1
  }
}
