import {
  OnGlobalQueueCompleted,
  OnQueueActive,
  Process,
  Processor
} from '@nestjs/bull'
import { Logger } from '@nestjs/common'
import { Job } from 'bull'

@Processor('queue')
export class QueueProcessor {
  private readonly logger = new Logger(QueueProcessor.name)

  @OnQueueActive()
  onActive(job: Job) {
    console.log(
      `Processing job ${job.id} of type ${job.name} with data ${job.data}...`
    )
  }
  async transcode(job: Job<unknown>) {
    let progress = 0
    for (let i = 0; i < 100; i++) {
      await console.log(job.data)
      progress += 10
      job.progress(progress)
    }
    return {}
  }

  @Process()
  handleTest(job: Job) {
    console.log('Start transcoding...')
    this.logger.debug('Start transcoding...')
    this.logger.debug(job.data)
    this.logger.debug('Transcoding completed')
  }
}
