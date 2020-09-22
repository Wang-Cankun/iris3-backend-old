import { InjectQueue } from '@nestjs/bull'
import { Controller, Get, Post } from '@nestjs/common'
import { Queue } from 'bull'

@Controller('queue')
export class QueueController {
  constructor(
    @InjectQueue('queue')
    private readonly jobQueue: Queue
  ) {}

  @Post('transcode')
  async transcode() {
    await this.jobQueue.add('transcode', {
      file: 'audio.mp3'
    })
  }

  @Get('/')
  index() {
    console.log('Add test event')
    this.jobQueue.add('testEvent', { data: 1, somedata: 2 })
  }
}
