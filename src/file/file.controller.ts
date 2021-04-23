import {
  Body,
  Controller,
  Post,
  UploadedFiles,
  UseInterceptors
} from '@nestjs/common'
import { FileFieldsInterceptor } from '@nestjs/platform-express'
import { Queue } from 'bull'
import { InjectQueue } from '@nestjs/bull'
import { FileService } from './file.service'

@Controller('file')
export class FileController {
  constructor(
    @InjectQueue('file')
    private readonly fileQueue: Queue,

    private readonly fileService: FileService
  ) {}

  @Post('upload')
  @UseInterceptors(
    FileFieldsInterceptor(
      [
        { name: 'expFile', maxCount: 3 },
        { name: 'labelFile', maxCount: 1 }
      ],
      { dest: 'tmp' }
    )
  )
  async uploadFile(@UploadedFiles() files, @Body() body) {
    console.log(files)
    const jobInfo = await this.fileQueue.add('upload', {
      file: files,
      body: body
    })
    //return { file, body, jobid }
    return jobInfo
  }
}
