import {
  Body,
  Controller,
  Get,
  Param,
  Post,
  UploadedFiles,
  UseInterceptors
} from '@nestjs/common'
import {
  AnyFilesInterceptor,
  FileFieldsInterceptor
} from '@nestjs/platform-express'
import { Queue } from 'bull'
import { InjectQueue } from '@nestjs/bull'
import { FileService } from './file.service'
import { ConfigService } from '@nestjs/config'
import { File } from './entities/file.entity'

@Controller('file')
export class FileController {
  constructor(
    private configService: ConfigService,

    @InjectQueue('file')
    private readonly fileQueue: Queue,

    private readonly fileService: FileService
  ) {}

  @Get(':id')
  async getJob(@Param('id') id: string) {
    return await this.fileQueue.getJob(id)
  }

  @Get('upload/:id')
  async findOne(@Param('id') id: string): Promise<File[]> {
    return await this.fileService.findAll(id)
  }

  @Post('upload')
  @UseInterceptors(AnyFilesInterceptor())
  async uploadFile(@UploadedFiles() files, @Body() body) {
    const uploadDataInfo = files.map((file) => ({
      ...file,
      jobid: body.jobid,
      index: body.index,
      species: body.species
    }))
    console.log(uploadDataInfo)
    if (body.jobid !== 'example') {
      for (const f of uploadDataInfo) {
        await this.fileService.create(f)
      }
    }

    return uploadDataInfo
  }
}
