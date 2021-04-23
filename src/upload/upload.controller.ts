import {
  Body,
  Controller,
  Get,
  Post,
  UploadedFile,
  UseInterceptors
} from '@nestjs/common'
import { FileInterceptor } from '@nestjs/platform-express'

@Controller('upload')
export class UploadController {
  @Get()
  findAll() {
    return 'test get'
  }
  @Post('data')
  @UseInterceptors(FileInterceptor('file'))
  uploadFile(@UploadedFile() file, @Body() body) {
    console.log(body)
    return body
  }
}
