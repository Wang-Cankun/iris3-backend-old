import { Controller, Get } from '@nestjs/common'

@Controller('upload')
export class UploadController {
  @Get()
  findAll() {
    return 'test get'
  }
}
