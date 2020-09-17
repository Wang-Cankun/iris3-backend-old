import { Controller, Get, Param } from '@nestjs/common'
import { DockerService } from './docker.service'

@Controller('docker')
export class DockerController {
  constructor(private readonly dockerService: DockerService) {}

  @Get('test')
  public test(@Param() params): any {
    return this.dockerService.test('flykun0620@gmail.com')
  }
}
