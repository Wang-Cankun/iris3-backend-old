import { Body, Controller, Get, Param, Post } from '@nestjs/common'
import { DockerService } from './docker.service'
import { ContainerName } from './dto/container-name.dto'
import Docker from 'dockerode'
import { RunContainerDto } from './dto/run-container.dto'
import { StopContainerDto } from './dto/stop-container.dto'
@Controller('docker')
export class DockerController {
  constructor(private readonly dockerService: DockerService) {}

  @Get('test')
  public test(@Param() params): any {
    console.log('test')
    return this.dockerService.test()
  }

  @Post('run')
  public async run(@Body() runContainerDto: RunContainerDto): Promise<Docker> {
    return await this.dockerService.runContainer(runContainerDto)
  }
  @Post('stop')
  public async stop(
    @Body() stopContainerDto: StopContainerDto
  ): Promise<Docker> {
    console.log(stopContainerDto)
    return await this.dockerService.stopContainer(stopContainerDto)
  }
}
