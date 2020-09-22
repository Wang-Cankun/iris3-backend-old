import { Module } from '@nestjs/common'
import { CommandService } from './command.service'
import { CommandController } from './command.controller'
import { DockerService } from 'src/docker/docker.service'

@Module({
  providers: [CommandService, DockerService],
  controllers: [CommandController]
})
export class CommandModule {}
