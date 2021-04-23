import { Injectable } from '@nestjs/common'
import { DockerService } from 'src/docker/docker.service'

@Injectable()
export class CommandService {
  constructor(private dockerservice: DockerService) {}

  async runPreprocessing() {
    return this.dockerservice.test()
  }
}
