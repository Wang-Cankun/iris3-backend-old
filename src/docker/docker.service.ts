import { Injectable } from '@nestjs/common'
import { compareSync } from 'bcrypt'
import * as Docker from 'dockerode'
import { ContainerName } from './dto/container-name.dto'
import { RunContainerDto } from './dto/run-container.dto'
import { StopContainerDto } from './dto/stop-container.dto'

@Injectable()
export class DockerService {
  private readonly dockerHost
  private readonly docker1
  private readonly CTRL_P = '\u0010\n'
  private readonly CTRL_Q = '\u0011\n'
  constructor() {
    this.dockerHost = new Docker()
    this.docker1 = new Docker({
      host: process.env.VM1_HOST,
      port: process.env.VM1_PORT
    })
  }

  /**
   * Exists the given stream and removes all listeners
   * @param stream The stream to exit
   */
  _exit(stream) {
    process.stdin.removeAllListeners()
    process.stdin.resume()
    stream.end()
  }

  /**
   * Ping the Docker API
   * @throws {DockerNotInstalledError} Gets thrown when could not ping Docker API
   */
  async ping() {
    try {
      return await this.docker1.ping()
    } catch (err) {
      throw new Error()
    }
  }

  async streamToString(stream) {
    const chunks = []
    return new Promise((resolve, reject) => {
      stream.on('data', (chunk) => chunks.push(chunk))
      stream.on('error', reject)
      stream.on('end', () => resolve(Buffer.concat(chunks).toString('utf8')))
    })
  }

  async test() {
    return await this.docker1.listContainers()
  }

  async listContainers(docker) {
    return await docker.listContainers()
  }

  async runContainer(runContainerDto: RunContainerDto): Promise<Docker> {
    const createOptions = {
      name: runContainerDto.port
    }
    const cmd = []
    const out = ''
    const result = this.docker1
      .run(runContainerDto.image, cmd, out, createOptions)
      .then((data) => {
        const output = data[0]
        const container = data[1]
        console.log(output)
        return container
      })
      .catch((err) => {
        console.log(err)
      })

    return result
  }
  async stopContainer(stopContainerDto: StopContainerDto): Promise<Docker> {
    const result = this.docker1
      .getContainer(stopContainerDto.id)
      .stop()
      .then((data) => {
        const output = data[0]
        const container = data[1]
        console.log(output)
        return container.remove()
      })
      .then((data) => {
        console.log('container removed', data)
      })
      .catch((err) => {
        console.log(err)
      })

    return result
  }
}
