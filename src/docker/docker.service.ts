import { Injectable } from '@nestjs/common'
import * as Docker from 'dockerode'

@Injectable()
export class DockerService {
  private readonly docker
  private readonly container
  private readonly CTRL_P = '\u0010\n'
  private readonly CTRL_Q = '\u0011\n'
  constructor() {
    this.docker = new Docker()
    this.container = this.docker.getContainer('iris3-workflow-env')
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
      return await this.docker.ping()
    } catch (err) {
      throw new Error()
    }
  }
  async testExec(param) {
    const logs = await this.container.logs({
      follow: true,
      stdout: true,
      stderr: true,
      details: false,
      tail: 50,
      timestamps: true
    })

    // const command = ['library(Seurat)']
    const options = {
      Cmd: ['R'],
      AttachStdout: true,
      AttachStderr: true,
      hijack: true,
      stdin: true,
      Detach: true,
      tty: false
    }

    const exec = await this.container.exec(options)

    return new Promise(async (resolve, reject) => {
      await exec.start(async (err, stream) => {
        if (err) return reject()
        let message = ''
        stream.on('data', (data) => (message += data.toString()))
        stream.on('end', () => resolve(message))
      })
    })
  }

  async test(param) {
    const options = {
      stream: true,
      stdin: true,
      stdout: true,
      stderr: true
    }
    console.log(await this.ping())
    return new Promise((resolve, reject) => {
      this.container.attach(options, async (err, stream) => {
        if (err) return reject(err)
        stream.write('R\ndate()\n')
        resolve('done')
      })
    })
  }
}
