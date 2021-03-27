import { HttpService, Injectable } from '@nestjs/common'

@Injectable()
export class PlumberService {
  constructor(private httpService: HttpService) {}
  async runCommand(endpoint, data) {
    const result = await this.httpService
      .post(`http://localhost:8000/${endpoint}`, data)
      .toPromise()
      .then((res) => res.data)
    return result
  }
  async runStaticImage(endpoint, data) {
    const result = await this.httpService
      .post(`http://localhost:8000/${endpoint}`, data, {
        responseType: 'arraybuffer'
      })
      .toPromise()
      .then(
        (response) =>
          'data:image/png;base64,' +
          Buffer.from(response.data, 'binary').toString('base64')
      )
    return result
  }
}
