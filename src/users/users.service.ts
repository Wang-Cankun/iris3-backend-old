import { Injectable } from '@nestjs/common'

export type User = any

@Injectable()
export class UsersService {
  private readonly users: User[]

  constructor() {
    this.users = [
      {
        userId: 1,
        username: 'john',
        info: 'Secret info 1',
        jobid: [12345, 45678],
        password: 'changeme'
      },
      {
        userId: 2,
        username: 'chris',
        info: 'Secret info 2',
        jobid: [12345, 45678],
        password: 'secret'
      },
      {
        userId: 3,
        username: 'maria',
        info: 'Secret info 3',
        jobid: [12345, 45678],
        password: 'guess'
      },
      {
        userId: 4,
        username: 'cankun.wang@osumc.edu',
        info: 'Cankun Wang information',
        jobid_public: [2020041684528, 2020080411720],
        jobid_private: [2019110564805, 2019102483326],
        password: '123456'
      }
    ]
  }

  async findOne(username: string): Promise<User | undefined> {
    return this.users.find((user) => user.username === username)
  }
  async queryUserInfo(username: string): Promise<User | undefined> {
    const { password, ...result } = this.users.find(
      (user) => user.username === username
    )
    return result
  }
}
