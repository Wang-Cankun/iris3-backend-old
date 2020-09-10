import { Injectable, HttpException, HttpStatus } from '@nestjs/common'
import { UsersService } from '../users/users.service'
import { JwtService } from '@nestjs/jwt'
import { User } from 'src/users/entities/user.entity'
import { Login } from './interfaces/login.interface'

@Injectable()
export class AuthService {
  constructor(
    private usersService: UsersService,
    private jwtService: JwtService
  ) {}

  async validateUser(login: Login): Promise<User> {
    const userDb = await this.usersService.findOneByEmail(login.email)
    const isPasswordMatching = await this.usersService.validatePassword(
      login.password,
      userDb.password
    )
    if (!!userDb && isPasswordMatching) {
      return userDb
    }
    throw new HttpException('Wrong password', HttpStatus.UNAUTHORIZED)
  }

  async login(user: User): Promise<unknown> {
    const payload = {
      id: user.id
    }
    return {
      id: user.id,
      access_token: this.jwtService.sign(payload)
    }
  }
}
