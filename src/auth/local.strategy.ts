import { Strategy } from 'passport-local'
import { PassportStrategy } from '@nestjs/passport'
import { Injectable, UnauthorizedException } from '@nestjs/common'
import { AuthService } from './auth.service'
import { User } from 'src/users/entities/user.entity'
import { UsersService } from 'src/users/users.service'

@Injectable()
export class LocalStrategy extends PassportStrategy(Strategy) {
  constructor(
    private authService: AuthService,
    private usersService: UsersService
  ) {
    super()
  }

  async validate(email: string, password: string): Promise<User> {
    const isValidUser = await this.authService.validateUser({
      email: email,
      password: password
    })
    if (!isValidUser) {
      throw new UnauthorizedException()
    }
    const user = await await this.usersService.findOneByEmail(email)
    return user
  }
}
