import {
  CanActivate,
  ExecutionContext,
  Injectable,
  HttpException,
  HttpStatus
} from '@nestjs/common'
import { Observable } from 'rxjs'
import { InjectRepository } from '@nestjs/typeorm'
import { Repository } from 'typeorm'
import { User } from '../entities/user.entity'

@Injectable()
export class UserRegisteredGuard implements CanActivate {
  constructor(
    @InjectRepository(User)
    private readonly userRepository: Repository<User>
  ) {}

  async validate(email: string): Promise<boolean> {
    const user = await this.userRepository.count({ email: email })
    if (user) {
      throw new HttpException(
        `Email ${email} already registered`,
        HttpStatus.CONFLICT
      )
    }
    return true
  }
  canActivate(
    context: ExecutionContext
  ): boolean | Promise<boolean> | Observable<boolean> {
    const user = context.switchToHttp().getRequest().body
    return this.validate(user.email)
  }
}
