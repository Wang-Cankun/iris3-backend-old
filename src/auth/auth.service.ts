import { Injectable, HttpException, HttpStatus } from '@nestjs/common'
import { UsersService } from '../users/users.service'
import { JwtService } from '@nestjs/jwt'
import { User } from 'src/users/entities/user.entity'
import { Login } from './interfaces/login.interface'
import { ChangePasswordDto } from './dto/change-password.dto'
import { Repository } from 'typeorm'
import { InjectRepository } from '@nestjs/typeorm'
import { CreateUserDto } from 'src/users/dto/create-user.dto'
import * as bcrypt from 'bcrypt'
import { Job } from 'src/users/entities/job.entity'
@Injectable()
export class AuthService {
  constructor(
    @InjectRepository(User)
    private readonly userRepository: Repository<User>,
    @InjectRepository(Job)
    private readonly jobRepository: Repository<Job>,
    private usersService: UsersService,
    private jwtService: JwtService
  ) {}

  async hashPassword(plainPassword: string): Promise<string> {
    // Set 10 as salt
    const hashedPassword = await bcrypt.hash(plainPassword, 10)
    return hashedPassword
  }

  async validatePassword(
    plainPassword: string,
    hashedPassword: string
  ): Promise<boolean> {
    // Set 10 as salt
    // const hashedPassword = await this.hashPassword(plainPassword)
    const isPasswordMatching = await bcrypt.compare(
      plainPassword,
      hashedPassword
    )
    return isPasswordMatching
  }

  async validateUser(login: Login): Promise<boolean> {
    const userDb = await this.usersService.findOneByEmail(login.email)
    const isPasswordMatching = await this.validatePassword(
      login.password,
      userDb.password
    )
    if (!!userDb && isPasswordMatching) {
      return true
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

  async createUser(createUserDto: CreateUserDto): Promise<User> {
    const job = await Promise.all(
      createUserDto.job.map((email) =>
        this.usersService.preloadJobByName(email)
      )
    )
    createUserDto.password = await this.hashPassword(createUserDto.password)
    const user = this.userRepository.create({
      ...createUserDto,
      job
    })
    return this.userRepository.save(user)
  }

  async updatePassword(changePasswordDto: ChangePasswordDto): Promise<User> {
    const userFromDb = await this.usersService.findOneByEmail(
      changePasswordDto.email
    )
    const isValidUser = await this.validateUser({
      email: changePasswordDto.email,
      password: changePasswordDto.currentPassword
    })
    if (isValidUser) {
      userFromDb.password = await this.hashPassword(
        changePasswordDto.newPassword
      )
    }
    return this.userRepository.save(userFromDb)
  }
}
