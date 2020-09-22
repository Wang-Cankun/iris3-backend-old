import { Injectable, HttpException, HttpStatus } from '@nestjs/common'
import { UsersService } from '../users/users.service'
import { JwtService } from '@nestjs/jwt'
import { User } from 'src/users/entities/user.entity'
import { ChangePasswordDto } from './dto/change-password.dto'
import { Repository } from 'typeorm'
import { InjectRepository } from '@nestjs/typeorm'
import { CreateUserDto } from 'src/users/dto/create-user.dto'
import * as bcrypt from 'bcrypt'
import * as jwt from 'jsonwebtoken'
import { Job } from 'src/users/entities/job.entity'
import { LoginDto } from './dto/login.dto'
import { EmailService } from 'src/email/email.service'
import { ResetPasswordDto } from './dto/reset-password.dto'
import { TokenDto } from './dto/token.dto'
@Injectable()
export class AuthService {
  constructor(
    @InjectRepository(User)
    private readonly userRepository: Repository<User>,
    @InjectRepository(Job)
    private readonly jobRepository: Repository<Job>,
    private usersService: UsersService,
    private emailService: EmailService,
    private jwtService: JwtService
  ) {}

  public async createToken(email: string): Promise<string> {
    const expiresIn = 3600
    const user = await this.usersService.findOneByEmail(email)
    const secret = user.password + '-' + user.createTime

    const token = jwt.sign({ email: email }, secret, {
      expiresIn: expiresIn // 1 hour
    })

    return token
  }

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

  async validateUser(loginDto: LoginDto): Promise<boolean> {
    const userDb = await this.usersService.findOneByEmail(loginDto.email)
    const isPasswordMatching = await this.validatePassword(
      loginDto.password,
      userDb.password
    )
    if (!!userDb && isPasswordMatching) {
      return true
    }
    throw new HttpException('Wrong password', HttpStatus.UNAUTHORIZED)
  }

  async login(user: User): Promise<any> {
    const payload = {
      email: user.email
    }
    return {
      email: user.email,
      access_token: this.jwtService.sign(payload)
    }
  }

  googleLogin(req) {
    if (!req.user) {
      return 'No user from google'
    }

    return {
      message: 'User information from google',
      user: req.user
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

  async setForgotPassword(resetPasswordDto: ResetPasswordDto): Promise<User> {
    const userFromDb = await this.usersService.findOneByEmail(
      resetPasswordDto.email
    )
    const payload = jwt.decode(resetPasswordDto.token)

    if (payload) {
      userFromDb.password = await this.hashPassword(resetPasswordDto.password)
    }
    return this.userRepository.save(userFromDb)
  }

  async removeAccount(loginDto: LoginDto): Promise<User> {
    const userFromDb = await this.usersService.findOneByEmail(loginDto.email)
    const isValidUser = await this.validateUser({
      email: loginDto.email,
      password: loginDto.password
    })
    if (isValidUser) {
      return this.userRepository.remove(userFromDb)
    }
    throw new HttpException('Wrong password', HttpStatus.UNAUTHORIZED)
  }

  async resetPassword(email: string): Promise<any> {
    const token = await this.createToken(email)
    const sendResetPasswordEmail = this.emailService.sendForgotPasswordEmail(
      email,
      token
    )
    return sendResetPasswordEmail
  }
}
