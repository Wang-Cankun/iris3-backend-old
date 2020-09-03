import {
  Injectable,
  NotFoundException,
  HttpException,
  HttpStatus
} from '@nestjs/common'
import { Repository, Connection } from 'typeorm'
import { InjectRepository } from '@nestjs/typeorm'
import { User } from './entities/user.entity'
import { Job } from './entities/job.entity'
import { CreateUserDto } from './dto/create-user.dto'
import { UpdateUserDto } from './dto/update-user.dto'
import * as nodemailer from 'nodemailer'

// export type User = any
export type Users = any
@Injectable()
export class UsersService {
  private readonly users: Users[]

  constructor(
    @InjectRepository(User)
    private readonly userRepository: Repository<User>,

    @InjectRepository(Job)
    private readonly jobRepository: Repository<Job>
  ) {}

  isValidEmail(email: string) {
    if (email) {
      const re = /^(([^<>()\[\]\\.,;:\s@"]+(\.[^<>()\[\]\\.,;:\s@"]+)*)|(".+"))@((\[[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}])|(([a-zA-Z\-0-9]+\.)+[a-zA-Z]{2,}))$/
      return re.test(email)
    } else return false
  }

  findAll() {
    return this.userRepository.find({
      relations: ['job']
    })
  }

  async findOneByEmail(email: string) {
    const user = await this.userRepository.findOne({
      where: {
        email: email
      }
    })
    if (!user) {
      return false
    }
    return user
  }

  async findOneById(id: string) {
    const user = await this.userRepository.findOne(id)
    if (!user) {
      throw new NotFoundException(`User id #${id} not found`)
    }
    return user
  }

  async findOne(email: string) {
    const user = await this.userRepository.findOne({
      where: {
        email: email
      }
    })
    if (!user) {
      throw new NotFoundException(`User email #${email} not found`)
    }
    return user
  }

  async create(createUserDto: CreateUserDto): Promise<User> {
    const job = await Promise.all(
      createUserDto.job.map((email) => this.preloadJobByName(email))
    )
    const userFromDb = await this.findOneByEmail(createUserDto.email)
    if (userFromDb)
      throw new HttpException(
        'REGISTRATION.USER_ALREADY_REGISTERED',
        HttpStatus.CONFLICT
      )
    const user = this.userRepository.create({
      ...createUserDto,
      job
    })
    return this.userRepository.save(user)
  }

  async update(id: string, updateUserDto: UpdateUserDto) {
    const job =
      updateUserDto.job &&
      (await Promise.all(
        updateUserDto.job.map((email) => this.preloadJobByName(email))
      ))

    const user = await this.userRepository.preload({
      id: +id,
      ...updateUserDto,
      job
    })
    if (!user) {
      throw new NotFoundException(`user #${id} not found`)
    }
    return this.userRepository.save(user)
  }

  async updatePassword(
    email: string,
    updateUserDto: UpdateUserDto
  ): Promise<User> {
    const userFromDb = await this.findOne(email)
    if (!userFromDb)
      throw new HttpException('LOGIN.USER_NOT_FOUND', HttpStatus.NOT_FOUND)
    userFromDb.password = updateUserDto.password
    return this.userRepository.save(userFromDb)
  }

  async updateProfile(
    email: string,
    updateUserDto: UpdateUserDto
  ): Promise<User> {
    const userFromDb = await this.findOne(email)
    userFromDb.firstName = updateUserDto.firstName
    userFromDb.lastName = updateUserDto.lastName
    userFromDb.institution = updateUserDto.institution
    userFromDb.newsletter = updateUserDto.newsletter
    return this.userRepository.save(userFromDb)
  }

  private async preloadJobByName(email: string): Promise<Job> {
    const existingJob = await this.jobRepository.findOne({ email })
    if (existingJob) {
      return existingJob
    }
    return this.jobRepository.create({ email })
  }

  async remove(id: string) {
    const coffee = await this.findOne(id)
    return this.userRepository.remove(coffee)
  }

  async queryUserInfo(username: string): Promise<User | undefined> {
    const { ...result } = await this.users.find(
      (user) => user.username === username
    )
    return result
  }

  sendEmailVerification(email: string): void {
    const content =
      `<table class='body' style='border-collapse: separate; mso-table-lspace: 0pt; mso-table-rspace: 0pt; width: 100%; background-color: #f6f6f6;' border='0' cellspacing='0' cellpadding='0'>
  <tbody>
  <tr>
  <td style='font-family: sans-serif; font-size: 14px; vertical-align: top;'>&nbsp;</td>
  <td class='container' style='font-family: sans-serif; font-size: 14px; vertical-align: top; display: block; margin: 0 auto; max-width: 700px; padding: 10px; width: 700px;'>
  <div class='content' style='box-sizing: border-box; display: block; margin: 0 auto; max-width: 700px; padding: 10px;'><!-- START CENTERED WHITE CONTAINER --> <span class='preheader' style='color: transparent; display: none; height: 0; max-height: 0; max-width: 0; opacity: 0; overflow: hidden; mso-hide: all; visibility: hidden; width: 0;'><g class='gr_ gr_65 gr-alert gr_spell gr_inline_cards gr_run_anim ContextualSpelling ins-del multiReplace' id='65' data-gr-id='65'></g></span>
  <table class='main' style='border-collapse: separate; mso-table-lspace: 0pt; mso-table-rspace: 0pt; width: 100%; background: #ffffff; border-radius: 3px;'><!-- START MAIN CONTENT AREA -->
  <tbody>
  <tr>
  <td class='wrapper' style='font-family: sans-serif; font-size: 14px; vertical-align: top; box-sizing: border-box; padding: 20px;'><br />
  <table style='border-collapse: separate; mso-table-lspace: 0pt; mso-table-rspace: 0pt; width: 100%;' border='0' cellspacing='0' cellpadding='0'>
  <tbody>
  <tr>
  <td style='font-family: sans-serif; font-size: 14px; vertical-align: top;'>
  <table border='0' width='100%' cellspacing='0' cellpadding='0'>
  <tbody>
  <tr>
  <td class='subhead'>Hello,</td>
  </tr>
  <tr>
  <td class='h1' style='padding: 5px 0 0 0;'><br />
  <div>
  <div><span>You have submitted a request to reset your password.<br /></span><span>Your email: ` +
      email +
      `</span></div>
  </div>
  </div>
  </td>
  </tr>
  </tbody>
  </table>
  <p style='font-family: sans-serif; font-size: 14px; font-weight: normal; margin: 0; margin-bottom: 15px;'></p>
  <span>&nbsp;</span></td>
  </tr>
  <tr><td style='font-family: sans-serif; font-size: 14px; vertical-align: top; background-color: #3498db; border-radius: 5px; text-align: center;'><a style='display: inline-block; color: #ffffff; background-color: #3498db; border: solid 1px #3498db; border-radius: 5px; box-sizing: border-box; cursor: pointer; text-decoration: none; font-size: 14px; font-weight: bold; margin: 0; padding: 12px 25px; text-transform: capitalize; border-color: #3498db;' href='https://bmbls.bmi.osumc.edu/iris3/' target='_blank' rel='noopener'>Click to reset your password</a></td></tr>
  </tbody>
  </table>
  </td>
  </tr>
  </tbody>
  </table>
  </td>
  </tr>
  <!-- END MAIN CONTENT AREA --></tbody>
  </table>
  <!-- START FOOTER -->
  <div class='footer' style='clear: both; margin-top: 10px; text-align: center; width: 100%;'>
  <table style='border-collapse: separate; mso-table-lspace: 0pt; mso-table-rspace: 0pt; width: 100%;' border='0' cellspacing='0' cellpadding='0'>
  <tbody>
  <tr>
  <td class='content-block' style='font-family: sans-serif; vertical-align: top; padding-bottom: 10px; padding-top: 10px; font-size: 12px; color: #999999; text-align: center;'><span>Copyright 2020 &copy; </span><a href='http://u.osu.edu/bmbl' target='_blank' rel='noopener'>BMBL</a><span>, </span><a href='https://medicine.osu.edu/departments/biomedical-informatics/' target='_blank' rel='noopener'>OSU</a><span>. All rights reserved. </span></td>
  </tr>
  <tr>
  <td class='content-block powered-by' style='font-family: sans-serif; vertical-align: top; padding-bottom: 10px; padding-top: 10px; font-size: 12px; color: #999999; text-align: center;'><a href='mailto:qin.ma\@osumc.edu' title='qin.ma\@osumc.edu'>Contact us: qin.ma\@osumc.edu</a><span> </span></td>
  </tr>
  </tbody>
  </table>
  </div>
   <!-- END CENTERED WHITE CONTAINER --></div>
  </td>
  <td style='font-family: sans-serif; font-size: 14px; vertical-align: top;'>&nbsp;</td>
  </tr>
  </tbody>
  </table>`
    const transporter = nodemailer.createTransport({
      service: 'Gmail',
      auth: {
        user: process.env.EMAIL_NAME,
        pass: process.env.EMAIL_PASSWORD
      }
    })

    const mailOptions = {
      from: 'IRIS3 <no-reply@bmbls.bmi.osumc.edu>',
      to: email,
      subject: 'Reset password from IRIS3',
      html: content
    }

    transporter.sendMail(mailOptions, function(error, info) {
      if (error) {
        console.log(error)
      } else {
        console.log('Email sent: ' + info.response)
      }
    })
  }
}
