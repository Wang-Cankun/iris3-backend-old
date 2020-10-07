import { Injectable, NotFoundException } from '@nestjs/common'
import { InjectRepository } from '@nestjs/typeorm'
import { User } from './entities/user.entity'
import { UpdateUserDto } from './dto/update-user.dto'
import { Repository } from 'typeorm'
import { Job } from 'src/job/entities/job.entity'

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

  isValidEmail(email: string): boolean {
    if (email) {
      const re = /^(([^<>()\[\]\\.,;:\s@"]+(\.[^<>()\[\]\\.,;:\s@"]+)*)|(".+"))@((\[[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}])|(([a-zA-Z\-0-9]+\.)+[a-zA-Z]{2,}))$/
      return re.test(email)
    } else return false
  }

  findAll(): Promise<User[]> {
    return this.userRepository.find({
      relations: ['job']
    })
  }

  async findOneByEmail(email: string): Promise<User> {
    const user = await this.userRepository.findOne({
      where: {
        email: email
      }
    })
    if (!user) {
      throw new NotFoundException(`Email ${email} not found`)
    }
    return user
  }

  async findOneById(id: string): Promise<User> {
    const user = await this.userRepository.findOne(id)
    if (!user) {
      throw new NotFoundException(`User ID #${id} not found`)
    }
    return user
  }

  async update(id: string, updateUserDto: UpdateUserDto): Promise<User> {
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
      throw new NotFoundException(`User ID #${id} not found`)
    }
    return this.userRepository.save(user)
  }

  async updateProfile(
    email: string,
    updateUserDto: UpdateUserDto
  ): Promise<User> {
    const userFromDb = await this.findOneByEmail(email)
    if (userFromDb.email !== updateUserDto.email) {
      userFromDb.email = updateUserDto.email
    }

    userFromDb.firstName = updateUserDto.firstName
    userFromDb.lastName = updateUserDto.lastName
    userFromDb.institution = updateUserDto.institution
    userFromDb.newsletter = updateUserDto.newsletter
    return this.userRepository.save(userFromDb)
  }

  async preloadJobByName(id: string): Promise<Job> {
    const existingJob = await this.jobRepository.findOne(id)
    if (existingJob) {
      return existingJob
    }
    return this.jobRepository.create(existingJob)
  }

  async test(updateUserDto: UpdateUserDto): Promise<any> {
    const user = await this.findOneByEmail(updateUserDto.email)

    console.log(user)
    return 1
  }
}
