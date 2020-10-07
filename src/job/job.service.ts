import { Injectable, NotFoundException } from '@nestjs/common'
import { InjectRepository } from '@nestjs/typeorm'
import { Repository } from 'typeorm'
import { CreateJobDto } from './dto/create-job.dto'
import { UpdateJobDto } from './dto/update-job.dto'
import { Job } from './entities/job.entity'

@Injectable()
export class JobService {
  constructor(
    @InjectRepository(Job)
    private readonly jobRepository: Repository<Job>
  ) {}

  async create(createJobDto: CreateJobDto): Promise<Job> {
    const job = this.jobRepository.create({
      ...createJobDto
    })
    return this.jobRepository.save(job)
  }

  findAll() {
    return this.jobRepository.find({
      relations: ['user']
    })
  }

  async findOne(id: number) {
    const job = await this.jobRepository.findOne(id)
    if (!job) {
      throw new NotFoundException(`Job ID #${id} not found`)
    }
    return job
  }

  async update(id: number, updateJobDto: UpdateJobDto) {
    const user = []
    const job = await this.jobRepository.preload({
      id: +id,
      ...updateJobDto,
      user
    })
    if (!user) {
      throw new NotFoundException(`Job ID #${id} not found`)
    }
    return this.jobRepository.save(job)
  }

  remove(id: number) {
    return `This action removes a #${id} job`
  }
  async preloadJobByName(id: string): Promise<Job> {
    const existingJob = await this.jobRepository.findOne(id)
    if (existingJob) {
      return existingJob
    }
    return this.jobRepository.create(existingJob)
  }
}
