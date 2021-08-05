import { Injectable, NotFoundException } from '@nestjs/common'
import { InjectRepository } from '@nestjs/typeorm'
import { Repository } from 'typeorm'
import { Project } from '../project/entities/project.entity'
import { CreateJobDto } from './dto/create-job.dto'
import { UpdateJobDto } from './dto/update-job.dto'
import { Job } from './entities/job.entity'

@Injectable()
export class JobService {
  constructor(
    @InjectRepository(Job)
    private readonly jobRepository: Repository<Job>,

    @InjectRepository(Project)
    private readonly projectRepository: Repository<Project>
  ) {}

  async create(createJobDto: CreateJobDto): Promise<Job> {
    const project = await this.preloadProjectById(createJobDto.projectUid)
    console.log(project)
    const job = this.jobRepository.create({
      ...createJobDto,
      project
    })
    return this.jobRepository.save(job)
  }

  findAll() {
    return this.jobRepository.find({
      relations: ['project']
    })
  }

  async listJobsById(projectUid: string) {
    const job = await this.jobRepository.find({
      where: { projectUid },
      relations: ['project']
    })
    if (!job) {
      throw new NotFoundException(`Project ID #${projectUid} not found`)
    }
    return job
  }

  async findOne(id: number) {
    const job = await this.jobRepository.findOne(id)
    if (!job) {
      throw new NotFoundException(`Job ID #${id} not found`)
    }
    return job
  }

  async update(id: number, updateJobDto: UpdateJobDto) {
    const job = await this.jobRepository.preload({
      id: +id,
      ...updateJobDto
    })

    return this.jobRepository.save(job)
  }

  remove(id: number) {
    return `This action removes a #${id} job`
  }
  async preloadProjectById(projectUid: string): Promise<Project> {
    const existingProject = await this.projectRepository.findOne({ projectUid })
    return existingProject
  }
}
