import { Injectable, NotFoundException } from '@nestjs/common'
import { InjectRepository } from '@nestjs/typeorm'
import { Repository } from 'typeorm'
import { Job } from '../job/entities/job.entity'
import { User } from '../users/entities/user.entity'
import { CreateProjectDto } from './dto/create-project.dto'
import { UpdateProjectDto } from './dto/update-project.dto'
import { Project } from './entities/project.entity'

@Injectable()
export class ProjectService {
  constructor(
    @InjectRepository(User)
    private readonly userRepository: Repository<User>,

    @InjectRepository(Project)
    private readonly projectRepository: Repository<Project>,

    @InjectRepository(Job)
    private readonly jobRepository: Repository<Job>
  ) {}

  async create(createProjectDto: CreateProjectDto): Promise<Project> {
    const user = await this.preloadUserById(createProjectDto.userId)

    const jobs = await Promise.all(
      createProjectDto.jobIds.map((jobId) => this.preloadJobById(jobId))
    )

    const project = this.projectRepository.create({
      ...createProjectDto,
      jobs,
      user
    })
    return this.projectRepository.save(project)
  }

  listJobsById(projectUid: string) {
    return this.projectRepository.find({
      where: { projectUid },
      relations: ['jobs']
    })
  }

  findAll() {
    return this.projectRepository.find({
      relations: ['user', 'jobs']
    })
  }

  async findOne(id: number) {
    const project = await this.projectRepository.findOne(id)
    if (!project) {
      throw new NotFoundException(`Job ID #${id} not found`)
    }
    return project
  }

  async update(id: number, updateProjectDto: UpdateProjectDto) {
    return this.projectRepository
  }

  remove(id: number) {
    return `This action removes a #${id} project`
  }
  async preloadJobById(jobId: string): Promise<Job> {
    const existingJob = await this.jobRepository.findOne({ jobId })
    if (existingJob) {
      return existingJob
    }
    return this.jobRepository.create(existingJob)
  }

  async preloadUserById(userId: string): Promise<User> {
    const existingUser = await this.userRepository.findOne({ userId })
    if (existingUser) {
      return existingUser
    }
    return this.userRepository.create(existingUser)
  }
}
