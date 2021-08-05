import { Module } from '@nestjs/common'
import { ProjectService } from './project.service'
import { ProjectController } from './project.controller'
import { TypeOrmModule } from '@nestjs/typeorm'
import { User } from 'src/users/entities/user.entity'
import { Project } from './entities/project.entity'
import { Job } from '../job/entities/job.entity'

@Module({
  imports: [TypeOrmModule.forFeature([Project, Job, User])],
  controllers: [ProjectController],
  providers: [ProjectService],
  exports: [ProjectService]
})
export class ProjectModule {}
