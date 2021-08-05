import { Module } from '@nestjs/common'
import { JobService } from './job.service'
import { JobController } from './job.controller'
import { TypeOrmModule } from '@nestjs/typeorm'
import { User } from 'src/users/entities/user.entity'
import { Job } from './entities/job.entity'
import { Project } from '../project/entities/project.entity'

@Module({
  imports: [TypeOrmModule.forFeature([User, Job, Project])],
  controllers: [JobController],
  providers: [JobService],
  exports: [JobService]
})
export class JobModule {}
