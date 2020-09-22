import { Module } from '@nestjs/common'
import { UsersService } from './users.service'
import { TypeOrmModule } from '@nestjs/typeorm'
import { User } from './entities/user.entity'
import { Job } from './entities/job.entity'
import { UsersController } from './users.controller'

@Module({
  imports: [TypeOrmModule.forFeature([User, Job])],
  providers: [UsersService],
  exports: [UsersService],
  controllers: [UsersController]
})
export class UsersModule {}
