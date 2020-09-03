import { Module } from '@nestjs/common'
import { UsersService } from './users.service'
import { TypeOrmModule } from '@nestjs/typeorm'
import { ConfigModule } from '@nestjs/config'
import { User } from './entities/user.entity'
import { Job } from './entities/job.entity'
import { UsersController } from './users.controller'

@Module({
  imports: [TypeOrmModule.forFeature([User, Job]), ConfigModule],
  providers: [UsersService],
  exports: [UsersService],
  controllers: [UsersController]
})
export class UsersModule {}
