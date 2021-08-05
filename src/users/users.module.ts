import { Module } from '@nestjs/common'
import { UsersService } from './users.service'
import { TypeOrmModule } from '@nestjs/typeorm'
import { User } from './entities/user.entity'
import { UsersController } from './users.controller'
import { Project } from '../project/entities/project.entity'

@Module({
  imports: [TypeOrmModule.forFeature([User, Project])],
  providers: [UsersService],
  exports: [UsersService],
  controllers: [UsersController]
})
export class UsersModule {}
