import { HttpModule, Module } from '@nestjs/common'
import { QueueService } from './queue.service'
import { QueueController } from './queue.controller'
import { BullModule } from '@nestjs/bull'
import { ConfigModule, ConfigService } from '@nestjs/config'
import { QueueProcessor } from './queue.processor'
import { JobService } from 'src/job/job.service'
import { JobModule } from 'src/job/job.module'
import { TypeOrmModule } from '@nestjs/typeorm'
import { Job } from 'src/job/entities/job.entity'
import { PlumberService } from '../plumber/plumber.service'
import { FileService } from '../file/file.service'
import { FileModule } from '../file/file.module'
import { File } from '../file/entities/file.entity'
import { Project } from '../project/entities/project.entity'

@Module({
  imports: [
    TypeOrmModule.forFeature([Job, File, Project]),
    HttpModule,
    BullModule.registerQueueAsync({
      name: 'task',
      imports: [ConfigModule],
      useFactory: async (configService: ConfigService) => ({
        redis: {
          host: configService.get('REDIS_HOST'),
          port: +configService.get('REDIS_PORT')
        }
      }),
      inject: [ConfigService]
    })
  ],
  providers: [
    FileService,
    QueueService,
    QueueProcessor,
    JobService,
    PlumberService
  ],
  exports: [QueueProcessor],
  controllers: [QueueController]
})
export class QueueModule {}
