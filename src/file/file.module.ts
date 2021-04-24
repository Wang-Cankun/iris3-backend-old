import { BullModule } from '@nestjs/bull'
import { Module } from '@nestjs/common'
import { ConfigModule, ConfigService } from '@nestjs/config'
import { MulterModule } from '@nestjs/platform-express'
import { TypeOrmModule } from '@nestjs/typeorm'
import { Job } from '../job/entities/job.entity'
import { File } from './entities/file.entity'
import { FileController } from './file.controller'
import { FileProcessor } from './file.processor'
import { FileService } from './file.service'

@Module({
  imports: [
    TypeOrmModule.forFeature([File, Job]),
    BullModule.registerQueueAsync({
      name: 'file',
      imports: [ConfigModule],
      useFactory: async (configService: ConfigService) => ({
        redis: {
          host: configService.get('REDIS_HOST'),
          port: +configService.get('REDIS_PORT')
        }
      }),
      inject: [ConfigService]
    }),
    MulterModule.registerAsync({
      imports: [ConfigModule],
      useFactory: async (configService: ConfigService) => ({
        dest: configService.get<string>('UPLOAD_PATH')
      }),
      inject: [ConfigService]
    })
  ],
  controllers: [FileController],
  providers: [FileService, FileProcessor]
})
export class FileModule {}
