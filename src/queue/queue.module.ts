import { Module } from '@nestjs/common'
import { QueueService } from './queue.service'
import { QueueController } from './queue.controller'
import { BullModule } from '@nestjs/bull'
import { ConfigModule, ConfigService } from '@nestjs/config'

@Module({
  imports: [
    BullModule.registerQueueAsync({
      name: 'queue',
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
  providers: [QueueService],
  controllers: [QueueController]
})
export class QueueModule {}
