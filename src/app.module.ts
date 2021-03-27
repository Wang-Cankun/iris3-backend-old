import { Module } from '@nestjs/common'
import { AppController } from './app.controller'
import { AppService } from './app.service'
import { UsersModule } from './users/users.module'
import { AuthModule } from './auth/auth.module'
import { UploadModule } from './upload/upload.module'
import { TypeOrmModule } from '@nestjs/typeorm'
import { ConfigModule } from '@nestjs/config'
import { EmailModule } from './email/email.module'
import { DockerModule } from './docker/docker.module'
import { CommandModule } from './command/command.module'
import { ServeStaticModule } from '@nestjs/serve-static'
import { join } from 'path'
import { BullModule } from '@nestjs/bull'
import { QueueModule } from './queue/queue.module'
import { EventsModule } from './events/events.module'
import * as Joi from '@hapi/joi'
import { MulterModule } from '@nestjs/platform-express'
import { JobModule } from './job/job.module'
import { PlumberService } from './plumber/plumber.service'
import { PlumberModule } from './plumber/plumber.module'

@Module({
  imports: [
    MulterModule.register({
      dest: '/tmp'
    }),
    ConfigModule.forRoot({
      isGlobal: true,

      envFilePath: '.env',
      validationSchema: Joi.object({
        NODE_ENV: Joi.string()
          .valid('development', 'production', 'test')
          .default('development'),
        PORT: Joi.number().default(9005),
        MYSQL_HOST: Joi.required(),
        MYSQL_PORT: Joi.number().default(3306),
        MYSQL_USER: Joi.required(),
        MYSQL_PASSWORD: Joi.required(),
        REDIS_PORT: Joi.number().default(6379),
        REDIS_HOST: Joi.required(),
        JWT_SECRET: Joi.required(),
        JWT_EXPIRATION_TIME: Joi.required(),
        FRONTEND_URL: Joi.required(),
        EMAIL_NAME: Joi.required(),
        EMAIL_PASSWORD: Joi.required(),
        GOOGLE_CLIENT_ID: Joi.required(),
        GOOGLE_SECRET: Joi.required()
      })
    }),
    TypeOrmModule.forRoot({
      type: 'mysql',
      host: process.env.MYSQL_HOST,
      port: +process.env.MYSQL_PORT,
      username: process.env.MYSQL_USER,
      password: process.env.MYSQL_PASSWORD,
      database: process.env.MYSQL_NAME,
      autoLoadEntities: true,
      synchronize: true
    }),
    AuthModule,
    UsersModule,
    UploadModule,
    EmailModule,
    DockerModule,
    CommandModule,
    QueueModule,
    EventsModule,
    JobModule,
    PlumberModule
  ],
  controllers: [AppController],
  providers: [AppService]
})
export class AppModule {}
