import { Module } from '@nestjs/common'
import { AppController } from './app.controller'
import { AppService } from './app.service'
import { UsersModule } from './users/users.module'
import { AuthModule } from './auth/auth.module'
import { FileModule } from './file/file.module'
import { TypeOrmModule } from '@nestjs/typeorm'
import { ConfigModule, ConfigService } from '@nestjs/config'
import { EmailModule } from './email/email.module'
import { DockerModule } from './docker/docker.module'
import { CommandModule } from './command/command.module'
import { QueueModule } from './queue/queue.module'
import * as Joi from '@hapi/joi'
import { MulterModule } from '@nestjs/platform-express'
import { JobModule } from './job/job.module'
import { PlumberModule } from './plumber/plumber.module'

@Module({
  imports: [
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
    MulterModule.registerAsync({
      imports: [ConfigModule],
      useFactory: async (configService: ConfigService) => ({
        dest: configService.get('UPLOAD_PATH')
      }),
      inject: [ConfigService]
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
    FileModule,
    AuthModule,
    UsersModule,
    EmailModule,
    DockerModule,
    CommandModule,
    QueueModule,
    JobModule,
    PlumberModule
  ],
  controllers: [AppController],
  providers: [AppService]
})
export class AppModule {}
