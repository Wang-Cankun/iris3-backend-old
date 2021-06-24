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
import { QueueModule } from './queue/queue.module'
import * as Joi from '@hapi/joi'
import { MulterModule } from '@nestjs/platform-express'
import { JobModule } from './job/job.module'
import { PlumberModule } from './plumber/plumber.module'
import { SharedModule } from './shared/shared.module'
import { EventsGateway } from './test.gateway'

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
        EMAIL_TOKEN: Joi.required(),
        GOOGLE_CLIENT_ID: Joi.required(),
        GOOGLE_SECRET: Joi.required()
      })
    }),
    MulterModule.registerAsync({
      imports: [ConfigModule],
      useFactory: async (configService: ConfigService) => ({
        dest: configService.get<string>('UPLOAD_PATH')
      }),
      inject: [ConfigService]
    }),

    TypeOrmModule.forRootAsync({
      imports: [ConfigModule],
      useFactory: (config: ConfigService) => ({
        type: 'mysql',
        host: config.get('MYSQL_HOST'),
        port: config.get('MYSQL_PORT'),
        username: config.get('MYSQL_USER'),
        password: config.get('MYSQL_PASSWORD'),
        database: config.get('MYSQL_NAME'),
        autoLoadEntities: true,
        synchronize: true
      }),
      inject: [ConfigService]
    }),
    FileModule,
    AuthModule,
    UsersModule,
    EmailModule,
    DockerModule,
    QueueModule,
    JobModule,
    PlumberModule,
    SharedModule
  ],
  controllers: [AppController],
  providers: [AppService, EventsGateway]
})
export class AppModule {}
