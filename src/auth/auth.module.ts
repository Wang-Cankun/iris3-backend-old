import { Module } from '@nestjs/common'
import { AuthService } from './auth.service'
import { LocalStrategy } from './local.strategy'
import { JwtStrategy } from './jwt.strategy'
import { GoogleStrategy } from './google.strategy'
import { UsersModule } from '../users/users.module'
import { PassportModule } from '@nestjs/passport'
import { JwtModule } from '@nestjs/jwt'
import { AuthController } from './auth.controller'
import { UsersService } from 'src/users/users.service'
import { TypeOrmModule } from '@nestjs/typeorm'
import { User } from 'src/users/entities/user.entity'
import { EmailService } from 'src/email/email.service'
import { ConfigModule, ConfigService } from '@nestjs/config'
import { Job } from 'src/job/entities/job.entity'
import { Project } from '../project/entities/project.entity'

@Module({
  imports: [
    TypeOrmModule.forFeature([User, Job, Project]),
    UsersModule,
    PassportModule,
    JwtModule.registerAsync({
      imports: [ConfigModule],
      useFactory: async (configService: ConfigService) => {
        return {
          secret: configService.get<string>('JWT_SECRET'),
          signOptions: {
            expiresIn: configService.get<string>('JWT_EXPIRATION_TIME')
          }
        }
      },
      inject: [ConfigService]
    }),
    ConfigModule
  ],
  providers: [
    AuthService,
    UsersService,
    EmailService,
    LocalStrategy,
    JwtStrategy,
    GoogleStrategy
  ],
  exports: [AuthService],
  controllers: [AuthController]
})
export class AuthModule {}
