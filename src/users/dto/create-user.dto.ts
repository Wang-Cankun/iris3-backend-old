import { IsString, IsBoolean } from 'class-validator'
import { ApiProperty } from '@nestjs/swagger'
import { isString, isBoolean } from 'util'
export class CreateUserDto {
  @ApiProperty({ description: 'The email of a user.' })
  @IsString()
  readonly email: string
  @IsString()
  password: string
  @IsString()
  readonly firstName: string
  @IsString()
  readonly lastName: string
  @IsString()
  readonly institution: string
  @IsBoolean()
  readonly newsletter: boolean
  @IsBoolean()
  readonly isActive: boolean
  @IsString({ each: true })
  readonly job: string[]
}
