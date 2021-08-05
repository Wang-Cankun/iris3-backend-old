import { IsString, IsBoolean } from 'class-validator'
import { ApiProperty } from '@nestjs/swagger'

export class CreateJobDto {
  @IsString()
  title: string

  @IsString()
  type: string

  @IsString()
  description: string

  @IsString()
  projectUid: string
}
