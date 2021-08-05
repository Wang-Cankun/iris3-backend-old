import { IsString, IsBoolean } from 'class-validator'
import { ApiProperty } from '@nestjs/swagger'

export class CreateProjectDto {
  @IsString()
  private: boolean

  @IsString()
  title: string

  @IsString()
  type: string

  @IsString()
  status: string

  @IsString({ each: true })
  tag: string[]

  @IsString()
  description: string

  @IsString()
  readonly userId: string

  @IsString({ each: true })
  readonly jobIds: string[]
}
