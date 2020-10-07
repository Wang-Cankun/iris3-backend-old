import { IsString, IsBoolean } from 'class-validator'
import { ApiProperty } from '@nestjs/swagger'

export class CreateJobDto {
  @ApiProperty({ description: 'The email of a user.' })
  @IsString()
  readonly email: string
  @IsString()
  readonly title: string
  @IsString()
  status: string
  @IsString()
  expFile: string
  @IsString()
  labelFile: string
  @IsString()
  species: string
  @IsString()
  description: string
}
