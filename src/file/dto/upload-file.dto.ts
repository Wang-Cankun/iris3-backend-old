import { IsString, IsBoolean, IsNumber } from 'class-validator'
export class UploadFileDto {
  @IsString()
  jobid: string

  @IsString()
  index: number

  @IsNumber()
  species: string

  @IsString()
  fieldname: string

  @IsString()
  originalname: string

  @IsString()
  encoding: string

  @IsString()
  mimetype: string

  @IsString()
  destination: string

  @IsString()
  filename: string

  @IsString()
  path: string

  @IsNumber()
  size: number
  @IsString()
  private: boolean
  @IsString()
  creator: string
  @IsString()
  tags: string
  @IsString()
  title: string
}
