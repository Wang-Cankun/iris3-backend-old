import { ApiProperty } from '@nestjs/swagger'

export class TokenDto {
  constructor(data: { expiresIn?: number; accessToken?: string }) {
    this.expiresIn = data.expiresIn
    this.accessToken = data.accessToken
  }

  @ApiProperty({ description: 'Expired time' }) expiresIn: number

  @ApiProperty() accessToken: string
}
