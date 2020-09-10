export class ChangePasswordDto {
  readonly email: string
  readonly newPassword: string
  readonly currentPassword: string
}
