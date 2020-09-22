import {
  Controller,
  Post,
  Get,
  Body,
  Request,
  Param,
  UseGuards,
  Patch,
  Req,
  Delete
} from '@nestjs/common'
import { ApiTags } from '@nestjs/swagger'
import { AuthService } from './auth.service'
import { UsersService } from 'src/users/users.service'
import { JwtAuthGuard } from './guards/jwt-auth.guard'
import { CreateUserDto } from 'src/users/dto/create-user.dto'
import { User } from 'src/users/entities/user.entity'
import { ChangePasswordDto } from './dto/change-password.dto'
import { UserRegisteredGuard } from 'src/users/guards/user-registered.guard'
import { UserExistGuard } from 'src/users/guards/user-exist.guard'
import { LoginDto } from './dto/login.dto'
import { EmailDto } from './dto/email.dto'
import { ResetPasswordDto } from './dto/reset-password.dto'
import { AuthGuard } from '@nestjs/passport'
import { LocalAuthGuard } from './guards/local-auth.guard'
import { JwtPayload } from './decorators/jwt-payload.decorator'
import { JwtPayloadDto } from './dto/jwt-payload.dto'

@ApiTags('Authentication')
@Controller('auth')
export class AuthController {
  constructor(
    private readonly authService: AuthService,
    private readonly usersService: UsersService
  ) {}

  @Get('test')
  public test(@Param() params): any {
    return this.authService.createToken('flykun0620@gmail.com')
  }

  @UseGuards(LocalAuthGuard)
  @Post('/login')
  async login(@Request() req) {
    return this.authService.login(req.user)
  }

  @UseGuards(JwtAuthGuard)
  @Post('/logout')
  async logout(@Request() req) {
    return 'logout to be add'
  }

  @Get('/google')
  @UseGuards(AuthGuard('google'))
  async googleAuth(@Req() req) {}

  @Get('/google/redirect')
  @UseGuards(AuthGuard('google'))
  googleAuthRedirect(@Req() req) {
    return this.authService.googleLogin(req)
  }

  @UseGuards(JwtAuthGuard)
  @Get('profile')
  getProfile(@JwtPayload() payload: JwtPayloadDto): JwtPayloadDto {
    return payload
  }

  @UseGuards(UserRegisteredGuard)
  @Post('/register')
  create(@Body() createUserDto: CreateUserDto): Promise<User> {
    return this.authService.createUser(createUserDto)
  }

  @UseGuards(UserExistGuard)
  @UseGuards(JwtAuthGuard)
  @Post('/password/change')
  updatePassword(@Body() changePasswordDto: ChangePasswordDto): Promise<User> {
    return this.authService.updatePassword(changePasswordDto)
  }

  @UseGuards(UserExistGuard)
  @Post('/password/forgot')
  setNewPassword(@Body() resetPasswordDto: ResetPasswordDto): Promise<User> {
    return this.authService.setForgotPassword(resetPasswordDto)
  }

  @UseGuards(UserExistGuard)
  @UseGuards(JwtAuthGuard)
  @Delete('/delete')
  remove(@Body() loginDto: LoginDto): Promise<User> {
    return this.authService.removeAccount(loginDto)
  }

  @UseGuards(UserExistGuard)
  @Post('/forgot')
  sendResetEmail(@Body() emailDto: EmailDto): any {
    return this.authService.resetPassword(emailDto.email)
  }
}
