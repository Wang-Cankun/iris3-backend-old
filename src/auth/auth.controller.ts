import {
  Controller,
  Post,
  HttpStatus,
  HttpCode,
  Get,
  Body,
  Request,
  Param,
  UseGuards,
  Patch,
  ParseIntPipe
} from '@nestjs/common'
import { AuthService } from './auth.service'
import { UsersService } from 'src/users/users.service'
import { JwtAuthGuard } from './jwt-auth.guard'
import { LocalAuthGuard } from './local-auth.guard'
import { CreateUserDto } from 'src/users/dto/create-user.dto'
import { User } from 'src/users/entities/user.entity'
import { ChangePasswordDto } from './dto/change-password.dto'
import { UserRegisteredGuard } from 'src/users/guards/user-registered.guard'
import { UserExistGuard } from 'src/users/guards/user-exist.guard'

@Controller('auth')
export class AuthController {
  constructor(
    private readonly authService: AuthService,
    private readonly usersService: UsersService
  ) {}

  @Get('test')
  public verifyEmail(@Param() params): any {
    return '123'
  }

  @UseGuards(LocalAuthGuard)
  @Post('/login')
  async login(@Request() req) {
    return this.authService.login(req.user)
  }

  @UseGuards(JwtAuthGuard)
  @Get('profile')
  getProfile(@Request() req) {
    return req.user
  }

  @UseGuards(UserRegisteredGuard)
  @Post('/register')
  create(@Body() createUserDto: CreateUserDto): Promise<User> {
    return this.authService.createUser(createUserDto)
  }

  @UseGuards(JwtAuthGuard)
  @UseGuards(UserExistGuard)
  @Patch('/password')
  updatePassword(@Body() changePasswordDto: ChangePasswordDto): Promise<User> {
    return this.authService.updatePassword(changePasswordDto)
  }
}
