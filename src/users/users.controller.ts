import {
  Body,
  Controller,
  Delete,
  Get,
  Param,
  Patch,
  Post,
  UseGuards,
  Request
} from '@nestjs/common'
import { UsersService } from './users.service'
import { UpdateUserDto } from './dto/update-user.dto'
import { ApiTags } from '@nestjs/swagger'
import { User } from './entities/user.entity'
import { UserExistGuard } from './guards/user-exist.guard'
import { JwtAuthGuard } from 'src/auth/guards/jwt-auth.guard'
import { JwtPayload } from 'src/auth/decorators/jwt-payload.decorator'

@ApiTags('users')
// @UseGuards(JwtAuthGuard)
@Controller('users')
export class UsersController {
  constructor(private readonly usersService: UsersService) {}

  // Delete this route later
  @Get()
  async findAll(): Promise<User[]> {
    // const { limit, offset } = paginationQuery;
    // await new Promise((resolve) => setTimeout(resolve, 500))
    return this.usersService.findAll()
  }

  @Get()
  @UseGuards(JwtAuthGuard)
  findOne(@JwtPayload('email') email: string): Promise<User> {
    return this.usersService.findOneByEmail(email)
  }

  @Get(':email')
  // @UseGuards(JwtAuthGuard)
  findOneEmail(@Param('email') email: string): Promise<User> {
    return this.usersService.findOneByEmail(email)
  }

  @Post('/test/')
  test(@Body() updateUserDto: UpdateUserDto): any {
    return this.usersService.test(updateUserDto)
  }

  @UseGuards(UserExistGuard)
  @Patch('/update/')
  updateProfile(
    @JwtPayload('email') email: string,
    @Body() updateUserDto: UpdateUserDto
  ): Promise<User> {
    return this.usersService.updateProfile(email, updateUserDto)
  }
}
