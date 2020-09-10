import {
  Body,
  Controller,
  Delete,
  Get,
  Param,
  Patch,
  Post,
  UseGuards,
  ParseIntPipe
} from '@nestjs/common'
import { UsersService } from './users.service'
import { CreateUserDto } from './dto/create-user.dto'
import { UpdateUserDto } from './dto/update-user.dto'
import { ApiTags } from '@nestjs/swagger'
import { JwtAuthGuard } from 'src/auth/jwt-auth.guard'
import { User } from './entities/user.entity'
import { UserExistGuard } from './guards/user-exist.guard'

@ApiTags('users')
@Controller('users')
export class UsersController {
  constructor(private readonly usersService: UsersService) {}

  @Get()
  async findAll(): Promise<User[]> {
    // const { limit, offset } = paginationQuery;
    // await new Promise((resolve) => setTimeout(resolve, 500))
    return this.usersService.findAll()
  }

  @UseGuards(JwtAuthGuard)
  @Get(':id')
  findOne(@Param('id', ParseIntPipe) id: string): Promise<User> {
    return this.usersService.findOneById(id)
  }

  @UseGuards(UserExistGuard)
  @Post()
  create(@Body() createUserDto: CreateUserDto): Promise<User> {
    return this.usersService.create(createUserDto)
  }

  @Post('/reset/:email')
  sendResetEmail(@Param('email') email: string): void {
    return this.usersService.sendEmailVerification(email)
  }

  @Post('/test/')
  test(@Body() updateUserDto: UpdateUserDto): any {
    return this.usersService.test(updateUserDto)
  }

  @UseGuards(JwtAuthGuard)
  @UseGuards(UserExistGuard)
  @Patch('/password/:id')
  updatePassword(
    @Param('id', ParseIntPipe) id: string,
    @Body() updateUserDto: UpdateUserDto
  ): Promise<User> {
    return this.usersService.updatePassword(id, updateUserDto)
  }

  @UseGuards(JwtAuthGuard)
  @UseGuards(UserExistGuard)
  @Patch(':id')
  updateProfile(
    @Param('id', ParseIntPipe) id: string,
    @Body() updateUserDto: UpdateUserDto
  ): Promise<User> {
    return this.usersService.updateProfile(id, updateUserDto)
  }

  @UseGuards(JwtAuthGuard)
  @UseGuards(UserExistGuard)
  @Delete(':email')
  remove(@Param('email') email: string): Promise<User> {
    return this.usersService.remove(email)
  }
}
