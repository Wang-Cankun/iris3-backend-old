import {
  Body,
  Controller,
  Delete,
  Get,
  Param,
  Patch,
  Post,
  UseGuards
} from '@nestjs/common'
import { UsersService } from './users.service'
import { CreateUserDto } from './dto/create-user.dto'
import { UpdateUserDto } from './dto/update-user.dto'
import { ApiTags } from '@nestjs/swagger'
import { JwtAuthGuard } from 'src/auth/jwt-auth.guard'

@ApiTags('users')
@Controller('users')
export class UsersController {
  constructor(private readonly usersService: UsersService) {}

  @Get()
  async findAll() {
    // const { limit, offset } = paginationQuery;
    // await new Promise((resolve) => setTimeout(resolve, 500))
    return this.usersService.findAll()
  }

  @UseGuards(JwtAuthGuard)
  @Get(':email')
  findOne(@Param('email') email: string) {
    return this.usersService.findOne(email)
  }

  @Post()
  create(@Body() createUserDto: CreateUserDto) {
    return this.usersService.create(createUserDto)
  }

  @Post('/reset/:email')
  sendResetEmail(@Param('email') email: string) {
    return this.usersService.sendEmailVerification(email)
  }

  @UseGuards(JwtAuthGuard)
  @Patch('/password/:id')
  updatePassword(
    @Param('id') id: string,
    @Body() updateUserDto: UpdateUserDto
  ) {
    return this.usersService.updatePassword(id, updateUserDto)
  }

  @UseGuards(JwtAuthGuard)
  @Patch(':id')
  updateProfile(@Param('id') id: string, @Body() updateUserDto: UpdateUserDto) {
    return this.usersService.updateProfile(id, updateUserDto)
  }

  @UseGuards(JwtAuthGuard)
  @Delete(':email')
  remove(@Param('email') email: string) {
    return this.usersService.remove(email)
  }
}
