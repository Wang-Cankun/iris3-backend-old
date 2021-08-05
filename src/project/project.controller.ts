import { Controller, Get, Post, Body, Put, Param, Delete } from '@nestjs/common'
import * as md5 from 'md5'
import { CreateProjectDto } from './dto/create-project.dto'
import { UpdateProjectDto } from './dto/update-project.dto'
import { ProjectService } from './project.service'

@Controller('project')
export class ProjectController {
  constructor(private readonly projectService: ProjectService) {}

  @Post()
  create(@Body() createProjectDto: CreateProjectDto) {
    return this.projectService.create(createProjectDto)
  }

  @Post('jobs')
  listJobs(@Body() projectUid: any) {
    console.log(projectUid)
    console.log(md5(JSON.stringify(projectUid)))
    return this.projectService.listJobsById(projectUid.projectUid)
  }

  @Get()
  findAll() {
    return this.projectService.findAll()
  }

  @Get(':id')
  findOne(@Param('id') id: string) {
    return this.projectService.findOne(+id)
  }

  @Put(':id')
  update(@Param('id') id: string, @Body() updateProjectDto: UpdateProjectDto) {
    return this.projectService.update(+id, updateProjectDto)
  }

  @Delete(':id')
  remove(@Param('id') id: string) {
    return this.projectService.remove(+id)
  }
}
