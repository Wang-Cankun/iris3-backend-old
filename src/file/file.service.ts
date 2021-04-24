import { Injectable } from '@nestjs/common'
import { InjectRepository } from '@nestjs/typeorm'
import { Repository } from 'typeorm'
import { UploadFileDto } from './dto/upload-file.dto'
import { File } from './entities/file.entity'

@Injectable()
export class FileService {
  constructor(
    @InjectRepository(File)
    private readonly fileRepo: Repository<File>
  ) {}

  async create(uploadFileDto: UploadFileDto): Promise<File> {
    const file = this.fileRepo.create({
      ...uploadFileDto
    })
    return this.fileRepo.save(file)
  }

  async findAll(jobid: string): Promise<File[]> {
    return await this.fileRepo.find({
      jobid: jobid
    })
  }
}
