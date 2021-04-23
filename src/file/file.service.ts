import { Injectable } from '@nestjs/common'
import { InjectRepository } from '@nestjs/typeorm'
import { Repository } from 'typeorm'
import { UploadFileDto } from './dto/upload-file.dto'
import { File } from './entities/file.entity'

@Injectable()
export class FileService {
  constructor(
    @InjectRepository(File)
    private readonly fileRepository: Repository<File>
  ) {}

  async create(uploadFileDto: UploadFileDto): Promise<File> {
    const file = this.fileRepository.create({
      ...uploadFileDto
    })
    return this.fileRepository.save(file)
  }

  findAll() {
    return this.fileRepository.find({
      relations: ['jobid']
    })
  }
}
