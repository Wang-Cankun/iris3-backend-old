import {
  Entity,
  PrimaryGeneratedColumn,
  Column,
  ManyToMany,
  CreateDateColumn
} from 'typeorm'
import { Job } from '../../job/entities/job.entity'

@Entity()
export class File {
  @PrimaryGeneratedColumn()
  id: number

  @Column()
  jobid: string

  @Column()
  private: boolean

  @Column()
  creator: string

  @Column()
  tags: string

  @Column({ default: 1 })
  index: number

  @Column({ default: 'active' })
  status: string

  @Column()
  species: string

  @CreateDateColumn()
  createTime: string

  @Column()
  fieldname: string

  @Column()
  originalname: string

  @Column()
  encoding: string

  @Column()
  mimetype: string

  @Column()
  destination: string

  @Column()
  filename: string
  @Column()
  title: string
  @Column()
  description: string

  @Column()
  path: string

  @Column()
  size: number

  @ManyToMany(
    (type) => Job,
    (job) => job.jobid
  )
  job: Job[]
}
