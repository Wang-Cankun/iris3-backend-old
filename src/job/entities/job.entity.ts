import { Project } from '../../project/entities/project.entity'

import {
  Entity,
  PrimaryGeneratedColumn,
  Column,
  CreateDateColumn,
  ManyToOne,
  Generated,
  JoinTable,
  JoinColumn,
  UpdateDateColumn
} from 'typeorm'

export enum STATUS {
  SUCCESS = 'success',
  ERROR = 'error',
  RUNNING = 'running',
  QUEUEING = 'queueing'
}

@Entity()
export class Job {
  @PrimaryGeneratedColumn()
  id: number

  @Column()
  @Generated('uuid')
  jobId: string

  @CreateDateColumn()
  createTime: string

  @UpdateDateColumn()
  updateTime: string

  @Column({
    type: 'enum',
    enum: STATUS,
    default: STATUS.QUEUEING
  })
  status: STATUS

  @Column()
  title: string

  @Column()
  type: string

  @Column()
  description: string

  @JoinColumn()
  @ManyToOne(
    (type) => Project,
    (project) => project.jobs,
    {
      cascade: true
    }
  )
  project: Project
}
