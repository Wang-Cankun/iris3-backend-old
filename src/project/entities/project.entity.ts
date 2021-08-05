import { User } from 'src/users/entities/user.entity'
import { Job } from '../../job/entities/job.entity'
import {
  Entity,
  PrimaryGeneratedColumn,
  Column,
  CreateDateColumn,
  OneToMany,
  ManyToOne,
  UpdateDateColumn,
  Generated,
  JoinTable,
  JoinColumn
} from 'typeorm'

@Entity()
export class Project {
  @PrimaryGeneratedColumn()
  id: number

  @Column()
  @Generated('uuid')
  projectUid: string

  @Column()
  private: boolean

  @Column()
  title: string

  @CreateDateColumn()
  createTime: string

  @UpdateDateColumn()
  updateTime: string

  @Column()
  status: string

  @Column('simple-array')
  tag: string[]

  @Column()
  description: string

  @OneToMany(
    (type) => Job,
    (job) => job.project
  )
  jobs: Job[]

  @ManyToOne(
    (type) => User,
    (user) => user.projects
  )
  @JoinColumn()
  user: User
}
