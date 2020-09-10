import {
  Entity,
  Column,
  PrimaryGeneratedColumn,
  JoinTable,
  ManyToMany,
  CreateDateColumn
} from 'typeorm'
import { Job } from './job.entity'

@Entity()
export class User {
  @PrimaryGeneratedColumn()
  id: number

  @Column({ unique: true })
  email: string

  // { select: false }
  @Column()
  password: string

  @Column()
  firstName: string

  @Column()
  lastName: string

  @Column()
  institution: string

  @Column()
  newsletter: boolean

  @CreateDateColumn()
  createTime: string

  @Column({ default: true })
  isActive: boolean

  @JoinTable()
  @ManyToMany(
    () => Job,
    (job) => job.user,
    {
      cascade: true
    }
  )
  job: Job[]
}
