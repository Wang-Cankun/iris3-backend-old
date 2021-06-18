import { User } from 'src/users/entities/user.entity'
import {
  Entity,
  PrimaryGeneratedColumn,
  Column,
  ManyToMany,
  CreateDateColumn
} from 'typeorm'

@Entity()
export class Job {
  @PrimaryGeneratedColumn()
  id: number

  @Column()
  jobid: string

  @Column()
  private: boolean

  @Column()
  creator: string

  @Column()
  email: string

  @Column()
  title: string

  @CreateDateColumn()
  createTime: string

  @Column()
  status: string

  @Column()
  expFile: string

  @Column()
  labelFile: string

  @Column()
  species: string

  @Column()
  tags: string

  @Column()
  description: string

  @ManyToMany(
    (type) => User,
    (user) => user.job
  )
  user: User[]
}
