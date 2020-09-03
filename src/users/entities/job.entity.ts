import { Entity, PrimaryGeneratedColumn, Column, ManyToMany } from 'typeorm'
import { User } from './user.entity'

@Entity()
export class Job {
  @PrimaryGeneratedColumn()
  id: number

  @Column()
  email: string

  @Column()
  name: string

  @Column()
  createTime: string

  @Column()
  status: string

  @Column()
  isPrivate: boolean

  @Column()
  size: string

  @ManyToMany(
    (type) => User,
    (user) => user.job
  )
  user: User[]
}
