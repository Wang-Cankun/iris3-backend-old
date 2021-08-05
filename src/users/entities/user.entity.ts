import { Project } from 'src/project/entities/project.entity'
import {
  Entity,
  Column,
  PrimaryGeneratedColumn,
  JoinTable,
  ManyToMany,
  CreateDateColumn,
  PrimaryColumn,
  OneToMany,
  Generated,
  JoinColumn
} from 'typeorm'

export enum UserRole {
  ADMIN = 'admin',
  EDITOR = 'editor',
  GUEST = 'guest'
}

@Entity()
export class User {
  @PrimaryGeneratedColumn()
  id: number

  @PrimaryColumn({ unique: true })
  email: string

  @PrimaryColumn()
  @Generated('uuid')
  userId: string

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

  @Column({
    type: 'enum',
    enum: UserRole,
    default: UserRole.GUEST
  })
  role: UserRole

  @OneToMany(
    (type) => Project,
    (project) => project.user,
    {
      cascade: true
    }
  )
  projects: Project[]
}
