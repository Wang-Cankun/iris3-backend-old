import { ContainerName } from './container-name.dto'

export class RunContainerDto {
  readonly image: ContainerName
  readonly port: string
}
