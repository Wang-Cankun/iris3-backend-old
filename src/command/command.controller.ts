import { Body, Controller, Post } from '@nestjs/common'
import { CommandService } from './command.service'

@Controller('command')
export class CommandController {
  constructor(private readonly commandService: CommandService) {}

  @Post('preprocess')
  public preprocess(@Body() params): any {
    console.log(params)
    // params = { nMitoGenes: '0.1', nPCs: '10', nResolution: '0.8' }
    return this.commandService.runPreprocessing(params)
  }
}
