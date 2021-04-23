import { Body, Controller, Post } from '@nestjs/common'
import { CommandService } from './command.service'
import {
  MessageBody,
  SubscribeMessage,
  WebSocketGateway,
  WebSocketServer,
  WsResponse
} from '@nestjs/websockets'
import { from, Observable } from 'rxjs'
import { map } from 'rxjs/operators'
import { Server } from 'socket.io'

@WebSocketGateway()
export class CommandController {
  constructor(private readonly commandService: CommandService) {}

  @WebSocketServer()
  server: Server

  @SubscribeMessage('preprocess')
  async preprocess(@Body() params) {
    console.log(params)
    this.commandService.runPreprocessing()
    return new Promise((resolve) => {
      setTimeout(() => {
        resolve(Date())
      }, 4500)
    })
  }

  @SubscribeMessage('submit')
  async handleEvent(@MessageBody() data: any) {
    return new Promise((resolve) => {
      setTimeout(() => {
        resolve(Date())
      }, 4000)
    })
  }
}
