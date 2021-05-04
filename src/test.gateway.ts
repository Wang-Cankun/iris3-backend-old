import { UseInterceptors } from '@nestjs/common'
import {
  MessageBody,
  SubscribeMessage,
  WebSocketGateway
} from '@nestjs/websockets'
import { from, Observable } from 'rxjs'
import { map } from 'rxjs/operators'

import { RedisPropagatorInterceptor } from './shared/redis-propagator/redis-propagator.interceptor'

@UseInterceptors(RedisPropagatorInterceptor)
@WebSocketGateway()
export class EventsGateway {
  @SubscribeMessage('events')
  public findAll(): Observable<any> {
    return from([1, 2, 30, 50]).pipe(
      map((item) => {
        return { event: 'events', data: item }
      })
    )
  }

  @SubscribeMessage('jobProgress')
  public returnStatus(@MessageBody() item: string): any {
    return item
  }

  @SubscribeMessage('')
  public hellp(@MessageBody() item: string): any {
    return { event: 'progress', data: 1 }
  }
}
