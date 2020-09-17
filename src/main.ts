import { NestFactory } from '@nestjs/core'
import { AppModule } from './app.module'
import { SwaggerModule, DocumentBuilder } from '@nestjs/swagger'

async function bootstrap() {
  const app = await NestFactory.create(AppModule, { cors: true })
  app.setGlobalPrefix('/iris3/api')
  const options = new DocumentBuilder()
    .setTitle('IRIS3 API')
    .setDescription('IRIS3 API application')
    .setVersion('1.0')
    .build()

  const document = SwaggerModule.createDocument(app, options)
  SwaggerModule.setup('/iris3/api/docs', app, document)
  await app.listen(9005)
}
bootstrap()
