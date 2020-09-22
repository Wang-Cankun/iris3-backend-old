import { NestFactory } from '@nestjs/core'
import { AppModule } from './app.module'
import { SwaggerModule, DocumentBuilder } from '@nestjs/swagger'
import * as express from 'express'
import { join } from 'path'
import * as compression from 'compression'
import { ConfigService } from '@nestjs/config'
async function bootstrap() {
  const app = await NestFactory.create(AppModule, { cors: true })

  // Get env configuration
  const configService = app.get(ConfigService)

  const port = configService.get('PORT')

  // Add global prefix on all routes
  app.setGlobalPrefix('/iris3/api')

  // Enable Express Compression
  app.use(compression())

  // Enbale OpenAI Swagger
  const swaggerOptions = new DocumentBuilder()
    .setTitle('IRIS3 API')
    .setDescription('IRIS3 API documentation')
    .setVersion('1.0')
    .build()
  const document = SwaggerModule.createDocument(app, swaggerOptions)
  SwaggerModule.setup('/iris3/api/docs', app, document)

  // Serve a static folder on ./tmp
  app.use(express.static(join(process.cwd(), '../tmp')))

  await app.listen(port)
  console.log(`Application is running on: ${await app.getUrl()}`)
}
bootstrap()
