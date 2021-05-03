import { NestFactory } from '@nestjs/core'
import { AppModule } from './app.module'
import { SwaggerModule, DocumentBuilder } from '@nestjs/swagger'
import * as compression from 'compression'
import { ConfigService } from '@nestjs/config'
import { initAdapters } from './adapters.init'

async function bootstrap(): Promise<void> {
  const app = await NestFactory.create(AppModule, { cors: true })

  // Get env configuration
  const configService = app.get(ConfigService)

  const port = configService.get('PORT')

  initAdapters(app)
  // Add global prefix on all routes
  app.setGlobalPrefix('/deepmaps/api')

  // Enable Express Compression
  app.use(compression())

  // Enbale OpenAI Swagger
  const swaggerOptions = new DocumentBuilder()
    .setTitle('deepmaps API')
    .setDescription('deepmaps API documentation')
    .setVersion('1.0')
    .build()
  const document = SwaggerModule.createDocument(app, swaggerOptions)
  SwaggerModule.setup('/deepmaps/api/docs', app, document)

  await app.listen(port)
  console.log(`Application is running on: ${await app.getUrl()}`)
}
bootstrap()
