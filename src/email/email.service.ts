import { Injectable } from '@nestjs/common'
import startJobEmail from './template/start-job'
import forgotPasswordEmail from './template/forgot-password'
import createUserEmail from './template/create-user'
import * as SibApiV3Sdk from 'sib-api-v3-sdk'
@Injectable()
export class EmailService {
  hostName: string = process.env.FRONTEND_HOST

  constructor() {
    SibApiV3Sdk.ApiClient.instance.authentications['api-key'].apiKey =
      process.env.EMAIL_TOKEN
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types
  sendForgotPasswordEmail(email, token) {
    const content = forgotPasswordEmail(email, token)
    const subject = 'DeepMAPS reset password'
    this.sendEmail(email, subject, content)
  }

  sendStartJobEmail(email: string, jobid: string): void {
    const content = startJobEmail(email, jobid)
    const subject = 'DeepMAPS job created'
    this.sendEmail(email, subject, content)
  }

  sendCreateUserEmail(email: string): void {
    const content = createUserEmail(email)
    const subject = 'Welcome to DeepMAPS'
    this.sendEmail(email, subject, content)
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types
  sendEmail(email, subject, content): void {
    new SibApiV3Sdk.TransactionalEmailsApi()
      .sendTransacEmail({
        subject: subject,
        sender: { email: 'no-reply@bmbls.bmi.osumc.edu', name: 'DeepMAPS' },
        bcc: [{ email: 'cankun.wang@osumc.edu', name: 'cankun' }],
        replyTo: {
          email: 'flykun0620@gmail.com',
          name: 'flykun0620@gmail.com'
        },
        to: [{ name: email, email: email }],
        htmlContent: content,
        params: { bodyMessage: '' }
      })
      .then(
        function(data) {
          console.log(data)
        },
        function(error) {
          console.error(error)
        }
      )
  }
}
