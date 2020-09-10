import { Injectable } from '@nestjs/common'

import * as nodemailer from 'nodemailer'
@Injectable()
export class EmailService {
  hostName: string = process.env.FRONTEND_HOST

  constructor() {}
  public sendForgotPasswordEmail(email: string, token: string): void {
    const resetPasswordLink =
      this.hostName + 'login/forgot/reset/?token=' + token + '&email=' + email
    const content =
      `<table class='body' style='border-collapse: separate; mso-table-lspace: 0pt; mso-table-rspace: 0pt; width: 100%; background-color: #f6f6f6;' border='0' cellspacing='0' cellpadding='0'>
  <tbody> 
  <tr>
  <td style='font-family: sans-serif; font-size: 14px; vertical-align: top;'>&nbsp;</td>
  <td class='container' style='font-family: sans-serif; font-size: 14px; vertical-align: top; display: block; margin: 0 auto; max-width: 700px; padding: 10px; width: 700px;'>
  <div class='content' style='box-sizing: border-box; display: block; margin: 0 auto; max-width: 700px; padding: 10px;'><!-- START CENTERED WHITE CONTAINER --> <span class='preheader' style='color: transparent; display: none; height: 0; max-height: 0; max-width: 0; opacity: 0; overflow: hidden; mso-hide: all; visibility: hidden; width: 0;'><g class='gr_ gr_65 gr-alert gr_spell gr_inline_cards gr_run_anim ContextualSpelling ins-del multiReplace' id='65' data-gr-id='65'></g></span>
  <table class='main' style='border-collapse: separate; mso-table-lspace: 0pt; mso-table-rspace: 0pt; width: 100%; background: #ffffff; border-radius: 3px;'><!-- START MAIN CONTENT AREA -->
  <tbody>
  <tr>
  <td class='wrapper' style='font-family: sans-serif; font-size: 14px; vertical-align: top; box-sizing: border-box; padding: 20px;'><br />
  <table style='border-collapse: separate; mso-table-lspace: 0pt; mso-table-rspace: 0pt; width: 100%;' border='0' cellspacing='0' cellpadding='0'>
  <tbody> 
  <tr>
  <td style='font-family: sans-serif; font-size: 14px; vertical-align: top;'>
  <table border='0' width='100%' cellspacing='0' cellpadding='0'>
  <tbody>
  <tr>
  <td class='subhead'>Hello,</td>
  </tr>
  <tr>
  <td class='h1' style='padding: 5px 0 0 0;'><br />
  <div>
  <div><span>You have submitted a request to reset your password.<br /></span><span>Your email: ` +
      email +
      `</span></div>
  </div>
  </div>
  </td>
  </tr>
  </tbody>
  </table>
  <p style='font-family: sans-serif; font-size: 14px; font-weight: normal; margin: 0; margin-bottom: 15px;'></p>
  <span>&nbsp;</span></td>
  </tr>
  <tr><td style='font-family: sans-serif; font-size: 14px; vertical-align: top; background-color: #3498db; border-radius: 5px; text-align: center;'><a style='display: inline-block; color: #ffffff; background-color: #3498db; border: solid 1px #3498db; border-radius: 5px; box-sizing: border-box; cursor: pointer; text-decoration: none; font-size: 14px; font-weight: bold; margin: 0; padding: 12px 25px; text-transform: capitalize; border-color: #3498db;' href='` +
      resetPasswordLink +
      `' target='_blank' rel='noopener'>Click to reset your password</a></td></tr>
  </tbody>
  </table>
  </td>
  </tr>
  </tbody>
  </table>
  </td>
  </tr>
  <!-- END MAIN CONTENT AREA --></tbody>
  </table>
  <!-- START FOOTER -->
  <div class='footer' style='clear: both; margin-top: 10px; text-align: center; width: 100%;'>
  <table style='border-collapse: separate; mso-table-lspace: 0pt; mso-table-rspace: 0pt; width: 100%;' border='0' cellspacing='0' cellpadding='0'>
  <tbody>
  <tr>
  <td class='content-block' style='font-family: sans-serif; vertical-align: top; padding-bottom: 10px; padding-top: 10px; font-size: 12px; color: #999999; text-align: center;'><span>Copyright 2020 &copy; </span><a href='http://u.osu.edu/bmbl' target='_blank' rel='noopener'>BMBL</a><span>, </span><a href='https://medicine.osu.edu/departments/biomedical-informatics/' target='_blank' rel='noopener'>OSU</a><span>. All rights reserved. </span></td>
  </tr>
  <tr>
  <td class='content-block powered-by' style='font-family: sans-serif; vertical-align: top; padding-bottom: 10px; padding-top: 10px; font-size: 12px; color: #999999; text-align: center;'><a href='mailto:qin.ma\@osumc.edu' title='qin.ma\@osumc.edu'>Contact us: qin.ma\@osumc.edu</a><span> </span></td>
  </tr>
  </tbody>
  </table>
  </div>
   <!-- END CENTERED WHITE CONTAINER --></div>
  </td>
  <td style='font-family: sans-serif; font-size: 14px; vertical-align: top;'>&nbsp;</td>
  </tr>
  </tbody>
  </table>`

    const transporter = nodemailer.createTransport({
      service: 'Gmail',
      auth: {
        user: process.env.EMAIL_NAME,
        pass: process.env.EMAIL_PASSWORD
      }
    })

    const mailOptions = {
      from: 'IRIS3 <no-reply@bmbls.bmi.osumc.edu>',
      to: email,
      subject: 'Reset password from IRIS3',
      html: content
    }

    transporter.sendMail(mailOptions, function(error, info) {
      if (error) {
        console.log(error)
      } else {
        console.log('Email sent: ' + info.response)
      }
    })
  }
}
