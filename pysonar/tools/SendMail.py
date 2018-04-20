from email.mime.text import MIMEText
from datetime import date
import smtplib



def send_email(data): 
    SMTP_server = 'smtp.gmail.com'
    SMTP_port = 587
    SMTP_username = 'ReportStatusOfAnalysis@gmail.com'
    SMTP_psw = 'HurraForDegSomFyller100'

    email_to = ['sindre.vatnehol@gmail.com','sindre.vatnehol@hi.no']
    email_from = 'sindre.vatnehol@gmail.com'
    email_subject = 'Report email:'
    email_space = ', '

    #data = 'Dette er en test' 


    msg=MIMEText(data)
    msg['Subject']=email_subject
    msg['To'] = email_space.join(email_to)
    msg['From'] = email_from
    mail = smtplib.SMTP(SMTP_server,SMTP_port) 
    mail.starttls()
    mail.login(SMTP_username,SMTP_psw) 
    mail.sendmail(email_from,email_to,msg.as_string())
    mail.quit()

    print('email send')

if __name__=='__main__': 
    send_email()

