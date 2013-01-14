"""
This module is for writing email through emacs.
Nothing more or less than that. It was decided 
after noticing that I'm so confused over using
Google Chrome to write stuffs. I really miss
all those fancy shortkey stuffs which makes me
work without having to put my hands on mouse.

Just for fun! 
 
               Youngung Jeong
               Materials Mechanics Laboratory
               GIFT-POSTECH,KOREA
"""


import smtplib #sending function
from email.mime.text import MIMEText

"""
msg = MIMEText(msg)
msg['Subject'] = 
"""

class compose:
    def __init__(self,mode='AUTO', 
                 msg=None, 
                 subj=None, 
                 From='youngung.jeong@gmail.com', 
                 To=None):


        ## Initial module 
        if msg!=None: self.msg = MIMEText(msg)
        else:         self.msg = MIMEText(self.__main__())
        if subj!=None: self.msg['Subject'] = subj
        else:          self.msg['Subject'] = self.__subject__()
        if From!=None: self.msg['From'] = msg
        else:          self.msg['From'] = self.__from__()
        if To!=None: self.msg['To'] = To
        else:        self.msg['To'] = self.__to__()


        ## Shoots on the screen
        print "*****************"
        print "information block"
        print "*****************"
        print self.msg.as_string()

        if mode=='auto' or mode=='AUTO':
            self.server()
            self.__login__(password=raw_input('type password >> '))
            self.__send__()
            self.__quit__()

    def server(self,):
        ## SMTP 
        self.s = smtplib.SMTP('smtp.gmail.com',587)
        self.s.ehlo()
        self.s.starttls()
        self.s.ehlo()

    def __login__(self,username='youngung.jeong',password=None):
        ## login --
        self.s.login(username, password)

    def __send__(self,):
        ###
        if type(self.msg['To']).__name__=='list':
            you =self.msg['To']
        elif type(self.msg['To']).__name__=='str':
            you = [self.msg['To']]
        self.s.sendmail(self.msg['To'], you, self.msg.as_string())

    def __quit__(self,):
        ###
        self.s.quit()


        


    def __from__(self,mode='fixed'):
        """ Returns the sender name"""
        if mode=='man': return raw_input("Type the sender's address>> ")
        elif mode=='fixed': return "youngung.jeong@gmail.com"



    def __to__(self,mode='man'):
        """ Returns the recipient  """
        if mode=='man':
            return raw_input('Please type the recipients email address >> ')
        else :  # different mode will be included in this module
            #plane modes : FILE, mans (multiple manual)
            pass


    def __subject__(self,):
        """ Subrject of the mail"""
        return raw_input("Please type the subject of the email >> ")


    def __main__(self,):
        """ Getting the message """
        return raw_input("type the message that you'd like to send >> ")
        pass

