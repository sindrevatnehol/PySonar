'''
TimeConverter.py

Last modify date: 24/11 2016

Usage: 
    TimeConverter(time_date): 
        - time_date - simrad time
        
    Output: 
        - full-date
        
        
Description: 
This function convert the time from the .raw file to actual utc time


    
    
Author: 
Dr. Sindre Vatnehol (PhD)
Institute of Marine Research, Norway

Mail: 
sindre.vatnehol@imr.no

Tlf.: 
+47 900 79 376


Project: 
REDUS (Reducing uncertainty in stock assessment)  
'''


import datetime

def TimeConverter(time_date):
    
    temp_date = (time_date[1]*2**32 + time_date[0])/10000
    fulldate = datetime.datetime.strptime('1601-01-01 00:00:00.000',"%Y-%m-%d %H:%M:%S.%f")
    
    
    fulldate = fulldate + datetime.timedelta(milliseconds=int(temp_date))
    return fulldate; 
    

if __name__ == '__main__':    
    TimeConverter()