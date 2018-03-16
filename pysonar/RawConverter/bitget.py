'''
ApplyTVG.py

Last modify date: 24/11 2016

Usage: 
    bitget(byteval,idx) 
        
Description: 
Similar as the matlab function with the same name


    
    
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

import math

    
    
def bitget(byteval,idx):
    if byteval == 0:
        return 0; 
    else: 
        if (int(math.log(byteval,2))+1)==idx:
            return 1;
        else:
            return 0;
    
        
if __name__ == '__main__': 
    bitget()