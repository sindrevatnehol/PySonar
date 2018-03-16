'''
StringArray(variabel).py

Last modify date: 24/11 2016

Usage: 
    StringArray(variabel)
        
        
Description: 
a function frequently used when converting .raw to .nc
made here to save space in another file.


    
    
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

def StringArray(variabel): 
    output = "" 
    for variabels in variabel:
        output = output+chr(variabels)
    return output; 

    
    
    
if __name__ == '__main__':    
    StringArray()
    
    
    