'''
ApplyTVG.py

Last modify date: 24/11 2016

Usage: 
    fread(fid, nelements, dtype):
        
Description: 
This modual does the same as the matlab function with the same name

    
    
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


import numpy as np


def fread(fid, nelements, dtype):
    
    
    if dtype is np.str:
        dt = np.uint8
    else: 
        dt = dtype 
        
    if (isinstance( nelements, int )== True)or (isinstance( nelements, np.int32 )== True):
        print(fid)
        data_array = np.fromfile(fid, dt, nelements)
    else: 
        data_array = np.fromfile(fid, dt, nelements[0])
        
    return data_array; 
    
    
if __name__ == '__main__':    
    fread(fid, nelements, dtype)
    