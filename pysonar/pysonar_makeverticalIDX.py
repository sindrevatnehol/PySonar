# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 11:39:22 2018

@author: sindrev
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 15:55:04 2018

@author: sindrev
"""
from tools import tools
import os
from netCDF4 import Dataset
#import matplotlib.pyplot as plt
import numpy as np

import scipy.io as scp
from math import radians, cos, sin, asin, sqrt

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371*1000 # Radius of earth in meters. Use 3956 for miles
    return c * r
    
    
    
    
    
    
    
    

def MakeVerticalIndex(ListOfFilesWithinTimeInterval,RemoveToCloseValues,
                      R_s,res,directory2Data,dirnc,beamgrp,start_time,
                      log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop,ping_idx): 

    
    
    #For bookkeeping
    transectID = log_start
    DataMatrix = []
    time0 = 0
    ping_counter=0
    r_mask = 0
    b_mask = 0
    ping_mask=0
    
    
    
    
    #Loop through all files within the time interval
    for filename_index in range(0,len(ListOfFilesWithinTimeInterval)):
        
        #Print the progression
        tools.printProgressBar(filename_index + 1, len(ListOfFilesWithinTimeInterval), prefix = 'Make Vertical:', suffix = 'Completed', length = 50)
        
        
        #Get the full path name
        filename = os.path.join(dirnc,ListOfFilesWithinTimeInterval[filename_index])
        
        
        #Load the nc file
        fileID = Dataset(filename,'r',format = 'NETCDF4')
        
        
        #Get the group with the vertical fan data
        if fileID.groups['Sonar'].groups['Beam_group1'].beam_mode == 'Vertical': 
            beamgrp = 'Beam_group1'
        elif fileID.groups['Sonar'].groups['Beam_group2'].beam_mode == 'Vertical':
            beamgrp = 'Beam_group2'
        
            
            
        #Get the vertical beam data
        #FIKS slik at den henter vertikal data
        variables = tools.GetVariablesFromNC(fileID,beamgrp,ListOfFilesWithinTimeInterval,ping_idx[filename_index])
        fileID.close()
            
            
#        try: 
        #The sonar data often includes corrputed value of 
        #the transmit power that destroys the analysis. 
        #This will fix this problum, but the sv values are not
        #correct. 
        #ADD this data should be labeled when making the work files
        #so the user can now that it is corrupted.
        if variables.transmitpower == 0: 
            variables.transmitpower = 4633
            


            
            
        #Compute the sv and TS 
        #ADD TS are not used here !!!
        sv, RangeOut= tools.ApplyTVG(10*np.log10(variables.BeamAmplitudeData),
                                variables.soundvelocity[0],
                                variables.sampleinterval,
                                variables.transmitpower,
                                variables.absorptioncoefficient[0],
                                variables.frequency,
                                variables.pulslength,
                                variables.gaintx,
                                variables.equivalentbeamangle,
                                variables.sacorrection,
                                variables.dirx)
            
        
        #Remove data too close to the vessel
        sv[np.where(RangeOut<=RemoveToCloseValues)] = np.nan
           
       
        sv[:,np.where(variables.dirx>=60)] = np.nan
#            sv[np.where(sv<-65)] = np.nan
           
           
        #Stack the vertical beam data
        if DataMatrix == []:
            sv_x,sv_y = sv.shape
            
            DataMatrix = 10**(sv/10)[:,:,np.newaxis]
            time0 = variables.time
            ping_counter = ping_counter+1
            
        else: 
            if sv_x>sv.shape[0]: 
                sv=np.append(sv,np.nan*np.ones((sv_x-sv.shape[0],sv_y)),axis=0)
            elif sv_x<sv.shape[0]: 
                DataMatrix=np.append(DataMatrix,np.nan*np.ones((sv.shape[0]-sv_x,sv_y,DataMatrix.shape[2])),axis=0)
                sv_x,sv_y = sv.shape
                
            DataMatrix = np.dstack((DataMatrix,10**(sv/10)[:,:,np.newaxis]))
            time0=np.hstack((time0,variables.time))
            ping_counter = ping_counter+1
            
#        except AttributeError:
#            print('* bad file')
        
            
            
            
    #Find the median
    Medianen = np.nanmedian(DataMatrix,axis=2)
    
    
    #Get index of everything above the median
    [r_mask0,b_mask0,ping_mask0] = np.where(DataMatrix>=4*np.repeat(Medianen[:,:,np.newaxis],len(DataMatrix[0,0,:]),axis=2))
    
    
    #
    if ping_mask0 != []: 
        ping_mask0 = time0[ping_mask0]


    #Stack the data
    r_mask = np.hstack((r_mask,r_mask0))
    b_mask = np.hstack((b_mask,b_mask0))
    ping_mask = np.hstack((ping_mask,ping_mask0))

    
    DataMatrix = []
    ping_counter = 0
        
        
    #Try saving the data
    try: 
        scp.savemat(directory2Data.dir_work+'/'+'Vertical_T'+str(transectID)+'.mat',mdict = {'r_mask':r_mask,'b_mask':b_mask,'ping_mask':ping_mask,
                    'start_time':start_time,'log_start':log_start,'stop_time':stop_time,'lat_start':lat_start,'lat_stop':lat_stop,'lon_start':lon_start,'lon_stop':lon_stop})
    except TypeError: 
        print(start_time,log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop)
    
    