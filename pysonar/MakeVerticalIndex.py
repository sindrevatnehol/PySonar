# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 15:55:04 2018

@author: sindrev
"""
import tools
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np



def MakeVerticalIndex(ListOfFilesWithinTimeInterval,RemoveToCloseValues,R_s,res,directory2Data,dirnc,beamgrp): 

    
    DataMatrix = []
    
    
    #Loop through all files within the time interval
    for filename_index in range(0,len(ListOfFilesWithinTimeInterval[:,0])):
        
        
        #Print the progression
        tools.printProgressBar(filename_index + 1, len(ListOfFilesWithinTimeInterval[:,0]), prefix = 'Make SearchMatrix:', suffix = 'Complete', length = 50)
        
        
        #Get the full path name
        filename = os.path.join(dirnc,ListOfFilesWithinTimeInterval[filename_index,1])
        
        
        #Load the nc file
        fileID = Dataset(filename,'r',format = 'NETCDF4')
        
        
        #Get the vertical beam data
        #FIKS slik at den henter vertikal data
        variables = tools.GetVariablesFromNC(fileID,beamgrp,ListOfFilesWithinTimeInterval,filename_index)
        
        
            
            
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
        sv, RangeOut= tools.ApplyTVG(variables.BeamAmplitudeData,
                                variables.soundvelocity,
                                variables.sampleinterval,
                                variables.transmitpower,
                                variables.absorptioncoefficient,
                                variables.frequency,
                                variables.pulslength,
                                gain,
                                variables.equivalentbeamangle,
                                variables.sacorrection,
                                variables.dirx)
            
        plt.figure(1)
        plt.clf()
        plt.imshow(sv,aspect='auto')
        plt.draw()
        plt.pause(0.001)
        
        
        #Set data below 40 degree tilt to nan
        sv[np.where(variables.dirx <40)] = np.nan

        
        #Stack the vertical beam data
        DataMatrix = np.dtack((DataMatrix,sv[:,:,np.newaxis]))
        
        
        
        
    #Out of loop
    
    BinaryIDX = np.where(DataMatrix >= (np.nanmedian(DataMatrix,axis=2)*1.1))
    
    
    #Remove everything below a median
    DataMatrix[np.where(DataMatrix <= np.nanmedian(DataMatrix,axis=2)*1.1)] = np.nan



    #Plot result
    
    
    
    #Cluster the data into schools
    
    
    
    
    #Save polygon
    
    
    #Make Luf20 file
    
    