# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 09:58:25 2018

@author: sindrev
"""


import scipy.io as scp
import os
import time, datetime
import numpy as np
from tools import tools
from xml.etree import ElementTree as ET
import random




from netCDF4 import Dataset





def SecondsBetweenTransect(startTime,unixtime):

    fulldate = datetime.datetime.strptime('1601-01-01 00:00:00.000',"%Y-%m-%d %H:%M:%S.%f")

    starten =  datetime.datetime.strptime(startTime,"%Y%m%d%H%M%S")

    fulldate = (starten-fulldate).total_seconds() #datetime.timedelta(milliseconds=int(startTime))

    seconds = fulldate-unixtime/10000000

    return seconds

    
    
    
    
    
    

    
    
    
    

def MakeHorizontalLuf20(CompleteListOfFiles,directory2Data,start_time,log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop, TimeIDX
                        ,nation,cruice_id,platform): 
    
    
    sonar_draft = 8
    Threshold = -75
    
    
    #Fiksed stuff
    integrator_dist = 1
    pel_ch_thickness = 1
    NumChannels = 300
#    ChannelSize = 350
    
    
    
                
    #Start information in the text
    root = ET.Element("echosounder_dataset")
    ET.SubElement(root,'report_time').text = str(datetime.datetime.fromtimestamp((time.time())).strftime('%Y-%m-%d %H:%M:%S'))
    ET.SubElement(root,'lsss_version').text = "PNMDformats v 0.1 - vertical"
    ET.SubElement(root,'nation').text = nation
    ET.SubElement(root,'platform').text = platform
    ET.SubElement(root,'cruise').text = cruice_id

    


    #Make the distance list hirarchy
    distance_list = ET.SubElement(root,'distance_list')
                
                
    
    
    #Get list of work files
    Listoffiles = os.listdir(directory2Data.dir_work)
            
            
            
    #Go through each log distance
    for i in range(len(start_time)): 
        
        
        #Make empty vector to start preparing for stacking
        SA = np.zeros((NumChannels,1))*np.nan
        
        
        #Convert the start and stop time to correct format
        stime = stop_time[i][:4]+'-'+start_time[i][4:6]+'-'+start_time[i][6:8]+' '+start_time[i][9:11]+':'+start_time[i][11:13]+':'+start_time[i][13:]
        etime = stop_time[i][:4]+'-'+stop_time[i][4:6]+'-'+stop_time[i][6:8]+' '+stop_time[i][9:11]+':'+stop_time[i][11:13]+':'+stop_time[i][13:]
        



        #Print the distance list information to xml
        distance = tools.addLogDistanceInfo(distance_list,lat_start[i],lon_start[i],
                                      stime,lat_stop[i],lon_stop[i],etime,integrator_dist,pel_ch_thickness,log_start[i])
                    
                    
                    
        
        #Print the frequecy level information to xlm
        sa_by_acocat = tools.addFrequencyLevelInfo(distance,26000,0,NumChannels,Threshold)
                            
       
                    
        
        
        for iii in range(CompleteListOfFiles.shape[0]): 
            
            ein = SecondsBetweenTransect(start_time[i].replace('T',''),int(CompleteListOfFiles[iii][0]))
            to = SecondsBetweenTransect(stop_time[i].replace('T',''),int(CompleteListOfFiles[iii][0]))
            
            
            
            
            
            
            if ein <=0 and to>=0: 
                tools.printProgressBar(iii,CompleteListOfFiles.shape[0],prefix = str(CompleteListOfFiles.shape[0]-iii)+' Status for HLuf20', suffix = 'Completed', length = 50)
                
                
                
            
                for file in Listoffiles: 
                    if 'Horizontal' in file: 
                        Mat = scp.loadmat(directory2Data.dir_work+'\\'+file)
        
                        WorkFile = Mat['WorkFile']
                        WorkIDX = Mat['WorkIDX'][0]
                        WorkTime = Mat['WorkTime'][0]
                        WorkPhi = Mat['WorkPhi'][0]
                        WorkBeam = Mat['WorkBeam'][0]


                        pingIDX = np.where(WorkTime==CompleteListOfFiles[iii][0])[0]


                        if len(pingIDX) >0: 
                            
                        
                            #Get the file name of the idx file
                            filname = directory2Data.dir_rawdata+'/'+np.unique(WorkFile[pingIDX])[0]
                
                            #Load file
                            fileID = Dataset(filname,'r',format = 'NETCDF4')
                
        
                            
                            #Get the group with the vertical fan data
                            if fileID.groups['Sonar'].groups['Beam_group1'].beam_mode == 'Horizontal': 
                                beamgrp = 'Beam_group1'
                            elif fileID.groups['Sonar'].groups['Beam_group2'].beam_mode == 'Horizontal':
                                beamgrp = 'Beam_group2'
                                
                            
        
        
                            #Get variables from .nc file
                            variables = tools.GetVariablesFromNC(fileID,beamgrp,False,int(np.unique(WorkIDX[pingIDX])[0]))
                            
                        
                            fileID.close()
                            
                            
                            
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
                
                            sv_ind = np.nan* sv
                            
                            
                            
                            for ipp in range(len(variables.diry)): 
                                WorkPhi[np.where(WorkPhi==variables.diry[ipp])] = int(ipp)
                            for ipp in range(len(np.unique(RangeOut))): 
                                WorkBeam[np.where(WorkBeam==np.unique(RangeOut)[ipp])] = int(ipp)
                            
                                
                            
                            Depth = np.dot(RangeOut[:,np.newaxis],np.sin(np.deg2rad(variables.dirx[:,np.newaxis])).T) +sonar_draft


                            
                            sv_ind[WorkBeam.astype(int),WorkPhi.astype(int)] = sv[WorkBeam.astype(int),WorkPhi.astype(int)]
        
                            #Threshold
                            sv_ind[np.where(sv_ind<Threshold)] = np.nan

                            #Go linear
                            sv_ind = 10**(sv_ind/10)

                            for depth_i in range(NumChannels):
                                a=np.where(np.round(Depth) == depth_i)
                                if len(a[0]) >0: 
                                    SA[depth_i,-1] = SA[depth_i,-1] + np.nansum(sv_ind[a])/len(a[0])


                            SA = np.hstack((SA,np.nan*np.zeros((NumChannels,1))))
                        
                            
    
                            
            elif ein<=0 and to <0: 
                SA2 = SA
                lenSA = SA2.shape[1]

                if lenSA >3: 
                    
                    A = np.nansum(SA2,axis=1)/lenSA*4*np.pi*1852*1852

                    


                    for ikk in range(len(A)): 
                    
                    
                        if A[ikk]>0: 
                            sa_value = ET.SubElement(sa_by_acocat,'sa')
                            sa_value.set('ch',str(ikk+1))
                            sa_value.text = str(A[ikk])
                            
                            

                            tools.indent(root)
                            tree = ET.ElementTree(root)
                            tree.write(directory2Data.dir_result+'\HLUF20.xml', xml_declaration=True, encoding='utf-8', method="xml")
                            
                break
            
                                 
        
            
            
        tools.indent(root)
        tree = ET.ElementTree(root)
        tree.write(directory2Data.dir_result+'\HLUF20.xml', xml_declaration=True, encoding='utf-8', method="xml")
    