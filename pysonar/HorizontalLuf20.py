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

    
    
    
    
    
    

    
def addFrequencyLevelInfo(distance,freq,tran,NumCH,TH): 
    '''Add information in the frequency level'''
    
    
    #Make new frequency with attributes
    frequency = ET.SubElement(distance,'frequency')
    frequency.set('freq',str(int(freq))  )
    frequency.set('transceiver',str(int(tran)))
    ET.SubElement(frequency,'quality').text = str(2)  
    ET.SubElement(frequency,'bubble_corr').text = str(0)  
    ET.SubElement(frequency,'threshold').text = str(TH)
    
    ET.SubElement(frequency,'num_pel_ch').text = str(NumCH)  
    
    ET.SubElement(frequency,'upper_interpret_depth').text = str(0)
    ET.SubElement(frequency,'lower_interpret_depth').text = str(0)
    ET.SubElement(frequency,'upper_interpret_depth').text = str(0)
    ET.SubElement(frequency,'lower_integrator_depth').text = str(0)

    

    ch_type = ET.SubElement(frequency,'ch_type')
    ch_type.set('type','P')
    sa_by_acocat = ET.SubElement(ch_type,'sa_by_acocat')
    sa_by_acocat.set('acocat',str(12))
    
    return sa_by_acocat
    
    
    
    
    
    
    
        
def indent(elem, level=0):
  '''
  Description: 
       Make the xml file more readable
  '''
  i = "\n" + level*"  "
  if len(elem):
    if not elem.text or not elem.text.strip():
      elem.text = i + "  "
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
    for elem in elem:
      indent(elem, level+1)
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
  else:
    if level and (not elem.tail or not elem.tail.strip()):
      elem.tail = i
      
      
      
      
    
    
      
      
def addLogDistanceInfo(distance_list,LatStart,LonStart,TimeStart,LatStop,LonStop,TimeStop,
                       integrator_dist,pel_ch_thickness,LogStart): 
    '''Add information in the log distance level'''
    
    
    #Get the unix time and convert it to time string
    correct_starttime =TimeStart# time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(TimeStart))
    correct_stoptime = TimeStop #time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(TimeStop))
    
    
    
    #Write new log distance with its attributes
    distance = ET.SubElement(distance_list,'distance')
    distance.set('log_start',str(LogStart))
    distance.set('start_time',str(correct_starttime))

    
    
    #Add variables to the log distance
    ET.SubElement(distance,'integrator_dist').text = str(integrator_dist)
    ET.SubElement(distance,'pel_ch_thickness').text = str(pel_ch_thickness)
    ET.SubElement(distance,'include_estimate').text = '1'
    ET.SubElement(distance,'lat_start').text = str(LatStart)
    ET.SubElement(distance,'lon_start').text = str(LonStart)
    ET.SubElement(distance,'lat_stop').text = str(LatStop)
    ET.SubElement(distance,'lon_stop').text = str(LonStop)
    ET.SubElement(distance,'stop_time').text = str(correct_stoptime)
    return distance
    
    
    
    
    
    
    
    
    

def MakeHorizontalLuf20(CompleteListOfFiles,directory2Data,start_time,log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop, TimeIDX ): 
    
    
    sonar_draft = 8
    
    
    
    #Fiksed stuff
    integration_dist = 1
    pel_ch_thickness = 1
    NumChannels = 300
    ChannelSize = 350
    nation = '578'
    platform = '4174'
    cruice_id = '2018105'
    
    
    
                
    #Start information in the text
    root = ET.Element("echosounder_dataset")
    ET.SubElement(root,'report_time').text = str(datetime.datetime.fromtimestamp((time.time())).strftime('%Y-%m-%d %H:%M:%S'))
    ET.SubElement(root,'lsss_version').text = "PNMDformats v 0.1 - vertical"
    ET.SubElement(root,'nation').text = nation
    ET.SubElement(root,'platform').text = platform
    ET.SubElement(root,'cruise').text = cruice_id

    
    #Make the distance list hirarchy
    distance_list = ET.SubElement(root,'distance_list')
                
                
    
    
    
    Listoffiles = os.listdir(directory2Data.dir_work)
            
            
            
                
    for i in range(len(start_time)): 
        print(start_time[i])
        
        
        SA = np.zeros((NumChannels,1))*np.nan
        
        
        stime = stop_time[i][:4]+'-'+start_time[i][4:6]+'-'+start_time[i][6:8]+' '+start_time[i][9:11]+':'+start_time[i][11:13]+':'+start_time[i][13:]
        etime = stop_time[i][:4]+'-'+stop_time[i][4:6]+'-'+stop_time[i][6:8]+' '+stop_time[i][9:11]+':'+stop_time[i][11:13]+':'+stop_time[i][13:]
        





        distance = addLogDistanceInfo(distance_list,lat_start[i],lon_start[i],
                                      stime,lat_stop[i],lon_stop[i],etime,1,10,log_start[i])
                    
                    
                    
        
        
        sa_by_acocat = addFrequencyLevelInfo(distance,26000,0,NumChannels,0)
                            
       
                    
        

        for iii in range(CompleteListOfFiles.shape[0]): 
            
            
            ein = SecondsBetweenTransect(start_time[i].replace('T',''),int(CompleteListOfFiles[iii][0]))
            to = SecondsBetweenTransect(stop_time[i].replace('T',''),int(CompleteListOfFiles[iii][0]))
            
            
            
            
            
            
            if ein <=0 and to>=0: 
                
                
                
            
                for file in Listoffiles: 
                    if 'Horizontal' in file: 
                        
                        Mat = scp.loadmat(directory2Data.dir_work+'\\'+file)
        
                
                        pingIDX = np.where(Mat['WorkFile'][:,0]==(CompleteListOfFiles[iii][0]))[0]



                        if len(pingIDX) >0: 
                            
                            
                        
                            #Get the file name of the idx file
                            filname = directory2Data.dir_rawdata+'/'+Mat['WorkFile'][int(Mat['WorkFile'][pingIDX,2]),1]
                
        
                            #Load file
                            fileID = Dataset(filname,'r',format = 'NETCDF4')
                
        
                            
                            #Get the group with the vertical fan data
                            if fileID.groups['Sonar'].groups['Beam_group1'].beam_mode == 'Horizontal': 
                                beamgrp = 'Beam_group1'
                            elif fileID.groups['Sonar'].groups['Beam_group2'].beam_mode == 'Horizontal':
                                beamgrp = 'Beam_group2'
                                
                            
        
        
                            #Get variables from .nc file
                            variables = tools.GetVariablesFromNC(fileID,beamgrp,False,int(Mat['WorkFile'][pingIDX,2]))
                            
                        
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
                
                            
                            
                            
                            BeamIDX = np.where(variables.diry == float(Mat['WorkFile'][pingIDX,3]))[0]
                            RIDX = np.where(RangeOut == float(Mat['WorkFile'][pingIDX,4]))[0]
                            Depth = RangeOut[RIDX]*np.sin(np.deg2rad(variables.dirx[BeamIDX])) +sonar_draft



                            sa = sv[RIDX,BeamIDX]
        


                            if sa <-65: 
                                sa = np.nan



                            sa = 10**(sa/10)



                            SA[int(np.round(Depth)-1),-1] = sa


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
                            
                            

                            indent(root)
                            tree = ET.ElementTree(root)
                            tree.write(directory2Data.dir_result+'HLUF20.xml', xml_declaration=True, encoding='utf-8', method="xml")
                            
                break
            
                                 
        
            
            
        indent(root)
        tree = ET.ElementTree(root)
        tree.write(directory2Data.dir_result+'HLUF20.xml', xml_declaration=True, encoding='utf-8', method="xml")
    