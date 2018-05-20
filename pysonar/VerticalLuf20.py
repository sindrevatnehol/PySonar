# -*- coding: utf-8 -*-
"""
Created on Sun May 13 15:32:06 2018

@author: sindrev
"""

 
import time, datetime, os, urllib
from xml.etree import ElementTree as ET
import scipy.io as scp
import numpy as np
import matplotlib.pyplot as plt
from tools import tools
from netCDF4 import Dataset


def TimeConverter(time0): 
    tid = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(int(time0)/10000000-int((datetime.datetime(1970,1,1,0,0) - datetime.datetime(1601,1,1,0,0)).total_seconds())))
    return tid
    
    
def GetBottomDepth(maxlat,minlat,maxlon,minlon):
    
     
    response = urllib.request.urlopen('http://coastwatch.pfeg.noaa.gov/erddap/griddap/usgsCeSrtm30v6.csv?topo[(' \
                                +str(maxlat)+'):1:('+str(minlat)+')][('+str(minlon)+'):1:('+str(maxlon)+')]').read().decode("utf-8")

    Depth = []
    for i in response.split('\n'): 
        try:  
            Depth = np.hstack((Depth,int(i.split(',')[2])))
        except: 
            k=1
    return np.min(Depth)

        
        
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
    
    
    
    
def MakeVerticalLuf20(CompleteListOfFiles_vertical,directory2Data): 
    
    #Some user input data
    integration_dist = 1
    pel_ch_thickness = 10
#    include_estimate = 1
    NumChannels = 35
    ChannelSize = 350
    nation = '578'
    platform = '4174'
    cruice_id = '2018105'
    
    
    
    #Convert stuff to integer
    Com = []
    for ik in range(len(CompleteListOfFiles_vertical[:,0])): 
        Com = np.hstack((Com,int(CompleteListOfFiles_vertical[ik,0])))

        
    #Something about the channels
    #This should be made more general
    Depth = np.linspace(0,ChannelSize,NumChannels)+ChannelSize/NumChannels/2
    Delta_Z = ChannelSize/NumChannels/2
    NASC = np.nan*np.linspace(0,ChannelSize,NumChannels).T[:,np.newaxis]
    
    #Something for bookkeeping
    start_Lat = []
    start_Lon = []
    LogDist = []
    LogDist_i = 0
    StartTime = []
    StopTime = []
    StartLat = []
    EndLat = []
    StartLon = []
    EndLon = []
    integration_idx = 0
    integration_idx_vec = 0

    NASC_out = []
    
    
    import random
    
    Work_list = os.listdir(directory2Data.dir_work)
    Work_I = np.arange(len(Work_list))
    random.shuffle(Work_I)
    #Loop through each work file
    for work_i in Work_I: 
        work = Work_list[work_i]
        NASC_out = []
        print(work)
        if not os.path.isfile(directory2Data.dir_result+"\2"+work.replace('mat','xml')): 
            
            
            #Start information in the text
            root = ET.Element("echosounder_dataset")
            ET.SubElement(root,'report_time').text = str(datetime.datetime.fromtimestamp((time.time())).strftime('%Y-%m-%d %H:%M:%S'))
            ET.SubElement(root,'lsss_version').text = "PNMDformats v 0.1 - vertical"
            ET.SubElement(root,'nation').text = nation
            ET.SubElement(root,'platform').text = platform
            ET.SubElement(root,'cruise').text = cruice_id
        
            
            #Make the distance list hirarchy
            distance_list = ET.SubElement(root,'distance_list')

    
            #Load mask of the filtered data
            workfile =  scp.loadmat(directory2Data.dir_work+'/'+work)
            pings = np.unique(workfile['ping_mask'])
            
            
            #Loop through each unique ping
            for ping_idx in range(len(pings[1:])): 
                tools.printProgressBar(ping_idx,len(pings),prefix = 'MakingLufOfTransect: ', suffix = 'Completed', length = 50)
    
    
                #Idx stuff
                idx = np.where(workfile['ping_mask']==pings[ping_idx+1])
                idx_file = np.where(Com == (pings[ping_idx+1]))
                
                #Get file name
                try: 
                    filname = directory2Data.dir_rawdata+'/'+CompleteListOfFiles_vertical[idx_file,1][0,0]
        
                    #Load file
                    fileID = Dataset(filname,'r',format = 'NETCDF4')
        
        
                    #Get the group with the vertical fan data
                    if fileID.groups['Sonar'].groups['Beam_group1'].beam_mode == 'Vertical': 
                        beamgrp = 'Beam_group1'
                    elif fileID.groups['Sonar'].groups['Beam_group2'].beam_mode == 'Vertical':
                        beamgrp = 'Beam_group2'
                        
        
                    #Get variables from .nc file
                    variables = tools.GetVariablesFromNC(fileID,beamgrp,False,CompleteListOfFiles_vertical[idx_file,2][0,0])
                    fileID.close()
        
                    
                    #Get distance information
                    if start_Lat == []: 
                        start_Lat = variables.Latitude
                        start_Lon = variables.Longitude
                        StartLat = variables.Latitude
                        StartLon = variables.Longitude
                        LogDist = 0
                        StartTime = TimeConverter(variables.time)
                        integration_idx = integration_idx+1
                        
                        
                        
                    else: 
                        stop_lat = variables.Latitude
                        stop_lon = variables.Longitude
        
                        int_dist = tools.haversine(start_Lon, start_Lat, stop_lon, stop_lat)/1852
                        integration_idx = integration_idx+1
        
                        if int_dist>integration_dist: 
                            start_Lat = variables.Latitude
                            start_Lon = variables.Longitude
                            
                            EndLon = np.hstack((EndLon,variables.Longitude))
                            EndLat = np.hstack((EndLat,variables.Latitude))
                            
                            EndTime = TimeConverter(variables.time)
                            
                            try: 
                                distance = addLogDistanceInfo(distance_list,StartLat[-1],StartLon[-1],
                                                              StartTime,EndLat[-1],EndLon[-1],EndTime,
                                       int_dist,pel_ch_thickness,LogDist[-1])
                            except TypeError: 
                                distance = addLogDistanceInfo(distance_list,StartLat[0],StartLon[0],
                                                              StartTime,EndLat[0],EndLon[0],EndTime,
                                       int_dist,pel_ch_thickness,LogDist)
                                
                            
                            sa_by_acocat = addFrequencyLevelInfo(distance,variables.frequency,0,NumChannels,0)
                                
                            
                            StartLat = np.hstack((StartLat,variables.Latitude))
                            StartLon = np.hstack((StartLon, variables.Longitude))
                            LogDist = np.hstack((LogDist,int_dist+LogDist_i))
                            LogDist_i = LogDist_i+int_dist
                            integration_idx_vec = np.hstack((integration_idx_vec,integration_idx))
                            
                        
                            #Get all depth within the sonar detection volume
                            delta_lat = 600/111111
                            delta_lon = 600/(111111*np.cos(np.deg2rad(StartLat[-1])))
                            
                            if StartLat[-1]>EndLat[-1]: 
                                delta_lat = -delta_lat
                            if StartLon[-1]>StartLon[-1]:
                                delta_lon = -delta_lon
                            
                            Depth_bottom = GetBottomDepth(StartLat[-1]-delta_lat,EndLat[-1]+delta_lat,
                                                          StartLon[-1]-delta_lon,EndLon[-1]+delta_lon)
#                            try: 
#                                Depth_bottom = GetBottomDepth(StartLat[-1],EndLat[-1],StartLon[-1],EndLon[-1])
#                            except: 
#                                Depth_bottom = GetBottomDepth(StartLat[-1],EndLat[-1]+0.0001,StartLon[-1],EndLon[-1]+0.0001)
#                                print('bad gps distance')
#                            Depth_bottom = -200
                            
                            A =  np.nansum(NASC,axis=1)[:,np.newaxis]/len(NASC[0,:])
                            A = A*(1852**2)*pel_ch_thickness
                            
                            if NASC_out == []:
                                NASC_out = A
                            else: 
                                NASC_out = np.hstack((NASC_out, A))
                                
                                
                            for ikk in range(len(A)): 
                                if Depth[ikk]>(abs(Depth_bottom)-100): 
                                    A[ikk] = 0
                                    NASC_out[ikk,-1] = 0
    
                                
                                if A[ikk]>0: 
                                    sa_value = ET.SubElement(sa_by_acocat,'sa')
                                    sa_value.set('ch',str(ikk+1))
                                    sa_value.text = str(A[ikk][0])
    
                                    
                            NASC = np.nan*np.linspace(0,ChannelSize,NumChannels).T[:,np.newaxis]
                  
            
                    
                    
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
        
                    #Go linear
                    sv = 10**(sv/10)
                    
                    #Make a temporarly amtrix
                    sv2 = sv#*np.nan
#                    sv2[workfile['r_mask'][idx],workfile['b_mask'][idx]]=sv[workfile['r_mask'][idx],workfile['b_mask'][idx]]
        
        
                    #Remove surface beam
                    sv2[:,np.where(variables.dirx < 3)] = np.nan
#                    sv[:,np.where(variables.dirx < 3)] = np.nan
        
                        
                    #Select only at a specific distance from vessel
                    X= np.cos(np.deg2rad(variables.dirx))*np.sin(np.deg2rad(variables.diry))
                    X = (X[:,np.newaxis]*RangeOut[:,np.newaxis].T).T
                    sv2[np.where(abs(X)>300)] = np.nan
                    sv2[np.where(abs(X)<250)] = np.nan
#                    sv[np.where(abs(X)>300)] = np.nan
#                    sv[np.where(abs(X)<250)] = np.nan


                    
        
                    #Get the depth of each pixel
                    Y= np.sin(np.deg2rad(variables.dirx))*np.sin(np.deg2rad(variables.diry))
                    Y = (abs(Y[:,np.newaxis]*RangeOut[:,np.newaxis].T)).T
        
        
                    #Integrate in depth channels
                    for i in range(len(Depth)): 
                        temp = sv2[np.where(abs(abs(Y)-i*ChannelSize/NumChannels +Delta_Z)<=Delta_Z)]
                        temp2 = np.count_nonzero(~np.isnan(temp))
#                        print(temp2,len(np.where(abs(abs(Y)-i*ChannelSize/NumChannels +Delta_Z)<=Delta_Z)[0]))
                        NASC[i,-1]= np.nansum(temp)/temp2
#                        NASC[i,-1]= np.nansum(sv2[np.where(abs(abs(Y)-i*ChannelSize/NumChannels +Delta_Z)<=Delta_Z)])/len(np.where(abs(abs(Y)-i*ChannelSize/NumChannels +Delta_Z)<=Delta_Z)[0])
#                        NASC[i,-1]= np.nansum(sv[np.where(abs(abs(Y)-i*ChannelSize/NumChannels +Delta_Z)<=Delta_Z)])/len(np.where(abs(abs(Y)-i*ChannelSize/NumChannels +Delta_Z)<=Delta_Z)[0])
                        
                        
                    #Stack NASC for each ping
                    NASC = np.hstack((NASC,np.nan*np.linspace(0,ChannelSize,NumChannels).T[:,np.newaxis]))
                        
                except IndexError: 
                    print('bad ping')
          
                    
            try: 
                plt.figure(2)
                plt.clf()
                plt.imshow(10*np.log10(NASC_out),aspect='auto',interpolation='none')
                plt.colorbar()
                plt.savefig(directory2Data.dir_result+'2/'+work.replace('mat','jpg'))
                plt.pause(0.1)
            except: 
                k=1
                                
            indent(root)
            tree = ET.ElementTree(root)
            tree.write(directory2Data.dir_result+"2/"+work.replace('mat','xml'), xml_declaration=True, encoding='utf-8', method="xml")
            
#            f = open(directory2Data.dir_src+'/Luf20prog/'+work.replace('mat','txt'),'w')
#            f.close()
            
    
#    indent(root)
# 
#    #write the file to xml format
#    tree = ET.ElementTree(root)
#    tree.write(directory2Data.dir_result+"/ListUserFile20.xml", xml_declaration=True, encoding='utf-8', method="xml")
    