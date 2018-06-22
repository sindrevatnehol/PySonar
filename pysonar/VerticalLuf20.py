# -*- coding: utf-8 -*-
"""
Created on Sun May 13 15:32:06 2018

@author: sindrev
"""

 
import time, datetime, os
import urllib.request as url
from xml.etree import ElementTree as ET
import scipy.io as scp
import numpy as np
#import matplotlib.pyplot as plt
from tools import tools
from netCDF4 import Dataset
import random





#def TimeConverter(time0): 
#    tid = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(int(time0)/10000000-int((datetime.datetime(1970,1,1,0,0) - datetime.datetime(1601,1,1,0,0)).total_seconds())))
#    return tid
    
    
    
    
    
    
def GetBottomDepth(maxlat,minlat,maxlon,minlon,delta_lat,delta_lon):
    try: 
        response = url.urlopen('http://coastwatch.pfeg.noaa.gov/erddap/griddap/usgsCeSrtm30v6.csv?topo[(' \
                                +str(maxlat+delta_lat)+'):1:('+str(minlat-delta_lat)+')][('+str(minlon-delta_lon)+'):1:('+str(maxlon+delta_lon)+')]').read().decode("utf-8")
    
        Depth = []
        for i in response.split('\n'): 
            try:  
                Depth = np.hstack((Depth,int(i.split(',')[2])))
            except: 
                dummy=1
    except: 
        Depth = 1000
        print('*    Bad Depth')
    return np.min(Depth)
        
    
      
      
      
      
    
    
    
    
def FastGetComplete(directory2Data,CompleteListOfFiles_vertical): 
    '''fast way to get teh complete list of files'''
    try:    
        Com =  scp.loadmat(directory2Data.dir_work+'/'+'Com.mat')['Com']
    except: 
        Com = []
        for ik in range(len(CompleteListOfFiles_vertical)): 
            tools.printProgressBar(ik,len(CompleteListOfFiles_vertical),prefix = 'Com: ', suffix = 'Completed', length = 50)
            Com = np.hstack((Com,int(CompleteListOfFiles_vertical[ik][0])))
        scp.savemat(directory2Data.dir_work+'/'+'Com.mat',mdict = {'Com':Com})
    Com =  scp.loadmat(directory2Data.dir_work+'/'+'Com.mat')['Com']
    
    return Com
    



    
def MakeVerticalLuf20(CompleteListOfFiles_vertical,directory2Data, nation,cruice_id,platform): 
    
    outloc = '/'
    
    #Some user input data
#    integration_dist = 1
    pel_ch_thickness = 10
    NumChannels = 35
    ChannelSize = 350
    
    
    
    #Something for bookkeeping
    NASC_out = []
    Depth_bottom = 0
    
    
    #Convert stuff to integer, 
    #This is to speed up the program
    Com = FastGetComplete(directory2Data,CompleteListOfFiles_vertical)


    #Something about the channels
    #This should be made more general
    Depth = np.linspace(0,ChannelSize,NumChannels)+ChannelSize/NumChannels/2
    Delta_Z = ChannelSize/NumChannels/2
    NASC = np.nan*np.linspace(0,ChannelSize,NumChannels).T[:,np.newaxis]
    
    
    Work_list = os.listdir(directory2Data.dir_work)
    Work_I = np.arange(len(Work_list))
    random.shuffle(Work_I)
    
    
    #Loop through each work file
    for work_i in Work_I: 
        work = Work_list[work_i]
        NASC_out = []
        print_sa = True
        if work != 'Com.mat':


            
            if not os.path.isfile(directory2Data.dir_result+outloc+'/Vertical/'+work.replace('mat','xml')): 
                
                
    
                
                #Load mask of the filtered data
                try: 
                    workfile =  scp.loadmat(directory2Data.dir_work+'/'+work)
                    
                    
                    #Start information in the text
                    root = ET.Element("echosounder_dataset")
                    ET.SubElement(root,'report_time').text = str(datetime.datetime.fromtimestamp((time.time())).strftime('%Y-%m-%d %H:%M:%S'))
                    ET.SubElement(root,'lsss_version').text = "PNMDformats v 0.1 - vertical"
                    ET.SubElement(root,'nation').text = nation
                    ET.SubElement(root,'platform').text = platform
                    ET.SubElement(root,'cruise').text = cruice_id
                
                    
                    #Make the distance list hirarchy
                    distance_list = ET.SubElement(root,'distance_list')
                
                    
                    
                    try: 
                        pings = np.unique(workfile['ping_mask'])
                        
                    
                    except: 
                        pings = []
                    
    
    
    
                    if not pings == []: 
                        try: 
                            tiime = workfile['start_time'][0]
                            stime = tiime[:4]+'-'+tiime[4:6]+'-'+tiime[6:8]+' '+tiime[9:11]+':'+tiime[11:13]+':'+tiime[13:]
                        except: 
                            stime = ''
                        try: 
                            tiime = workfile['stop_time'][0]
                            etime = tiime[:4]+'-'+tiime[4:6]+'-'+tiime[6:8]+' '+tiime[9:11]+':'+tiime[11:13]+':'+tiime[13:]
                        except: 
                            etime = ''
                        try: 
                            latstart=workfile['lat_start'][0]
                        except: 
                            latstart = ''
                        try: 
                            latstop = workfile['lat_stop'][0]
                        except: 
                            latstop = ''
                        try: 
                            lonstart = workfile['lon_start'][0]
                        except: 
                            lonstart = ''
                        try: 
                            lonstop = workfile['lon_stop'][0]
                        except: 
                            lonstop = ''
    
                        try: 
                            logstart = workfile['log_start'][0]
                        except: 
                            logstart = ''
                        
                            
                            
                        distance = tools.addLogDistanceInfo(distance_list,latstart,lonstart,
                                                      stime,latstop,lonstop,etime,1,10,logstart)
                        
                        
                        
                        
                        
                        #Loop through each unique ping
                        for ping_idx in range(len(pings[1:])): 
                            tools.printProgressBar(ping_idx,len(pings),prefix = 'LUF20 of '+work+':', suffix = 'Completed', length = 50)
                
                
                            #Idx stuff
                            idx = np.where(workfile['ping_mask']==pings[ping_idx+1])
                            idx_file = np.where(Com == (pings[ping_idx+1]))
                            
                            
                            #Get the file name of the idx file
                            filname = directory2Data.dir_rawdata+'/'+CompleteListOfFiles_vertical[idx_file[1][0],1]
                
        
                            #Load file
                            fileID = Dataset(filname,'r',format = 'NETCDF4')
                
        
                            
                            
                            #Get the group with the vertical fan data
                            if fileID.groups['Sonar'].groups['Beam_group1'].beam_mode == 'Vertical': 
                                beamgrp = 'Beam_group1'
                            elif fileID.groups['Sonar'].groups['Beam_group2'].beam_mode == 'Vertical':
                                beamgrp = 'Beam_group2'
                                
                            
                            #Get the complete list of files
                            Comp = CompleteListOfFiles_vertical[idx_file[1][0],2]
        
        
        
                            #Get variables from .nc file
                            variables = tools.GetVariablesFromNC(fileID,beamgrp,False,Comp)
                            
                            
                            #Close file
                            fileID.close()
                            
                            
                            if print_sa == True: 
                                sa_by_acocat = tools.addFrequencyLevelInfo(distance,variables.frequency,0,NumChannels,0)
                                print_sa = False
                                
                                
                                
                            
                            if Depth_bottom == 0: 
                                
                                
                                delta_lat = 600/111111
                                try: 
                                    delta_lon = 600/(111111*np.cos(np.deg2rad(float(workfile['lat_start'][0]))))
                                except: 
                                    delta_lon = 600/111111
                            
    #                            try: 
                                Depth_bottom = GetBottomDepth(float(workfile['lat_start'][0]),float(workfile['lat_stop'][0]),
                                                              float(workfile['lon_start'][0]),float(workfile['lon_stop'][0]),delta_lat,delta_lon)
    #                            except: 
    #                                Depth_bottom = 1000
    #                                print('bad bottom ')
                                    
                                    
                                
                            
                            #Get file name
                            try: 
                                
                                
                                
                                
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
                    
                                
                                
        #                        sv[np.where(RangeOut<=RemoveToCloseValues)] = np.nan
                                   
                                sv[np.where(RangeOut<=30)] = np.nan
                                sv[:,np.where(variables.dirx>=60)] = np.nan
                                sv[np.where(sv<-70)] = np.nan
                        
                                
                                #Go linear
    #                            sv = 10**(sv/10)
                                
                                #Make a temporarly amtrix
                                sv2 = sv#*np.nan
                                
                                sv2 = sv*np.nan
                                sv2[workfile['r_mask'][idx],workfile['b_mask'][idx]]=sv[workfile['r_mask'][idx],workfile['b_mask'][idx]]
                    
                    
    #                            sv2[np.where(sv2<-60)] = np.nan
                                sv2 = 10**(sv2/10)
                                sv2[np.isnan(sv2)] = 0
                                #Remove surface beam
        #                        sv2[:,np.where(variables.dirx < 3)] = np.nan
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
                                    NASC[i,-1]= np.nansum(temp)/temp2
    
    
                                    
                            except IndexError: 
                                print('bad ping')
                    
                            
                        A =  np.nansum(NASC,axis=1)[:,np.newaxis]/len(NASC[0,:])
                        A = A*(1852**2)*pel_ch_thickness*4*np.pi
                            
        #                print(A)
                        if NASC_out == []:
                            NASC_out = A
                        else: 
                            NASC_out = np.hstack((NASC_out, A))
                        
                        
                        for ikk in range(len(A)): 
                            if Depth[ikk]>(abs(Depth_bottom)-150): 
                                A[ikk] = 0
                                NASC_out[ikk,-1] = 0
                        
                        
                            if A[ikk]>0: 
                                sa_value = ET.SubElement(sa_by_acocat,'sa')
                                sa_value.set('ch',str(ikk+1))
                                sa_value.text = str(A[ikk][0])
        
                                
                                
                except: 
                     print('bad file')                   
                tools.indent(root)
                tree = ET.ElementTree(root)
                tree.write(directory2Data.dir_result+outloc+'/Vertical/'+work.replace('mat','xml'), xml_declaration=True, encoding='utf-8', method="xml")
                
