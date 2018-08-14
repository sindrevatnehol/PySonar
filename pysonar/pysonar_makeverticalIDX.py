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
    
    
        
        
        
        
        
        
        
        
        
    
def GetBottomDepth(maxlat,minlat,maxlon,minlon,delta_lat,delta_lon):
    import urllib.request as url
    
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
    


    
    
        
        
        
        
        
        
        
def MakeReport(CompleteListOfFiles,directory2Data, nation,cruice_id,platform): 
    import random
    import datetime, time
    
    outloc = '/'
    
    #Some user input data
#    integration_dist = 1
    pel_ch_thickness = 10
    NumChannels = 35
    PhantomEchoDist = 275
    PhantomEchoWidth = 50
    TH = -70
    ChannelSize = NumChannels*pel_ch_thickness
    
    Depth_bottom = 0
    
    #Something for bookkeeping
    NASC_out = []
    
    #Convert stuff to integer, 
    #This is to speed up the program
    Com = FastGetComplete(directory2Data,CompleteListOfFiles)


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
        if work != 'Com.mat':

            try: 
                if not os.path.isfile(directory2Data.dir_result+outloc+'/Vertical/Report_'+work): 
                    
                
    
                
                    #Load mask of the filtered data
                    workfile =  scp.loadmat(directory2Data.dir_work+'/'+work)
                    
                    
                    
                    try: 
                        pings = np.unique(workfile['ping_mask'])
                    except: 
                        pings = []
                    
        
                        
                        
                    
                    #Loop through each unique ping
                    for ping_idx in range(len(pings[1:])): 
                        
                        tools.printProgressBar(ping_idx,len(pings),prefix = 'LUF20 of '+work+':', suffix = 'Completed', length = 50)
            
            
                        #Idx stuff
                        idx = np.where(workfile['ping_mask']==pings[ping_idx+1])
                        idx_file = np.where(Com == (pings[ping_idx+1]))
                        
                        
                        File = CompleteListOfFiles['FileList'][np.where(CompleteListOfFiles['ping_time']==pings[ping_idx+1])][0]
                        IDX = CompleteListOfFiles['IDX'][np.where(CompleteListOfFiles['ping_time']==pings[ping_idx+1])][0]
        
        
                                
                        #Get the file name of the idx file
                        filname = directory2Data.dir_rawdata+'/'+CompleteListOfFiles['FileList'][idx_file[0]]
                        filname = directory2Data.dir_rawdata+'/'+File
                        
                        #Load file
                        fileID = Dataset(filname,'r',format = 'NETCDF4')
            
                        
                        #Get the group with the vertical fan data
                        if fileID.groups['Sonar'].groups['Beam_group1'].beam_mode == 'Vertical': 
                            beamgrp = 'Beam_group1'
                        elif fileID.groups['Sonar'].groups['Beam_group2'].beam_mode == 'Vertical':
                            beamgrp = 'Beam_group2'
        
                            
        
        
                        #Get variables from .nc file
                        variables = tools.GetVariablesFromNC(fileID,beamgrp,False,IDX)
                        
                        
                        #Close file
                        fileID.close()
                        
                        
        #                        if print_sa == True: 
        #                            sa_by_acocat = tools.addFrequencyLevelInfo(distance,variables.frequency,0,NumChannels,0)
        #                            print_sa = False
                            
                            
                            
                        
                        if Depth_bottom == 0: 
                            
                                
                            delta_lat = 2000/111111
                            try: 
                                delta_lon = 2000/(111111*np.cos(np.deg2rad(float(workfile['lat_start'][0]))))
                            except: 
                                delta_lon = 2000/111111
                        
                            Depth_bottom = GetBottomDepth(float(workfile['lat_start'][0]),float(workfile['lat_stop'][0]),
                                                          float(workfile['lon_start'][0]),float(workfile['lon_stop'][0]),delta_lat,delta_lon)
                            
                        
                        #Get file name
        #                            try: 
                            
                            
                        
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
        #                   
                        sv[np.where(RangeOut<=30)] = np.nan
                        sv[:,np.where(variables.dirx>=60)] = np.nan
                        sv[np.where(sv<TH)] = np.nan
                           
                           
                        #Get the depth of each pixel
                        Y= np.sin(np.deg2rad(variables.dirx))*np.sin(np.deg2rad(variables.diry))
                        Y = (abs(Y[:,np.newaxis]*RangeOut[:,np.newaxis].T)).T
                        sv[np.where(Y>=abs(Depth_bottom+150))] = np.nan
                           
                           
                        X= np.cos(np.deg2rad(variables.dirx))*np.sin(np.deg2rad(variables.diry))
                        X = (X[:,np.newaxis]*RangeOut[:,np.newaxis].T).T
                        sv[np.where(abs(X)>(PhantomEchoDist+PhantomEchoWidth/2))] = np.nan
                        sv[np.where(abs(X)<(PhantomEchoDist-PhantomEchoWidth/2))] = np.nan
                        
                        
                        
                        sv2 = sv*np.nan
                        sv2[workfile['r_mask'][idx],workfile['b_mask'][idx]]=sv[workfile['r_mask'][idx],workfile['b_mask'][idx]]
        #    
        #    
                        sv2 = 10**(sv2/10)
                        sv2[np.isnan(sv2)] = 0
                            
                            
                        
                        #Remove surface beam
        #                        sv2[:,np.where(variables.dirx < 3)] = np.nan
        #                    sv[:,np.where(variables.dirx < 3)] = np.nan
            
        #                    
                        #Select only at a specific distance from vessel
        #                    sv[np.where(abs(X)>300)] = np.nan
        #                    sv[np.where(abs(X)<250)] = np.nan
        
        #
        #    
        #    
        #    
                        #Integrate in depth channels
                        for i in range(len(Depth)): 
                            temp = sv2[np.where(abs(abs(Y)-i*ChannelSize/NumChannels +Delta_Z)<=Delta_Z)]
                            temp2 = np.count_nonzero(~np.isnan(temp))
                            NASC[i,-1]= np.nansum(temp)/temp2
        
        #                    
        ##                        except IndexError: 
        ##                            print('bad ping')
        #    
        #            
                A =  np.nansum(NASC,axis=1)[:,np.newaxis]/len(NASC[0,:])
                A = A*(1852**2)*pel_ch_thickness*4*np.pi
        
                
                
                Liste = {}
                Liste['Report_time'] = str(datetime.datetime.fromtimestamp((time.time())).strftime('%Y-%m-%d %H:%M:%S'))
                
                #General overview
                Liste['Instrument'] = {}
                Liste['Calibration'] = {}
                Liste['DataAcquisition'] = {}
                Liste['DataProcessingMethod'] = {}
                Liste['Cruice'] = {}
        
        
        
                #Instrument section
                
                #Må gå inn i xml fil
                Liste['Instrument']['Frequency'] = variables.frequency
                Liste['Instrument']['TransducerLocation'] = 'AA'
                Liste['Instrument']['TransducerManufacturer'] = 'Simrad'
                Liste['Instrument']['TransducerModel'] = 'SU90'
                Liste['Instrument']['TransducerSerial'] = 'Unknown'
                Liste['Instrument']['TransducerBeamType'] = 'M2'
                Liste['Instrument']['TransducerDepth'] = ''
                Liste['Instrument']['TransducerOrientation'] = ''
                Liste['Instrument']['TransducerPSI'] = ''
                Liste['Instrument']['TransducerBeamAngleMajor'] = ''
                Liste['Instrument']['TransducerBeamAngleMinor'] = ''
                Liste['Instrument']['TransceiverManufacturer'] = ''
                Liste['Instrument']['TransceiverModel'] = ''
                Liste['Instrument']['TransceiverSerial'] = ''
                Liste['Instrument']['TransducerOrientation'] = ''
        
            
                #Calibration Section
                Liste['Calibration']['Date'] = ''
                Liste['Calibration']['AcquisitionMethod'] = ''
                Liste['Calibration']['ProcessingMethod'] = ''
                Liste['Calibration']['AccuracyEstimate'] = ''
                
        
                #Calibration Section
                Liste['DataAcquisition']['SoftwareName'] = ''
                Liste['DataAcquisition']['SoftwareVersion'] = ''
                Liste['DataAcquisition']['StoredDataFormat'] = ''
                Liste['DataAcquisition']['StoredDataFormatVersion'] = ''
                Liste['DataAcquisition']['ConvertedDataFormat'] = ''
                Liste['DataAcquisition']['ConvertedDataFormatVersion'] = ''
                Liste['DataAcquisition']['PingDutyCycle'] = ''
        
        
                #Data Processing
                Liste['DataProcessingMethod']['SoftwareName'] = 'pysonar'
                Liste['DataProcessingMethod']['SoftwareVersion'] = '0.1'
                Liste['DataProcessingMethod']['TriwaveCorrection '] = 'NA'
                Liste['DataProcessingMethod']['ChannelID'] = ''
                Liste['DataProcessingMethod']['Bandwidth'] = ''
                Liste['DataProcessingMethod']['Frequency'] = variables.frequency
                Liste['DataProcessingMethod']['TransceiverPower'] = variables.transmitpower
                Liste['DataProcessingMethod']['TransmitPulseLength'] = variables.pulslength
                Liste['DataProcessingMethod']['OnAxisGain'] = variables.gaintx
                Liste['DataProcessingMethod']['OnAxisGainUnit'] = 'dB'
                Liste['DataProcessingMethod']['SaCorrection'] = variables.sacorrection
                Liste['DataProcessingMethod']['Absorption'] = variables.absorptioncoefficient[0]
                Liste['DataProcessingMethod']['AbsorptionDescription'] = 'Nominal'
                Liste['DataProcessingMethod']['SoundSpeed'] = variables.soundvelocity[0]
                Liste['DataProcessingMethod']['SoundSpeedDescription'] = 'Nominal'
                Liste['DataProcessingMethod']['TransducerPSI'] = ''
        
                
        
                #Cruice
                Liste['Cruice']['Survey'] = ''
                Liste['Cruice']['Country'] = 'NO'
                Liste['Cruice']['Nation'] = nation
                Liste['Cruice']['Platform'] = platform
                Liste['Cruice']['StartDate'] = ''
                Liste['Cruice']['StopDate'] = ''
                Liste['Cruice']['Organisation'] = ''
                Liste['Cruice']['LocalID'] = cruice_id
        
                print(Depth_bottom)
                print(directory2Data.dir_result+outloc+'/Vertical/Report_'+work)
        
        
        #        try: 
                tiime = workfile['start_time'][0]
                stime = tiime[:4]+'-'+tiime[4:6]+'-'+tiime[6:8]+' '+tiime[9:11]+':'+tiime[11:13]+':'+tiime[13:]
        #        except: 
        #            stime = ''
        #        try: 
                tiime = workfile['stop_time'][0]
                etime = tiime[:4]+'-'+tiime[4:6]+'-'+tiime[6:8]+' '+tiime[9:11]+':'+tiime[11:13]+':'+tiime[13:]
        #        except: 
        #            etime = ''
        
        
                Liste['Cruice']['Log'] = {}
    
    
                Liste['Cruice']['Log']['Distance'] = workfile['log_start'][0]
                Liste['Cruice']['Log']['TimeStart'] = stime
                Liste['Cruice']['Log']['TimeStop'] = etime
                Liste['Cruice']['Log']['Latitude_start'] = workfile['lat_start'][0]
                Liste['Cruice']['Log']['Longitude_start'] = workfile['lon_start'][0]
                Liste['Cruice']['Log']['Latitude_stop'] = workfile['lat_stop'][0]
                Liste['Cruice']['Log']['Longitude_stop'] = workfile['lon_stop'][0]
                Liste['Cruice']['Log']['Bottom_depth'] = Depth_bottom
                Liste['Cruice']['Log']['Origin'] = ''
                Liste['Cruice']['Log']['Validity'] = ''
        
        
                Liste['Cruice']['Log']['Sample'] = {}
                Liste['Cruice']['Log']['Sample']['ChannelThickness'] = pel_ch_thickness
                Liste['Cruice']['Log']['Sample']['Acocat'] = 'Unknown'
                Liste['Cruice']['Log']['Sample']['SvThreshold'] = str(TH)
                Liste['Cruice']['Log']['Sample']['PingAxisInterval'] = 1  #samme som integration thickness
                Liste['Cruice']['Log']['Sample']['PingAxisIntervalType'] = 'distance'
                Liste['Cruice']['Log']['Sample']['PingAxisIntervalUnit'] = 'nmi'
                Liste['Cruice']['Log']['Sample']['DataValue'] = A
                Liste['Cruice']['Log']['Sample']['DataType'] = 'C'
                Liste['Cruice']['Log']['Sample']['DataUnit'] = 'm2nm-2'
                Liste['Cruice']['Log']['Sample']['PhantomEchoDistance'] = PhantomEchoDist
                Liste['Cruice']['Log']['Sample']['PhantomEchoWidth'] = PhantomEchoWidth
                
        
                import scipy.io as sc
                sc.savemat(directory2Data.dir_result+outloc+'/Vertical/Report_'+work,mdict=Liste)
                
    
                Depth_bottom = 0
            except: 
                
                from SendMail import send_email
                send_email('bad stuff for Report_'+work)
        #Må generaliseres










#            
##                print(A)
#        if NASC_out == []:
#            NASC_out = A
#        else: 
#            NASC_out = np.hstack((NASC_out, A))
#        
#        
#        for ikk in range(len(A)): 
#            if Depth[ikk]>(abs(Depth_bottom)-150): 
#                A[ikk] = 0
#                NASC_out[ikk,-1] = 0
#        
#        
#            if A[ikk]>0: 
##                        sa_value = ET.SubElement(sa_by_acocat,'sa')
#                sa_value.set('ch',str(ikk+1))
#                sa_value.text = str(A[ikk][0])
#
#                
#                            
#                except: 
#                     print('bad file')                   
#                tools.indent(root)
#                tree = ET.ElementTree(root)
#                tree.write(directory2Data.dir_result+outloc+'/Vertical/'+work.replace('mat','xml'), xml_declaration=True, encoding='utf-8', method="xml")