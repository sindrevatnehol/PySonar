# -*- coding: utf-8 -*-
"""
Created on Mon May  7 10:39:10 2018

@author: sindrev
"""

import os, scipy
import numpy as np
from shutil import copyfile
import Raw2NetcdfConverter




def printProgressBar(iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = ' '):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    
        

        
        
        
        
        
        
        
def WelcomeScreen(projectName): 
    
    print(projectName)
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('    >=>      /\      >=>          >=>           ')
    print('  >=>       /  \           >=>         >=>      ')
    print('           / >=>\                               ')
    print('   >=>    /      \                 >=>          ')
    print('         /   >=>  \      >=>                >=> ')
    print('        / >=>      \              >=>    >=>    ')
    print('_______/________>=>_\___________________________')
    print('                              Sindre N. Vatnehol')
    print('                    Institute of Marine Research')
    print('                                          Norway')
    print('                     Mail: sindre.vatnehol@hi.no')

    
       
    
    
    
    

def UnpackBeam(BeamAmplitudeData):
    '''Unpack the beam data from the nc structure'''
    
    BeamAmplitude = np.zeros((len(BeamAmplitudeData[0]),64))
    
    for i in range(len(BeamAmplitude[0,:])):
        BeamAmplitude[:,i] = BeamAmplitudeData[i]

    return BeamAmplitude

    
    
            
            
    
class FolderStructure(object):
    '''
    Defines the folder structure used in the project.
    the PySonar... wil probably be changed
    '''
    def __init__(self, dir_cruice,equipment):
        self.dir_cruice = dir_cruice+'/'
        self.dir_acoustic_data = dir_cruice+'/ACOUSTIC_DATA/'
        self.dir_su90 = dir_cruice+'/ACOUSTIC_DATA/'+equipment+'/'
        self.dir_PySonar = dir_cruice+'/ACOUSTIC_DATA/'+equipment+'/PySonar'
        self.dir_nc = dir_cruice+'/ACOUSTIC_DATA/'+equipment+'/PySonar/netcdf'
        self.dir_work = dir_cruice+'/ACOUSTIC_DATA/'+equipment+'/PySonar/WorkFiles'
        self.dir_search = dir_cruice+'/ACOUSTIC_DATA/'+equipment+'/PySonar/Search'
        self.dir_result = dir_cruice+'/ACOUSTIC_DATA/'+equipment+'/PySonar/Result'
        self.dir_rawdata = dir_cruice+'/ACOUSTIC_DATA/'+equipment+'/RAWDATA'
        self.dir_originalrawdata = dir_cruice+'/ACOUSTIC_DATA/'+equipment+'/ORIGINALRAWDATA'
        self.dir_src = dir_cruice+'/ACOUSTIC_DATA/'+equipment+'/PySonar/src'
        self.dir_NCconvertProgress = dir_cruice+'/ACOUSTIC_DATA/'+equipment+'/PySonar/src/NCconvertProgress'
            
        
        
        
    
        
def MakeNewFolders(directory2Data):
    '''Process to make new folder if they do not exist'''
    if not os.path.exists(directory2Data.dir_cruice):
        os.makedirs(directory2Data.dir_cruice)
    if not os.path.exists(directory2Data.dir_acoustic_data):
        os.makedirs(directory2Data.dir_acoustic_data)
    if not os.path.exists(directory2Data.dir_su90):
        os.makedirs(directory2Data.dir_su90)
    if not os.path.exists(directory2Data.dir_work):
        os.makedirs(directory2Data.dir_work)
    if not os.path.exists(directory2Data.dir_PySonar):
        os.makedirs(directory2Data.dir_PySonar)
    if not os.path.exists(directory2Data.dir_search):
        os.makedirs(directory2Data.dir_search)
    if not os.path.exists(directory2Data.dir_result):
        os.makedirs(directory2Data.dir_result)
    if not os.path.exists(directory2Data.dir_nc):
        os.makedirs(directory2Data.dir_nc)
    if not os.path.exists(directory2Data.dir_originalrawdata):
        os.makedirs(directory2Data.dir_originalrawdata)
    if not os.path.exists(directory2Data.dir_rawdata):
        os.makedirs(directory2Data.dir_rawdata)
    if not os.path.exists(directory2Data.dir_src):
        os.makedirs(directory2Data.dir_src)
    if not os.path.exists(directory2Data.dir_work):
        os.makedirs(directory2Data.dir_work)
    if not os.path.exists(directory2Data.dir_NCconvertProgress):
        os.makedirs(directory2Data.dir_NCconvertProgress)
        
               
        
        
        


def DataConverter(CruiceIndex,WorkDirectory,current_dir,maxPingInFile,
                  MaxNumberOfFilesInNC,directory2Data, reconvert) :
    ''' 
    Protocoll to convert the data
    '''

    
    #Vessel and platform information
    vessel_name = CruiceIndex.getAttribute('vessel')
    platform_type = CruiceIndex.getAttribute('platform_type')
    platform_code = CruiceIndex.getAttribute('platform_code')
    
    
    
    
    #Loop through each fishery sonar type
    #Denne vil endres til å bli mer generell
    for i in ['SU90','SX90','SH90']:         
            
        
        os.chdir(current_dir)
        Raw2NetcdfConverter.convert.Raw2NetcdfConverter(directory2Data.dir_originalrawdata,vessel_name,platform_type,
                    platform_code,maxPingInFile,MaxNumberOfFilesInNC,
                    directory2Data.dir_rawdata,reconvert)
            
            
        #Convert the .raw data to nc. 
#        if not os.path.exists(directory2Data.dir_NCconvertProgress + '/finishedNC.txt'):
#            os.chdir(current_dir)
#            Raw2NetcdfConverter.convert.Raw2NetcdfConverter(directory2Data.dir_originalrawdata,vessel_name,platform_type,
#                        platform_code,maxPingInFile,MaxNumberOfFilesInNC,
#                        directory2Data.dir_rawdata)
#            
#            
#            #When convertion is finnished, make a txt file to indicate this
#            f = open(directory2Data.dir_NCconvertProgress + '/finishedNC.txt','w')
#            f.write('0')
#            f.close()
        
       
            
            
        
def TransectTimeIDX(CruiceIndex): 
    TimeIDX =[]
    TransCode = CruiceIndex.getElementsByTagName("transect")
    for t in TransCode: 
        if TimeIDX == []: 
            TimeIDX = np.array([CruiceIndex.getAttribute("code")+'_'+t.getAttribute('num'),
                                t.getAttribute('startTime'),t.getAttribute('endTime')
                                ])[:,np.newaxis].T
        else: 
            TimeIDX  = np.vstack((TimeIDX,np.array([CruiceIndex.getAttribute("code")+'_'+t.getAttribute('num'),
                                t.getAttribute('startTime'),t.getAttribute('endTime')
                                ])[:,np.newaxis].T))      
    return TimeIDX
    
    
    
    
    
        
        
        
def OrginizeData(CruiceIndex,WorkDirectory,OS): 
    
    
    #Get the vessel, cruice and equipment information
    vesselCode = CruiceIndex.getAttribute("vesselCode")
    cruiceCode = CruiceIndex.getAttribute("code")
    equipment =  CruiceIndex.getAttribute("Equipment")
    
    
    #Loop through each fishery sonar type
    for i in ['SU90','SX90','SH90']: 
        
        
        if i == equipment: 
            
            #Make a correct sonar folder structure
            directory2Data =FolderStructure(WorkDirectory+'/'+cruiceCode[:4]+'/S'+cruiceCode+'_P'+vesselCode,i)       
            
            MakeNewFolders(directory2Data)
        
            
#            #Get list of files that has not been copied to the correct structure
#            print('    -Get files from server  ',end='\r')
#            ListFromServer = os.listdir(OS+CruiceIndex.getAttribute('CruicePath'))
#            ListOfFilesNotcopied = list(set(ListFromServer)-set(os.listdir(directory2Data.dir_originalrawdata)))
#            
#            
#            #Go through each file
#            for i in np.arange(len(ListOfFilesNotcopied)): 
#                
#                #Print a progressbar
#                printProgressBar(i, len(ListOfFilesNotcopied), prefix = 'Copy files', suffix = 'Files left: '+
#                                 str(len(ListOfFilesNotcopied)-i) , decimals = 1, length = 50)
#                
#                #Copy the files
#                try:  
#                    try:
#                        copyfile(OS+CruiceIndex.getAttribute('CruicePath')+'/'+ListOfFilesNotcopied[i], 
#                             directory2Data.dir_originalrawdata+'/'+ListOfFilesNotcopied[i])
#                    except PermissionError: 
#                        print('', end = '\r')
#                except IsADirectoryError:  
#                    print('', end = '\r')
#                    
    
    return(directory2Data)
   

    
    

def GainAdjustment(pulslength,frequency,inngain):
    '''Find the correct gain to be used when scrutinising the sonar data. 
    #NB: the pulselength must be in ms and the frequency in kHz
    #
    #The correction for the pulse length is similar as to the one given in (Macauley et al. 2016)
    #The correction for frequency is different due to new avaliable information. 
    '''
    
    a = 6.99915703
    b = -6.24304318
    PulslengthGain = a +b/pulslength
    
    
    
    a = 48.4116
    b = -10.9486
    c = 0.6446
    d = -0.0108
    FrequencyGain = a+b*frequency +c*frequency**2+d*frequency**3
    
    outgain = FrequencyGain+PulslengthGain+inngain
    
    return outgain, FrequencyGain, PulslengthGain
    
    
    
    
    
    
def ApplyTVG(AmplitudeData,soundvelocity,sampleinterval,transmitpower,
             absorptioncoefficient,frequency,pulslength,gain,
             equivalentbeamangle,sacorrection,tiltAngle): 
    
    
    
    #There is a time delay in the sonar. This corresponds to 3 m
    rangeCorr = 3; 
    
    
    effective_pulselength = 0.7*pulslength
    
    #The distance between each sample in meter
    samplespace =soundvelocity*sampleinterval/2
    
    
    #Declare the range vector to be outputed
    Range = np.arange(0,AmplitudeData[:].shape[0])*samplespace-rangeCorr
                      
                      
                      
    #Remove bad values  
    Range[np.where(Range<=0)] = 1E-200
               

    #Make the TVG function (20 log) as a vector
    tvg20 = np.repeat(20*np.log10(Range[:, np.newaxis]) + 2*absorptioncoefficient*Range[:, np.newaxis], AmplitudeData[:].shape[1], axis=1)
                         
                      
    tvg20[np.where(tvg20<0)] = 0
    
    #Correction function for the loss off energy when steering the beam
    tiltcorr =40*np.log10(abs(np.cos(np.deg2rad(-tiltAngle))))
    
    #compute the wavelength of the signal
    wavelength = soundvelocity/frequency

    
    
    #compute the correction function used when computing the SV
    svconstSV =10*np.log10(transmitpower*wavelength**2*soundvelocity*effective_pulselength/(32*np.pi**2))+tiltcorr+gain+10*np.log10(equivalentbeamangle)+sacorrection
    
        
    #Compute the volume backscattering coefficient
    sv= AmplitudeData+tvg20-np.repeat(svconstSV[:,np.newaxis],AmplitudeData.shape[0],axis=1).transpose()

    #Return the SV the range and the TS from the data.
    return sv, Range; 


    
    
    
def ConvertToechogram(Wdist,sv_mat):
    '''
    Purpose
    Search through the data in a banana shape to find the sA
    Not likely to be changed unless the methodology changes
    
    
    It requieres large memory
    '''
    
    
    Wdist = Wdist.toarray()
    
    #Sett weight to 1 for all relevant distances
    Wdist[np.where(Wdist>0)] = 1
    

    #find the index of the files without 0, This for speed improvemment
    idxWeight = np.where(Wdist>0)
                
    #Find the sa value in linear domain. 
    sA= np.nansum(sv_mat[idxWeight])
    
    return sA
    
    
    
    
def GetDistanceMatrix(DistanceMatrix,RangeMatrix,BeamDirectionMatrix,svMatrix,theta_tilt,BananaTool):   
    '''
    This function computes the distance matrix. 
    For values outside the tubed school are sett to 0. 
    
    It requieres large memory
    
    '''
    
    
    
    #make the distance matrixes as a dictionary for easyer access
    Wdist_port = dict()
    Wdist_stb= dict()
    
    
    
    
    #convert angles into radians
    theta_tilt = theta_tilt * np.pi/180
    
    
    
    
    #Distance traveled by the vessel, centered at origo
    '''This will be made more general in the future. 
    E.g. from number of sampels in time and range, mean ping time, and vessel speed, make both
    the travelDistance and the x/y/z pixel matrixes
    '''

    travelDistance = DistanceMatrix[0,0,:]-np.nanmin(DistanceMatrix)
    
    travelDistance = np.max(travelDistance)/len(travelDistance)*np.arange(len(travelDistance))
    
    travelDistance=travelDistance-np.max(travelDistance)/2

    #Compute the location of all points in cartesian coordinates
    x_pixel = RangeMatrix*np.sin(theta_tilt)*np.cos(BeamDirectionMatrix*np.pi/180)
    
    x_pixel = x_pixel+DistanceMatrix
    
    x_pixel = x_pixel-np.min(x_pixel)
    
    x_pixel = x_pixel-np.max(x_pixel)/2    
    
    y_pixel = RangeMatrix*np.sin(theta_tilt)*np.sin(BeamDirectionMatrix*np.pi/180)
    
    z_pixel = RangeMatrix*np.cos(theta_tilt)
    
    
    correction = np.repeat(np.repeat((RangeMatrix[:,0,0]*np.sin(2*np.pi/64))[:,np.newaxis],64,axis=1)[:,:,np.newaxis],len(RangeMatrix[0,0,:]),axis=2)
    
    index=0
    #Run through all pings and find the distance of each pixel 
    for r_ind_stbex in range(int(BananaTool[1]),int(BananaTool[2]),int(BananaTool[3])):
        printProgressBar(r_ind_stbex+1,int(BananaTool[2]), prefix = 'Make Distance:', suffix = 'Completed     ', length = 50)
        
        
        #Get the minimum range to sonar trajectory
        r_0 = RangeMatrix[r_ind_stbex,0,0]



        #Model for school trajectory, assuming a stationary school      
        x_traj = np.zeros(travelDistance.shape) #travelDistance          
        y_traj_port = r_0*np.ones(x_traj.shape) *np.sin(theta_tilt)    
        y_traj_stb = -r_0*np.ones(x_traj.shape) *np.sin(theta_tilt)
        z_traj = np.sqrt(x_traj**2+y_traj_port**2)/np.tan(theta_tilt)
        
        
               
        #Empty range matrix for bookkeeping
        R_port = np.zeros(svMatrix.shape)
        R_stb = np.zeros(svMatrix.shape)
        R2_port = np.zeros(svMatrix.shape)
        R2_stb = np.zeros(svMatrix.shape)
        
        
        
        
        #Compute the range of all pixels from the ideal line        
        for i in range(len(x_traj)): 
            R_port[:,:,i] = np.sqrt((x_pixel[:,:,i]-x_traj[i])**2+(y_pixel[:,:,i]
                                -y_traj_port[i])**2+(z_pixel[:,:,i]-z_traj[i])**2)
            R_stb[:,:,i] = np.sqrt((x_pixel[:,:,i]-x_traj[i])**2+(y_pixel[:,:,i]
                                -y_traj_stb[i])**2+(z_pixel[:,:,i]-z_traj[i])**2)
        
            
           
        indeksen = np.repeat(np.repeat(np.reshape(np.amin(np.amin(R_port,axis=1),
                              axis=0),(1,1,-1)),R_port.shape[0],axis=0),64,axis=1)
        indeksen_stb = np.repeat(np.repeat(np.reshape(np.amin(np.amin(R_stb,axis=1),
                              axis=0),(1,1,-1)),R_stb.shape[0],axis=0),64,axis=1)
            
        
        R2_port[R_port<=(indeksen+correction+int(BananaTool[0]))]=1
        R2_stb[R_stb<=(indeksen_stb+correction+int(BananaTool[0]))]=1
        
               
               
        #Save the distance matrixe to memory/alternatively it can be stoored in a seperate file
        Wdist_port[index] =  scipy.sparse.csr_matrix(np.reshape(R2_port,(-1,1)))
        Wdist_stb[index] =  scipy.sparse.csr_matrix(np.reshape(R2_stb,(-1,1)))
        
        index+=1
        

    return Wdist_port,Wdist_stb

    
    
    
    
    
    
    
    
    
    
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
    
    

    
def ComputeDistance(travelDist,lat,lon): 
    '''
    This function computes the distance traveled within each transect. 
    For distance the hamorsine function is used. 
    
    Parts of the code is awkward due to bugs in the files. Will be changed later. 
    '''
    
    #We need some points to evaluate the distanc traveled. 
    if len(lat)==1: 
        travelDist = 0
    else: 
        travelDist = np.hstack((travelDist,haversine(lon[-1], lat[-1], lon[-2], lat[-2])))
    DistanceTraveled = np.cumsum(travelDist)
#        travelDist=np.hstack((travelDist,np.nan))
        
        
        
#    elif len(lat)>=20:
#        #find the delta in longitude and latitude
#        delta_lat = lat[-1]-lat[-2]
#        delta_lon = lon[-1]-lon[-2]
#
#
#
#        #Computing the distance using hammersine fucntion
#        a = (np.sin(delta_lat/2*np.pi/180))**2 + np.cos(lat[-2]*np.pi/180) * np.cos(lat[-1]*np.pi/180) * (np.sin(delta_lon/2*np.pi/180)**2)
#        c = 2 * np.arctan2(np.sqrt(a),np.sqrt(1-a))
#        
#        
#        #Compute the distance to meter and stack the variable
#        travelDist = np.hstack((travelDist,6371000*c))
#        
#        
#        
#        #Correct for some bug in the software. 
#        #if the travelDistance is not a number compute the meidan distance
#        if len(travelDist)>=4:
#            if travelDist[-2]==np.nan:
#                travelDist[-2]== (travelDist[-1]+travelDist[-3])/2
#    
#
#
#    #travelDist = distance between each point
#    #DistanceTraveled is the total length traveled for each ping
#    DistanceTraveled = travelDist
#
#    DistanceTraveled[np.where(np.isnan(DistanceTraveled))] = 0
#    DistanceTraveled = np.cumsum(DistanceTraveled)
        
    
    return DistanceTraveled, travelDist;


    
    
    
class GetVariablesFromNC(object):
    
    def __init__(self,fileID,beamgrp,Files,fileIDX):
        if type(beamgrp) == np.str: 
            bmgrp = beamgrp
        else:
            bmgrp = beamgrp[fileIDX]
        
        if Files == False:
            fIDX = int(fileIDX)
        else: 
            fIDX = int(Files[fileIDX,2])
        
        #Get beam sonar configuration info   
        try: 
            pingtime = fileID.groups['Sonar'].groups[bmgrp].variables['ping_time'][fIDX]
            self.time = pingtime
            self.frequency = fileID.groups['Sonar'].groups[bmgrp].variables['transmit_frequency_start'][fIDX]
            self.transmitpower = fileID.groups['Sonar'].groups[bmgrp].variables['transmit_power'][fIDX]
            self.pulslength = fileID.groups['Sonar'].groups[bmgrp].variables['transmit_duration_nominal'][fIDX]
            self.gaintx = fileID.groups['Sonar'].groups[bmgrp].variables['transducer_gain'][fIDX]
            self.gainrx = fileID.groups['Sonar'].groups[bmgrp].variables['receiver_sensitivity'][fIDX]
            self.sampleinterval = fileID.groups['Sonar'].groups[bmgrp].variables['sample_interval'][fIDX]
            self.equivalentbeamangle=fileID.groups['Sonar'].groups[bmgrp].variables['equivalent_beam_angle'][fIDX,:]
            self.sacorrection = 0
    
            #Get environment data
            self.soundvelocity = fileID.groups['Environment'].variables['sound_speed_indicative'][:]
            self.absorptioncoefficient = fileID.groups['Environment'].variables['absorption_indicative'][:]
    
          
            
            #Get beam configuration
            beam_direction_x=fileID.groups['Sonar'].groups[bmgrp].variables['beam_direction_x'][fIDX,:]
            beam_direction_y=fileID.groups['Sonar'].groups[bmgrp].variables['beam_direction_y'][fIDX,:]
            beam_direction_z=fileID.groups['Sonar'].groups[bmgrp].variables['beam_direction_z'][fIDX,:]


            self.dirx = np.arcsin(beam_direction_z)/np.pi*180
            self.diry = np.arctan2(beam_direction_y,beam_direction_x)*180/np.pi
            
    
            #Unpack and write beam data  
            BeamIM = fileID.groups['Sonar'].groups[bmgrp].variables['backscatter_i']
            BeamReal = fileID.groups['Sonar'].groups[bmgrp].variables['backscatter_r']
            BeamAmplitudeDataIM=UnpackBeam(BeamIM[fIDX,:])
            BeamAmplitudeDataReal=UnpackBeam(BeamReal[fIDX,:])
            
            self.BeamAmplitudeData =(BeamAmplitudeDataIM**2 + BeamAmplitudeDataReal**2)
            
            
            NMEA_time= fileID.groups['Platform'].variables[fileID.groups['Platform'].variables['longitude'].dimensions[0]][:]/100
            Latitude = fileID.groups['Platform'].variables['latitude'][:]
            Longitude = fileID.groups['Platform'].variables['longitude'][:]
            
            NMEA_idx = np.where(abs(NMEA_time-pingtime)==np.min(abs(NMEA_time-pingtime)))
            #Get NMEA data
            self.Longitude = Longitude[NMEA_idx]
            self.Latitude = Latitude[NMEA_idx]
            self.NMEA_time = NMEA_time[NMEA_idx]
            

        except IndexError: 
            k=1



