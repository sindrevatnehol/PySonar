# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 18:48:05 2018

@author: sindrev
"""


from netCDF4 import Dataset
import scipy
import numpy as np
from tools import tools




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
    samplespace =soundvelocity[0]*sampleinterval[0]/2
    
    
    #Declare the range vector to be outputed
    Range = np.arange(0,AmplitudeData[:].shape[0])*samplespace-rangeCorr
                      
                      
                      
    #Remove bad values  
    Range[np.where(Range<0)] = 1E-200
               
    #Make the TVG function (20 log) as a vector
    tvg20 = np.repeat(20*np.log10(Range[:, np.newaxis]) + 2*absorptioncoefficient[0]*Range[:, np.newaxis], AmplitudeData[:].shape[1], axis=1)
                         
                      
    
    #Correction function for the loss off energy when steering the beam
    tiltcorr =40*np.log10(abs(np.cos(np.deg2rad(-tiltAngle))))
    
    
    
    #compute the wavelength of the signal
    wavelength = soundvelocity/frequency

    
    
    #compute the correction function used when computing the SV
    svconstSV =10*np.log10(transmitpower*wavelength**2*soundvelocity*effective_pulselength/(32*np.pi**2))+tiltcorr+gain+equivalentbeamangle+sacorrection
    
    
    #Compute the volume backscattering coefficient
    sv= AmplitudeData+tvg20-np.repeat(svconstSV,AmplitudeData[:].shape[0],axis=1).transpose()

    
    
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
        tools.printProgressBar(r_ind_stbex+1,int(BananaTool[2]), prefix = 'Make Distance:', suffix = 'Complete', length = 50)
        
        
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




def NMEAconverter(GPS):
    lat = float(GPS)
    lat = np.floor(lat) + (np.floor((lat-np.floor(lat))*100))/60 + (((lat-np.floor(lat))*100-np.floor((lat-np.floor(lat))*100))/100)
    
    
    
def ComputeDistance(travelDist,lat,lon): 
    '''
    This function computes the distance traveled within each transect. 
    For distance the hamorsine function is used. 
    
    Parts of the code is awkward due to bugs in the files. Will be changed later. 
    '''
    
    
    #We need some points to evaluate the distanc traveled. 
    if len(lat)==1: 
        travelDist=np.hstack((travelDist,np.nan))
        
        
        
    elif len(lat)>=20:
        print(lat)
        print(lon) 
        #find the delta in longitude and latitude
        delta_lat = lat[-1]-lat[-2]
        delta_lon = lon[-1]-lon[-2]



        #Computing the distance using hammersine fucntion
        a = (np.sin(delta_lat/2*np.pi/180))**2 + np.cos(lat[-2]*np.pi/180) * np.cos(lat[-1]*np.pi/180) * (np.sin(delta_lon/2*np.pi/180)**2)
        c = 2 * np.arctan2(np.sqrt(a),np.sqrt(1-a))
        
        
        #Compute the distance to meter and stack the variable
        travelDist = np.hstack((travelDist,6371000*c))
        
        
        
        #Correct for some bug in the software. 
        #if the travelDistance is not a number compute the meidan distance
        if len(travelDist)>=4:
            if travelDist[-2]==np.nan:
                travelDist[-2]== (travelDist[-1]+travelDist[-3])/2
    


    #travelDist = distance between each point
    #DistanceTraveled is the total length traveled for each ping
    DistanceTraveled = travelDist

    DistanceTraveled[np.where(np.isnan(DistanceTraveled))] = 0
    DistanceTraveled = np.cumsum(DistanceTraveled)
        
    
    return DistanceTraveled, travelDist;


    
    
def MakeSearch(ListOfFilesWithinTimeInterval,RemoveToCloseValues,R_s,res,directory2Data,dirnc): 

    PingCount = True
    MakeWdistStuff= True
    MaximumDetectionRange = 10000
    NominalCalibraitonGain = []  

    lat = np.array([])
    lon = np.array([])    
    travelDist = np.array([])
    TimeStamp = np.array([])
    svMatrix = np.array([])
    DistanceMatrix = np.array([])
    BeamDirectionMatrix = np.array([])
    RangeMatrix = np.array([])
    TiltVec = np.array([])

    
    
    #Report to the user that files are being loaded
            
    #Loop through all files within the time interval
    for filename_index in range(0,len(ListOfFilesWithinTimeInterval)): 
        tools.printProgressBar(filename_index + 1, len(ListOfFilesWithinTimeInterval), prefix = 'Make SearchMatrix:', suffix = 'Complete', length = 50)
        
        #Name of the current file in the short list
        CurrentFileName = ListOfFilesWithinTimeInterval[filename_index].replace(' ','')
        import os 
        os.chdir(dirnc)
        #Open the NETCDF file
        try: 
            fileID = Dataset(CurrentFileName,'r',format='NETCDF4')
        except OSError:
            fileID = Dataset(CurrentFileName+'.nc','r',format='NETCDF4')
        
        
        
        #Get all the relevant data from NETCDF file
        #ADD With a new netcdf format, this protocoll will 
        #be changed. It should then identify the correct location
        #within the file
        BeamAmplitudeData = fileID.variables['BeamAmplitudeData'][:]
        soundvelocity = fileID.variables['soundvelocity'][:]
        sampleinterval = fileID.variables['sampleinterval'][:]
        absorptioncoefficient = fileID.variables['absorptioncoefficient'][:]
        dirx = fileID.variables['dirx'][:]
        diry = fileID.variables['diry'][:]
        dirz = fileID.variables['dirz'][:]
        frequency = fileID.variables['frequency'][:]
        transmitpower = fileID.variables['transmitpower'][:]
        pulslength = fileID.variables['pulslength'][:]
        gaintx = fileID.variables['gaintx'][:]
        gainrx = fileID.variables['gainrx'][:]
        equivalentbeamangle = fileID.variables['equivalentbeamangle'][:]
        sacorrection = fileID.variables['sacorrection'][:]
        Longitude = fileID.variables['Longitude'][:]
        Latitude = fileID.variables['Latitude'][:]
        fileID.close()
                      
        
        
        
        #Convert NMEA gps data to degrees
        #ADD This will be replaced with a proper function 
        #using the telegram itselfe. Will be done with
        #the new netcdf location!!!
        lat = np.hstack((lat,NMEAconverter(str(Latitude[0][0]))))
        lon = np.hstack((lon,NMEAconverter(str(Longitude[0][0]))))   
        
        
        #Correct for sonar alignment
        #The rotation helps to investigate ghost schools
        #ADD this will be included in the new netcdf format
        #Therefore it should be deleted once that is made!!!
        diry = diry+dirz+0
            
        
        
        
        #Correction of angles
        #ADD Not necessary and will be deleted once 
        #properly tested  !!!
        diry[np.where(diry<-180)]=360-diry[np.where(diry<-180)]
        diry[np.where(diry>180)]=diry[np.where(diry>180)]-360
        
             
        
        #Find the distance the vessel has traveled
        #ADD check if the function can be simplified, or if
        #a similar function is avaliable to be downloaded
        DistanceTraveled, travelDist = ComputeDistance(travelDist,lat,lon)
        
        
           
        #A bug fix
        if len(lat)>5: 
            DT,td = ComputeDistance(travelDist,
                                    np.array(lat[:1],lat[-1]),
                                    np.array(lon[:1],lon[-1]))
            DistanceTraveled = np.linspace(0,np.max(DT),len(lat))
            
        
        
#            print(np.nanmax(DistanceTraveled))
        #Get the calibration gain and add it to the data
        #ADD calibration cain is not jet avaliable. 
        #the function may be chainged once tis is ready#
        if not NominalCalibraitonGain: 
            gain,FrequencyGain,PulslengthGain = GainAdjustment(pulslength[0][0]*1E3,frequency[0][0]/1E3,gaintx[0]+gainrx[0]) 
            gain = gaintx[0]+gainrx[0]+PulslengthGain
        else: 
            gain,FrequencyGain,PulslengthGain = GainAdjustment(pulslength[0][0]*1E3,frequency[0][0]/1E3,NominalCalibraitonGain) 
            
            
            
            
        #The sonar data often includes corrputed value of 
        #the transmit power that destroys the analysis. 
        #This will fix this problum, but the sv values are not
        #correct. 
        #ADD this data should be labeled when making the work files
        #so the user can now that it is corrupted.
        if transmitpower == 0: 
            transmitpower = 4633
            
            
            
            
        #Compute the sv and TS 
        #ADD TS are not used here !!!
        sv, RangeOut= ApplyTVG(BeamAmplitudeData,
                                soundvelocity,
                                sampleinterval,
                                transmitpower,
                                absorptioncoefficient,
                                frequency,
                                pulslength,
                                gain,
                                equivalentbeamangle,
                                sacorrection,
                                dirx)
            
        
        
        
        #Remove data too close to the vessel
        sv[np.where(RangeOut<=RemoveToCloseValues)] = np.nan

        
           
        
        #ADD Temporary filter that will be deleted once checked !!!
        if np.nanmean(sv)>0: 
            sv = np.nan*np.ones(sv.shape)
           
        
                   
                   
        #Protocoll to identify if the batch of data is 
        #finished loaded
        if PingCount == True:
            NumberOfPingsInBatch = 5000
            if DistanceTraveled.max()>=MaximumDetectionRange:
                NumberOfPingsInBatch=len(DistanceTraveled)
                PingCount = False  


                
        
        #Some vertical pings are labeled as horizontal
        #This solution is a temporarly fix that check the 
        #direction of the beams. 
        #ADD it does not always work this part of the sonar
        #data are quite unstable. Find a better fix if 
        #possible !!!
        if dirx[0][0]>3:       
            
                
            #Get the timestamp of the ping
            #Needed when generating work file
            TimeStamp = np.hstack((TimeStamp,CurrentFileName[5:-6]))
            TiltVec = np.hstack((TiltVec,dirx[0][0]))
            
            
            if dirx[0][0]<np.nanmedian(TiltVec):
                print('Tilt ble endret',end='\r')
                sv = np.nan*np.ones(sv.shape)
            
            sv[:,np.where(abs(diry)>165)]=np.nan 
		#Remov noise from vessel vake  
            #sv[np.where(1)]
            
            #Implement the inteferance removal filter
            #ADD it is still under development and must be 
            #properly tested !!!
            #Uncommet to remove spikes from the data
            #ADD a more proper spike detection and removal must be made in the future !!!
#                                if SpikeIDX.shape[0] == 0: 
#                                    SpikeIDX = np.nanmean(sv,axis=1)[:,np.newaxis]
#                                elif SpikeIDX.shape[1] <=5: 
#                                    SpikeIDX = np.hstack((SpikeIDX,np.nanmean(sv,axis=1)[:,np.newaxis]))
#                                    
#                                else:
#                                    SpikeIDX = np.hstack((SpikeIDX[:,1:],np.nanmean(sv,axis=1)[:,np.newaxis]))
#                                   
#                                spike = np.where(( np.nanmean(sv,axis=1))>( np.nanmean(SpikeIDX,axis=1)+6))
#                                sv[spike,:] = np.nan
#                                    
                
            
            
            #Start the protocol for making the search matrix
#                if filename_index == 0 or filename_index == 1: 
                
                
                
                
                #Set the variable size to define the 
                #trajectory lines
            
            #print(np.nanmax(RangeOut))#MaximumDetectionRange) 
                
                
            #If the first file of the transect
            if len(svMatrix) == 0: 
                
                
                BananaTool = [R_s,0,len(RangeOut),res]
                MaximumDetectionRange = 2*np.nanmax(RangeOut)
                
                
                #Start making the buffer matrix
                svMatrix = sv
                DistanceMatrix = np.ones((len(sv),64))*DistanceTraveled[-1]
                BeamDirectionMatrix = np.repeat(diry.T,len(sv),axis=0)
                RangeMatrix=np.repeat(RangeOut[:,np.newaxis],64,axis=1) 
                
                
            
                
            #If the buffer size is large enough 
            elif filename_index >= NumberOfPingsInBatch:
                
                
                
                
                #if first ping after buffersize is large enough. 
                #ADD make the if sentence more simpler by
                #Only have the makeWdiststuff !!!
#                    if ((filename_index == NumberOfPingsInBatch) or (filename_index == (NumberOfPingsInBatch+1)))and MakeWdistStuff == True: 
                if MakeWdistStuff == True: 
                    
                    #Start making the distance matrix for 
                    #the transect
                    #ADD needs a way to correct this
                    #if the vessel speed or direction has 
                    #changed    !!!
                    print('Generating Distance Matrix: ',end='\r')
                    Wdist_port,Wdist_stb = GetDistanceMatrix(DistanceMatrix,
                                                             RangeMatrix,
                                                             BeamDirectionMatrix,
                                                             svMatrix,
                                                             int(dirx[0][0]+90)
                                                             ,BananaTool)
                    
                    
                    
                    #Start making the ghost school matrix
                    print('Generating Distance Matrix for ghost: ',end='\r')
                    Wdist_portGhost,Wdist_stbGhost = GetDistanceMatrix(DistanceMatrix,RangeMatrix,BeamDirectionMatrix+90,svMatrix,int(dirx[0][0]+90),BananaTool)




                    #This will make the program to stopp 
                    #increasing the buffer 
    
    

                    
                    
                #Correct the size of the sv and range
                #ADD check if it is necesarry or if it is 
                #already fixed above !!!
                AddNaN2Matrix = len(svMatrix[:,0])-len(RangeOut)
#                    sv = np.vstack((sv,np.nan*np.ones((AddNaN2Matrix,64))))
                
                if AddNaN2Matrix >0:
                    sv = np.vstack((sv,np.nan*np.ones((AddNaN2Matrix,64))))
                    RangeOut = np.vstack((RangeOut[:,np.newaxis],np.nan*np.ones((AddNaN2Matrix,1))))
                    
                else: 
                    RangeOut = RangeOut[:,np.newaxis]
                    sv = sv[:len(svMatrix[:,0])]



                
                
                #Give a progress report for the user
                #ADD should also include transect number
                print('        Progress:  '+ 
                      str(((filename_index+1)/len(ListOfFilesWithinTimeInterval))*100)+' %',end='\r')
                
                
                
                
                
                #Add new data and delete the redundent ping
                #to the buffer matrix
                svMatrix = np.dstack((svMatrix,sv))[:,:,1:]

                #ADD sjekk om dette er nødvendig !!!
                DistanceMatrix = np.dstack((DistanceMatrix, np.ones((len(sv),64))*DistanceTraveled[-1]))[:,:,1:]
#                BeamDirectionMatrix = np.dstack((BeamDirectionMatrix, np.repeat(diry.T,len(sv),axis=0)))[:,:,1:]
#                RangeMatrix = np.dstack((RangeMatrix, np.repeat(RangeOut,64,axis=1)))[:,:,1:]

                

                #print(np.nanmax(DistanceTraveled))                               

                #Transform into a 1D array, and make it linear
                sv_mat = 10**(np.reshape(svMatrix,(-1,1))/10)
                

                #sV_port = (Parallel(n_jobs=1)(delayed(ConvertToechogram)(Wdist_port[indeks],sv_mat) for indeks in range(len(range(int(BananaTool[1]),int(BananaTool[2]),int(BananaTool[3]))))))
             
                #sV_port2 = np.asarray(sV_port)
                sV_stb = []
                sV_port = []
                sV_portGhost = []
                sV_stbGhost = []
                for indeks in range(len(range(int(BananaTool[1]),int(BananaTool[2]),int(BananaTool[3])))): 
                    sV_port = np.hstack((sV_port,ConvertToechogram(Wdist_port[indeks],sv_mat)))
                    sV_stb = np.hstack((sV_stb,ConvertToechogram(Wdist_stb[indeks],sv_mat)))
                    sV_portGhost = np.hstack((sV_portGhost,ConvertToechogram(Wdist_portGhost[indeks],sv_mat)))
                    sV_stbGhost = np.hstack((sV_stbGhost,ConvertToechogram(Wdist_stbGhost[indeks],sv_mat)))
               # print(sV_port-sV_port2)
                #sV_port = np.asarray(sV_port)
                #sV_stb = (Parallel(n_jobs=1)(delayed(ConvertToechogram)(Wdist_stb[indeks],sv_mat) for indeks in range(len(range(int(BananaTool[1]),int(BananaTool[2]),int(BananaTool[3]))))))
               # sV_stb = np.asarray(sV_stb)
  
                
                
                #Making the ghost school search matrix
                #on both sides. 
                #ADD check if this can replace the call of 
                #function in the run.py  !!!
               # sV_portGhost = np.asarray(Parallel(n_jobs=1)(delayed(ConvertToechogram)(Wdist_portGhost[indeks],sv_mat) for indeks in range(len(range(int(BananaTool[1]),int(BananaTool[2]),int(BananaTool[3]))))))
                #sV_stbGhost = np.asarray(Parallel(n_jobs=1)(delayed(ConvertToechogram)(Wdist_stbGhost[indeks],sv_mat) for indeks in range(len(range(int(BananaTool[1]),int(BananaTool[2]),int(BananaTool[3]))))))
               
                


                #Add dimension on the data, a bug fix
                #ADD also add the ghost school   !!!
                sV_port = sV_port[:,np.newaxis]
                sV_stb = sV_stb[:,np.newaxis]
                sV_portGhost = sV_portGhost[:,np.newaxis]
                sV_stbGhost = sV_stbGhost[:,np.newaxis]
                

                                

                #Start making the the search matrix
                #ADD replace the stacking rutine with 
                #somehting nicer.  !!!
#                    if (filename_index == NumberOfPingsInBatch) or (filename_index == (NumberOfPingsInBatch+1)):     
                if MakeWdistStuff == True: 
                    
                    
                    #Make the first value in search matrix
                    #ADD in future make it possible to 
                    #analyze incomplete buffer.  !!!
                    #ADD also add the chost school  !!!
                    SVres_port = sV_port
                    SVres_stb = sV_stb
                    SVres_portGhost = sV_portGhost
                    SVres_stbGhost = sV_stbGhost
                    MakeWdistStuff = False
                    
                    
                    
                else: 
                    
                    #add new value to the search matrix
                    #ADD also add the chost school  !!
                    sV_port = np.vstack((sV_port,np.nan*np.ones(len(SVres_port[:,0])-len(sV_port[:,0]))[:,np.newaxis]))
                    sV_stb = np.vstack((sV_stb,np.nan*np.ones(len(SVres_stb[:,0])-len(sV_stb[:,0]))[:,np.newaxis]))
                    sV_portGhost = np.vstack((sV_portGhost,np.nan*np.ones(len(SVres_portGhost[:,0])-len(sV_portGhost[:,0]))[:,np.newaxis]))
                    sV_stbGhost = np.vstack((sV_stbGhost,np.nan*np.ones(len(SVres_stbGhost[:,0])-len(sV_stbGhost[:,0]))[:,np.newaxis]))
                        
                        
                    SVres_port = np.hstack((SVres_port,sV_port))
                    SVres_stb = np.hstack((SVres_stb,sV_stb))
                    SVres_portGhost = np.hstack((SVres_portGhost,sV_portGhost))
                    SVres_stbGhost = np.hstack((SVres_stbGhost,sV_stbGhost))
               
                      
                   # print(DirectoryToRESULT+'/SearchMatrix'+str(transectCode)+'.mat')
                   # scipy.io.savemat(DirectoryToRESULT+'/SearchMatrix'+str(transectCode)+'.mat',
                   #                    mdict={'SVres_port': (SVres_port),
                   #                   'SVres_stb':(SVres_stb),
                   #                   'SVres_portGhost':(SVres_portGhost),
                   #                   'SVres_stbGhost':(SVres_stbGhost), 
                   #                   'DistanceTraveled':DistanceTraveled,
                   #                   'ListOfFilesWithinTimeInterval':ListOfFilesWithinTimeInterval})
            #If the buffer size is insufficient, keep stacking
            else:
               
                AddNaN = len(svMatrix[:,0])-len(RangeOut)
                
                if AddNaN >0:
                    sv = np.vstack((sv,np.nan*np.ones((AddNaN,64))))
                    RangeOut = np.vstack((RangeOut[:,np.newaxis],np.nan*np.ones((AddNaN,1))))
                else: 
                    sv=sv[:svMatrix.shape[0],:]
                    RangeOut = RangeOut[:svMatrix.shape[0],np.newaxis]

                
                svMatrix = np.dstack((svMatrix,sv))
                DistanceMatrix = np.dstack((DistanceMatrix, np.ones((len(sv),64))*DistanceTraveled[-1]))
                BeamDirectionMatrix = np.dstack((BeamDirectionMatrix, np.repeat(diry.T,len(sv),axis=0)))
                RangeMatrix = np.dstack((RangeMatrix, np.repeat(RangeOut[:,np.newaxis],64,axis=1)))
                
                

    #Save the data
    try: 
        scipy.io.savemat(directory2Data, 
                     mdict={'SVres_port': (SVres_port),

                     'SVres_stb':(SVres_stb),'R_s':R_s,'res':res,
                     'SVres_portGhost':(SVres_portGhost),
                     'SVres_stbGhost':(SVres_stbGhost),
                     'DistanceTraveled':DistanceTraveled,
                     'ListOfFilesWithinTimeInterval':ListOfFilesWithinTimeInterval})
    except UnboundLocalError: 
        print('empty files',end='\r')

