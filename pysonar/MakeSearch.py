# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 18:48:05 2018

@author: sindrev
"""


from netCDF4 import Dataset
import scipy, os
import numpy as np
from tools import tools





    
def MakeSearch(ListOfFilesWithinTimeInterval,RemoveToCloseValues,R_s,res,directory2Data,dirnc,beamgrp): 

    
    
    #Something for bookkeeping and progress
    #This should be cleaned
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

    
    
    
    
    #Loop through all files within the time interval
    for filename_index in range(0,len(ListOfFilesWithinTimeInterval[:,0])):
        
        #Print the progression
        tools.printProgressBar(filename_index + 1, len(ListOfFilesWithinTimeInterval[:,0]), prefix = 'Make SearchMatrix:', suffix = 'Complete', length = 50)
        
        
        #Get the full path name
        filename = os.path.join(dirnc,ListOfFilesWithinTimeInterval[filename_index,1])
        
        
        #Load the nc file
        fileID = Dataset(filename,'r',format = 'NETCDF4')
        
        
        #Get data from the nc file
        variables = tools.GetVariablesFromNC(fileID,beamgrp,ListOfFilesWithinTimeInterval,filename_index)
        
        
        #Stack the nmea position data
        NMEA_idx = np.where(abs(variables.NMEA_time[:]-int(ListOfFilesWithinTimeInterval[filename_index,0])) == 
                            np.min(abs(variables.NMEA_time[:]-int(ListOfFilesWithinTimeInterval[filename_index,0]))))
        
        lat = np.hstack((lat,variables.Latitude[NMEA_idx]))
        lon = np.hstack((lat,variables.Longitude[NMEA_idx]))
        
        
        #close the nc file
        fileID.close()
                      
        
        
        #Find the distance the vessel has traveled
        #ADD check if the function can be simplified, or if
        #a similar function is avaliable to be downloaded
        DistanceTraveled, travelDist = tools.ComputeDistance(travelDist,lat,lon)
        
        
           
        #A bug fix
        if len(lat)>5: 
            DT,td = tools.ComputeDistance(travelDist,
                                    np.array(lat[:1],lat[-1]),
                                    np.array(lon[:1],lon[-1]))
            DistanceTraveled = np.linspace(0,np.max(DT),len(lat))
            
        
        
        #Get the calibration gain and add it to the data
        #ADD calibration cain is not jet avaliable. 
        #the function may be chainged once tis is ready#
        if not NominalCalibraitonGain: 
            gain,FrequencyGain,PulslengthGain = tools.GainAdjustment(variables.pulslength*1E3,variables.frequency/1E3,variables.gaintx+variables.gainrx) 
            gain = variables.gaintx+variables.gainrx+PulslengthGain
        else: 
            gain,FrequencyGain,PulslengthGain = tools.GainAdjustment(variables.pulslength*1E3,variables.frequency/1E3,NominalCalibraitonGain) 
            
            
            
            
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
#        if dirx[0][0]>3:       
            
            
        #Get the timestamp of the ping
        #Needed when generating work file
        TimeStamp = np.hstack((TimeStamp,filename))
        TiltVec = np.hstack((TiltVec,variables.dirx[0]))
        
        
        if variables.dirx[0]<np.nanmedian(TiltVec):
            print('Tilt ble endret',end='\r')
            sv = np.nan*np.ones(sv.shape)
        
        sv[:,np.where(abs(variables.diry)>165)]=np.nan 
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
            BeamDirectionMatrix = np.repeat(variables.diry.T,len(sv),axis=0)
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
                Wdist_port,Wdist_stb = tools.GetDistanceMatrix(DistanceMatrix,
                                                         RangeMatrix,
                                                         BeamDirectionMatrix,
                                                         svMatrix,
                                                         int(variables.variables.dirx[0]+90)
                                                         ,BananaTool)
                
                
                
                #Start making the ghost school matrix
                print('Generating Distance Matrix for ghost: ',end='\r')
                Wdist_portGhost,Wdist_stbGhost = tools.GetDistanceMatrix(DistanceMatrix,RangeMatrix,BeamDirectionMatrix+90,svMatrix,int(variables.dirx[0]+90),BananaTool)




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
                sV_port = np.hstack((sV_port,tools.ConvertToechogram(Wdist_port[indeks],sv_mat)))
                sV_stb = np.hstack((sV_stb,tools.ConvertToechogram(Wdist_stb[indeks],sv_mat)))
                sV_portGhost = np.hstack((sV_portGhost,tools.ConvertToechogram(Wdist_portGhost[indeks],sv_mat)))
                sV_stbGhost = np.hstack((sV_stbGhost,tools.ConvertToechogram(Wdist_stbGhost[indeks],sv_mat)))
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
            BeamDirectionMatrix = np.dstack((BeamDirectionMatrix, np.repeat(variables.diry.T,len(sv),axis=0)))
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

