''' 
MakeWork.py

Last modify date: 17.10.2017


Description: 
This script is made only for the development of the algorithm
It will be replaced with the MakePingWork.py script


Currently this will make the work matrix based on the information from the 
search matrix

Than it will plot all the figures for publication 

    
    
    
Author: 
Dr. Sindre Vatnehol (PhD)
Institute of Marine Research, Norway

Mail: 
sindre.vatnehol@imr.no

Tlf.: 
+47 900 79 376


Project: 
REDUS (Reducing uncertainty in stock assessment)  
www.redus.no
'''


import scipy.io as sc
import scipy as scpy
import numpy as np
#from sklearn.cluster import DBSCAN
import os
from tools import tools
from glob import glob
from MakeIndex import MakeIndex
from netCDF4 import Dataset


    
    
    
    
    
    
    
def FindLocationOfSchool(CatMat,BananaTool,DistanceTraveled):  
    
    '''Find the location of pixels'''
    if len(BananaTool)==1:
        BananaTool = BananaTool[0]
    
    

    #Get the scrutinized search matrix
    r_ind ,ping_ind = np.where(CatMat>0)


    y__school_traj = np.arange(0,1000,BananaTool[0][0])[r_ind] #y__school_traj.T 
    

    x__school_traj = DistanceTraveled[0,ping_ind]


    return x__school_traj, y__school_traj
                                        
    
    
def GetSpikesFromEchogram(SVres,threshold,S): 
    '''Purpose: 
        Threshold the data, and get index of the leftover
        ''' 
    
    #Make the cathegorization matrix
    CathegorizationMatrix = np.zeros(SVres.shape)
    

    median = np.nanmedian(SVres,axis=1)
    median_filter = (np.repeat(median[:,np.newaxis],SVres.shape[1],axis=1))
    
    #Find all pixels likly to have fish inside it (larger then background noise)
    CathegorizationMatrix[np.where(SVres>(median_filter*threshold))] = 1
    
    
    x_i,y_i = np.where(CathegorizationMatrix>0)

    Y_vec = -x_i[:,np.newaxis]
    X_vec = y_i[:,np.newaxis]
                                        
    #Return index file                                  
    return CathegorizationMatrix,np.hstack((X_vec,Y_vec)),
    





#def GetSpikesFromEchogramGhost(SVres,threshold): 
#    '''Purpose: 
#        Threshold the data, and get index of the leftover
#        ''' 
#    
#    #Make the cathegorization matrix
#    CathegorizationMatrix = np.zeros(SVres.shape)
#    
#    
#    #Find all pixels likly to have fish inside it (larger then background noise)
#    CathegorizationMatrix[np.where(SVres>=threshold)] = 1
#    
#    x_i,y_i = np.where(CathegorizationMatrix>0)
#
#    Y_vec = -x_i[:,np.newaxis]
#    X_vec = y_i[:,np.newaxis]
#    
#      
#    #Return index file                                  
#    return CathegorizationMatrix,np.hstack((X_vec,Y_vec)), 
    








    
    
#The main function    
#MakeWork(True, directory2Data,'','','',1,3)
def MakeWork(makeNewWork,data_directory,SearchMatrixName,
             prefix,compensate,clusterTH,Threshold):
    
    
    
    
    #Sjekk om dett er cluster TH eller Threshold
    CatThr = 3

    
    
    
    
    #Get all mat files
    os.chdir(data_directory.dir_search)
    FileNames = glob('*.mat')
    
    
    
    #Loop through each search matrix (transect)
    for SearchMatrixName in FileNames: 
        
        
        if not os.path.isfile(data_directory.dir_work+'/Horizontal_'+SearchMatrixName) == True:
            
            
            
            #Load search matrix
            SearchMatrix = sc.loadmat(data_directory.dir_search+'/'+SearchMatrixName)
        
            
            
            #Get the data from the files
            SVres_port = SearchMatrix['SVres_port']
            SVres_stb = SearchMatrix['SVres_stb']
            R_s = SearchMatrix['R_s']
            res = SearchMatrix['res']
            SVres_portGhost = SearchMatrix['SVres_portGhost']
            SVres_stbGhost = SearchMatrix['SVres_stbGhost']
            DistanceTraveled = SearchMatrix['DistanceTraveled']
            ListOfFilesWithinTimeInterval = SearchMatrix['ListOfFilesWithinTimeInterval']
        
    
    
    
    
            
            # Finding the location of the schools in the search matrix
            # ADD may be replaced with an other function
            Cat_port, X_port = GetSpikesFromEchogram(SVres_port,CatThr,
                                        np.hstack((np.reshape(SVres_port,(-1,)),
                                        np.reshape(SVres_stb,(-1,)))))     
            Cat_stb, X_stb = GetSpikesFromEchogram(SVres_stb,CatThr,
                                        np.hstack((np.reshape(SVres_port,(-1,)),
                                        np.reshape(SVres_stb,(-1,)))))
        
            
            
            
            
            
            #Get the location of each school in xy coordinates
            Bananatool = [R_s,0,0,res]
            x__school_traj_port,y__school_traj_port = FindLocationOfSchool(Cat_port, Bananatool, DistanceTraveled)
            x__school_traj_stb,y__school_traj_stb = FindLocationOfSchool(Cat_stb, Bananatool, DistanceTraveled)
                
            
            
                
            
            
            #Start writing work files.
            #Go through each ping
            WorkFileOutput = np.array([])
            for ii in range(len(DistanceTraveled.T)): 
                
                
                
                #Progress for user
                tools.printProgressBar(ii + 1, len(ListOfFilesWithinTimeInterval[:,0]), prefix = 'Scanning files:', suffix = len(DistanceTraveled.T)-ii, length = 50)
                
                
                
                #If first ping
                if ii == 0: 
                    
                    
                    #Open the file and get variables
                    fileID = Dataset(data_directory.dir_rawdata+'/'+ ListOfFilesWithinTimeInterval[ii][1],'r',format = 'NETCDF4')
                    variables = tools.GetVariablesFromNC(fileID,'Beam_group1',ListOfFilesWithinTimeInterval,ii)
                    fileID.close()
                    
                    
                    
                    #Compute the range
                    rangeCorr = 3; 
                    samplespace =variables.soundvelocity*variables.sampleinterval/2
                    Range = np.arange(0,len(variables.BeamAmplitudeData))*samplespace-rangeCorr
    
    
                        
                    #Compute the pixel location with direction and range
                    X =np.reshape(np.dot(Range[:,np.newaxis], np.cos(np.deg2rad(variables.diry))[:,np.newaxis].T),(-1))
                    Y =np.reshape(np.dot(Range[:,np.newaxis], np.sin(np.deg2rad(variables.diry))[:,np.newaxis].T),(-1))
                    Phi = np.reshape(np.dot(np.zeros(Range[:,np.newaxis].shape),(variables.diry)[:,np.newaxis].T),(-1))
                    Range = np.reshape(np.dot(Range[:,np.newaxis],np.zeros((variables.diry)[:,np.newaxis].T.shape)),(-1))
        
    
    
                
                    
                #Skipp those pings that don't have data
                if np.nansum([(x__school_traj_port<=np.nanmax(X)+DistanceTraveled[0,ii])&(x__school_traj_port>=np.nanmin(X)+DistanceTraveled[0,ii])]) >0: 
                    
                    
                    #go through each pixel and find if it is relevant data
                    for iik in range(len(X)): 
                        r_stb = np.sqrt((x__school_traj_stb-(X[iik]+DistanceTraveled[0,ii]))**2+(y__school_traj_stb-Y[iik])**2)
                        
                        
                        if np.min(r_stb)<=R_s[0]: 
                            if len(WorkFileOutput) == 0: 
                                WorkFileOutput=np.hstack((int(ListOfFilesWithinTimeInterval[ii][0]),Phi[iik],Range[iik]))
                            else: 
                                WorkFileOutput = np.vstack((WorkFileOutput,np.hstack((int(ListOfFilesWithinTimeInterval[ii][0]),Phi[iik],Range[iik]))))
                            
                            
                        r_port = np.sqrt((x__school_traj_port-(X[iik]+DistanceTraveled[0,ii]))**2+(y__school_traj_port-Y[iik])**2)
                        if np.min(r_port)<=R_s[0]: 
                            
                            if len(WorkFileOutput) == 0: 
                                WorkFileOutput=np.hstack((int(ListOfFilesWithinTimeInterval[ii][0]),Phi[iik],Range[iik]))
                            else: 
                                WorkFileOutput = np.vstack((WorkFileOutput,np.hstack((int(ListOfFilesWithinTimeInterval[ii][0]),Phi[iik],Range[iik]))))
                            
                            
            sc.savemat(data_directory.dir_work+'/Horizontal_'+SearchMatrixName,mdict={'WorkFile':WorkFileOutput})
            
            
            
            
            
            
            
            try: 
                import matplotlib.pyplot as plt
        
                plt.figure(1)
                plt.clf()
                plt.imshow((Cat_port),aspect = 'auto')
                plt.colorbar()
                plt.draw()
                plt.savefig(data_directory.dir_search+'/'+SearchMatrixName.replace('.mat','2.jpg'))
                
                
                
                plt.figure(2)
                plt.clf()
                plt.plot(x__school_traj_port,y__school_traj_port,'r.')
                plt.plot(x__school_traj_stb,-y__school_traj_stb,'b.')
                plt.draw()
                plt.savefig(data_directory.dir_search+'/'+SearchMatrixName.replace('.mat','3.jpg'))
            except: 
                k=1

        



#
#    
#def processinnput(Work,x__school_traj,X_traj,y__school_traj,Y_traj,case,corr): 
#    
#    if len(case['BananaTool'])==1:
#        BananaTool = case['BananaTool'][0]
#    else: 
#        BananaTool = case['BananaTool']
#
#
#    for ii in range(len(x__school_traj)):   
#        R = (x__school_traj[ii]-X_traj)**2+(y__school_traj[ii]-Y_traj)**2
#        if ii == 0: 
#            r_test = R
#        else: 
#            r_test = np.minimum(r_test,R)
#    
#    x_i,y_i,z_i = np.where(np.sqrt(r_test)<=(int(BananaTool[0])+corr))
#    Work[x_i,y_i,z_i]=1
#
#
#
#    return Work
#        
#    
#    
#    
#    
#def MakeIndex(case,Work,x__school_traj_port,y__school_traj_port,X_traj,Y_traj,
#              x__school_traj_stb,y__school_traj_stb,filename): 
#    '''Make index file '''
#    
#    
#    #Correction function when the distance of beams increases
#    
#    corr =  case['MakeRange']*np.sin(2*np.pi/64)
#
#    print('    Start (1/2)')
#    Work = processinnput(Work,x__school_traj_port,X_traj,y__school_traj_port,Y_traj,case,corr)
#
#    
#    print('    Start (2/2)')
#    Work = processinnput(Work,x__school_traj_stb,X_traj,y__school_traj_stb,Y_traj,case,corr)
#
#    #Save the work file for future use
#    sc.savemat(filename,mdict={'Work':Work})
#    
    
    
#
#def ComputePearsonR(SA_profos,SA_profos2,SA_AutoWithoutGhost):
#    #Compute the Pearson r between the two models
#
#    sap = np.zeros(len(SA_profos))
#    sap2 = np.zeros(len(SA_profos))
#    saa = np.zeros(len(SA_profos))
#    for ik in range(len(SA_profos)): 
#        if ik == 0: 
#            sap[ik] = SA_profos[ik]
#            sap2[ik] = SA_profos2[ik]
#            saa[ik] = SA_AutoWithoutGhost[ik]
#        else:   
#            sap[ik] = SA_profos[ik]-SA_profos[(ik-1)]
#            sap2[ik] = SA_profos2[ik]-SA_profos2[(ik-1)]
#            saa[ik] = SA_AutoWithoutGhost[ik]-SA_AutoWithoutGhost[(ik-1)]
#            
#    PR_Before = scpy.stats.pearsonr(sap,saa)
#    PR_After = scpy.stats.pearsonr(sap2,saa)
#    
#    return PR_Before, PR_After
#    
    
    
        

#def WaveRemover(dim,pitch,roll,limit):
#    #This function identify when there is a chanse that, due to vessel movement, 
#    #there is a chanse of surface reveberation.
#    #It will then remove these pings
#    #ADD could be included in the MakeSearchMatrix.py  !!!
#    
#    WaveRemover = np.ones(dim)
#    WaveRemover[:,:,np.where(abs(pitch-np.mean(pitch))>limit)] = np.nan
#    WaveRemover[:,:,np.where((roll-np.mean(roll))>limit)] = np.nan
#
#    return WaveRemover


    
    
    

#def ProtocolForCluster(X_port,X_port_ghost,X_stb,X_stb_ghost,Bananatool,
#                      MakeDist,Cat_port,min_samples):
#    '''
#    This is the protocol for making the cluster
#    
#    In future this will change so it don't merge small schools closly togetter in range direction
#    '''
#    
#    
#    #Find only data without ghost
#    #Bruk heller x__xchool_trajS
#    X_port_G = np.zeros((1,2))
#    for i in np.arange(len(X_port[:,0])-1,0,-1):
#        for ii in range(len(X_port_ghost[:,0])-1,0,-1): 
#            if (X_port_ghost[ii,0]==X_port[i,0]): 
#                if (X_port_ghost[ii,1]==X_port[i,1]): 
#                    X_port_G = np.vstack((X_port_G,X_port_ghost[ii,:]))
#    X_port_G = X_port_G[1:,:]
#
#
#    X_stb_G = np.zeros((1,2))
#    for i in np.arange(len(X_stb[:,0])-1,0,-1):
#        for ii in range(len(X_stb_ghost[:,0])-1,0,-1): 
#            if (X_stb_ghost[ii,0]==X_stb[i,0]): 
#                if (X_stb_ghost[ii,1]==X_stb[i,1]): 
#                    X_stb_G = np.vstack((X_stb_G,X_stb_ghost[ii,:]))
#    X_stb_G = X_stb_G[1:,:]
#
#
#
#    '''Make cluster out of the data'''
#    
#    
#    X_port_G = X_port_G.astype(int)
#    X_stb_G = X_stb_G.astype(int)
#    
#    
#    if len(Bananatool)==1:
#        BananaTool = Bananatool[0]
#    else: 
#        BananaTool = Bananatool
#        
#        
#        
#        
#    
#    #Convert range data to m
#    #The constant 0.18 is the sample distance
#    X_port_G[:,1] = abs(X_port_G[:,1])*int(BananaTool[3])*0.18
#    X_stb_G[:,1] = X_stb_G[:,1]*int(BananaTool[3])*0.18
#    
#
#
#
#    #Help to find the correct travel distance of the pixel
#    Distance = np.unique(MakeDist)
#    StartFrom = int((len(Distance)-Cat_port.shape[1])/2)+1
#
#
#
#
#    #Convert travele distane to m
#    X_port_G[:,0] = Distance[StartFrom+X_port_G[:,0]]-Distance[StartFrom]
#    X_stb_G[:,0] = Distance[StartFrom+X_stb_G[:,0]]-Distance[StartFrom]*1.
#
#
#
#    
#    #Stack the two sides togheter
#    #ADD skift til den andre posisjonskoordinatencd 
#    X_db = np.vstack((X_port_G,X_stb_G))
#
#    
#    
#    
#    #Do the density algorithm, It may be replaced with another if suitable
#    distanceBetweenPings = [x - np.unique(MakeDist)[i - 1] for i, x in enumerate(np.unique(MakeDist))][1:]
#    dB_ = DBSCAN(np.nanmean(distanceBetweenPings)*2, min_samples).fit(X_db)
#        
#    return X_db,dB_
#            
#     
    
    