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
import numpy as np
#from sklearn.cluster import DBSCAN
import os, random
from tools import tools
from glob import glob
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
    
    
    try: 
        import matplotlib.pyplot as plt
    except: 
        dummy = 1
        
        
    
    #Sjekk om dett er cluster TH eller Threshold
    CatThr = 3

    
    
    
    
    #Get all mat files
    os.chdir(data_directory.dir_search)
    FileNames = glob('*.mat')
    
    random.shuffle(FileNames)
    
    #Loop through each search matrix (transect)
    for SearchMatrixName in FileNames: 
        
        
        ticken = 0
            
            
        
        if not os.path.isfile(data_directory.dir_work+'/Horizontal_'+str(ticken)+'_'+SearchMatrixName) == True:
            
            
            
            #Load search matrix
            SearchMatrix = sc.loadmat(data_directory.dir_search+'/'+SearchMatrixName)
        
            
            
            #Get the data from the files
            SVres_port = SearchMatrix['SVres_port']
            SVres_stb = SearchMatrix['SVres_stb']
            R_s = SearchMatrix['R_s']
            res = SearchMatrix['res']
#            SVres_portGhost = SearchMatrix['SVres_portGhost']
#            SVres_stbGhost = SearchMatrix['SVres_stbGhost']
            DistanceTraveled = SearchMatrix['DistanceTraveled']
            ListOfFilesWithinTimeInterval = SearchMatrix['ListOfFilesWithinTimeInterval']
        
            try: 
                plt.figure(1)
                plt.clf()
                plt.subplot(2,1,1)
                plt.imshow(10*np.log10(SVres_port),aspect = 'auto')
                plt.subplot(2,1,2)
                plt.imshow(10*np.log10(SVres_stb),aspect = 'auto')
                plt.savefig(data_directory.dir_search+'/'+SearchMatrixName.replace('.mat','search.jpg'))
            except: 
                dummy = 1
    
            
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
                
            y__school_traj_stb = -y__school_traj_stb
            
            
            try: 
                plt.figure(1)
                plt.clf()
                plt.plot(x__school_traj_port,y__school_traj_port,'k.')
                plt.plot(x__school_traj_stb,y__school_traj_stb,'k.')
                plt.savefig(data_directory.dir_search+'/'+SearchMatrixName.replace('.mat','filt.jpg'))
            except: 
                dummy = 1
                
            
            
            #Start writing work files.
            #Go through each ping
            WorkPhi = np.array([])
            WorkBeam = np.array([])
            WorkFile = np.array([])
            WorkTime = np.array([])
            WorkIDX = np.array([])
            
            
            
            
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
                    Phi = np.reshape(np.dot(np.ones(Range[:,np.newaxis].shape),(variables.diry)[:,np.newaxis].T),(-1))
                    Range = np.reshape(np.dot(Range[:,np.newaxis],np.ones((variables.diry)[:,np.newaxis].T.shape)),(-1))
        
                    
            
    
                idx = np.where(abs(x__school_traj_port-DistanceTraveled[0,ii]) <=700)[0]
                for iik in range(len(idx)): 
                    r_port = np.sqrt((x__school_traj_port[idx[iik]]-(X+DistanceTraveled[0,ii]))**2+(y__school_traj_port[idx[iik]]-Y)**2)
                    WorkPhi = np.hstack((WorkPhi,Phi[np.where(r_port<=R_s[0])[0]]))
                    WorkBeam = np.hstack((WorkBeam,Range[np.where(r_port<=R_s[0])[0]]))
                    WorkTime = np.hstack((WorkTime,np.repeat(int(ListOfFilesWithinTimeInterval[ii][0]),len(Phi[np.where(r_port<=R_s[0])[0]]))))
                    WorkFile = np.hstack((WorkFile,np.repeat(ListOfFilesWithinTimeInterval[ii][1],len(Phi[np.where(r_port<=R_s[0])[0]]))))
                    WorkIDX = np.hstack((WorkIDX,np.repeat(int(ListOfFilesWithinTimeInterval[ii][2]),len(Phi[np.where(r_port<=R_s[0])[0]]))))
                    
                    
                    
#                    plt.draw()
#                    plt.pause(0.001)
                    
#                    
#                    for iip in np.where(r_port<=R_s[0])[0]: 
##                            if np.min(r_port)<=R_s[0]: 
#                        if len(WorkFileOutput) == 0: 
#                            WorkFileOutput=np.hstack((int(ListOfFilesWithinTimeInterval[ii][0]),
#                                                      ListOfFilesWithinTimeInterval[ii][1],
#                                                      int(ListOfFilesWithinTimeInterval[ii][2]),
#                                                      Phi[iip],Range[iip]))
#                            
#                            
#                        else: 
#                            WorkFileOutput = np.vstack((WorkFileOutput
#                                                        ,np.hstack((int(ListOfFilesWithinTimeInterval[ii][0]),
#                                                                    ListOfFilesWithinTimeInterval[ii][1],
#                                                                    int(ListOfFilesWithinTimeInterval[ii][2])
#                                                                    ,Phi[iip],Range[iip]))))
                                    
                                 
                idx = np.where(abs(x__school_traj_stb-DistanceTraveled[0,ii]) <=700)[0]
                for iik in range(len(idx)): 
                    
                    r_stb = np.sqrt((x__school_traj_stb[idx[iik]]-(X+DistanceTraveled[0,ii]))**2+(y__school_traj_stb[idx[iik]]-Y)**2)
                    WorkPhi = np.hstack((WorkPhi,Phi[np.where(r_stb<=R_s[0])[0]]))
                    WorkBeam = np.hstack((WorkBeam,Range[np.where(r_stb<=R_s[0])[0]]))
                    WorkPhi = np.hstack((WorkPhi,Phi[np.where(r_stb<=R_s[0])[0]]))
                    WorkBeam = np.hstack((WorkBeam,Range[np.where(r_stb<=R_s[0])[0]]))
                    WorkTime = np.hstack((WorkTime,np.repeat(int(ListOfFilesWithinTimeInterval[ii][0]),len(Phi[np.where(r_stb<=R_s[0])[0]]))))
                    WorkFile = np.hstack((WorkFile,np.repeat(ListOfFilesWithinTimeInterval[ii][1],len(Phi[np.where(r_stb<=R_s[0])[0]]))))
                    WorkIDX = np.hstack((WorkIDX,np.repeat(int(ListOfFilesWithinTimeInterval[ii][2]),len(Phi[np.where(r_stb<=R_s[0])[0]]))))
                    

#                print(WorkFileOutput.shape[0])
                #Make a new file so it does not become too large
                if len(WorkPhi)>2000000: 
                    sc.savemat(data_directory.dir_work+'/Horizontal_'+str(ticken)+'_'+SearchMatrixName,mdict={'WorkPhi':WorkPhi,'WorkBeam':WorkBeam,'WorkFile':WorkFile,'WorkIDX':WorkIDX,'WorkTime':WorkTime})
                    ticken = ticken+1
                    
                    WorkPhi = np.array([])
                    WorkBeam = np.array([])
                    WorkFile = np.array([])
                    WorkTime = np.array([])
                    WorkIDX = np.array([])
            sc.savemat(data_directory.dir_work+'/Horizontal_'+str(ticken)+'_'+SearchMatrixName,mdict={'WorkPhi':WorkPhi,'WorkBeam':WorkBeam,'WorkFile':WorkFile,'WorkIDX':WorkIDX,'WorkTime':WorkTime})
            
            
                    
                    
#                #Skip those pings that don't have relevant data
#                if np.nansum([(x__school_traj_port<=np.nanmax(X)+DistanceTraveled[0,ii])&(x__school_traj_port>=np.nanmin(X)+DistanceTraveled[0,ii])]) >0: 
#                    
#                    
#                    
#                    
#                    #go through each pixel and find if it is relevant data
#                    for iik in range(len(X)): 
#                        
#                        print(len(X)-iik,end = '\r')
#                        
#                        
#                        #Find the distance between pixels and the search matrix output
#                        r_stb = np.sqrt((x__school_traj_stb-(X[iik]+DistanceTraveled[0,ii]))**2+(y__school_traj_stb-Y[iik])**2)
#                        r_port = np.sqrt((x__school_traj_port-(X[iik]+DistanceTraveled[0,ii]))**2+(y__school_traj_port-Y[iik])**2)
#                        
#                        
#                        
#                        if np.min(r_stb)<=R_s[0]: 
#                            if len(WorkFileOutput) == 0: 
#                                WorkFileOutput=np.hstack((int(ListOfFilesWithinTimeInterval[ii][0]),
#                                                          ListOfFilesWithinTimeInterval[ii][1],
#                                                          int(ListOfFilesWithinTimeInterval[ii][2]),
#                                                          Phi[iik],Range[iik]))
#                                
#                                
#                            else: 
#                                WorkFileOutput = np.vstack((WorkFileOutput
#                                                            ,np.hstack((int(ListOfFilesWithinTimeInterval[ii][0]),
#                                                                        ListOfFilesWithinTimeInterval[ii][1],
#                                                                        int(ListOfFilesWithinTimeInterval[ii][2])
#                                                                        ,Phi[iik],Range[iik]))))
#                            
#                            
#                        if np.min(r_port)<=R_s[0]: 
#                            if len(WorkFileOutput) == 0: 
#                                WorkFileOutput=np.hstack((int(ListOfFilesWithinTimeInterval[ii][0]),
#                                                          ListOfFilesWithinTimeInterval[ii][1],
#                                                          int(ListOfFilesWithinTimeInterval[ii][2]),
#                                                          Phi[iik],Range[iik]))
#                                
#                                
#                            else: 
#                                WorkFileOutput = np.vstack((WorkFileOutput
#                                                            ,np.hstack((int(ListOfFilesWithinTimeInterval[ii][0]),
#                                                                        ListOfFilesWithinTimeInterval[ii][1],
#                                                                        int(ListOfFilesWithinTimeInterval[ii][2])
#                                                                        ,Phi[iik],Range[iik]))))
#                                
#                                
#                    
                                
            
            
            
            
            
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
                dummy=1

        

                
    