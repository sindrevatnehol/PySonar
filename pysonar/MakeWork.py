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
from sklearn.cluster import DBSCAN
import os
from glob import glob
#from MakeFigures import plotFigure
#from PySonTools import MakeIntegrationLine,GetSpikesFromEchogram,FindLocationOfSchool
from MakeIndex import MakeIndex





def ComputePearsonR(SA_profos,SA_profos2,SA_AutoWithoutGhost):
    #Compute the Pearson r between the two models

    sap = np.zeros(len(SA_profos))
    sap2 = np.zeros(len(SA_profos))
    saa = np.zeros(len(SA_profos))
    for ik in range(len(SA_profos)): 
        if ik == 0: 
            sap[ik] = SA_profos[ik]
            sap2[ik] = SA_profos2[ik]
            saa[ik] = SA_AutoWithoutGhost[ik]
        else:   
            sap[ik] = SA_profos[ik]-SA_profos[(ik-1)]
            sap2[ik] = SA_profos2[ik]-SA_profos2[(ik-1)]
            saa[ik] = SA_AutoWithoutGhost[ik]-SA_AutoWithoutGhost[(ik-1)]
            
    PR_Before = scpy.stats.pearsonr(sap,saa)
    PR_After = scpy.stats.pearsonr(sap2,saa)
    
    return PR_Before, PR_After
    
    
    
    
    
        

def WaveRemover(dim,pitch,roll,limit):
    #This function identify when there is a chanse that, due to vessel movement, 
    #there is a chanse of surface reveberation.
    #It will then remove these pings
    #ADD could be included in the MakeSearchMatrix.py  !!!
    
    WaveRemover = np.ones(dim)
    WaveRemover[:,:,np.where(abs(pitch-np.mean(pitch))>limit)] = np.nan
    WaveRemover[:,:,np.where((roll-np.mean(roll))>limit)] = np.nan

    return WaveRemover


    
    
    

def ProtocolForCluster(X_port,X_port_ghost,X_stb,X_stb_ghost,Bananatool,
                      MakeDist,Cat_port,min_samples):
    '''
    This is the protocol for making the cluster
    
    In future this will change so it don't merge small schools closly togetter in range direction
    '''
    
    
    #Find only data without ghost
    #Bruk heller x__xchool_trajS
    X_port_G = np.zeros((1,2))
    for i in np.arange(len(X_port[:,0])-1,0,-1):
        for ii in range(len(X_port_ghost[:,0])-1,0,-1): 
            if (X_port_ghost[ii,0]==X_port[i,0]): 
                if (X_port_ghost[ii,1]==X_port[i,1]): 
                    X_port_G = np.vstack((X_port_G,X_port_ghost[ii,:]))
    X_port_G = X_port_G[1:,:]


    X_stb_G = np.zeros((1,2))
    for i in np.arange(len(X_stb[:,0])-1,0,-1):
        for ii in range(len(X_stb_ghost[:,0])-1,0,-1): 
            if (X_stb_ghost[ii,0]==X_stb[i,0]): 
                if (X_stb_ghost[ii,1]==X_stb[i,1]): 
                    X_stb_G = np.vstack((X_stb_G,X_stb_ghost[ii,:]))
    X_stb_G = X_stb_G[1:,:]



    '''Make cluster out of the data'''
    
    
    X_port_G = X_port_G.astype(int)
    X_stb_G = X_stb_G.astype(int)
    
    
    if len(Bananatool)==1:
        BananaTool = Bananatool[0]
    else: 
        BananaTool = Bananatool
        
        
        
        
    
    #Convert range data to m
    #The constant 0.18 is the sample distance
    X_port_G[:,1] = abs(X_port_G[:,1])*int(BananaTool[3])*0.18
    X_stb_G[:,1] = X_stb_G[:,1]*int(BananaTool[3])*0.18
    



    #Help to find the correct travel distance of the pixel
    Distance = np.unique(MakeDist)
    StartFrom = int((len(Distance)-Cat_port.shape[1])/2)+1




    #Convert travele distane to m
    X_port_G[:,0] = Distance[StartFrom+X_port_G[:,0]]-Distance[StartFrom]
    X_stb_G[:,0] = Distance[StartFrom+X_stb_G[:,0]]-Distance[StartFrom]*1.



    
    #Stack the two sides togheter
    #ADD skift til den andre posisjonskoordinatencd 
    X_db = np.vstack((X_port_G,X_stb_G))

    
    
    
    #Do the density algorithm, It may be replaced with another if suitable
    distanceBetweenPings = [x - np.unique(MakeDist)[i - 1] for i, x in enumerate(np.unique(MakeDist))][1:]
    dB_ = DBSCAN(np.nanmean(distanceBetweenPings)*2, min_samples).fit(X_db)
        
    return X_db,dB_
            
     
    
    
    
def GetSpikesFromEchogram(SVres,threshold,S): 
    '''Purpose: 
        Threshold the data, and get index of the leftover
        ''' 
    
    #Make the cathegorization matrix
    CathegorizationMatrix = np.zeros(SVres.shape)
    

#    median = np.nanmedian(10**(SVres/10),axis=1)
    median = np.nanmedian(SVres,axis=1)
    
#    median_filter = 10*np.log10(np.repeat(median[:,np.newaxis],SVres.shape[1],axis=1))
    median_filter = (np.repeat(median[:,np.newaxis],SVres.shape[1],axis=1))
    
    #Find all pixels likly to have fish inside it (larger then background noise)
#    CathegorizationMatrix[np.where(SVres>(10*np.log10(PP)+threshold))] = 1
    CathegorizationMatrix[np.where(SVres>(median_filter*threshold))] = 1
    
    
    x_i,y_i = np.where(CathegorizationMatrix>0)

    Y_vec = -x_i[:,np.newaxis]
    X_vec = y_i[:,np.newaxis]
                                        
    #Return index file                                  
    return CathegorizationMatrix,np.hstack((X_vec,Y_vec)),
    


def GetSpikesFromEchogramGhost(SVres,threshold): 
    '''Purpose: 
        Threshold the data, and get index of the leftover
        ''' 
    
    #Make the cathegorization matrix
    CathegorizationMatrix = np.zeros(SVres.shape)
    
    
    #Find all pixels likly to have fish inside it (larger then background noise)
    CathegorizationMatrix[np.where(SVres>=threshold)] = 1
    
    x_i,y_i = np.where(CathegorizationMatrix>0)

    Y_vec = -x_i[:,np.newaxis]
    X_vec = y_i[:,np.newaxis]
    
      
    #Return index file                                  
    return CathegorizationMatrix,np.hstack((X_vec,Y_vec)), 
    


    
    
#The main function    
#MakeWork(True, directory2Data,'','','',1,3)
def MakeWork(makeNewWork, data_directory,SearchMatrixName,
             prefix,compensate,clusterTH,Threshold):
    
    
    #Sjekk om dett er cluster TH eller Threshold
    CatThr = 1

    
    #Loop through each search matrix (transect)
    for SearchMatrixName in os.listdir(data_directory.dir_search): 
        
        
        
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
    

#        print(ListOfFilesWithinTimeInterval)
        
        # Finding the location of the schools in the search matrix
        # ADD may be replaced with an other function
        Cat_port, X_port = GetSpikesFromEchogram(SVres_port,CatThr,
                                    np.hstack((np.reshape(SVres_port,(-1,)),
                                    np.reshape(SVres_stb,(-1,)))))     
        Cat_stb, X_stb = GetSpikesFromEchogram(SVres_stb,CatThr,
                                    np.hstack((np.reshape(SVres_port,(-1,)),
                                    np.reshape(SVres_stb,(-1,)))))
    
        
        
        try: 
            import matplotlib.pyplot as plt
    
            plt.figure(1)
            plt.clf()
            plt.imshow(SVres_port,aspect = 'auto')
            plt.colorbar()
            plt.draw()
            plt.savefig(data_directory.dir_search+'/'+SearchMatrixName.replace('.mat','2.jpg'))
        except: 
            k=1
            
        print(SVres_port)
        print(X_stb)
        
#        x__school_traj_port,y__school_traj_port = FindLocationOfSchool.FindLocationOfSchool(Cat_port,
#                                                                                            Bananatool,
#                                                                                            NumberOfPingsInBatch/2,
#                                                                                            RangeOut, 
#                                                                                            DistanceMatrix)
#        x__school_traj_stb,y__school_traj_stb = FindLocationOfSchool.FindLocationOfSchool(Cat_stb,
#                                                                                          Bananatool,
#                                                                                          NumberOfPingsInBatch/2,
#                                                                                          RangeOut,
#                                                                                          DistanceMatrix)
        
        
        
                                                                                
    #Get the current directory
#    current_directory = os.getcwd()
    
    
    
    #Change the directory to the location of the data
#    os.chdir(data_directory)
    
    
    
    #Get all of the subfolders within folder
#    ListOfCases = glob('*/')
    
    
    
    #Some bokkeeping for results
#    Compar_auto = []
#    Compar_manual = []
#    Compar_manualcorr = []



    #Go through each cruice
    for subdir in ListOfCases: 
        #Print msg to the user
        print('Preparing cruice:   ' +subdir)
        print('Start loading files')
        
        
        
        #Set the path to the search matrix
        #ADD this may be changed later  !!!
        dir2file = data_directory+'/'+subdir+'/ACOUSTIC_DATA/SU90/survey/netcdf/'
        
        
        
        #Change directory to the location of the search matrix
        os.chdir(dir2file)
        
        
        
        #Load GhostSchool data 
        caseGhost = sc.loadmat(dir2file+'figure90DegRot.mat')

        
        
        #Load school data 
        case = sc.loadmat(dir2file+SearchMatrixName)
        
        
        
        
        #Load PROFOS work data 
#        caseProfos = sc.loadmat(dir2file+'Profosnormal.mat')
        
        
        
        
        #Get the proper variables
        SVres_port = case['SVres_port']
        SVres_stb = case['SVres_stb']
        CatThr = Threshold #case['CatThreshold']        
        Bananatool = case['BananaTool']
        NumberOfPingsInBatch = case['NumberOfPingsInBatch'][0][0]
        RangeOut = case['RangeOut']
        DistanceMatrix = case['DistanceMatrix']
        GhostSV_port = caseGhost['SVres_port']
        GhostSV_stb = caseGhost['SVres_stb']
        MakeRange = case['MakeRange']
        MakeBeam = case['MakeBeam']
        MakeDist = case['MakeDist']
        MakeSV = case['MakeSV']
        MakeProfosWork = caseProfos['MakeProfosWork']
        MakeProfosWork2 = caseProfos['MakeProfosWork2']
        pitch = caseProfos['pitch'][0]
        roll = caseProfos['roll'][0]
        RemoveToCloseValues = case['RemoveToCloseValues']
        Work3DMatrix = case['Work3DMatrix']
        RangeMatrix = case['RangeMatrix']
        DistanceTraveled = case['DistanceTraveled']
        BeamDirectionMatrix = case['BeamDirectionMatrix']
#        TiltVec = case['Tilt']
#        sV_port = case['sV_port']




        #Make some structural changes of the matrix, so it is more compatable
        RangeOut = RangeOut[:,:MakeDist.shape[0]]
        RangeMatrix =np.repeat(np.repeat(MakeRange[:,0,0][:,np.newaxis],len(MakeDist[0,:,0]),axis=1)[:,:,np.newaxis],len(MakeDist[0,0,:]),axis=2)
        DistanceMatrix = MakeDist
        
        
        
        
        #Small bug correction
        #ADD This can be fixed earlier in the script !!!
        min_stuff = np.min((RangeMatrix.shape[2],BeamDirectionMatrix.shape[2]))
        RangeMatrix  =RangeMatrix[:,:,:min_stuff]
        BeamDirectionMatrix  =BeamDirectionMatrix[:,:,:min_stuff]
        MakeDist = MakeDist[:,:,:min_stuff]
        Work3DMatrix = Work3DMatrix[:,:,:min_stuff]
        


        print('    Files are loaded')
        
        
        
        # Finding the location of the schools in the search matrix
        # ADD may be replaced with an other function
        Cat_port, X_port = GetSpikesFromEchogram.GetSpikesFromEchogram(SVres_port,
                                                                       CatThr,
                                                                       np.hstack((np.reshape(SVres_port,(-1,)),
                                                                                  np.reshape(SVres_stb,(-1,)))))     
        Cat_stb, X_stb = GetSpikesFromEchogram.GetSpikesFromEchogram(SVres_stb,
                                                                     CatThr,
                                                                     np.hstack((np.reshape(SVres_port,(-1,)),
                                                                                np.reshape(SVres_stb,(-1,)))))
    
        
    
        #AutoSchoolDetection, find the loaction in the geo-reference coordinate
        #system 
        x__school_traj_port,y__school_traj_port = FindLocationOfSchool.FindLocationOfSchool(Cat_port,
                                                                                            Bananatool,
                                                                                            NumberOfPingsInBatch/2,
                                                                                            RangeOut, 
                                                                                            DistanceMatrix)
        x__school_traj_stb,y__school_traj_stb = FindLocationOfSchool.FindLocationOfSchool(Cat_stb,
                                                                                          Bananatool,
                                                                                          NumberOfPingsInBatch/2,
                                                                                          RangeOut,
                                                                                          DistanceMatrix)
        
        
        
        #Make matrix indexing the chost schools
        CaseGhost_port = np.zeros((SVres_stb.shape))
        CaseGhost_port[np.where(-np.maximum(GhostSV_port-SVres_port,GhostSV_stb-SVres_port)>0)] = 1
        CaseGhost_stb = np.zeros((SVres_stb.shape))
        CaseGhost_stb[np.where(-np.maximum(GhostSV_port-SVres_stb,GhostSV_stb-SVres_stb)>0)] = 1
                                          
                                         

                      
        # Finding the location of the ghost school in search matrix
        Cat_port_ghost, X_port_ghost = GetSpikesFromEchogram.GetSpikesFromEchogramGhost(CaseGhost_port,clusterTH)
        Cat_stb_ghost, X_stb_ghost = GetSpikesFromEchogram.GetSpikesFromEchogramGhost(CaseGhost_stb,clusterTH)
    
        
        
        # ADD this is no longer part of the manuscript and can therefore be removed !!!
#        #Find the location of ghost schools
#        Cat_port_ghost_max, X_port_ghost_max = GetSpikesFromEchogram.GetSpikesFromEchogramGhost(-np.maximum(GhostSV_port-SVres_port,
#                               GhostSV_stb-SVres_port),0)
#        Cat_stb_ghost_max, X_stb_ghost_max = GetSpikesFromEchogram.GetSpikesFromEchogramGhost(-np.maximum(GhostSV_port-SVres_stb,
#                              GhostSV_stb-SVres_stb),0)
        

        
        #Ghost school detection, find geografic location of relative vessel start location
        # ADD  NumberOfPingsInBatch/4  is redundent !!!
        x__school_traj_port_Ghost,y__school_traj_port_Ghost = FindLocationOfSchool.FindLocationOfSchool(Cat_port_ghost,Bananatool,NumberOfPingsInBatch/2,RangeOut,DistanceMatrix)
        x__school_traj_stb_Ghost,y__school_traj_stb_Ghost = FindLocationOfSchool.FindLocationOfSchool(Cat_stb_ghost,Bananatool,NumberOfPingsInBatch/2,RangeOut,DistanceMatrix)

        

        #Location of all data in the geo-reference coordinate system
        #ADD include vessel heading and speed change as well  !!!
        X_traj = MakeRange*np.cos(MakeBeam*np.pi/180+compensate)+MakeDist
        Y_traj = MakeRange*np.sin(MakeBeam*np.pi/180+compensate)     
        
        
        
        
        print('    Location of pixels are computed')
        
        
        
        
        #Remove those data that are likly ghost
        x__new__port = []
        y__new__port = []
        for i in np.arange(len(x__school_traj_port[:])-1,0,-1):
            for ii in range(len(x__school_traj_port_Ghost[:])-1,0,-1): 
                if (x__school_traj_port_Ghost[ii]==x__school_traj_port[i]): 
                    if (y__school_traj_port_Ghost[ii]==y__school_traj_port[i]): 
                        x__new__port = np.hstack((x__new__port,x__school_traj_port_Ghost[ii]))
                        y__new__port = np.hstack((y__new__port,y__school_traj_port_Ghost[ii]))
        
        x__new__stb = []
        y__new__stb = []
        for i in np.arange(len(x__school_traj_stb[:])-1,0,-1):
            for ii in range(len(x__school_traj_stb_Ghost[:])-1,0,-1): 
                if (x__school_traj_stb_Ghost[ii]==x__school_traj_stb[i]): 
                    if (y__school_traj_stb_Ghost[ii]==y__school_traj_stb[i]): 
                        x__new__stb = np.hstack((x__new__stb,x__school_traj_stb_Ghost[ii]))
                        y__new__stb = np.hstack((y__new__stb,y__school_traj_stb_Ghost[ii]))
                        
    
        
        #Make work matrixes
        Work = np.zeros(X_traj.shape)
        WorkGhost = np.zeros(X_traj.shape)
        
        
        
        #Make the work output from the automatic stuff. 
        if makeNewWork == True: 
            print('Making work index: ')
            MakeIndex.MakeIndex(case,Work,x__school_traj_port,y__school_traj_port,
                                X_traj,Y_traj,x__school_traj_stb,-y__school_traj_stb,
                                dir2file+prefix+'AutoWork.mat')
            
        
            print('Making Ghost work index: ')
            MakeIndex.MakeIndex(case,WorkGhost,x__new__port,y__new__port,X_traj,
                                Y_traj,x__new__stb,-y__new__stb,
                                dir2file+prefix+'GhostSchoolWork.mat')
        else: 
            
            
            
            #If there is no work file in folder, make one
            if not os.path.isfile(dir2file+prefix+'AutoWork.mat'): 
                print('Making work index: ')
                MakeIndex.MakeIndex(case,Work,x__school_traj_port,
                                    y__school_traj_port,X_traj,Y_traj,
                                    x__school_traj_stb,-y__school_traj_stb,
                                    dir2file+prefix+'AutoWork.mat')
                
                
            if not os.path.isfile(dir2file+prefix+'GhostSchoolWork.mat'): 
                print('Making Ghost work index: ')
                MakeIndex.MakeIndex(case,WorkGhost,x__new__port,y__new__port,
                                    X_traj,Y_traj,x__new__stb,-y__new__stb,
                                    dir2file+prefix+'GhostSchoolWork.mat')


#                
#        #Load work files
#        Work = sc.loadmat(dir2file+prefix+'AutoWork.mat')['Work']
#        WorkGhost = sc.loadmat(dir2file+prefix+'GhostSchoolWork.mat')['Work']
#
#        
#        #Remove pixels to close to these vessel   
#        Work[np.where(MakeRange<RemoveToCloseValues)] = 0
#             
#
#        #Find the index of where there is registered fish in the work files. 
#        workIDX = np.where(Work>0)
#        
#        
#        
#
#
#        
#        #Small fix for case 13, check if it can be removed        
##        MakeSV[np.where(MakeRange<=60)] = np.nan
#        
#        
#        
#        
#        #Make the integrator line for both profos and the automatic algorithm. 
#        if makeNewWork == True: 
#            print('Integrating with ghost')
#            MakeIntegrationLine.MakeIntegrationLine(MakeProfosWork,
#                                                    MakeDist,
#                                                    MakeSV,workIDX,
#                                                    dir2file+prefix+'Integrator.mat')
#            print('Integrating with ghost v2')
#            MakeIntegrationLine.MakeIntegrationLine(MakeProfosWork2,
#                                                    MakeDist,MakeSV,
#                                                    workIDX,
#                                                    dir2file+prefix+'Integrator2.mat')
#        else: 
#            if not os.path.isfile(dir2file+prefix+'Integrator.mat'): 
#                print('Integrating with ghost')
#                MakeIntegrationLine.MakeIntegrationLine(MakeProfosWork,MakeDist,
#                                                        MakeSV,workIDX,
#                                                        dir2file+prefix+'Integrator.mat')
#            elif not os.path.isfile(dir2file+prefix+'Integrator2.mat'):
#                print('Integrating with ghost v2')
#                MakeIntegrationLine.MakeIntegrationLine(MakeProfosWork2,
#                                                        MakeDist,MakeSV,
#                                                        workIDX,
#                                                        dir2file+prefix+'Integrator2.mat')
#
#                
#
#        #remove data from ghost schools
#        Work = Work*WorkGhost
#                                     
#             
#        #Find the index of where there is registered fish in the work files. 
#        workIDX = np.where(Work>0)
#         
#        
#        
#        
#        
#        #Remove pings where there is a larger possibillity of surface reveberation
#        #ADD make the Tilt and BW_vertical avaliable in the loaded files  !!!
#        #WaveRem = WaveRemover(MakeSV.shape,pitch,roll,(Tilt-BW_vert/2))
##        WaveRem = WaveRemover(MakeSV.shape,pitch,roll,TiltVec-4)
#        WaveRem = WaveRemover(MakeSV.shape,pitch,roll,3)
#        
#        
#        
#        #Small fix for case 13, check if it can be removed        
#        MakeSV[np.where(MakeRange<=80)] = np.nan
#
#
#
#        if makeNewWork == True: 
#            print('Integrating without ghost')
#            MakeIntegrationLine.MakeIntegrationLine(MakeProfosWork,MakeDist,MakeSV*WaveRem,workIDX,dir2file+prefix+'IntegratorWithoutGhost.mat')
#        else: 
#            if not os.path.isfile(dir2file+prefix+'IntegratorWithoutGhost.mat'): 
#                print('Integrating without ghost')
#                MakeIntegrationLine.MakeIntegrationLine(MakeProfosWork,MakeDist,MakeSV*WaveRem,workIDX,dir2file+prefix+'IntegratorWithoutGhost.mat')
#
#                
#    
#
#
#        #Load the data from the integration 
#        Integrator = sc.loadmat(dir2file+prefix+'Integrator.mat')
#        Integrator2 = sc.loadmat(dir2file+prefix+'Integrator2.mat')
#        IntegratorWithoutGhost = sc.loadmat(dir2file+prefix+'IntegratorWithoutGhost.mat')
#        
#        
#        #Get values from the loaded files
#        SA_Auto = Integrator['SA_Auto'][0][:-1]
#        SA_profos = Integrator['SA_profos'][0][:-1]
#        SA_profos2 = Integrator2['SA_profos'][0][:-1]
#        SA_dist = Integrator['SA_dist'][0][:-1]
#        SA_AutoWithoutGhost = IntegratorWithoutGhost['SA_Auto'][0][:-1]
#
#
#
#
#
#        #Stack the data for further analysiz
#        Compar_auto = np.hstack((Compar_auto,np.nansum(SA_AutoWithoutGhost)/MakeProfosWork.shape[2]/(np.max(MakeDist)/1850)))
#        Compar_manual = np.hstack((Compar_manual,np.nansum(SA_profos)/MakeProfosWork.shape[2]/(np.max(MakeDist)/1850)))
#        Compar_manualcorr = np.hstack((Compar_manualcorr,np.nansum(SA_profos2)/MakeProfosWork.shape[2]/(np.max(MakeDist)/1850)))
#        
#        
#        
#        
#        #Get the location of the schools from PROFOS
#        X_profos = (MakeDist+ Work3DMatrix* 
#                    RangeMatrix*np.cos(BeamDirectionMatrix*np.pi/180+compensate))
#        
#        
#        
#        #Remove those data whithout relevant data
#        Work3DMatrix[np.where(X_profos<0)] = 0
#        Work3DMatrix[np.where(Work3DMatrix ==0)] = np.nan
#
#
#
#        #The cartesion position of the data
#        X_profos =np.reshape(X_profos,(-1) )
#        X_auto = np.reshape(DistanceMatrix+ Work* RangeMatrix*np.cos(BeamDirectionMatrix*np.pi/180+compensate),(-1) )        
#                                  
#        
#        #Remove those data whithout relevant data
#        Work[np.where(X_auto<0)] = 0
#        Work[np.where(Work ==0)] = np.nan
#        
#        
#        
#        #Make result directory if it does not exist
#        if os.path.exists(current_directory+'/Result')==False: 
#            os.mkdir(current_directory+'/Result')
#        if os.path.exists(current_directory+'/Result/'+subdir)==False: 
#            os.mkdir(current_directory+'/Result/'+subdir)
#        
#            
#        
#        
#        #Start making figures and movies
#        
##        
##        #Make the 3D figure of all the PROFOS data
##        plotFigure.plotFigure1(MakeProfosWork,DistanceMatrix,
##                    BeamDirectionMatrix+compensate*180/np.pi,RangeMatrix,
##                    MakeSV,current_directory+'/Result/'+subdir+prefix+'Figure1.jpg') 
###        
##        
##        #Make the movie of the 3D profos data. 
##        #Does not work when runned in terminal
###        plotFigure.plotMovie1(MakeProfosWork,DistanceMatrix,
###                    BeamDirectionMatrix+compensate*180/np.pi,RangeMatrix,
###                    MakeSV,current_directory+'/Result/'+subdir+prefix+'Movie1.mp4') 
##        
##        
##        
##        #Make the figure of the search matrix
##        plotFigure.plotFigure2(SVres_port,SVres_stb,DistanceMatrix,
##                    RangeMatrix,current_directory+'/Result/'+subdir+prefix+'Figure3.jpg',subdir[:-1])
##        
##        
###        #Make the figure of the noise level
###        plotFigure.plotFigure3(sV_port,SVres_port,
###                    current_directory+'/Result/'+subdir+prefix+'Figure4.jpg')
##    
##        
##        #Make the figure of the scrutinize searrch matrix
##        plotFigure.plotFigure4(Cat_port,Cat_stb,DistanceMatrix,RangeMatrix,
##                    DistanceTraveled,NumberOfPingsInBatch,
##                    current_directory+'/Result/'+subdir+prefix+'Figure4A.jpg','A')
##        
##        
#        
#        
#        
##        #Make the figure of the ghost schools
##        plotFigure.plotFigure4(Cat_port_ghost_max,Cat_stb_ghost_max,DistanceMatrix,RangeMatrix,
##                    DistanceTraveled,NumberOfPingsInBatch,current_directory+'/Result/'+subdir+prefix+'Figure4B.jpg','B')
##
##        
##        
##        
##        
##        
##        #Make the figure of schools without ghost
##        plotFigure.plotFigure4(Cat_port_ghost_max*Cat_port,Cat_stb_ghost_max*Cat_stb,DistanceMatrix,RangeMatrix,
##                    DistanceTraveled,NumberOfPingsInBatch,current_directory+'/Result/'+subdir+prefix+'Figure4C.jpg','C')
##        
#        
#        
#        
#        
##        
##        #Compute the cluster information
##        X_db,dB_ = ProtocolForCluster(X_port,X_port_ghost,X_stb,X_stb_ghost,Bananatool,MakeDist,Cat_port,min_samples)
##        
##        
##        #Plot the cluster stuff
##        plotFigure.plotFigure6(dB_,X_db,current_directory+'/Result/'+subdir+prefix+'Figure5.jpg',subdir[:-1])
##
#
#        
#        
#        
#        
#        
#        #Make the figure of the integration lines
#        plotFigure.plotFigure5(SA_dist,SA_profos,SA_profos2,SA_Auto,SA_AutoWithoutGhost,subdir,
#                               current_directory+'/Result/'+subdir+prefix+'Figure6.jpg')
#
#
#        
#        
#        
#        #Print the similarity to the user
#        #ADD Print it to an text file
#        PR_Before, PR_After = ComputePearsonR(SA_profos,SA_profos2,SA_AutoWithoutGhost)        
#
#        print('Similarity: ')
#        print(PR_Before)
#        print(PR_After)
#        
#        print(len(SA_profos))
#        PR_Before, PR_After = ComputePearsonR(SA_profos[:-5],SA_profos2[:-5],SA_AutoWithoutGhost[:-5])        
#
#        print('Similarity: ')
#        print(PR_Before)
#        print(PR_After)
#       
#       
#
#
#        #Set the size of the figure
#        figureWidth = 70*0.03937*2
#        figureHeight = 70*0.03937*1.7
#    
#    
##        fig = plt.figure(8,figsize=(figureWidth,figureHeight))
##        plt.clf()
##        
##        ax = plt.gca()
##        for i in range(len(Compar_auto)): 
##            ax.text(i+1, Compar_auto[i], 'A', size=18, color = 'blue',ha='center', va='center')
##            ax.text(i+1, Compar_manual[i], 'M', size=18, color = 'red',ha='center', va='center')
##            
##            
##        plt.xlabel('Case number',fontsize=14)
##        h=plt.ylabel('$\sum s_v  / nmi $',fontsize=14)
##        
###        h.set_rotation(0)
##
##        plt.xlim((0,len(Compar_auto)+1))
##        plt.ylim((-np.min(Compar_auto),np.max(Compar_auto)*1.1))
##        
##        plt.gcf().subplots_adjust(left = 0.21)
##    
##        plt.savefig(current_directory+'/Result/'+'Figure10.jpg',dpi=600)
##        plt.draw()
#        
#        
#        
#        fig = plt.figure(9,figsize=(figureWidth,figureHeight))
#        plt.clf()
#        
#        ax = plt.gca()
#        
#        plt.plot(np.arange(len(Compar_auto))+1+0.2,Compar_auto,'kx',label = 'Auto')
#        plt.plot(np.arange(len(Compar_manual))+1,Compar_manual,'ko',label = 'Semi-auto')
#        plt.plot(np.arange(len(Compar_manualcorr))+1-0.2,Compar_manualcorr,'k*', label = 'Semi-auto, corr.')
##        for i in range(len(Compar_auto)): 
##            ax.text(i+1, Compar_auto[i], 'A', size=18, color = 'blue',ha='center', va='center')
##            ax.text(i+1, Compar_manual[i], 'M', size=18, color = 'red',ha='center', va='center')
##            ax.text(i+1, Compar_manualcorr[i], 'Mc', size=18, color = 'green',ha='center', va='center')
#            
#            
#        plt.xlabel('Case number',fontsize=14)
#        h=plt.ylabel('$\sum s_v  / nmi $',fontsize=14)
#        
#        h.set_rotation(90)
#
#        plt.xlim((0,len(Compar_auto)+1))
#        plt.ylim((-np.min(Compar_auto),np.max(Compar_auto)*1.1))
#
#
#
#    
#        #Add figure legends
#        plt.legend(frameon=False,fontsize=12,loc = 2,numpoints=1)
#        
#        
#        plt.gcf().subplots_adjust(left = 0.21)
#    
#        plt.savefig(current_directory+'/Result/'+'Figure7.jpg',dpi=1000)
#        plt.draw()
#        
#        
#        sc.savemat(current_directory+'/Result/'+'MatforOverallfigure.mat',mdict={'Compar_auto':Compar_auto,'Compar_manual':Compar_manual,'Compar_manualcorr':Compar_manualcorr})
#        
        
        
#if __name__ == '__main__':   
#    
#    makeNewWork = True
#    
#    #for clustering
##    eps = 45                #Not Used anymore
##    
##    min_samples = 1         #Not Used anymore
#    
#    data_directory = ''
#    
#    data_directory = (os.getcwd()+'\Data')
#    
#    SearchMatrixName = 'figurenormal.mat'
#    
#    prefix = ''
#    
#    compensate = 0          #Not Used anymore
#    
#    clusterTH = 1           #This is for finding the ghost school
#    
#    MakeWork(makeNewWorkdata_directory,SearchMatrixName,prefix,compensate,clusterTH)
#    
#    
#    
#    
        