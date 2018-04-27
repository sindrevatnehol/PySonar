# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 13:47:24 2017

@author: sindrev
"""




import xml.dom.minidom as minidom
import os
##from shutil import copyfile

#import MakeSearchMatrix, makePingWork, MakeSearchMatrixVert 
import numpy as np
#from SendMail import send_email
import platform
from tools import tools
from tools import SendMail
import scipy.io as sc
from MakeSearch import MakeSearch
#from MakeVerticalIndex import MakeVerticalIndex
#from MakeWork import MakeWork
from netCDF4 import Dataset



def GetShortListOfFiles(CompleteListOfFiles,startTime,endTime): 
    ShortListOfFiles = []
    for i in range(len(CompleteListOfFiles[:,0])): 
        start = SecondsBetweenTransect(startTime,float(CompleteListOfFiles[i,0]))
        end = SecondsBetweenTransect(endTime,float(CompleteListOfFiles[i,0]))
            
        if start<0: 
            if end >0: 
                if ShortListOfFiles == []: 
                    ShortListOfFiles = CompleteListOfFiles[i,:]
                else: 
                    ShortListOfFiles = np.vstack((ShortListOfFiles,CompleteListOfFiles[i,:]))
                                    
    return ShortListOfFiles
                            


        
def ComListOfFiles(directory2Data,beam_mode): 
    #This function makes an idx of all files
    
    
    #Open the IDX files
    try: 
        fd = open(directory2Data.dir_src+'/'+'pingtime.csv','r')    
        fd2 = open(directory2Data.dir_src+'/'+'FileList.csv','r')
        fd3 = open(directory2Data.dir_src+'/'+'IDX.csv','r')
        fd4 = open(directory2Data.dir_src+'/'+'BeamIDX.csv','r')
    
        
        A = np.asarray(fd.read().split(','))[:-1].astype(np.float)
        B = np.asarray(fd2.read().split(','))[:-1]
        C = np.asarray(fd3.read().split(','))[:-1].astype(np.int)
        D = np.asarray(fd4.read().split(','))[:-1]
        CompleteListOfFiles = np.array((A.T,B.T,C.T)).T

        

    # If the files don't exist, create them
    except FileNotFoundError: 
        
        
        #Sort all the .nc files
        ListOfFiles = np.sort(os.listdir(directory2Data.dir_rawdata))
    
    
    
    
        #open and rewrite IDX files
        fd = open(directory2Data.dir_src+'/'+'pingtime.csv','w')    
        fd2 = open(directory2Data.dir_src+'/'+'FileList.csv','w')
        fd3 = open(directory2Data.dir_src+'/'+'IDX.csv','w')
        fd4 = open(directory2Data.dir_src+'/'+'BeamIDX.csv','w')
        
            
            
        #This helps to prevent opening and closing the same file multiple times
        oldFileName = ''
    
    
        
        #Loop through each file
        CompleteListOfFiles = []
        for i in range(len(ListOfFiles)): 
            
            #Print progress
            tools.printProgressBar(i+1,len(ListOfFiles),prefix = 'Stacking: ', suffix = 'Completed', length = 50)
    
            
            
            #If this is a new file, close the last and open this
            if ListOfFiles[i]!=oldFileName:
                try: 
                    fid_nc.close()
                except: 
                    k=1
                    
                fid_nc = Dataset(os.path.join(directory2Data.dir_rawdata, ListOfFiles[i]),'r')
                oldFileName = ListOfFiles[i]
            

            #Get the correct beam group
            if fid_nc.groups['Sonar'].groups['Beam_group1'].beam_mode == beam_mode: 
                beamgrp = 'Beam_group1'
                
            elif fid_nc.groups['Sonar'].groups['Beam_group2'].beam_mode == beam_mode:
                beamgrp = 'Beam_group2'
            
            
                
            #Get the beam data
            beam_data = fid_nc.groups['Sonar'].groups[beamgrp]
    
    
            #Loop through each ping
            for ii in range(len(beam_data.variables['ping_time'])): 
                fd.write(str(beam_data.variables['ping_time'][ii])+',')
                fd2.write(ListOfFiles[i]+',')
                fd3.write(str(ii)+',')
                fd4.write(beamgrp+',')
    
        #Close files
        fd.close()
        fd2.close()
        fd3.close()
        fd4.close()
        
        
        #Reopen the files
        fd=open(directory2Data.dir_src+'/'+'pingtime.csv','r')
        fd2=open(directory2Data.dir_src+'/'+'FileList.csv','r')
        fd3=open(directory2Data.dir_src+'/'+'IDX.csv','r')
        fd4 = open(directory2Data.dir_src+'/'+'BeamIDX.csv','r')
        
        
        #Get the IDX information
        A = np.asarray(fd.read().split(','))[:-1].astype(np.float)
        B = np.asarray(fd2.read().split(','))[:-1]
        C = np.asarray(fd3.read().split(','))[:-1].astype(np.int)
        D = np.asarray(fd4.read().split(','))[:-1]
        CompleteListOfFiles = np.array((A.T,B.T,C.T)).T



        #Close the files
        fd.close()
        fd2.close()
        fd3.close()
        fd4.close()
        
        
        
    return CompleteListOfFiles, D
        
    
    

    
    
def SecondsBetweenTransect(startTime,unixtime): 

    import datetime
    
    fulldate = datetime.datetime.strptime('1601-01-01 00:00:00.000',"%Y-%m-%d %H:%M:%S.%f")

    starten =  datetime.datetime.strptime(startTime,"%Y%m%d%H%M%S")


    fulldate = (starten-fulldate).total_seconds() #datetime.timedelta(milliseconds=int(startTime))

    seconds = fulldate-unixtime/10000000
                           
    return seconds


    
    
def GetTransectTimes(filename,code): 
        
    from lxml import etree
    doc2 = etree.parse(filename)
    distance_list = doc2.find('distance_list')
    old_stop = ''
    Start = 0
    End = []
    transect_IDX = []
    transect_i = 1
    new_start = False
    
    for i in distance_list:
        start_time = i.get('start_time').replace(' ','T').replace('-','').replace(':','')
        stop_time = i.find('stop_time').text.replace(' ','T').replace('-','').replace(':','')

        if Start == 0:
            Start = start_time
            transect_IDX = np.hstack((transect_IDX,code + '_'+str(transect_i)))
            transect_i = transect_i+1
        elif new_start == True: 
            Start = np.hstack((Start,start_time))
            new_start = False
            transect_IDX = np.hstack((transect_IDX,code + '_'+str(transect_i)))
            transect_i = transect_i+1
            
        if old_stop!='':
            if old_stop != start_time: 
                End = np.hstack((End,stop_time))
                new_start = True
            
        
        old_stop = stop_time
        
      
    #add last time
    End = np.hstack((End,stop_time))
    TimeIDX  = np.vstack((transect_IDX.T,Start.T,End.T)).T    
    return TimeIDX
    
    
    
    
    
def main(TS = 0): 
    
    
    #Clear the terminal window and display a welcome screen.
    os.system('cls' if os.name == 'nt' else 'clear')
    tools.WelcomeScreen('PySonar.py')
    
    
    
    threshold = 2
    RemoveToCloseValues = 30

    #Size of the school trajectory, (Do not change from 3 unless needed)
    R_s = 3
    res = R_s/0.18
    

    maxPingInFile = 2000
    MaxNumberOfFilesInNC  = 100
    recompute = True
#    current_dir = os.getcwd()
    
    SonarEquipment = ['SU90','SX90','SH90']
    
    
    #Something to adapt for data organization at the IMR
    if platform.system() == 'Linux': 
        OS = '//data'
    else: 
        OS = '//ces.imr.no'
        
        
        
        
    #Get the work directory to where all the files is stoored
    WorkDirectory = OS+ '/mea/2018_Redus'
    
    
    
    #Timeseries = GetTimeSeriesInfoFromNMDAPI()
    
    #for Cruice in Timeseries: 
    
    
    
    #Parse the structure file, and get information of each survey
    doc = minidom.parse(os.getcwd()+'/src/DataStructure.xml')
    
    
    
    #Prepare to go through all cruices
    CruiceCode = doc.getElementsByTagName("CruiceCode") 
    
    
    
    #Loop through each cruice
    for CruiceIndex in CruiceCode:
        print('Start on survey: '+CruiceIndex.getAttribute('code'))
        
#        SendMail.send_email('Start on survey: '+CruiceIndex.getAttribute('code'))
        
        
        #Organise and copy the data from tapeserver to scratch
        print('    *Organise data      ',end='\r')
        
        try: 
            directory2Data =tools.OrginizeData(CruiceIndex,WorkDirectory,OS)
        except: 
            SendMail.send_email('Failed to organize '+CruiceIndex.getAttribute('code'))
        
        
            
        #Get transecttimes from Luf20
        for dirpath,_,filenames in os.walk(directory2Data.dir_src + '/EKLUF20/'):
            for f in filenames:
                filename =  os.path.abspath(os.path.join(dirpath, f))
                TimeIDX = GetTransectTimes(filename,str(CruiceIndex.getAttribute('code')))
        
                
                
                    
        CompleteListOfFiles, beamgrp = ComListOfFiles(directory2Data,'Horizontal')
              
        
                    
#        tools.DataConverter(CruiceIndex,WorkDirectory,os.getcwd(),maxPingInFile,
#                  MaxNumberOfFilesInNC,directory2Data)
    
    
        
        
        
        
        #Loop through each transect in a randomized maner
        NumTransect = np.arange(len(TimeIDX))
        for Transect in NumTransect: 
            print('    -Start on transect: '+ TimeIDX[Transect,0])
        
            
            #Go through each equipment
            for eqip in SonarEquipment: 
                if eqip == CruiceIndex.getAttribute("Equipment"): 
                    
                    
                    #Get list of files in one transect
                    ShortListOfFiles = GetShortListOfFiles(CompleteListOfFiles,TimeIDX[Transect,1].replace('T',''),TimeIDX[Transect,2].replace('T',''))
                    

                    
                    #Get fish from vertical fan of beams
#                    MakeVerticalIndex(ShortListOfFiles,RemoveToCloseValues,R_s,res,directory2Data,directory2Data.dir_rawdata,beamgrp)                    
                    
                    
                    
                    #Make the search matrix
#                    try: 
                    if os.path.isfile(directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat') == False: 
                        MakeSearch(ShortListOfFiles,RemoveToCloseValues,R_s,res,directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat',directory2Data.dir_rawdata,beamgrp)
                    elif recompute == True: 
                         MakeSearch(ShortListOfFiles,RemoveToCloseValues,R_s,res,directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat',directory2Data.dir_rawdata,beamgrp)
#                    except: 
#                        SendMail.send_email('failed to make search matrix for transect '+ TimeIDX[Transect,0])
                        

                    #Do vertical
                    




                        
                    #Make the work stuff    
#                    if os.path.isfile(directory2Data.dir_work+'/'+TimeIDX[Transect,0]+'.txt') == False: 
#                        if os.path.isfile(directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat') == True: 
#                            print('    *Make Work',end='\r')
#                            MakeWork(True, directory2Data,'',
#             '','',1,3)
#                            f = open(directory2Data.dir_work+'/'+TimeIDX[Transect,0]+'.txt','w')
#                            f.write('0')
#                            f.close()
#                    elif recompute == True: 
#                        if os.path.isfile(directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat') == True: 
#                            print('    *Make Work',end='\r')
#                            MakeWork(True, directory2Data,'',
#             '','',1,3)
#                            f = open(directory2Data.dir_work+'/'+TimeIDX[Transect,0]+'.txt','w')
#                            f.write('0')
#                            f.close()
#                            
#                            
                        
#                    
#        print('    *Make LUF20',end='\r')
#                    
        
        
        
                    
                    
                    
        #Get List of files 
    
            
 
                    
        
        
        #Go to transect    
        
        
            #Get idX of list of files
        
        
            #Make Search matrix
            
            
            #Make work files
        
            
            
        #Make LUF20
        
        
        
        
        
        
#        
#        #Bookeeping      
#        transectCode = []
#        startTime = []
#        endTime = []
#        
#
#        #Get vessel information
#        
#
#        print('Start on cruice ID: ' + str(cruiceCode))        
#        
#        
#        #Get the directory of the data (This should be standard)
#        dir_temp = OS+'/cruise_data/'+cruiceCode[:4]+'/S'+cruiceCode+'_P'+vesselCode+'/ACOUSTIC_DATA'
#        
#
#        #Make new folder if the do not exist 
#        ProtocolToMakeNewFolders(WorkDirectory,cruiceCode,vesselCode)
#
#
#        
#        #Set the directory path
#        DirectoryToNC = WorkDirectory+'/'+cruiceCode[:4]+'/S'+cruiceCode+'_P'+vesselCode+'/ACOUSTIC_DATA/SU90/netcdf'
#        DirectoryToWORK = WorkDirectory+'/'+cruiceCode[:4]+'/S'+cruiceCode+'_P'+vesselCode+'/ACOUSTIC_DATA/SU90/WorkFiles'
#        DirectoryToRESULT = WorkDirectory+'/'+cruiceCode[:4]+'/S'+cruiceCode+'_P'+vesselCode+'/ACOUSTIC_DATA/SU90/Result'
##        DirectoryToSRC = WorkDirectory+'/'+cruiceCode[:4]+'/S'+cruiceCode+'_P'+vesselCode+'/ACOUSTIC_DATA/SU90/src'
#        DirectoryToNCProg = WorkDirectory+'/'+cruiceCode[:4]+'/S'+cruiceCode+'_P'+vesselCode+'/ACOUSTIC_DATA/SU90/src/NCconvertProgress'
#        
#
#
#        print('    *Start converting files',end='\r')
##        protocolForNetCDFmaker(dir_temp,DirectoryToNC,CruiceTruePath,DirectoryToNCProg)
#        print('    -Files Converted       ')
##    
#                
#        
#                
#
#            
#                
#        elements =(np.arange(len(transectCode)))
#        
#        
#        
#        if Random == True: 
#            np.random.shuffle(elements)
#        
#            
#            
#        for i in elements: 
#            print('    -Start analysing transect: ' +str(transectCode[i]))
##            MakeSearchMatrixVert.MakeSearchMatrixVert(DirectoryToNC,DirectoryToWORK,DirectoryToRESULT,
##                                          transectCode[i],startTime[i],endTime[i],RemoveToCloseValues,R_s,res)
#
#
#
#
#            if os.path.isfile((DirectoryToWORK+'/Vertical/vert_'+str(transectCode[i])+'.mat')) == False: 
##                try: 
#                print('    *Start making search matrix vertical', end = '\r')
#                MakeSearchMatrixVert.MakeSearchMatrixVert(DirectoryToNC,DirectoryToWORK,DirectoryToRESULT,
#                                              transectCode[i],startTime[i],endTime[i],RemoveToCloseValues,R_s,res)
#                fid = open(DirectoryToRESULT+'/SearchMatrixVert'+str(transectCode[i])+'.txt','w')
#                fid.close()
#                print('    *Search matrix vertical is finnished', end = '\r')
#                send_email('Search Matrix Vertical'+transectCode[i]+ ' is finnished')
#    
##                except Exception as e:
##                    
##                    try: 
##                        exc_type,exc_value,exc_traceback = sys.exc_info() 
##                        send_email(transectCode[i]+'   error occured:    ') # +exc_type+ '    ' + exc_value + '    ' + exc_traceback )
##                    except TimeoutError: 
##                        send_email(transectCode[i]+'    TimeOut: ')
#            elif ReCompute == True: 
#                try: 
#                    print('    *Start making vertical search matrix', end = '\r')
#                    MakeSearchMatrixVert.MakeSearchMatrixVert(DirectoryToNC,DirectoryToWORK,DirectoryToRESULT,
#                                                  transectCode[i],startTime[i],endTime[i],RemoveToCloseValues,R_s,res)
#                    fid = open(DirectoryToRESULT+'/SearchMatrixVert'+str(transectCode[i])+'.txt','w')
#                    fid.close()
#                    print('    *Search matrix vertical is finnished', end = '\r')
#                    send_email('Search Matrix vertical '+transectCode[i]+ ' is finnished')
#    
#                except Exception as e:
#                    
#                    try: 
#                        exc_type,exc_value,exc_traceback = sys.exc_info() 
#                        send_email(transectCode[i]+'   error occured:    ') # +exc_type+ '    ' + exc_value + '    ' + exc_traceback )
#                    except TimeoutError: 
#                        send_email(transectCode[i]+'    TimeOut: ')
                        
#            else: 
#                print('')
#                print('Viewing')
#                print('')
#                MakeSearchMatrixVert.MakeSearchMatrixVertRun(DirectoryToNC,DirectoryToWORK,DirectoryToRESULT,
#                                                  transectCode[i],startTime[i],endTime[i],RemoveToCloseValues,R_s,res)
                
                
#                
#                
#                        
#            if not os.path.isfile((DirectoryToRESULT+'/SearchMatrix'+str(transectCode[i])+'.mat')): 
#                try: 
#                    print('    *Start making search matrix', end = '\r')
#                    MakeSearchMatrix.MakeSearchMatrix(DirectoryToNC,DirectoryToWORK,DirectoryToRESULT,
#                                                  transectCode[i],startTime[i],endTime[i],RemoveToCloseValues,R_s,res)
#                    print('    *Search matrix is finnished', end = '\r')
#                    send_email('Search Matrix '+transectCode[i]+ ' is finnished')
#    
#                except Exception as e:
#                    
#                    try: 
#                        exc_type,exc_value,exc_traceback = sys.exc_info() 
#                        send_email(transectCode[i]+'   error occured:    ') # +exc_type+ '    ' + exc_value + '    ' + exc_traceback )
#                    except TimeoutError: 
#                        send_email(transectCode[i]+'    TimeOut: ')
#            elif ReCompute == True: 
#                try: 
#                    print('    *Start making search matrix', end = '\r')
#                    MakeSearchMatrix.MakeSearchMatrix(DirectoryToNC,DirectoryToWORK,DirectoryToRESULT,
#                                                  transectCode[i],startTime[i],endTime[i],RemoveToCloseValues,R_s,res)
#                    print('    *Search matrix is finnished', end = '\r')
#                    send_email('Search Matrix '+transectCode[i]+ ' is finnished')
#    
#                except Exception as e:
#                    
#                    try: 
#                        exc_type,exc_value,exc_traceback = sys.exc_info() 
#                        send_email(transectCode[i]+'   error occured:    ') # +exc_type+ '    ' + exc_value + '    ' + exc_traceback )
#                    except TimeoutError: 
#                        send_email(transectCode[i]+'    TimeOut: ')
#            

#            #MakeWorkStuffi
#            if not os.path.isfile(DirectoryToWORK+'/'+transectCode[i]+'.txt'): 
#                print('    *Start making work files',end='\r')
#                makePingWork.makePingWork(DirectoryToNC,DirectoryToWORK,DirectoryToRESULT,
#                                           transectCode[i],threshold,res)
#                print('    *Stop making work files', end = '\r')
#                try: 
#                    print('    *Start making work files',end='\r')
#                    makePingWork.makePingWork(DirectoryToNC,DirectoryToWORK,DirectoryToRESULT,
#                                               transectCode[i],threshold,res)
#                    print('    *Stop making work files', end = '\r')
#                
#                
#                    fid = open(DirectoryToWORK+'/'+transectCode[i]+'.txt','w')
#                    fid.close()
#                    send_email('Work stuff from '+transectCode[i]+' is finnished')
#                except FileNotFoundError : 
#                    print('fileNotFound')
#                except Exception: 
#                    exc_type,exc_value,exc_traceback = sys.exc_info() 
#                    send_email(transectCode[i]+'   Feil i Work')
#            elif ReCompute == True: 
#                try: 
#                    print('    *Start making work files',end='\r')
#                    makePingWork.makePingWork(DirectoryToNC,DirectoryToWORK,DirectoryToRESULT,
#                                               transectCode[i],threshold,res)
#                    print('    *Stop making work files', end = '\r')
#                
#                
#                    fid = open(DirectoryToWORK+'/'+transectCode[i]+'.txt','w')
#                    fid.close()
#                    send_email('Work stuff from '+transectCode[i]+' is finnished')
#                except FileNotFoundError : 
#                    print('fileNotFound')
#                except Exception: 
#                    exc_type,exc_value,exc_traceback = sys.exc_info() 
#                    send_email(transectCode[i]+'   Feil i Work')

#        print('    *Cruice is finnished')
if __name__ == '__main__':    
    
#    first_arg = sys.argv[1]
#
#    print(first_arg)
#    second_arg = sys.argv[2]

    
    main(TS = 0)
    
