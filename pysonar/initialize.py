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
import scipy.io as sc
from MakeSearch import MakeSearch



def ListOfFiles(dir2file,dir2nc): 
    if os.path.exists(dir2file+'/'+'listOfFiles.mat')==False: 
        print('    *Making ListOfFiles',end='\r')
        
        ListOfFiles = os.listdir(dir2nc)
        
        mat = []
#        f=open(DirectoryToRESULT+'/'+'listOfFiles.txt' ,'w')
        for i in range(len(ListOfFiles)): 
            tools.printProgressBar(i + 1, len(ListOfFiles), prefix = 'Progress:', suffix = 'Complete', length = 50)
#            f.write(str(ListOfFiles[i])+'\n')
            fileIndex = ListOfFiles[i]

            try: 
                Test = (fileIndex[5:-4]+(20-len(fileIndex[5:-4]))*'0')
                mat=np.hstack((mat,Test))#int(fileIndex[5:-4]+(20-len(fileIndex[5:-4]))*'0')))
            except ValueError: 
                print('Bad file name', end='\r')   
            
#        f.close()
        sc.savemat(dir2file+'/'+'listOfFiles.mat',mdict={'mat':mat})
    else: 
        mat = sc.loadmat(dir2file+'/'+'listOfFiles.mat')['mat']

    return mat
    
    
def AlternativTransectTime(filename,code): 

        
    from lxml import etree
    doc2 = etree.parse(filename)
    distance_list = doc2.find('distance_list')
    old_stop = ''
    Start = []
    End = []
    transect_IDX = []
    transect_i = 1
    new_start = False
    for i in distance_list:
        
        start_time = i.get('start_time').replace(' ','T').replace('-','').replace(':','')
        stop_time = i.find('stop_time').text.replace(' ','T').replace('-','').replace(':','')
        
        if Start == []:
            Start = np.hstack((Start,start_time))
            transect_IDX = np.hstack((transect_IDX,code + '_'+str(transect_i)))
            transect_i = transect_i+1
        elif new_start == True: 
            Start = np.hstack((Start,start_time))
            new_start = False
            transect_IDX = np.hstack((transect_IDX,code + '_'+str(transect_i)))
            transect_i = transect_i+1
            
        if old_stop!=start_time: 
            End = np.hstack((End,stop_time))
            new_start = True
        old_stop = stop_time
      
    #add last time
    End = np.hstack((End,stop_time))
    
    
    TimeIDX  = np.vstack((transect_IDX.T,Start.T,End.T)).T    
    
    return TimeIDX
    
    
def main(): 
    
    #Get the instruction file
    #The instruction file includes the location of each survey.
    os.system('cls' if os.name == 'nt' else 'clear')
    tools.WelcomeScreen('PySonar.py')
    
    
    
    threshold = 2
    RemoveToCloseValues = 30

    #Size of the school trajectory, (Do not change from 3 unless needed)
    R_s = 3
    res = R_s/0.18
    

    maxPingInFile = 2000
    MaxNumberOfFilesInNC  = 100
    recompute = False
    current_dir = os.getcwd()
    
    
    
    #Something to adapt for data organization at the IMR
    if platform.system() == 'Linux': 
        OS = '//data'
    else: 
        OS = '//ces.imr.no'
        
        
        
        
    #Get the work directory to where all the files is stoored
    WorkDirectory = OS+ '/mea/2018_Redus'
    
    
    
    #Parse the structure file, and get information of each survey
    doc = minidom.parse(os.getcwd()+'/src/DataStructure.xml')
    
    
    
    #Prepare to go through all cruices
    CruiceCode = doc.getElementsByTagName("CruiceCode") 
    
    
    
    #Loop through each cruice
    for CruiceIndex in CruiceCode:
        print('Start on survey: '+CruiceIndex.getAttribute('code'))
        
        
        
        
        #Organise and copy the data from tapeserver to scratch
        print('    *Orginizing data      ',end='\r')
        directory2Data =tools.OrginizeData(CruiceIndex,WorkDirectory,OS)

        
        
        
        #Convert data to netcdf
        print('    *Convert data         ',end='\r')
#        os.chdir(current_dir)
#        tools.DataConverter(CruiceIndex,WorkDirectory,current_dir,
#                            maxPingInFile,MaxNumberOfFilesInNC,directory2Data)
        
        
        
        
        #Get list of transect with start and stop times
        TimeIDX=tools.TransectTimeIDX(CruiceIndex)
        
        if TimeIDX == []: 
#            import glob
#            print(glob.glob(directory2Data.dir_src + '/EKLUF20'))
            
            filename = 'C:/Users/sindrev/Desktop/ListUserFile20__L7284.0-1959.0.txt'
            TimeIDX = AlternativTransectTime(filename,str(CruiceIndex.getAttribute('code')))
        
        print(TimeIDX)
#            
#        
#        
#        
#        #Get list of files in cruice
#        ListOfFilesInFolder = ListOfFiles(directory2Data.dir_src,directory2Data.dir_nc)
        
        
        
        
        #Loop through each transect in a randomized maner
#        NumTransect = np.arange(len(TimeIDX))
#        np.random.shuffle(NumTransect)
#        for Transect in NumTransect: 
#            print('    -Start on transect: '+ TimeIDX[Transect,0])
#        
#            
#            #Go through each equipment
#            for i in ['SU90','SX90','SH90']: 
#                if i == CruiceIndex.getAttribute("Equipment"): 
#                    
#                    
#                    #Get list of files
#                    #Need to fix the equipment stuff.
#                    #Try ListOfFilesInFolder[0][:4]
#                    ListOfFilesWithinTimeInterval = [i+'-'+str(i)+'.nc' for i in np.sort(ListOfFilesInFolder) if TimeIDX[Transect,1]*1E6 <= i <=TimeIDX[Transect,2]*1E6]
#                    
#
#
#                    #Make the search matrix
#                    if os.path.isfile(directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat') == False: 
#                        MakeSearch(ListOfFilesWithinTimeInterval,RemoveToCloseValues,R_s,res,directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat')
#                        print('    *Make New Search',end='\r')
#                    elif recompute == True: 
#                        print('    *Make New Search',end='\r')
#                        MakeSearch(ListOfFilesWithinTimeInterval,RemoveToCloseValues,R_s,res,directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat')
                        
                        
#                        
#                    #Make the work stuff    
#                    if os.path.isfile(directory2Data.dir_work+'/'+TimeIDX[Transect,0]+'.txt') == False: 
#                        if os.path.isfile(directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat') == True: 
#                            print('    *Make Work',end='\r')
#                            MakeWork()
#                            f = open(directory2Data.dir_work+'/'+TimeIDX[Transect,0]+'.txt','w')
#                            f.write('0')
#                            f.close()
#                    elif recompute == True: 
#                        if os.path.isfile(directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat') == True: 
#                            print('    *Make Work',end='\r')
#                            MakeWork()
#                            f = open(directory2Data.dir_work+'/'+TimeIDX[Transect,0]+'.txt','w')
#                            f.write('0')
#                            f.close()
#                            
                            
                        
                    
        print('    *Make LUF20',end='\r')
                    
                    
                    
                    
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
    main()
    
