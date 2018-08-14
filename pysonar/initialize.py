
import xml.dom.minidom as minidom
import os,datetime,platform, glob, random
from lxml import etree
import numpy as np
from SendMail import send_email
from tools import tools
from shutil import copyfile
import scipy.io as sc
from MakeSearch import MakeSearch
#from Index import MakeVerticalIndex
from MakeVerticalIndex import haversine
from VerticalLuf20 import MakeVerticalLuf20
from MakeWork import MakeWork
from HorizontalLuf20 import MakeHorizontalLuf20
from netCDF4 import Dataset
from numba import jit




    

#@jit
def GetShortListOfFiles(CompleteListOfFiles,startTime,endTime):
    ShortListOfFiles = []

    for i in range(len(CompleteListOfFiles[:][:])):
        
        tools.printProgressBar(i + 1, len(CompleteListOfFiles[:,0]), prefix = 'Make short list:', suffix = 'Completed       ', length = 50)
        
#        print(CompleteListOfFiles[:][i])
        
        start = SecondsBetweenTransect(startTime,int(CompleteListOfFiles[:][i][0]))
        end = SecondsBetweenTransect(endTime,int(CompleteListOfFiles[:][i][0]))
        
        
        if start<0:
            if end >0:
                if ShortListOfFiles == []:
                    ShortListOfFiles = CompleteListOfFiles[i,:]
                else:
                    ShortListOfFiles = np.vstack((ShortListOfFiles,CompleteListOfFiles[i,:]))

    return ShortListOfFiles
    
    
def GetTransectTimes(filename,code):

    
    doc2 = etree.parse(filename)
    distance_list = doc2.find('distance_list')
    old_stop = ''
    Start = 0
    End = []
    transect_IDX = []
    transect_i = 1
    new_start = False

    start_time = []
    log_start = []
    stop_time = []
    lat_start = []
    lon_start = []
    lat_stop = []
    lon_stop = []
    
    
    for i in distance_list:
        
        start_time = np.hstack((start_time,i.get('start_time').replace(' ','T').replace('-','').replace(':','')))
        log_start = np.hstack((log_start,i.get('log_start')))
        stop_time = np.hstack((stop_time,i.find('stop_time').text.replace(' ','T').replace('-','').replace(':','')))
        lat_start = np.hstack((lat_start,i.find('lat_start').text))
        lon_start = np.hstack((lon_start,i.find('lon_start').text))
        lat_stop = np.hstack((lat_stop,i.find('lat_stop').text))
        lon_stop = np.hstack((lon_stop,i.find('lon_stop').text))

    return start_time,log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop
    
    
    
def GetLuf20Info(directory2Data,CruiceIndex):
    
    for dirpath,_,filenames in os.walk(directory2Data.dir_src + '/EKLUF20/'):
        for f in filenames:
            filename =  os.path.abspath(os.path.join(dirpath, f))
            TimeIDX = GetTransectTimesIDX(filename,str(CruiceIndex.getAttribute('code')))
            start_time,log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop = GetTransectTimes(filename,str(CruiceIndex.getAttribute('code')))
        if filenames == []: 
            
            
            for path, subdirs, files in os.walk(directory2Data.dir_acoustic_data):
                for name in files:
                    if 'ListUserFile20__L'in os.path.join(path,name): 
                        try: 
                            copyfile(os.path.join(path,name),directory2Data.dir_src+'\EKLUF20\\'+name)
                        except shutil.SameFileError: 
                            k=1
                        break
                    
            for dirpath,_,filenames in os.walk(directory2Data.dir_src + '/EKLUF20/'):
                for f in filenames:
                    filename =  os.path.abspath(os.path.join(dirpath, f))
                    TimeIDX = GetTransectTimesIDX(filename,str(CruiceIndex.getAttribute('code')))
                    start_time,log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop = GetTransectTimes(filename,str(CruiceIndex.getAttribute('code')))
                
        
    return  start_time,log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop, TimeIDX
    
    
    
    
def ComListOfFiles(directory2Data,beam_mode):
    #This function makes an idx of all files


    #Open the IDX files
    try:
        fd = open(directory2Data.dir_src+'/'+beam_mode+'pingtime.csv','r')
        fd2 = open(directory2Data.dir_src+'/'+beam_mode+'FileList.csv','r')
        fd3 = open(directory2Data.dir_src+'/'+beam_mode+'IDX.csv','r')
        fd4 = open(directory2Data.dir_src+'/'+beam_mode+'BeamIDX.csv','r')
#        fd5 = open(directory2Data.dir_src+'/'+beam_mode+'LogDist.csv','r')


        A = np.asarray(fd.read().split(','))[:-1]
        B = np.asarray(fd2.read().split(','))[:-1]
        C = np.asarray(fd3.read().split(','))[:-1]
        D = np.asarray(fd4.read().split(','))[:-1]
        CompleteListOfFiles = np.array((A.T,B.T,C.T)).T



    # If the files don't exist, create them
    except FileNotFoundError:


        #Sort all the .nc files
        ListOfFiles = np.sort(os.listdir(directory2Data.dir_rawdata))




        #open and rewrite IDX files
        fd = open(directory2Data.dir_src+'/'+beam_mode+'pingtime.csv','w')
        fd2 = open(directory2Data.dir_src+'/'+beam_mode+'FileList.csv','w')
        fd3 = open(directory2Data.dir_src+'/'+beam_mode+'IDX.csv','w')
        fd4 = open(directory2Data.dir_src+'/'+beam_mode+'BeamIDX.csv','w')
#        fd5 = open(directory2Data.dir_src+'/'+beam_mode+'LogDist.csv','w')



        #This helps to prevent opening and closing the same file multiple times
        oldFileName = ''

        
        #Loop through each file
        CompleteListOfFiles = []
        for i in range(len(ListOfFiles)):

            #Print progress
            tools.printProgressBar(i+1,len(ListOfFiles),prefix = 'Stacking: ', suffix = 'Completed', length = 50)


            try: 
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
                    timeID = np.where(abs(beam_data.variables['ping_time'][ii]-fid_nc.groups['Platform'].variables['time1'][:]/100)==
                                          np.nanmin(abs(beam_data.variables['ping_time'][ii]-fid_nc.groups['Platform'].variables['time1'][:]/100)))
                    
                    
                    fd.write(str(beam_data.variables['ping_time'][ii])+',')
                    fd2.write(ListOfFiles[i]+',')
                    fd3.write(str(ii)+',')
                    fd4.write(beamgrp+',')
            except: 
                print('   *bad ping*')
                
        #Close files
        fd.close()
        fd2.close()
        fd3.close()
        fd4.close()
        
        
        #Reopen the files
        fd=open(directory2Data.dir_src+'/'+beam_mode+'pingtime.csv','r')
        fd2=open(directory2Data.dir_src+'/'+beam_mode+'FileList.csv','r')
        fd3=open(directory2Data.dir_src+'/'+beam_mode+'IDX.csv','r')
        fd4 = open(directory2Data.dir_src+'/'+beam_mode+'BeamIDX.csv','r')
#        fd5 = open(directory2Data.dir_src+'/'+beam_mode+'LogDist.csv','r')


        #Get the IDX information
        A = np.asarray(fd.read().split(','))[:-1]
        B = np.asarray(fd2.read().split(','))[:-1]
        C = np.asarray(fd3.read().split(','))[:-1]
        D = np.asarray(fd4.read().split(','))[:-1]
#        logdist = np.asarray(fd5.read().split(','))[:-1]
        CompleteListOfFiles = np.array((A.T,B.T,C.T)).T



        #Close the files
        fd.close()
        fd2.close()
        fd3.close()
        fd4.close()
        fd5.close()


    
    return CompleteListOfFiles, D
    
    
def SecondsBetweenTransect(startTime,unixtime):

    fulldate = datetime.datetime.strptime('1601-01-01 00:00:00.000',"%Y-%m-%d %H:%M:%S.%f")

    starten =  datetime.datetime.strptime(startTime,"%Y%m%d%H%M%S")

    fulldate = (starten-fulldate).total_seconds() #datetime.timedelta(milliseconds=int(startTime))

    seconds = fulldate-unixtime/10000000

    return seconds
    
    
    
def GetTransectTimesIDX(filename,code):

    
    doc2 = etree.parse(filename)
    distance_list = doc2.find('distance_list')
    old_stop = ''
    Start = 0
    End = []
    transect_IDX = []
    transect_i = 1
    new_start = False

    start_time = []
    log_start = []
    stop_time = []
    lat_start = []
    lon_start = []
    lat_stop = []
    lon_stop = []
    
    
    for i in distance_list:
        
        start_time = i.get('start_time').replace(' ','T').replace('-','').replace(':','')
        log_start = i.get('log_start')
        stop_time = i.find('stop_time').text.replace(' ','T').replace('-','').replace(':','')
        lat_start = i.find('lat_start').text
        lon_start = i.find('lon_start').text
        lat_stop = i.find('lat_stop').text
        lon_stop = i.find('lon_start').text
#
#    return start_time,log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop
#        print(start_time)
#        print(log_start,lon_stop)
#
        if transect_i <10: 
            trnsID = '00'+str(transect_i)
        elif transect_i <100: 
            trnsID = '0'+str(transect_i)
        else:
            trnsID = str(transect_i)
        
            
        if old_stop!=start_time: 
            End = np.hstack((End,old_stop))
            Start = np.hstack((Start,start_time))
            transect_IDX = np.hstack((transect_IDX,code + '_'+trnsID))
            transect_i = transect_i+1
            
            
            
        if Start == 0:
            Start = start_time
            transect_IDX = np.hstack((transect_IDX,code + '_'+trnsID))
            transect_i = transect_i+1
    
            
        old_stop = stop_time

    #add last time
    End = np.hstack((End,stop_time))
    
    
    TimeIDX  = np.vstack((transect_IDX.T,Start[1:].T,End[1:].T)).T
    return TimeIDX
    
    
def main(threshold, RemoveToCloseValues, R_s, recompute, reconvert,GO_horizontal, GO_vertical, Stack):

    
    
    
    #Clear the terminal window and display a welcome screen.
    os.system('cls' if os.name == 'nt' else 'clear')
    tools.WelcomeScreen('PySonar.py')


    
    
    #Some user innputs that will be moved
    res = R_s/0.18
    maxPingInFile = 2000
    MaxNumberOfFilesInNC  = 100
    
    

    
    
    
    
    #Sonar names
    SonarEquipment = ['SU90','SX90','SH90']






    #Something to adapt for data organization at the IMR
    #For other organisations, this should be changed
    if platform.system() == 'Linux':
        OS = '//data'
    else:
        OS = '//ces.imr.no'



        
        

    #Get the work directory to where all the files is stoored
    #This should be merged with the one above
    WorkDirectory = OS+ '/mea/2018_Redus'
    WorkDirectory = 'F:'





    #Parse the structure file, and get information of each survey
    doc = minidom.parse(os.getcwd()+'/src/DataStructure.xml')



    
    
    #Prepare to go through all cruices
    CruiceCode = doc.getElementsByTagName("CruiceCode")



    
    #Loop through each cruice
    for CruiceIndex in CruiceCode:
        
        print(CruiceIndex)
        
        
        os.system('cls' if os.name == 'nt' else 'clear')
        tools.WelcomeScreen('PySonar.py')
        print('Start on survey: '+CruiceIndex.getAttribute('code'))

        
        
        
        #Get cruice info
        nation = CruiceIndex.getAttribute('code')
        cruice_id =CruiceIndex.getAttribute('code')
        vplatform = CruiceIndex.getAttribute('VCode')
        
        
        
        
        
        #Uncomment to enable
#        send_email('Start on survey: '+CruiceIndex.getAttribute('code'))

        

        

        #Organise and copy the data from tapeserver to scratch
        print('    *Organise data                ',end='\r')
#        try:
        directory2Data =tools.OrginizeData(CruiceIndex,WorkDirectory,OS)
#        except:
#            send_email('Failed to organize '+CruiceIndex.getAttribute('code'))

            
            
            
            
            

        #This will delete the finnished file, and all raw data are converted
        #into netcdf. 
        if reconvert == True: 
            try: 
                os.remove(directory2Data.dir_src+'/NCconvertProgress/'+os.listdir(directory2Data.dir_src+'/NCconvertProgress')[0])
            except IndexError: 
                k=1
                
                
                
                
                
                
#        tools.DataConverter(CruiceIndex,WorkDirectory,os.getcwd(),maxPingInFile,
#                  MaxNumberOfFilesInNC,directory2Data, reconvert)
                
#        #Convert raw file to netcdf                   
#        print('    *Convert .raw data                     ',end='\r') 
#        try: 
#            tools.DataConverter(CruiceIndex,WorkDirectory,os.getcwd(),maxPingInFile,
#                      MaxNumberOfFilesInNC,directory2Data, reconvert)
#            print('    *.raw data donverted                      ',end='\r') 
#        except: 
#            send_email('Failed to convert files')
#        
        
        
        print('Continue on survey: '+CruiceIndex.getAttribute('code'))
            
            
        

        #Delete old index files
        if reconvert == True: 
            for i in os.listdir(directory2Data.dir_src) : 
                if '.csv' in i: 
                    os.remove(directory2Data.dir_src+'/'+i)

                    
                        
                        
                
                        
            
        if Stack == True: 
                        
                        
                        
            #Make new index files, and/or load it
            print('    *Load IDX                       ',end='\r') 
            try: 
                CompleteListOfFiles_horizontal,beamgrp_hor = ComListOfFiles(directory2Data,'Horizontal')
            except: 
                send_email('Make List of files for horizontal failed!!!')
        
            try: 
                CompleteListOfFiles_vertical,beamgrp_vert = ComListOfFiles(directory2Data,'Vertical')
            except: 
                send_email('Make List of files for vertical failed!!!')
               
            CompleteListOfFiles_horizontal,beamgrp_hor = ComListOfFiles(directory2Data,'Horizontal')
            CompleteListOfFiles_vertical,beamgrp_vert = ComListOfFiles(directory2Data,'Vertical')
                
                
            print('    *IDX loaded                      ',end='\r') 
                
                
                
                
                
                
            
            
            
        #Get transect-times from Luf20
        #A system to copy the files from server must be implemented
        start_time,log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop, TimeIDX = GetLuf20Info(directory2Data,CruiceIndex)
        print('LUF20 scanned',end='\r')
        print('Scanning LUF20',end='\r')
                
                
                
                
                
                
                    
                    
                
                
                
                
                
        #Loop through each transect in a randomized maner
        if GO_horizontal == True: 
            print('Start making search matrix',end='\r')
            
            
            
            
            #Get list of all transect and randomize it
            NumTransect = np.arange(TimeIDX.shape[0])
            random.shuffle(NumTransect)
            
            
            
            
            
            #Loop through each transect
            for Transect in NumTransect:
                
    
                
                
                
                #Make the search matrix
                if os.path.isfile(directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat') == False:
                    print('    -Start on transect: '+ TimeIDX[Transect,0],end='\r')
#                    send_email('Start on horizontal transect: '+ TimeIDX[Transect,0])
            
                    
                    
                
                    #Get list of files within time intervals
                    ShortListOfFiles_horizontal = GetShortListOfFiles(CompleteListOfFiles_horizontal,
                                                                      TimeIDX[Transect,1].replace('T','')
                                                                      ,TimeIDX[Transect,2].replace('T',''))
                    
                    
                    
                    #Make the search matrix
                    MakeSearch(ShortListOfFiles_horizontal,
                               RemoveToCloseValues,
                               R_s,
                               res,
                               directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat',
                               directory2Data.dir_rawdata,
                               beamgrp_hor)
                    
                    
                    
                    
                elif recompute == True:
                    print('    -Start on transect: '+ TimeIDX[Transect,0],end='\r')
                    send_email('Start on horizontal transect: '+ TimeIDX[Transect,0])
            
                
                    
                    
                    #Get list of files within time intervals
                    ShortListOfFiles_horizontal = GetShortListOfFiles(CompleteListOfFiles_horizontal,
                                                                      TimeIDX[Transect,1].replace('T','')
                                                                      ,TimeIDX[Transect,2].replace('T',''))
                    
                    
                    
                    
                    
                    #Make search matrix
                    MakeSearch(ShortListOfFiles_horizontal,
                               RemoveToCloseValues,
                               R_s,
                               res,
                               directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat',
                               directory2Data.dir_rawdata,
                               beamgrp_hor)
                    
                    
                    
                    
                    
                
                
            #Make Work files
           
            #Get list of all transect and randomize it
            NumTransect = np.arange(TimeIDX.shape[0])
                        
            
            
            #Loop through each transect
            for Transect in NumTransect:
                
                
                
                
                print('Start on horizontal transect: '+ str(Transect),end='\r')
    #            send_email('Start on horizontal transect: '+ str(Transect))
                
                
                
                
                if os.path.isfile(directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat') == True:
                    MakeWork(True, directory2Data,'','','',1,3)
    
                
                    
                    
                    
                
                
            #Start making the luf20 files
    #        send_email('Start Making Horizontal LUF20: '+ str(Transect))
            MakeHorizontalLuf20(CompleteListOfFiles_horizontal,directory2Data,start_time,
                                log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop,
                                TimeIDX, nation,cruice_id,vplatform)
                    
                    
            
            
            
            
            
            
                
        #Loop through each transect in a randomized maner
#        if GO_vertical == True: 
#            
#            
#            
#            
#            
#            #Get transact list and randomise it
#            NumTransect = np.arange(len(start_time))
#            random.shuffle(NumTransect)
#            
#            
#            
#            
#            
#            
#            #Loop through each transect
#            for Transect in NumTransect:
#                print('    -Start on log distance: '+ log_start[Transect]=#,end='\r')
#    
#                
#                
#                
#                
#                
#                
#                #Go through each equipment
#                for eqip in SonarEquipment:
#                    if eqip == CruiceIndex.getAttribute("Equipment"):
#                        if not os.path.isfile(directory2Data.dir_work+'/'+'Vertical_T'+str(log_start[Transect])+'.mat'): 
#                    
#                            
#                            print('Start making short list')
#                            #Get list of files in one transect
#                            ShortListOfFiles = GetShortListOfFiles(CompleteListOfFiles_vertical,
#                                                                   start_time[Transect].replace('T','')
#                                                                   ,stop_time[Transect].replace('T',''))
#    
#                            
#                            
#                            
#                            if not ShortListOfFiles == []:
#                                MakeVerticalIndex(ShortListOfFiles
#                                                  ,RemoveToCloseValues
#                                                  ,R_s,res
#                                                  ,directory2Data
#                                                  ,directory2Data.dir_rawdata
#                                                  ,beamgrp_vert,
#                                                  start_time[Transect],
#                                                  log_start[Transect],
#                                                  stop_time[Transect],
#                                                  lat_start[Transect],
#                                                  lat_stop[Transect],
#                                                  lon_start[Transect],
#                                                  lon_stop[Transect])
#                                
#                                
#
#                            else: 
#                                sc.savemat(directory2Data.dir_work+'/'+'Vertical_T'+str(log_start[Transect])+'.mat',mdict={'empty':0})
#                                
#                                
#                                
#                        elif recompute == True: 
#                    
#                            
#                            
#                            #Get list of files in one transect
#                            ShortListOfFiles = GetShortListOfFiles(CompleteListOfFiles_vertical,
#                                                                   start_time[Transect].replace('T','')
#                                                                   ,stop_time[Transect].replace('T',''))
#    
#                            
#                            
#                            if not ShortListOfFiles == []:
#                                MakeVerticalIndex(ShortListOfFiles
#                                                  ,RemoveToCloseValues
#                                                  ,R_s,res
#                                                  ,directory2Data
#                                                  ,directory2Data.dir_rawdata
#                                                  ,beamgrp_vert,
#                                                  start_time[Transect],
#                                                  log_start[Transect],
#                                                  stop_time[Transect],
#                                                  lat_start[Transect],
#                                                  lat_stop[Transect],
#                                                  lon_start[Transect],
#                                                  lon_stop[Transect])
#                                
#                                
#
#                            else: 
#                                sc.savemat(directory2Data.dir_work+'/'+'Vertical_T'+str(log_start[Transect])+'.mat',mdict={'empty':0})
#    
#
#                                
#                            
#            print('     -finnished making log distance ',end='\r')
#                
#                                
#                                
#                                
#                        
#            #Make Report files               
#            MakeVerticalLuf20(CompleteListOfFiles_vertical,directory2Data, nation,cruice_id,vplatform)
#            tools.mergexml(directory2Data.dir_result+'\Vertical',directory2Data.dir_result)
    
        
        
        


if __name__ == '__main__':

    #TODO: enable so these can be selected from commandline
    
    threshold = 2
    RemoveToCloseValues = 30
    R_s = 15
    
    
    
    
    #Some user innputt that will be selectable as inputs
    recompute = False
    reconvert = False
    GO_horizontal = True
    GO_vertical = False
    Stack = True
    
    


    main(threshold, RemoveToCloseValues, R_s, recompute, reconvert,GO_horizontal, GO_vertical, Stack)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             