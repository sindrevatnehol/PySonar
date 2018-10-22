# -*- coding: utf-8 -*-
"""
Created on Wed Jul  4 08:23:00 2018

@author: sindrev
"""


    
def SecondsBetweenTransect(startTime,unixtime):
    import datetime

    fulldate = datetime.datetime.strptime('1601-01-01 00:00:00.000',"%Y-%m-%d %H:%M:%S.%f")

    starten =  datetime.datetime.strptime(startTime,"%Y%m%d%H%M%S")

    fulldate = (starten-fulldate).total_seconds() #datetime.timedelta(milliseconds=int(startTime))

    seconds = fulldate-unixtime/10000000

    return seconds

    
    
    
    
    
def GetTransectTimesIDX(filename,code):
    from lxml import etree
    import numpy as np
    
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
    
def GetTransectTimes(filename,code):

    from lxml import etree
    import numpy as np
    
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
    
    
    
def GetLuf20Info(directory2Data,code):
    import os
    import shutil
    
    print('Startet')
    
    
    for dirpath,_,filenames in os.walk(directory2Data.dir_src + '/EKLUF20/'):
        for f in filenames:
            filename =  os.path.abspath(os.path.join(dirpath, f))
            TimeIDX = GetTransectTimesIDX(filename,str(code))
            start_time,log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop = GetTransectTimes(filename,str(code))
            
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
                    TimeIDX = GetTransectTimesIDX(filename,str(code))
                    start_time,log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop = GetTransectTimes(filename,str(code))
                
        
    return  start_time,log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop, TimeIDX
    
    
    

    
    
    
    
def doHorizontalProcess(directory2Data,CruiceIndex): 
    import numpy as np
    import random, os
    from tools import tools
    
    
    start_time,log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop, TimeIDX = GetLuf20Info(directory2Data,CruiceIndex)
    
    
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
            ShortListOfFiles_horizontal = tools.GetShortListOfFiles(CompleteListOfFiles_horizontal,
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
            ShortListOfFiles_horizontal = tools.GetShortListOfFiles(CompleteListOfFiles_horizontal,
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
            
            
    
    
    
    
    
            
    
    
    


def doVerticalProcess(directory2Data,idx_list,LUF20_info_list,liste,RemoveToCloseValues,R_s,res, randomize= True, recompute = False): 
    
    import numpy as np
    import random, os
    from pysonar_makeverticalIDX import MakeVerticalIndex
    from VerticalLuf20 import MakeVerticalLuf20
    
    
    
    
    start_time= LUF20_info_list['start_time']
    log_start = LUF20_info_list['log_start']
    stop_time = LUF20_info_list['stop_time']
    lat_start = LUF20_info_list['lat_start']
    lat_stop = LUF20_info_list['lat_stop']
    lon_start = LUF20_info_list['lon_start']
    lon_stop = LUF20_info_list['lon_stop']
    
    
    #Get transact list and randomise it
    NumTransect = np.arange(len(start_time))
    
    
    
    if randomize == True: 
        random.shuffle(NumTransect)
    
    
    
#    idx_list['ping_time'][np.where(idx_list['ping_time'])==0] = np.nan


    
                
        
    #Loop through each transect
    for Transect in NumTransect:
        
        
        
    
        
        start = SecondsBetweenTransect(start_time[Transect].replace('T',''),idx_list['ping_time'])
        end = SecondsBetweenTransect(stop_time[Transect].replace('T',''),idx_list['ping_time'])
        


        try: 
            short_idx = np.arange(np.min(np.where(start<=0)),np.max(np.where(end>=0)))
        except: 
            short_idx = []



#short_idx

        if len(short_idx)>0: 
            print('    -Start on log distance: '+ log_start[Transect],end='\r')
            
            
            
            if not os.path.isfile(directory2Data.dir_verticalwork+'/'+'Vertical_T'+str(log_start[Transect])+'.mat'):
                MakeVerticalIndex(idx_list['FileList'][short_idx]
                                              ,RemoveToCloseValues
                                              ,R_s,res
                                              ,directory2Data
                                              ,directory2Data.dir_rawdata
                                              ,0,
                                              start_time[Transect],
                                              log_start[Transect],
                                              stop_time[Transect],
                                              lat_start[Transect],
                                              lat_stop[Transect],
                                              lon_start[Transect],
                                              lon_stop[Transect],
                                              idx_list['IDX'][short_idx])
            elif recompute == True: 
                MakeVerticalIndex(idx_list['FileList'][short_idx]
                                              ,RemoveToCloseValues
                                              ,R_s,res
                                              ,directory2Data
                                              ,directory2Data.dir_rawdata
                                              ,0,
                                              start_time[Transect],
                                              log_start[Transect],
                                              stop_time[Transect],
                                              lat_start[Transect],
                                              lat_stop[Transect],
                                              lon_start[Transect],
                                              lon_stop[Transect],
                                              idx_list['IDX'][short_idx])

            
#                print(np.where(abs(end) == np.min(abs(end))))
#                print(np.where(start<=0))
#                print(np.where(end>=0))
                
                        
                        
    