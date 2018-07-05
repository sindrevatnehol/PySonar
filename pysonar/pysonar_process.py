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

    
    

    
    


def doVerticalProcess(directory2Data,idx_list,LUF20_info_list,liste,RemoveToCloseValues,R_s,res, randomize= True, recompute = False): 
    
    import numpy as np
    import random, os
    from pysonar_makeverticalIDX import MakeVerticalIndex
    
    
    
    start_time=LUF20_info_list['start_time']
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
            
            
            
            if not os.path.isfile(directory2Data.dir_work+'/'+'Vertical_T'+str(log_start[Transect])+'.mat'):
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
                
                        
                        
                
    #Make Report files               
#    tools.MakeVerticalLuf20(CompleteListOfFiles_vertical,directory2Data, nation,cruice_id,vplatform)
#    tools.mergexml(directory2Data.dir_result+'\Vertical',directory2Data.dir_result)
    


    