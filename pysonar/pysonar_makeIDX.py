# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 19:06:45 2018

@author: sindrev
"""

 


def getLuf20Info(directory2Data,code):
    import os
    import numpy as np
    from lxml import etree
    
    
    
    
    def GetTransectTimes(filename,code):

        
        #Bookkeeping
        
        start_time= []
        log_start= []
        stop_time= []
        lat_start= []
        lon_start = []
        lat_stop= []
        lon_stop= []
        End= []
        transect_IDX = []
        transect_i = 1
        old_stop = ''
        Start = 0
        
        
        #Load LUF20 file
        doc2 = etree.parse(filename)
        distance_list = doc2.find('distance_list')
        
        
        for i in distance_list:
            
            
            start_time0 = i.get('start_time').replace(' ','T').replace('-','').replace(':','')
            stop_time0 = i.find('stop_time').text.replace(' ','T').replace('-','').replace(':','')
            start_time = np.hstack((start_time,i.get('start_time').replace(' ','T').replace('-','').replace(':','')))
            log_start = np.hstack((log_start,i.get('log_start')))
            stop_time = np.hstack((stop_time,i.find('stop_time').text.replace(' ','T').replace('-','').replace(':','')))
            lat_start = np.hstack((lat_start,i.find('lat_start').text))
            lon_start = np.hstack((lon_start,i.find('lon_start').text))
            lat_stop = np.hstack((lat_stop,i.find('lat_stop').text))
            lon_stop = np.hstack((lon_stop,i.find('lon_stop').text))
            
            if transect_i <10: 
                trnsID = '00'+str(transect_i)
            elif transect_i <100: 
                trnsID = '0'+str(transect_i)
            else:
                trnsID = str(transect_i)
            
                
            if old_stop!=start_time0: 
                End = np.hstack((End,old_stop))
                Start = np.hstack((Start,start_time0))
                transect_IDX = np.hstack((transect_IDX,code + '_'+trnsID))
                transect_i = transect_i+1
    
                
            if Start == 0:
                Start = start_time0
                transect_IDX = np.hstack((transect_IDX,code + '_'+trnsID))
                transect_i = transect_i+1
        
                
            old_stop = stop_time0
    
        #add last time
        End = np.hstack((End,stop_time0))
        
        
        TimeIDX  = np.vstack((transect_IDX.T,Start[1:].T,End[1:].T)).T
                              
                              
        return start_time,log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop, TimeIDX


    
        
    
    
    
    for dirpath,_,filenames in os.walk(directory2Data.dir_src + '/EKLUF20/'):
        for f in filenames:
            filename =  os.path.abspath(os.path.join(dirpath, f))
            start_time,log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop,TimeIDX = GetTransectTimes(filename,code)

            
            
            
        
                    
    liste = {}
    liste['start_time'] = start_time
    liste['log_start'] = log_start
    liste['stop_time'] = stop_time             
    liste['lat_start'] = lat_start
    liste['lat_stop'] = lat_stop
    liste['lon_start'] = lon_start
    liste['lon_stop'] = lon_stop
    liste['TimeIDX'] = TimeIDX
    return  liste
    
    
    
    
    
    
    

def readIDX(directory2Data,beam_mode):
    
    from netCDF4 import Dataset
    import numpy as np
    
    print('Read IDX')
    f = Dataset(directory2Data.dir_src+'/IDX_'+beam_mode+'.nc','r')
    
    
    Time = f.groups['IDXinfo'].variables['ping_time'][:]

        
        
    

    liste = {}

    liste['ping_time'] = f.groups['IDXinfo'].variables['ping_time'][:]
    liste['FileList'] = f.groups['IDXinfo'].variables['FileList'][:]
    liste['IDX'] = f.groups['IDXinfo'].variables['IDX'][:]
    liste['BeamIDX'] = f.groups['IDXinfo'].variables['BeamIDX'][:]
    
    f.close()
    
    print('IDX is read')
    return liste

    
    
    
    
    
    
def makeIDX(directory2Data,beam_mode): 
    
    import numpy as np
    import os
    from tools import tools
    from netCDF4 import Dataset


    
    
    
    #Sort all the .nc files
    print('Get list of files')
    ListOfFiles = np.sort(os.listdir(directory2Data.dir_rawdata))
    print('Got list of files')

    
    
    
    
    
    ping_time = []
    FileList = []
    IDX = []
    BeamIDX = []

    tot_index = 0

                
    f = Dataset(directory2Data.dir_src+'/IDX_'+beam_mode+'.nc','w', format='NETCDF4')
    
    tempgrp = f.createGroup('IDXinfo')
    tempgrp.createDimension('dim', None)
    nchar_t = tempgrp.createVLType(np.str,'nchar')

    ping = tempgrp.createVariable('ping_time', np.int64, ('dim',), chunksizes = (1024,))
    FileL = tempgrp.createVariable('FileList',nchar_t,('dim',), chunksizes = (512,))
    beamidx = tempgrp.createVariable('BeamIDX',nchar_t,('dim'), chunksizes = (512,))
    idx = tempgrp.createVariable('IDX',int,'dim',chunksizes = (512,))

    
    
    #This helps to prevent opening and closing the same file multiple times
    oldFileName = ''

    #Loop through each file
    for i in range(len(ListOfFiles)):

        
        
        #Print progress
        tools.printProgressBar(i+1,len(ListOfFiles),prefix = 'Stacking: ', suffix = 'Completed', length = 50)

        
        

#        try: 
            #If this is a new file, close the last and open this
        if ListOfFiles[i]!=oldFileName:
            try:
                fid_nc.close()
            except:
                dummy=1

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
        for ii in range(len(beam_data.variables['ping_time'][:])):
            ping_time = np.hstack((ping_time,beam_data.variables['ping_time'][ii]))
            
#                FileList = np.hstack((FileList,ListOfFiles[i]))
#                IDX = np.hstack((IDX,ii))
#                BeamIDX = np.hstack((BeamIDX,beamgrp))

            
            FileL[tot_index] = str(ListOfFiles[i])#FileList
            idx[tot_index] = ii
            beamidx[tot_index] = str(beamgrp)

            ping_out = int(beam_data.variables['ping_time'][ii])
            ping[tot_index] = ping_out#ping_time
            

    
            tot_index = tot_index +1
                    
        
#        except: 
#            dummy = 1

        
    f.close()
