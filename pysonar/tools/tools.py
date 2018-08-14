# -*- coding: utf-8 -*-
"""
Created on Mon May  7 10:39:10 2018

@author: sindrev
"""

import os, scipy, shutil
import numpy as np
from xml.etree import ElementTree
import glob, datetime, time
from shutil import copyfile
from xml.etree import ElementTree as ET
import Raw2NetcdfConverter
#from numba import jit 
from lxml import etree
from netCDF4 import Dataset








def printProgressBar(iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = ' '):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    
        

        
        
        
    
    
    
        
        
        
        
def WelcomeScreen(projectName): 
    
    print(projectName)
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('    >=>      /\      >=>          >=>           ')
    print('  >=>       /  \           >=>         >=>      ')
    print('           / >=>\                               ')
    print('   >=>    /      \                 >=>          ')
    print('         /   >=>  \      >=>                >=> ')
    print('        / >=>      \              >=>    >=>    ')
    print('_______/________>=>_\___________________________')
    print('                              Sindre N. Vatnehol')
    print('                    Institute of Marine Research')
    print('                                          Norway')
    print('                     Mail: sindre.vatnehol@hi.no')
    print('')

    
      
    
    
    
    
    
    
    

def UnpackBeam(BeamAmplitudeData):
    '''Unpack the beam data from the nc structure'''
    
    BeamAmplitude = np.zeros((len(BeamAmplitudeData[0]),64))
    
    for i in range(len(BeamAmplitude[0,:])):
        BeamAmplitude[:,i] = BeamAmplitudeData[i]

    return BeamAmplitude

    
    
            
            
    
    
    
class FolderStructure(object):
    '''
    Defines the folder structure used in the project.
    the PySonar
    '''
    def __init__(self, dir_cruice,equipment):
        self.dir_cruice = dir_cruice+'/'
        self.dir_acoustic_data = dir_cruice+'/ACOUSTIC/'
        self.dir_su90 = dir_cruice+'/ACOUSTIC/'+equipment+'/'
        self.dir_rawdata = dir_cruice+'/ACOUSTIC/'+equipment+'/'+equipment+'_RAWDATA'
        self.dir_originalrawdata = dir_cruice+'/ACOUSTIC/'+equipment+'/'+equipment+'_ORIGINALRAWDATA'
        self.dir_PySonar = dir_cruice+'/ACOUSTIC/PySonar/'+equipment
        self.dir_nc = dir_cruice+'/ACOUSTIC/PySonar/'+equipment+'/netcdf'
        self.dir_work = dir_cruice+'/ACOUSTIC/PySonar/'+equipment+'/WorkFiles'
        self.dir_search = dir_cruice+'/ACOUSTIC/PySonar/'+equipment+'/Search'
        self.dir_result = dir_cruice+'/ACOUSTIC/PySonar/'+equipment+'/Report'
        self.dir_src = dir_cruice+'/ACOUSTIC/PySonar/'+equipment+'/src'
        self.dir_LSSS_report = dir_cruice+'/ACOUSTIC/LSSS/REPORTS'
        self.dir_PROFOS_report = dir_cruice+'/ACOUSTIC/LSSS/REPORTS/PROFOS'
        self.dir_PROMUS_report = dir_cruice+'/ACOUSTIC/LSSS/REPORTS/PROMUS'
            
        
        
        
    




def mergexml(directory,res_dir): 
    '''Protocol to merge smaller xml files into one larger'''
    root = ET.Element("echosounder_dataset")
    
    
    xml_files = glob.glob(directory +"/*.xml")
    
    first = 0
    
    ET.SubElement(root,'report_time').text = str(datetime.datetime.fromtimestamp((time.time())).strftime('%Y-%m-%d %H:%M:%S'))
    ET.SubElement(root,'lsss_version').text = "PNMDformats v 0.1 - vertical"
    

    for xml_file in xml_files:
        data = ElementTree.parse(xml_file).getroot()
        for result in data:#.iter('echosounder_dataset'):
            if first == 0: 
                
                if result.tag == 'nation':
                    nation = result.text
                    ET.SubElement(root,'nation').text = nation
                if result.tag == 'platform':
                    platform = result.text
                    ET.SubElement(root,'platform').text = platform
                if result.tag == 'cruice_id':
                    cruice_id = result.text
                    ET.SubElement(root,'cruice_id').text = cruice_id

            if result.tag == 'distance_list': 
                
                
                if first == 0: 
                    distance_list = ET.SubElement(root,'distance_list')
                    
                    
                for dist in result: 
                    
                    distance = ET.SubElement(distance_list,'distance')
                    
                    distance.set('log_start',dist.attrib['log_start'])
                    distance.set('start_time',dist.attrib['start_time'])
                    
                    for var in dist: 
                        if var.tag == 'integrator_dist':
                            ET.SubElement(distance,'integrator_dist').text = var.text
                        if var.tag == 'pel_ch_thickness':
                            ET.SubElement(distance,'pel_ch_thickness').text = var.text
                        if var.tag == 'include_estimate':
                            ET.SubElement(distance,'include_estimate').text = var.text
                        if var.tag == 'lat_start':
                            ET.SubElement(distance,'lat_start').text = var.text
                        if var.tag == 'lon_start':
                            ET.SubElement(distance,'lon_start').text = var.text
                        if var.tag == 'lat_stop':
                            ET.SubElement(distance,'lat_stop').text = var.text
                        if var.tag == 'lon_stop':
                            ET.SubElement(distance,'lon_stop').text = var.text
                        if var.tag == 'stop_time':
                            ET.SubElement(distance,'stop_time').text = var.text
                        if var.tag == 'frequency':
                            
                            frequency = ET.SubElement(distance,'frequency')
                            frequency.set('freq',var.attrib['freq'] )
                            frequency.set('transceiver',var.attrib['transceiver'])
                            
                            for freq in var: 
                                if freq.tag == 'quality': 
                                    ET.SubElement(frequency,'quality').text = freq.text
                                if freq.tag == 'bubble_corr': 
                                    ET.SubElement(frequency,'bubble_corr').text = freq.text
                                if freq.tag == 'threshold': 
                                    ET.SubElement(frequency,'threshold').text = '-65'
                                if freq.tag == 'num_pel_ch': 
                                    ET.SubElement(frequency,'num_pel_ch').text = freq.text
                                if freq.tag == 'upper_interpret_depth': 
                                    ET.SubElement(frequency,'upper_interpret_depth').text = freq.text
                                if freq.tag == 'lower_interpret_depth': 
                                    ET.SubElement(frequency,'lower_interpret_depth').text = freq.text
                                if freq.tag == 'upper_interpret_depth': 
                                    ET.SubElement(frequency,'upper_interpret_depth').text = freq.text
                                if freq.tag == 'lower_integrator_depth': 
                                    ET.SubElement(frequency,'lower_integrator_depth').text = freq.text
                                if freq.tag == 'ch_type': 
                                    ch_type = ET.SubElement(frequency,'ch_type')
                                    ch_type.set('type','P')
                                    for chn in freq: 
                                        sa_by = ET.SubElement(ch_type,'sa_by_acocat')
                                        sa_by.set('acocat', chn.attrib['acocat'])
                                        for sa in chn: 
                                            
                                            sa_value = ET.SubElement(sa_by,'sa')
                                            sa_value.set('ch',sa.attrib['ch'])
                                            sa_value.text= sa.text
        first = 1
          
        
    indent(root)
    tree = ET.ElementTree(root)
    tree.write(res_dir+'/ListUserFile20_SU90_vertical.xml')
                
    
    
    
    
    
    
        
def MakeNewFolders(directory2Data):
    '''Process to make new folder if they do not exist'''
    
    
    if not os.path.exists(directory2Data.dir_cruice):
        os.makedirs(directory2Data.dir_cruice)
    if not os.path.exists(directory2Data.dir_acoustic_data):
        os.makedirs(directory2Data.dir_acoustic_data)
    if not os.path.exists(directory2Data.dir_su90):
        os.makedirs(directory2Data.dir_su90)
    if not os.path.exists(directory2Data.dir_work):
        os.makedirs(directory2Data.dir_work)
    if not os.path.exists(directory2Data.dir_PySonar):
        os.makedirs(directory2Data.dir_PySonar)
    if not os.path.exists(directory2Data.dir_search):
        os.makedirs(directory2Data.dir_search)
    if not os.path.exists(directory2Data.dir_result):
        os.makedirs(directory2Data.dir_result)
    if not os.path.exists(directory2Data.dir_result+'\Horizontal'):
        os.makedirs(directory2Data.dir_result+'\Horizontal')
    if not os.path.exists(directory2Data.dir_result+'\Vertical'):
        os.makedirs(directory2Data.dir_result+'\Vertical')
    if not os.path.exists(directory2Data.dir_nc):
        os.makedirs(directory2Data.dir_nc)
    if not os.path.exists(directory2Data.dir_originalrawdata):
        os.makedirs(directory2Data.dir_originalrawdata)
    if not os.path.exists(directory2Data.dir_rawdata):
        os.makedirs(directory2Data.dir_rawdata)
    if not os.path.exists(directory2Data.dir_src):
        os.makedirs(directory2Data.dir_src)
    if not os.path.exists(directory2Data.dir_work):
        os.makedirs(directory2Data.dir_work)
    if not os.path.exists(directory2Data.dir_LSSS_report):
        os.makedirs(directory2Data.dir_LSSS_report)
    if not os.path.exists(directory2Data.dir_PROFOS_report):
        os.makedirs(directory2Data.dir_PROFOS_report)
        
               

    




        
def indent(elem, level=0):
  '''
  Description: 
       Make the xml file more readable
  '''
  i = "\n" + level*"  "
  if len(elem):
    if not elem.text or not elem.text.strip():
      elem.text = i + "  "
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
    for elem in elem:
      indent(elem, level+1)
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
  else:
    if level and (not elem.tail or not elem.tail.strip()):
      elem.tail = i
      
      
      
      
      
      
      
      
      
    
def addFrequencyLevelInfo(distance,freq,tran,NumCH,TH): 
    '''Add information in the frequency level in luf20 file'''
    
    
    #Make new frequency with attributes
    frequency = ET.SubElement(distance,'frequency')
    frequency.set('freq',str(int(freq))  )
    frequency.set('transceiver',str(int(tran)))
    ET.SubElement(frequency,'quality').text = str(2)  
    ET.SubElement(frequency,'bubble_corr').text = str(0)  
    ET.SubElement(frequency,'threshold').text = str(TH)
    
    ET.SubElement(frequency,'num_pel_ch').text = str(NumCH)  
    
    ET.SubElement(frequency,'upper_interpret_depth').text = str(0)
    ET.SubElement(frequency,'lower_interpret_depth').text = str(0)
    ET.SubElement(frequency,'upper_interpret_depth').text = str(0)
    ET.SubElement(frequency,'lower_integrator_depth').text = str(0)

    

    ch_type = ET.SubElement(frequency,'ch_type')
    ch_type.set('type','P')
    sa_by_acocat = ET.SubElement(ch_type,'sa_by_acocat')
    sa_by_acocat.set('acocat',str(12))
    
    return sa_by_acocat
    
    
    
    
    
    
def GetTransectTimes(filename,code):

    
    doc2 = etree.parse(filename)
    distance_list = doc2.find('distance_list')
#    old_stop = ''
#    Start = 0
#    End = []
#    transect_IDX = []
#    transect_i = 1
#    new_start = False

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
            printProgressBar(i+1,len(ListOfFiles),prefix = 'Stacking: ', suffix = 'Completed', length = 50)


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
#                    timeID = np.where(abs(beam_data.variables['ping_time'][ii]-fid_nc.groups['Platform'].variables['time1'][:]/100)==
#                                          np.nanmin(abs(beam_data.variables['ping_time'][ii]-fid_nc.groups['Platform'].variables['time1'][:]/100)))
                    
                    
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
#        fd5.close()
        
        
    return CompleteListOfFiles, D



def GetTransectTimesIDX(filename,code):

    
    doc2 = etree.parse(filename)
    distance_list = doc2.find('distance_list')
    old_stop = ''
    Start = 0
    End = []
    transect_IDX = []
    transect_i = 1
#    new_start = False

    start_time = []
#    log_start = []
    stop_time = []
#    lat_start = []
#    lon_start = []
#    lat_stop = []
#    lon_stop = []
    
    
    for i in distance_list:
        
        start_time = i.get('start_time').replace(' ','T').replace('-','').replace(':','')
#        log_start = i.get('log_start')
        stop_time = i.find('stop_time').text.replace(' ','T').replace('-','').replace(':','')
#        lat_start = i.find('lat_start').text
#        lon_start = i.find('lon_start').text
#        lat_stop = i.find('lat_stop').text
#        lon_stop = i.find('lon_start').text
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
      
    
    
      
      
def addLogDistanceInfo(distance_list,LatStart,LonStart,TimeStart,LatStop,LonStop,TimeStop,
                       integrator_dist,pel_ch_thickness,LogStart): 
    '''Add information in the log distance level in luf20 xml'''
    
    
    #Get the unix time and convert it to time string
    correct_starttime =TimeStart# time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(TimeStart))
    correct_stoptime = TimeStop #time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(TimeStop))
    
    
    
    #Write new log distance with its attributes
    distance = ET.SubElement(distance_list,'distance')
    distance.set('log_start',str(LogStart))
    distance.set('start_time',str(correct_starttime))

    
    
    #Add variables to the log distance
    ET.SubElement(distance,'integrator_dist').text = str(integrator_dist)
    ET.SubElement(distance,'pel_ch_thickness').text = str(pel_ch_thickness)
    ET.SubElement(distance,'include_estimate').text = '1'
    ET.SubElement(distance,'lat_start').text = str(LatStart)
    ET.SubElement(distance,'lon_start').text = str(LonStart)
    ET.SubElement(distance,'lat_stop').text = str(LatStop)
    ET.SubElement(distance,'lon_stop').text = str(LonStop)
    ET.SubElement(distance,'stop_time').text = str(correct_stoptime)
    return distance
    
    
    
    
    
    
    


def DataConverter(CruiceIndex,WorkDirectory,current_dir,maxPingInFile,
                  MaxNumberOfFilesInNC,directory2Data, reconvert) :
    ''' 
    Protocoll to convert the data
    '''

    
    #Vessel and platform information
    #This information is taken from the recepy. 
    #It is important that the platform name and 
    vessel_name = CruiceIndex.getAttribute('vessel')
    platform_type = CruiceIndex.getAttribute('platform_type')
    platform_code = CruiceIndex.getAttribute('platform_code')
    
    
    
    
    #Loop through each fishery sonar type
    #Denne vil endres til å bli mer generell
    for i in ['SU90','SX90','SH90']:         
            
        
        os.chdir(current_dir)
        Raw2NetcdfConverter.convert.Raw2NetcdfConverter(directory2Data.dir_originalrawdata,vessel_name,platform_type,
                    platform_code,maxPingInFile,MaxNumberOfFilesInNC,
                    directory2Data.dir_rawdata,reconvert)
            
            
        
       
        
        
        
        
        
            
            
        
def TransectTimeIDX(CruiceIndex): 
    TimeIDX =[]
    TransCode = CruiceIndex.getElementsByTagName("transect")
    for t in TransCode: 
        if TimeIDX == []: 
            TimeIDX = np.array([CruiceIndex.getAttribute("code")+'_'+t.getAttribute('num'),
                                t.getAttribute('startTime'),t.getAttribute('endTime')
                                ])[:,np.newaxis].T
        else: 
            TimeIDX  = np.vstack((TimeIDX,np.array([CruiceIndex.getAttribute("code")+'_'+t.getAttribute('num'),
                                t.getAttribute('startTime'),t.getAttribute('endTime')
                                ])[:,np.newaxis].T))      
    return TimeIDX
    
    
    
    
    #@jit
def GetShortListOfFiles(CompleteListOfFiles,startTime,endTime):
    ShortListOfFiles = []

    for i in range(len(CompleteListOfFiles[:][:])):
        
        printProgressBar(i + 1, len(CompleteListOfFiles[:,0]), prefix = 'Make short list:', suffix = 'Completed       ', length = 50)
        
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
    
    
    
def SecondsBetweenTransect(startTime,unixtime):

    fulldate = datetime.datetime.strptime('1601-01-01 00:00:00.000',"%Y-%m-%d %H:%M:%S.%f")

    starten =  datetime.datetime.strptime(startTime,"%Y%m%d%H%M%S")

    fulldate = (starten-fulldate).total_seconds() #datetime.timedelta(milliseconds=int(startTime))

    seconds = fulldate-unixtime/10000000

    return seconds

    
        
        
        
def OrginizeData(CruiceIndex,WorkDirectory,OS): 
    '''Protocol to organize the data on the local server.'''
    
    
    #Get the vessel, cruice and equipment information
    vesselCode = CruiceIndex.getAttribute("vesselCode")
    cruiceCode = CruiceIndex.getAttribute("code")
    equipment =  CruiceIndex.getAttribute("Equipment")
    
    
    #Loop through each fishery sonar type
    for i in ['SU90','SX90','SH90']: 
        
        
        if i == equipment: 
            
            #Make a correct sonar folder structure
            directory2Data =FolderStructure(WorkDirectory+'/'+cruiceCode[:4]+'/S'+cruiceCode+'_P'+vesselCode,i)       
            
            
            
            #Make the folder according to standard
            print('Make Structure               ')
            MakeNewFolders(directory2Data)
        
            
            
            #Get list of files that has not been copied to the correct structure
            print('    -Get files from server  ',end='\r')
            
            
            #Check if file exist, this is a bugfix when working localy
            if os.path.isdir(OS+CruiceIndex.getAttribute('CruicePath')) == True: 
            
                
                
                #Get list of files on server and list that has not been copied
                ListFromServer = os.listdir(OS+CruiceIndex.getAttribute('CruicePath'))
                ListOfFilesNotcopied = list(set(ListFromServer)-set(os.listdir(directory2Data.dir_originalrawdata)))

                #Go through each file that has not been copyed
                for i in np.arange(len(ListOfFilesNotcopied)): 
                    
                    
                    
                    #Print a progressbar for the user
                    printProgressBar(i, len(ListOfFilesNotcopied), prefix = 'Copy files', suffix = 'Files left: '+
                                     str(len(ListOfFilesNotcopied)-i) , decimals = 1, length = 50)
                    
                    
                    if '.raw' in ListOfFilesNotcopied[i]: 
                        #Copy the files
                        try:
                            copyfile(OS+CruiceIndex.getAttribute('CruicePath')+'/'+ListOfFilesNotcopied[i], 
                                 directory2Data.dir_originalrawdata+'/'+ListOfFilesNotcopied[i])
                        except: 
                            print('      * Bad file', end = '\r')
                    else: 
                        print('Skipped: '+ListOfFilesNotcopied[i])
                        
                        
                        
    return(directory2Data)
   

    
    

def GainAdjustment(pulslength,frequency,inngain):
    '''Find the correct gain to be used when scrutinising the sonar data. 
    #NB: the pulselength must be in ms and the frequency in kHz
    #
    #The correction for the pulse length is similar as to the one given in (Macauley et al. 2016)
    #The correction for frequency is different due to new avaliable information. 
    '''
    
    a = 6.99915703
    b = -6.24304318
    PulslengthGain = a +b/pulslength
    
    
    
    a = 48.4116
    b = -10.9486
    c = 0.6446
    d = -0.0108
    FrequencyGain = a+b*frequency +c*frequency**2+d*frequency**3
    
    outgain = FrequencyGain+PulslengthGain+inngain
    
    return outgain, FrequencyGain, PulslengthGain
    
    
    
    
    
    
def ApplyTVG(AmplitudeData,soundvelocity,sampleinterval,transmitpower,
             absorptioncoefficient,frequency,pulslength,gain,
             equivalentbeamangle,sacorrection,tiltAngle): 
    
    
    
    #There is a time delay in the sonar. This corresponds to 3 m
    rangeCorr = 3; 
    
    
    effective_pulselength = 0.7*pulslength
    
    #The distance between each sample in meter
    samplespace =soundvelocity*sampleinterval/2
    
    
    #Declare the range vector to be outputed
    Range = np.arange(0,AmplitudeData[:].shape[0])*samplespace-rangeCorr
                      
                      
                      
    #Remove bad values  
    Range[np.where(Range<=0)] = 1E-200
               

    #Make the TVG function (20 log) as a vector
    tvg20 = np.repeat(20*np.log10(Range[:, np.newaxis]) + 2*absorptioncoefficient*Range[:, np.newaxis], AmplitudeData[:].shape[1], axis=1)
                         
                      
    tvg20[np.where(tvg20<0)] = 0
    
    #Correction function for the loss off energy when steering the beam
    tiltcorr =40*np.log10(abs(np.cos(np.deg2rad(-tiltAngle))))
    
    #compute the wavelength of the signal
    wavelength = soundvelocity/frequency

    
    
    #compute the correction function used when computing the SV
    svconstSV =10*np.log10(transmitpower*wavelength**2*soundvelocity*effective_pulselength/(32*np.pi**2))+tiltcorr+gain+10*np.log10(equivalentbeamangle)+sacorrection
    
        
    #Compute the volume backscattering coefficient
    sv= AmplitudeData+tvg20-np.repeat(svconstSV[:,np.newaxis],AmplitudeData.shape[0],axis=1).transpose()

    #Return the SV the range and the TS from the data.
    return sv, Range; 


    
    
    
#@jit
def ConvertToechogram(Wdist,sv_mat):
    '''
    Purpose
    Search through the data in a banana shape to find the sA
    Not likely to be changed unless the methodology changes
    
    
    It requieres large memory
    '''
    
    
    Wdist = Wdist.toarray()
    
    #Sett weight to 1 for all relevant distances
    Wdist[np.where(Wdist>0)] = 1
    

    #find the index of the files without 0, This for speed improvemment
    idxWeight = np.where(Wdist>0)
                
    #Find the sa value in linear domain. 
    sA= np.nansum(sv_mat[idxWeight])
    
    return sA
    
    
#@jit
def ConvertToechogram2(Wdist1,Wdist2,Wdist3,Wdist4,sv_mat,rangen):
    sA = np.array([])
    sA2 = np.array([])
    sA3 = np.array([])
    sA4 = np.array([])
    
    for i in rangen: 
        idxWeight = np.where(Wdist1[i].toarray()>0)
        idxWeight2 = np.where(Wdist2[i].toarray()>0)
        idxWeight3 = np.where(Wdist3[i].toarray()>0)
        idxWeight4 = np.where(Wdist4[i].toarray()>0)
        
        #Find the sa value in linear domain. 
        sA= np.hstack((sA,np.nansum(sv_mat[idxWeight])))
        sA2= np.hstack((sA2,np.nansum(sv_mat[idxWeight2])))
        sA3= np.hstack((sA3,np.nansum(sv_mat[idxWeight3])))
        sA4= np.hstack((sA4,np.nansum(sv_mat[idxWeight4])))
    
    return sA, sA2, sA3, sA4
    
    
    
    
    
    
    
#@jit
def GetDistanceMatrix(DistanceMatrix,RangeMatrix,BeamDirectionMatrix,svMatrix,theta_tilt,BananaTool):   
    '''
    This function computes the distance matrix. 
    For values outside the tubed school are sett to 0. 
    
    It requieres large memory
    
    '''
    
    
    
    #make the distance matrixes as a dictionary for easyer access
    Wdist_port = dict()
    Wdist_stb= dict()
    
    
    
    
    #convert angles into radians
    theta_tilt = theta_tilt * np.pi/180
    
    
    
    
    #Distance traveled by the vessel, centered at origo
    '''This will be made more general in the future. 
    E.g. from number of sampels in time and range, mean ping time, and vessel speed, make both
    the travelDistance and the x/y/z pixel matrixes
    '''

    travelDistance = DistanceMatrix[0,0,:]-np.nanmin(DistanceMatrix)
    
    travelDistance = np.max(travelDistance)/len(travelDistance)*np.arange(len(travelDistance))
    
    travelDistance=travelDistance-np.max(travelDistance)/2

    #Compute the location of all points in cartesian coordinates
    x_pixel = RangeMatrix*np.sin(theta_tilt)*np.cos(BeamDirectionMatrix*np.pi/180)
    
    x_pixel = x_pixel+DistanceMatrix
    
    x_pixel = x_pixel-np.min(x_pixel)
    
    x_pixel = x_pixel-np.max(x_pixel)/2    
    
    y_pixel = RangeMatrix*np.sin(theta_tilt)*np.sin(BeamDirectionMatrix*np.pi/180)
    
    z_pixel = RangeMatrix*np.cos(theta_tilt)
    
    
    correction = np.repeat(np.repeat((RangeMatrix[:,0,0]*np.sin(2*np.pi/64))[:,np.newaxis],64,axis=1)[:,:,np.newaxis],len(RangeMatrix[0,0,:]),axis=2)
    
    index=0
    #Run through all pings and find the distance of each pixel 
    for r_ind_stbex in range(int(BananaTool[1]),int(BananaTool[2]),int(BananaTool[3])):
        printProgressBar(r_ind_stbex+1,int(BananaTool[2]), prefix = 'Make Distance:', suffix = 'Completed     ', length = 50)
        
        
        #Get the minimum range to sonar trajectory
        r_0 = RangeMatrix[r_ind_stbex,0,0]



        #Model for school trajectory, assuming a stationary school      
        x_traj = np.zeros(travelDistance.shape) #travelDistance          
        y_traj_port = r_0*np.ones(x_traj.shape) *np.sin(theta_tilt)    
        y_traj_stb = -r_0*np.ones(x_traj.shape) *np.sin(theta_tilt)
        z_traj = np.sqrt(x_traj**2+y_traj_port**2)/np.tan(theta_tilt)
        
        
               
        #Empty range matrix for bookkeeping
        R_port = np.zeros(svMatrix.shape)
        R_stb = np.zeros(svMatrix.shape)
        R2_port = np.zeros(svMatrix.shape)
        R2_stb = np.zeros(svMatrix.shape)
        
        
        
        
        #Compute the range of all pixels from the ideal line        
        for i in range(len(x_traj)): 
            R_port[:,:,i] = np.sqrt((x_pixel[:,:,i]-x_traj[i])**2+(y_pixel[:,:,i]
                                -y_traj_port[i])**2+(z_pixel[:,:,i]-z_traj[i])**2)
            R_stb[:,:,i] = np.sqrt((x_pixel[:,:,i]-x_traj[i])**2+(y_pixel[:,:,i]
                                -y_traj_stb[i])**2+(z_pixel[:,:,i]-z_traj[i])**2)
        
            
           
        indeksen = np.repeat(np.repeat(np.reshape(np.amin(np.amin(R_port,axis=1),
                              axis=0),(1,1,-1)),R_port.shape[0],axis=0),64,axis=1)
        indeksen_stb = np.repeat(np.repeat(np.reshape(np.amin(np.amin(R_stb,axis=1),
                              axis=0),(1,1,-1)),R_stb.shape[0],axis=0),64,axis=1)
            
        
        R2_port[R_port<=(indeksen+correction+int(BananaTool[0]))]=1
        R2_stb[R_stb<=(indeksen_stb+correction+int(BananaTool[0]))]=1
        
               
               
        #Save the distance matrixe to memory/alternatively it can be stoored in a seperate file
        Wdist_port[index] =  scipy.sparse.csr_matrix(np.reshape(R2_port,(-1,1)))
        Wdist_stb[index] =  scipy.sparse.csr_matrix(np.reshape(R2_stb,(-1,1)))
        
        index+=1
        

    return Wdist_port,Wdist_stb

    
    
    
    
    
    
    
    
    
    
from math import radians, cos, sin, asin, sqrt

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371*1000 # Radius of earth in meters. Use 3956 for miles
    return c * r
    
    

    
def ComputeDistance(travelDist,lat,lon): 
    '''
    This function computes the distance traveled within each transect. 
    For distance the hamorsine function is used. 
    
    Parts of the code is awkward due to bugs in the files. Will be changed later. 
    '''
    
    #We need some points to evaluate the distanc traveled. 
    if len(lat)==1: 
        travelDist = 0
    else: 
        travelDist = np.hstack((travelDist,haversine(lon[-1], lat[-1], lon[-2], lat[-2])))
    DistanceTraveled = np.cumsum(travelDist)
#        travelDist=np.hstack((travelDist,np.nan))
        
        
        
#    elif len(lat)>=20:
#        #find the delta in longitude and latitude
#        delta_lat = lat[-1]-lat[-2]
#        delta_lon = lon[-1]-lon[-2]
#
#
#
#        #Computing the distance using hammersine fucntion
#        a = (np.sin(delta_lat/2*np.pi/180))**2 + np.cos(lat[-2]*np.pi/180) * np.cos(lat[-1]*np.pi/180) * (np.sin(delta_lon/2*np.pi/180)**2)
#        c = 2 * np.arctan2(np.sqrt(a),np.sqrt(1-a))
#        
#        
#        #Compute the distance to meter and stack the variable
#        travelDist = np.hstack((travelDist,6371000*c))
#        
#        
#        
#        #Correct for some bug in the software. 
#        #if the travelDistance is not a number compute the meidan distance
#        if len(travelDist)>=4:
#            if travelDist[-2]==np.nan:
#                travelDist[-2]== (travelDist[-1]+travelDist[-3])/2
#    
#
#
#    #travelDist = distance between each point
#    #DistanceTraveled is the total length traveled for each ping
#    DistanceTraveled = travelDist
#
#    DistanceTraveled[np.where(np.isnan(DistanceTraveled))] = 0
#    DistanceTraveled = np.cumsum(DistanceTraveled)
        
    
    return DistanceTraveled, travelDist;


    
    
    
class GetVariablesFromNC(object):
    
    def __init__(self,fileID,beamgrp,Files,fileIDX):
        
        
        if type(beamgrp) == np.str: 
            bmgrp = beamgrp
        else:
            bmgrp = beamgrp[fileIDX]
        
        fIDX = fileIDX
           
#        if Files == False:
#            fIDX = int(fileIDX)
#        else: 
#            fIDX = int(Files[fileIDX,2])
        
        #Get beam sonar configuration info   
        try: 
            pingtime = fileID.groups['Sonar'].groups[bmgrp].variables['ping_time'][fIDX]
            self.time = pingtime
            self.frequency = fileID.groups['Sonar'].groups[bmgrp].variables['transmit_frequency_start'][fIDX]
            self.transmitpower = fileID.groups['Sonar'].groups[bmgrp].variables['transmit_power'][fIDX]
            self.pulslength = fileID.groups['Sonar'].groups[bmgrp].variables['transmit_duration_nominal'][fIDX]
            self.gaintx = fileID.groups['Sonar'].groups[bmgrp].variables['transducer_gain'][fIDX]
            self.gainrx = fileID.groups['Sonar'].groups[bmgrp].variables['receiver_sensitivity'][fIDX]
            self.sampleinterval = fileID.groups['Sonar'].groups[bmgrp].variables['sample_interval'][fIDX]
            self.equivalentbeamangle=fileID.groups['Sonar'].groups[bmgrp].variables['equivalent_beam_angle'][fIDX,:]
            self.sacorrection = 0
    
            #Get environment data
            self.soundvelocity = fileID.groups['Environment'].variables['sound_speed_indicative'][:]
            self.absorptioncoefficient = fileID.groups['Environment'].variables['absorption_indicative'][:]
    
          
            
            #Get beam configuration
            beam_direction_x=fileID.groups['Sonar'].groups[bmgrp].variables['beam_direction_x'][fIDX,:]
            beam_direction_y=fileID.groups['Sonar'].groups[bmgrp].variables['beam_direction_y'][fIDX,:]
            beam_direction_z=fileID.groups['Sonar'].groups[bmgrp].variables['beam_direction_z'][fIDX,:]


            self.dirx = np.arcsin(beam_direction_z)/np.pi*180
            self.diry = np.arctan2(beam_direction_y,beam_direction_x)*180/np.pi
            
    
            #Unpack and write beam data  
            BeamIM = fileID.groups['Sonar'].groups[bmgrp].variables['backscatter_i']
            BeamReal = fileID.groups['Sonar'].groups[bmgrp].variables['backscatter_r']
            BeamAmplitudeDataIM=UnpackBeam(BeamIM[fIDX,:])
            BeamAmplitudeDataReal=UnpackBeam(BeamReal[fIDX,:])
            
            self.BeamAmplitudeData =(BeamAmplitudeDataIM**2 + BeamAmplitudeDataReal**2)
            
            
            NMEA_time= fileID.groups['Platform'].variables[fileID.groups['Platform'].variables['longitude'].dimensions[0]][:]/100
            Latitude = fileID.groups['Platform'].variables['latitude'][:]
            Longitude = fileID.groups['Platform'].variables['longitude'][:]
            
            NMEA_idx = np.where(abs(NMEA_time-pingtime)==np.min(abs(NMEA_time-pingtime)))
            #Get NMEA data
            self.Longitude = Longitude[NMEA_idx]
            self.Latitude = Latitude[NMEA_idx]
            self.NMEA_time = NMEA_time[NMEA_idx]
            

        except IndexError: 
            k=1
        except KeyError: 
            print('Bad variable innput in file')



