# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 14:53:31 2018

@author: sindrev
"""



        
        


    
def getListOfSurveys(directory):
    '''
    
    Description: This function convert the xml instruction file 
    to a list function. Alternatively, only the list can be made. 
    '''
    
    
    #import packages
    import xml.dom.minidom as minidom
    
    
    
    
    #Parse the structure file, and get information of each survey
    doc = minidom.parse(directory+'/src/DataStructure.xml')



    
    
    #Prepare to go through all cruices
    CruiceCode = doc.getElementsByTagName("CruiceCode")

    

    liste = {}
    for CruiceIndex in CruiceCode:
        liste[CruiceIndex.getAttribute('code')]= {}
        liste[CruiceIndex.getAttribute('code')]['cruice_id'] =CruiceIndex.getAttribute('code')
        liste[CruiceIndex.getAttribute('code')]['nation'] =CruiceIndex.getAttribute('nation')
        liste[CruiceIndex.getAttribute('code')]['platform_code'] =CruiceIndex.getAttribute('VCode')
        liste[CruiceIndex.getAttribute('code')]['platform_name'] =CruiceIndex.getAttribute('vessel')
        liste[CruiceIndex.getAttribute('code')]['platform_type'] =CruiceIndex.getAttribute('platform_type')
        liste[CruiceIndex.getAttribute('code')]['platform_code'] =CruiceIndex.getAttribute('platform_code')
        liste[CruiceIndex.getAttribute('code')]['SonarEquipment'] =CruiceIndex.getAttribute('Equipment')
        liste[CruiceIndex.getAttribute('code')]['vesselCode'] =CruiceIndex.getAttribute('vesselCode')
        liste[CruiceIndex.getAttribute('code')]['SonarDataPath'] =CruiceIndex.getAttribute('CruicePath')
    
    return liste
    
    

    
    
    
    
                        
                        
                        
    


        
        
        
    
    
                
                
def getOrganisedData(liste,WorkDirectory,OS, copy = False, reconvert = False): 
    '''Protocol to organize the data on the local server.'''
    
    import os
    from tools import tools
    import numpy as np
    from shutil import copyfile



    
    
        
        
        
        
    
            
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
            
            
            
            
            
            
    def CopyFileFromTapeServer(OS,liste,directory2Data): 
    
        #Get list of files that has not been copied to the correct structure
        print('    -Copy files from tape-server  ',end='\r')
        
        
        #Check if file exist, this is a bugfix when working localy
        if os.path.isdir(OS+liste['SonarDataPath']) == True: 
        
            print('    -Search files on tape-server  ',end='\r')       
            
            #Get list of files on server and list that has not been copied
            ListFromServer = os.listdir(OS+liste['SonarDataPath'])
            
            
            print('    -Get files not copied  ',end='\r')       
            ListOfFilesNotcopied = list(set(ListFromServer)-set(os.listdir(directory2Data.dir_originalrawdata)))
            
            
            print('    -Start copy files from tape-server  ',end='\r')
            #Go through each file that has not been copyed
            for i in np.arange(len(ListOfFilesNotcopied)): 
                
                
                
                #Print a progressbar for the user
                tools.printProgressBar(i, len(ListOfFilesNotcopied), prefix = 'Copy files', suffix = 'Files left: '+
                                 str(len(ListOfFilesNotcopied)-i) , decimals = 1, length = 50)
                
                
                if '.raw' in ListOfFilesNotcopied[i]: 
                    #Copy the files
                    try:
                        copyfile(OS+liste['SonarDataPath']+'/'+ListOfFilesNotcopied[i], 
                             directory2Data.dir_originalrawdata+'/'+ListOfFilesNotcopied[i])
                    except: 
                        print('      * Bad file', end = '\r')
                else: 
                    dummy = 1
                    
                    
        
                    
                    
    #Get the vessel, cruice and equipment information
    vesselCode = liste['vesselCode']
    cruiceCode = liste['cruice_id']
    
    
    #Loop through each fishery sonar type
    for i in ['SU90','SX90','SH90']: 
        if i == liste['SonarEquipment']: 
            
            
            #Make a correct sonar folder structure
            directory2Data =FolderStructure(WorkDirectory+'/'+cruiceCode[:4]+'/S'+cruiceCode+'_P'+vesselCode,i)       
            
            
            #Delete everything and start recompute
            if reconvert == True: 
                for i in os.listdir(directory2Data.dir_src) : 
                    if '.csv' in i: 
                        os.remove(directory2Data.dir_src+'/'+i)
                    
                    
            #Make the folder according to standard
            MakeNewFolders(directory2Data)
        
            if copy == True: 
                #Copy files from tape server
                CopyFileFromTapeServer(OS,liste,directory2Data)
            
            
                        
    return(directory2Data)
   
    
            
        
        
        