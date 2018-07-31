# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 14:23:28 2018

@author: sindrev
"""


import os,platform, random
#import numpy as np
from SendMail import send_email
from tools import tools
#import scipy.io as sc
#from MakeSearch import MakeSearch
from pysonar_makeverticalIDX import MakeReport
#from MakeWork import MakeWork
#from HorizontalLuf20 import MakeHorizontalLuf20
from pysonar_organiser import getOrganisedData, getListOfSurveys
import Raw2NetcdfConverter
from pysonar_makeIDX import makeIDX, readIDX, getLuf20Info
from pysonar_process import doVerticalProcess
from pysonar_lufgenerator import reader


    
    


def main(threshold, RemoveToCloseValues, R_s, recompute, reconvert,GO_horizontal, GO_vertical, Stack, convert):

    from pysonar_makeverticalIDX import MakeVerticalIndex
    
    
    
    #Clear the terminal window and display a welcome screen.
    os.system('cls' if os.name == 'nt' else 'clear')
    tools.WelcomeScreen('PySonar.py')


    
    
    #Some user innputs that will be moved
    res = R_s/0.18
    maxPingInFile = 2000
    MaxNumberOfFilesInNC  = 100
    
    

    





    #Something to adapt for data organization at the IMR
    #For other organisations, this should be changed
    if platform.system() == 'Linux':
        OS = '//data'
    else:
        OS = '//ces.imr.no'



        
        

    #Get the work directory to where all the files is stoored
    #This should be merged with the one above
    WorkDirectory = OS+ '/mea/2018_Redus'
#    WorkDirectory = 'F:'



    #Get the list of surveys from xml file
    liste = getListOfSurveys(directory = os.getcwd())

    
    
    #Loop through all surveys
    for lista in liste: 
#        lista = '2018105'
        
        
        #something for the user
        os.system('cls' if os.name == 'nt' else 'clear')
        tools.WelcomeScreen('PySonar.py')
        print('Start on survey: '+lista)
#        send_email('Start on survey: ' +lista)
        
        
        #Get the structure of all data
#        try: 
        directory2Data = getOrganisedData(liste[lista],WorkDirectory,OS)
#        except: 
            
#            send_email('Get directory2data failed for' +lista)
                
        
        #Convert .raw to .nc 
        #This is from an secondary package
        if convert == True: 
            try: 
                Raw2NetcdfConverter.convert.Raw2NetcdfConverter(directory2Data.dir_originalrawdata,liste[lista]['platform_name'],
                                liste[lista]['platform_type'],liste[lista]['platform_code'],
                                maxPingInFile,MaxNumberOfFilesInNC,
                                directory2Data.dir_rawdata,reconvert)
            except: 
                
                    send_email('Convertion of files failed for ' +lista)
        
        
            os.system('cls' if os.name == 'nt' else 'clear')
            tools.WelcomeScreen('PySonar.py')
            print('Continue on survey: '+lista)
        
        
        
        
        
        
                        
#        try: 
#            idx_list_horizontal = readIDX(directory2Data,'Horizontal')
#        except: 
#            makeIDX(directory2Data,'Horizontal')
#            idx_list_horizontal = readIDX(directory2Data,'Horizontal')
        try: 
            idx_list_vertical = readIDX(directory2Data,'Vertical')
        except: 
            makeIDX(directory2Data,'Vertical')
            idx_list_vertical = readIDX(directory2Data,'Vertical')
            
        
        
        try: 
            LUF20_info_list = getLuf20Info(directory2Data,lista)
        except: 
            
            send_email('Read EKLUF20 failed for ' +lista)
            
        
        if GO_vertical == True: 
            print('Start process vertical')
            doVerticalProcess(directory2Data,idx_list_vertical,LUF20_info_list,liste[lista],RemoveToCloseValues,R_s,res,randomize=True)
            nation = liste[lista]['nation']
            cruice_id = liste[lista]['cruice_id']
            vplatform = liste[lista]['vesselCode']
            
            
                
            #Make Report files               
            MakeReport(idx_list_vertical,directory2Data, nation,cruice_id,vplatform)
            
            
            #Make the LUF report files
            reader(directory2Data, LUF20 = True, LUFICES = False)
        

        
#                    
#                
#                
#                
#                
#                
#        #Loop through each transect in a randomized maner
#        if GO_horizontal == True: 
#            print('Start making search matrix',end='\r')
#            
#            
#            
#            
#            #Get list of all transect and randomize it
#            NumTransect = np.arange(TimeIDX.shape[0])
#            random.shuffle(NumTransect)
#            
#            
#            
#            
#            
#            #Loop through each transect
#            for Transect in NumTransect:
#                
#    
#                
#                
#                
#                #Make the search matrix
#                if os.path.isfile(directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat') == False:
#                    print('    -Start on transect: '+ TimeIDX[Transect,0],end='\r')
##                    send_email('Start on horizontal transect: '+ TimeIDX[Transect,0])
#            
#                    
#                    
#                
#                    #Get list of files within time intervals
#                    ShortListOfFiles_horizontal = tools.GetShortListOfFiles(CompleteListOfFiles_horizontal,
#                                                                      TimeIDX[Transect,1].replace('T','')
#                                                                      ,TimeIDX[Transect,2].replace('T',''))
#                    
#                    
#                    
#                    #Make the search matrix
#                    MakeSearch(ShortListOfFiles_horizontal,
#                               RemoveToCloseValues,
#                               R_s,
#                               res,
#                               directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat',
#                               directory2Data.dir_rawdata,
#                               beamgrp_hor)
#                    
#                    
#                    
#                    
#                elif recompute == True:
#                    print('    -Start on transect: '+ TimeIDX[Transect,0],end='\r')
#                    send_email('Start on horizontal transect: '+ TimeIDX[Transect,0])
#            
#                
#                    
#                    
#                    #Get list of files within time intervals
#                    ShortListOfFiles_horizontal = tools.GetShortListOfFiles(CompleteListOfFiles_horizontal,
#                                                                      TimeIDX[Transect,1].replace('T','')
#                                                                      ,TimeIDX[Transect,2].replace('T',''))
#                    
#                    
#                    
#                    
#                    
#                    #Make search matrix
#                    MakeSearch(ShortListOfFiles_horizontal,
#                               RemoveToCloseValues,
#                               R_s,
#                               res,
#                               directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat',
#                               directory2Data.dir_rawdata,
#                               beamgrp_hor)
#                    
#                    
#                    
#                    
#                    
#                
#                
#            #Make Work files
#           
#            #Get list of all transect and randomize it
#            NumTransect = np.arange(TimeIDX.shape[0])
#                        
#            
#            
#            #Loop through each transect
#            for Transect in NumTransect:
#                
#                
#                
#                
#                print('Start on horizontal transect: '+ str(Transect),end='\r')
#    #            send_email('Start on horizontal transect: '+ str(Transect))
#                
#                
#                
#                
#                if os.path.isfile(directory2Data.dir_search+'/'+TimeIDX[Transect,0]+'.mat') == True:
#                    MakeWork(True, directory2Data,'','','',1,3)
#    
#                
#                    
#                    
#                    
#                
#                
#            #Start making the luf20 files
#    #        send_email('Start Making Horizontal LUF20: '+ str(Transect))
#            MakeHorizontalLuf20(CompleteListOfFiles_horizontal,directory2Data,start_time,
#                                log_start,stop_time,lat_start,lat_stop,lon_start,lon_stop,
#                                TimeIDX, nation,cruice_id,vplatform)
#                    
#                    
#            
#            


if __name__ == '__main__':

    #TODO: enable so these can be selected from commandline
    
    threshold = 2
    RemoveToCloseValues = 30
    R_s = 15
    
    
    
    
    #Some user innputt that will be selectable as inputs
    recompute = False
    reconvert = False
#    convert = True
    GO_horizontal = False
    GO_vertical = True
    Stack = True
    
    


    main(threshold, RemoveToCloseValues, R_s, recompute, reconvert,GO_horizontal, GO_vertical, Stack, convert = False)
