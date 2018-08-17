# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 14:23:28 2018

@author: sindrev
"""


import os,platform
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
from pysonar_process import doVerticalProcess,doHorizontalProcess
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
        lista = '2016844'
        lista = '2017836'
        
        
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
        
        
        print('Start reading LUF20')
#        try: 
        LUF20_info_list = getLuf20Info(directory2Data,lista)
#        except: 
            
#            send_email('Read EKLUF20 failed for ' +lista)
        print('Luf20 has been read')
        
        if GO_vertical == True: 
            
            
            try: 
                idx_list_vertical = readIDX(directory2Data,'Vertical')
            except: 
                makeIDX(directory2Data,'Vertical')
                idx_list_vertical = readIDX(directory2Data,'Vertical')
            
            
            print('Start process vertical')
#            doVerticalProcess(directory2Data,idx_list_vertical,LUF20_info_list,liste[lista],RemoveToCloseValues,R_s,res,randomize=True)
            
            nation = liste[lista]['nation']
            cruice_id = liste[lista]['cruice_id']
            vplatform = liste[lista]['vesselCode']
            
            print('Start making report per log')
            #Make Report files               
            MakeReport(idx_list_vertical,directory2Data, nation,cruice_id,vplatform)
            
            
            print('Merge report')
            #Make the LUF report files
            reader(directory2Data, LUF20 = True, LUFICES = False)
        

        
                    
                
                
                
                
            
            
            
            
            
            
            
                
        #Loop through each transect in a randomized maner
        if GO_horizontal == True: 
            
            
            try: 
                idx_list_vertical = readIDX(directory2Data,'Horizontal')
            except: 
                makeIDX(directory2Data,'Horizontal')
                idx_list_vertical = readIDX(directory2Data,'Horizontal')
                
            
            doHorizontalProcess(directory2Data,lista)
            


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
