'''
ReadRawData.py

Last modify date: 24/11 2016

Usage: 
    ReadRawData()
        
        
        
Description: 
This is a function copied from a matlab script originated from SIMRAD


    
    
Author: 
Dr. Sindre Vatnehol (PhD)
Institute of Marine Research, Norway

Mail: 
sindre.vatnehol@imr.no

Tlf.: 
+47 900 79 376


Project: 
REDUS (Reducing uncertainty in stock assessment)  
'''

import numpy as np

from RawConverter import fread,StringArray,TimeConverter,bitget

#class ReadRawData:
def ReadRawData(fid,headerlength):
    
    class structtype(): 
        pass
    
    #Make the Filedata structurable 
    FileData = structtype()    
    
    #Define all the dictionary needed
    FileData.ChannelID=dict()
    FileData.beamtype=dict()
    FileData.frequency = dict()
    FileData.gain = dict()
    FileData.equivalentbeamangle = dict()
    FileData.beamwidthalongship = dict()
    FileData.beamwidthathwartship =dict()
    FileData.anglesensitivityalongship = dict()
    FileData.anglesensitivityathwartship = dict()
    FileData.angleoffsetalongship = dict()
    FileData.angleoffsetathwartship = dict()
    FileData.posx = dict()
    FileData.posy = dict()
    FileData.posz = dict()
    FileData.dirx = dict()
    FileData.diry = dict()
    FileData.dirz = dict()
    FileData.pulselengthtable = dict()
    FileData.spare2=dict()
    FileData.gaintable = dict()
    FileData.spare3=dict()
    FileData.sacorrectiontable = dict()
    FileData.spare4=dict()
    FileData.NMEA = dict()
    FileData.NMEA_time = dict()
    FileData.NMEA_info = dict()
    FileData.PingData=dict()
    
    
    #reade some line that is not saved
    fread.fread(fid,1,np.int32)
    
    
    #Write the datagram type
    FileData.datagramtype=StringArray.StringArray(fread.fread(fid,4,np.str) )

    #Datetime
    FileData.datetime=TimeConverter.TimeConverter(fread.fread(fid,2,np.uint32))
    
    #Write data about the sounder
    FileData.surveyname=StringArray.StringArray(fread.fread(fid,128,np.str))
    FileData.transectname=StringArray.StringArray(fread.fread(fid,128,np.str))
    FileData.soundername=StringArray.StringArray(fread.fread(fid,128,np.str))
    FileData.version=StringArray.StringArray(fread.fread(fid,30,np.str))
    FileData.multiplexing = fread.fread(fid,1,np.int16)
    FileData.timebias = fread.fread(fid,1,np.int32)
    FileData.soundvelocityaverage = fread.fread(fid,1,np.float32)
    FileData.soundvelocitytransducer = fread.fread(fid,1,np.float32)
    FileData.mruoffsetx = fread.fread(fid,1,np.float32)
    FileData.mruoffsety = fread.fread(fid,1,np.float32)
    FileData.mruoffsetz = fread.fread(fid,1,np.float32)
    FileData.mrualphax = fread.fread(fid,1,np.float32)
    FileData.mrualphay = fread.fread(fid,1,np.float32)
    FileData.mrualphaz = fread.fread(fid,1,np.float32)
    FileData.gpsoffsetx = fread.fread(fid,1,np.float32)
    FileData.gpsoffsety = fread.fread(fid,1,np.float32)
    FileData.gpsoffsetz = fread.fread(fid,1,np.float32)
    FileData.spare = StringArray.StringArray(fread.fread(fid,48,np.str))
    FileData.transducercount = fread.fread(fid,1,np.int32)
    
  
    
    for i in range(0,FileData.transducercount[0]):
        FileData.ChannelID[i]=StringArray.StringArray(fread.fread(fid,128,np.str))
        FileData.beamtype[i] = fread.fread(fid,1,np.int32)
        FileData.frequency[i] = fread.fread(fid,1,np.float32)
        FileData.gain[i] = fread.fread(fid,1,np.float32)
        FileData.equivalentbeamangle[i] = fread.fread(fid,1,np.float32)
        FileData.beamwidthalongship[i] = fread.fread(fid,1,np.float32)
        FileData.beamwidthathwartship[i] = fread.fread(fid,1,np.float32)
        FileData.anglesensitivityalongship[i] = fread.fread(fid,1,np.float32)
        FileData.anglesensitivityathwartship[i] = fread.fread(fid,1,np.float32)
        FileData.angleoffsetalongship[i] = fread.fread(fid,1,np.float32)
        FileData.angleoffsetathwartship[i] = fread.fread(fid,1,np.float32)
        FileData.posx[i] = fread.fread(fid,1,np.float32)
        FileData.posy[i] = fread.fread(fid,1,np.float32)
        FileData.posz[i] = fread.fread(fid,1,np.float32)
        FileData.dirx[i] = fread.fread(fid,1,np.float32)
        FileData.diry[i] = fread.fread(fid,1,np.float32)
        FileData.dirz[i] = fread.fread(fid,1,np.float32)
        FileData.pulselengthtable[i] = fread.fread(fid,5,np.float32)
        FileData.spare2[i]=StringArray.StringArray(fread.fread(fid,8,np.str))
        FileData.gaintable[i] = fread.fread(fid,5,np.float32)
        FileData.spare3[i]=StringArray.StringArray(fread.fread(fid,8,np.str))
        FileData.sacorrectiontable[i] = fread.fread(fid,5,np.float32)
        FileData.spare4[i]=StringArray.StringArray(fread.fread(fid,52,np.str))
    

    #Something loaded byt never used
    fread.fread(fid,1,np.int32)
 
    
    #Some bookkeeping
    ping_count = 0
    NMEA_count = 0
    
    
    #Continiue to read to the end of the file
    while(1):
        bytes2read = fread.fread(fid,1,np.int32)
        datagramtype=StringArray.StringArray(fread.fread(fid,4,np.str))
        
        if not bytes2read:
            break
        
        #get time
        NMEA= fread.fread(fid,1,np.uint64)
       
    
        #Writn con1 information
        if datagramtype == 'CON1':
            FileData.text = StringArray.StringArray(fread.fread(fid,(bytes2read-headerlength),np.str))
        
            
        #Write NMEA information
        elif datagramtype == 'NME0':
            FileData.NMEA_time[NMEA_count]=NMEA
            FileData.NMEA_info[NMEA_count]=StringArray.StringArray(fread.fread(fid,(bytes2read-headerlength),np.str))
            NMEA_count = NMEA_count+1
            
        elif datagramtype == 'TAG0':
            FileData.text = StringArray.StringArray(fread.fread(fid,(bytes2read-headerlength),np.str))
            
            
        elif datagramtype == 'RAW0':
            print('missing this function. It should be echo sounder data')
            
        elif datagramtype == 'RAW1':
            
            
            BeamCount = fread.fread(fid,1,np.int16)[0]
            
                        
                        
            if BeamCount == 1: 
                ping_count = ping_count+1
                
                FileData.PingData[ping_count-1] = structtype()
                    #Beam number
                FileData.PingData[ping_count-1].channel=dict()
                    #Datatypes
                FileData.PingData[ping_count-1].datatype=dict()
                    #number of complex values per sample
                FileData.PingData[ping_count-1].ncomplexpersample = dict()
                    #Gain TX, used in conversion to Sc
                FileData.PingData[ping_count-1].gaintx=np.array([])
                    #Frequency
                FileData.PingData[ping_count-1].frequency=np.array([])
                    #Transmit power - Used in conversion to Sv
                FileData.PingData[ping_count-1].transmitpower=np.array([])
                    #Signal pulse length
                FileData.PingData[ping_count-1].pulslength=np.array([])
                    #Signal band width
                FileData.PingData[ping_count-1].bandwidth=dict()
                FileData.PingData[ping_count-1].sampleinterval=dict()
                FileData.PingData[ping_count-1].soundvelocity=np.array([])
                FileData.PingData[ping_count-1].absorptioncoefficient=dict()
                FileData.PingData[ping_count-1].heave=dict()
                FileData.PingData[ping_count-1].roll=dict()
                FileData.PingData[ping_count-1].pitch=dict()
                FileData.PingData[ping_count-1].temperature=dict()
                FileData.PingData[ping_count-1].heading=dict()
                FileData.PingData[ping_count-1].transmitmode = dict()
                FileData.PingData[ping_count-1].pulseform = dict()
                FileData.PingData[ping_count-1].dirx = np.array([])
                FileData.PingData[ping_count-1].diry = np.array([])
                FileData.PingData[ping_count-1].dirz = np.array([])
                FileData.PingData[ping_count-1].gainrx = np.array([])
                FileData.PingData[ping_count-1].sacorrection = np.array([])
                FileData.PingData[ping_count-1].equivalentbeamangle =np.array([])
                FileData.PingData[ping_count-1].beamwidthalongshiprx = np.array([])
                FileData.PingData[ping_count-1].beamwidthathwartshiprx = np.array([])
                
                
                FileData.PingData[ping_count-1].anglesensitivityalongship = dict()
                FileData.PingData[ping_count-1].anglesensitivityathwartship = dict()
                FileData.PingData[ping_count-1].angleoffsetalongship = dict()
                FileData.PingData[ping_count-1].angleoffsetathwartship = dict()
                FileData.PingData[ping_count-1].spare1=dict()
                FileData.PingData[ping_count-1].noisefilter =dict()
                FileData.PingData[ping_count-1].beamwidthmode = dict()
                FileData.PingData[ping_count-1].beammode =dict()
                FileData.PingData[ping_count-1].beamwidthhorizontaltx = dict()
                FileData.PingData[ping_count-1].beamwidthverticaltx = dict()
                FileData.PingData[ping_count-1].offset =dict()
                FileData.PingData[ping_count-1].count =dict()
                FileData.PingData[ping_count-1].BeamAmplitudeData=np.array([])
                FileData.PingData[ping_count-1].BeamAmplitudeData_imaginary=np.array([])
                FileData.PingData[ping_count-1].power=np.array([])
                
                
                FileData.NMEA[ping_count-1] = NMEA
                

            #read information
            t=fread.fread(fid,2,np.int8)     
            
            #write inforamtion
            FileData.PingData[ping_count-1].datatype[BeamCount-1] =t[0]
            FileData.PingData[ping_count-1].ncomplexpersample[BeamCount-1]  = t[1]
            
            #read information
            t=fread.fread(fid,13,np.float32)
            
            #Write information
            FileData.PingData[ping_count-1].gaintx=np.hstack((FileData.PingData[ping_count-1].gaintx,t[0])) #[BeamCount-1] =t[0]
            FileData.PingData[ping_count-1].frequency=np.hstack((FileData.PingData[ping_count-1].frequency,t[1]))
            FileData.PingData[ping_count-1].transmitpower=np.hstack((FileData.PingData[ping_count-1].transmitpower,t[2]))
            FileData.PingData[ping_count-1].pulslength=np.hstack((FileData.PingData[ping_count-1].pulslength,t[3]))
            FileData.PingData[ping_count-1].bandwidth[BeamCount-1] =t[4]
            FileData.PingData[ping_count-1].sampleinterval[BeamCount-1] =t[5]
            FileData.PingData[ping_count-1].soundvelocity =np.hstack((FileData.PingData[ping_count-1].soundvelocity,t[6]))
            FileData.PingData[ping_count-1].absorptioncoefficient[BeamCount-1] =t[7]
            FileData.PingData[ping_count-1].heave[BeamCount-1] =t[8]
            FileData.PingData[ping_count-1].roll[BeamCount-1] =t[9]
            FileData.PingData[ping_count-1].pitch[BeamCount-1] =t[10]
            FileData.PingData[ping_count-1].temperature[BeamCount-1] =t[11]
            FileData.PingData[ping_count-1].heading[BeamCount-1] =t[12]

            #read information
            t=fread.fread(fid,2,np.int16)
            
            #write information
            FileData.PingData[ping_count-1].transmitmode[BeamCount-1]  = t[0]
            FileData.PingData[ping_count-1].pulseform[BeamCount-1]  = t[1]

            #read information
            t=fread.fread(fid,12,np.float32)
           
            
            #Write information
            FileData.PingData[ping_count-1].dirx  = np.hstack((FileData.PingData[ping_count-1].dirx,t[0]))
            FileData.PingData[ping_count-1].diry= np.hstack((FileData.PingData[ping_count-1].diry,t[1]))
            FileData.PingData[ping_count-1].dirz = np.hstack((FileData.PingData[ping_count-1].dirz,t[2]))
            FileData.PingData[ping_count-1].gainrx=np.hstack((FileData.PingData[ping_count-1].gainrx,t[3]))
            FileData.PingData[ping_count-1].sacorrection=np.hstack((FileData.PingData[ping_count-1].sacorrection,t[4]))
            FileData.PingData[ping_count-1].equivalentbeamangle=np.hstack((FileData.PingData[ping_count-1].equivalentbeamangle,t[5]))
            FileData.PingData[ping_count-1].beamwidthalongshiprx=np.hstack((FileData.PingData[ping_count-1].beamwidthalongshiprx,t[6]))
            FileData.PingData[ping_count-1].beamwidthathwartshiprx=np.hstack((FileData.PingData[ping_count-1].beamwidthathwartshiprx,t[7]))
            FileData.PingData[ping_count-1].anglesensitivityalongship[BeamCount-1]  = t[8]
            FileData.PingData[ping_count-1].anglesensitivityathwartship[BeamCount-1]  = t[9]
            FileData.PingData[ping_count-1].angleoffsetalongship[BeamCount-1]  = t[10]
            FileData.PingData[ping_count-1].angleoffsetathwartship[BeamCount-1]  = t[11]
            FileData.PingData[ping_count-1].spare1[BeamCount-1] =StringArray.StringArray(fread.fread(fid,2,np.str))
            
            #read information
            t=fread.fread(fid,3,np.int16)
            
            #write information
            FileData.PingData[ping_count-1].noisefilter[BeamCount-1]  = t[0]
            FileData.PingData[ping_count-1].beamwidthmode[BeamCount-1]  = t[1]
            FileData.PingData[ping_count-1].beammode[BeamCount-1]  = t[2]

            #read informaiton
            t=fread.fread(fid,2,np.float32)
            
            #write information
            
            FileData.PingData[ping_count-1].beamwidthhorizontaltx[BeamCount-1]  = t[0]
            FileData.PingData[ping_count-1].beamwidthverticaltx[BeamCount-1]  = t[1]

            #read information
            t=fread.fread(fid,2,np.int32)            
            
            #write information
            FileData.PingData[ping_count-1].offset[BeamCount-1]  = t[0]
            FileData.PingData[ping_count-1].count[BeamCount-1]  = t[1]


            #read the acousitcal data from raw file
            if bitget.bitget(FileData.PingData[ping_count-1].datatype[BeamCount-1] ,4)==1: 
                

                #read the data
                data = fread.fread(fid,2*(FileData.PingData[ping_count-1].count[BeamCount-1]),
                                          np.float32)
                #print(data.shape)
                try:
                    data = data.reshape(FileData.PingData[ping_count-1].count[BeamCount-1],2).transpose() 
                except ValueError: 
                    data = np.nan*np.ones((FileData.PingData[ping_count-1].count[BeamCount-1],2)).transpose()

                #sort the data into real and imaginary
                realstuff = data[0,0:t[1]]
                imaginarystuff = data[1,0:t[1]]
              
                
                #If the data is the first beam
                if len(data)>len(FileData.PingData[ping_count-1].BeamAmplitudeData):
                    FileData.PingData[ping_count-1].BeamAmplitudeData=realstuff[:, np.newaxis]
                    FileData.PingData[ping_count-1].BeamAmplitudeData_imaginary=imaginarystuff[:, np.newaxis]

                #For other beams
                else:
                    FileData.PingData[ping_count-1].BeamAmplitudeData=np.hstack((FileData.PingData[ping_count-1].BeamAmplitudeData,realstuff[:, np.newaxis]))
                    FileData.PingData[ping_count-1].BeamAmplitudeData_imaginary=np.hstack((FileData.PingData[ping_count-1].BeamAmplitudeData_imaginary,imaginarystuff[:, np.newaxis]))

                
        #somehting read but never used        
        fread.fread(fid,1,np.int32)
        
    #return the data
    return FileData
    
    
