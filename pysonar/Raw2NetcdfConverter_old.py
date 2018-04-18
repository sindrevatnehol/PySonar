'''
Description: 
This Modul convert the .raw files to .nc files. 
In future, the function should be abe to store several pings within one file


to install pynmea2 type following in command line
    pip install pynmea2    


    
Author: 
Dr. Sindre Vatnehol (PhD)
Institute of Marine Research, Norway

Mail: 
sindre.vatnehol@imr.no

Tlf.: 
+47 900 79 376


Project: 
REDUS (Reducing uncertainty in stock assessment)  
www.redus.no
'''



def TimeConverter(time_date):
    import datetime   
    fulldate = datetime.datetime.strptime('1601-01-01 00:00:00.000',"%Y-%m-%d %H:%M:%S.%f")
    try:  
        temp_date = (time_date[1]*2**32 + time_date[0])/10000
        fulldate = datetime.datetime.strptime('1601-01-01 00:00:00.000',"%Y-%m-%d %H:%M:%S.%f")
    
    
        fulldate = fulldate + datetime.timedelta(milliseconds=int(temp_date))

    except IndexError: 
        fulldate = 0
    return fulldate; 
    



def Raw2NetcdfConverter2(directory,directoryOutput,DirectoryToNCProg): 
    
    import os, glob, pynmea2
    from RawConverter import ReadRawData
    from netCDF4 import Dataset
    import numpy as np
    
    #Give teh header length        
    headerlength = 12
    
    
    #NumberOfBeams. Rewrite this to something more general
    numberOfBeams = 64 
    
    LastFile = False
    lastFile = '' 
    
    #Check status file exist
    if not os.path.exists(DirectoryToNCProg+'/ConvertedFiles.txt'): 
        LastFile = True
    else: 
        f_id = open(DirectoryToNCProg+'/ConvertedFiles.txt','r')
        lastFile = f_id.read()
        f_id.close()
        
        
        
        
                    
    #Change directory to the raw folder
    os.chdir(directory)
    
    
    #Loop through all raw files inside folder
    for filename in sorted(glob.glob("*.raw")):
        if filename == lastFile:
            LastFile = True
            
            
        if LastFile == True: 
            os.chdir(directory)

            #open .raw file
            fid = open(filename,mode='rb')
            
            
            #Read the .raw file
            FileData = ReadRawData.ReadRawData(fid,headerlength)
            
            print('File '+filename+' is loaded             ' , end='\r')
            
            #Change directory to netcdf folder
            os.chdir(directoryOutput)
            
            
            
            #Some old stuff. May be deleted when properly tested
            oldTime = '0'
            
            
            
            
            NMEATime = []
            Lat = []
            Lon = []
            for i in range(len(FileData.NMEA_info)):
                try: 
                    msg = pynmea2.parse(FileData.NMEA_info[i])
                    try: 
                        try: 
                            np.float(msg.lon)/100+float(msg.lat)/100+int(str(msg.timestamp).replace(':','').replace('.',''))
                            Lon = np.hstack((Lon,np.float(msg.lon)/100))
                            Lat = np.hstack((Lat,np.float(msg.lat)/100))
                            NMEATime = np.hstack((NMEATime,int(str(msg.timestamp).replace(':','').replace('.',''))))
                        except ValueError: 
                            k=1
                    except AttributeError: 
                        k=1
                except pynmea2.nmea.ParseError: 
                    k=1
                    
                            
                    
                    
            #Loop through every ping in raw fil
            for increment in range(0,len(FileData.PingData)): 
                
                
                #If there is information inside the ping continue
                if not FileData.NMEA[increment] ==[]:
                    #Get the time of the ping from nmea
                    Time = FileData.NMEA[increment]


                    import datetime
                    fulldate = datetime.datetime.strptime('1601-01-01 00:00:00.000',"%Y-%m-%d %H:%M:%S.%f")
                    time = fulldate + datetime.timedelta(milliseconds=int(Time[0]/10000))
                    #a = datetime.datetime.fromtimestamp(time).strftime('%Y-%m-%f %H:%M:%S.%f')
                    #time= TimeConverter(Time)
                    time = str(time).replace("-","").replace(" ","").replace(":","").replace(".","")
                    
                    Time = time
 
                    IdxNMEAtime =  np.where((np.abs(NMEATime-int(str(Time)[8:])))==np.min(np.abs(NMEATime-int(str(Time)[8:]))))
                            
            
                    '''
                    Temporarly fix for fixing a bug in time indexing
                    It should be redundant now. 
                    To be deleted when properly tested
                    '''
                    if Time == oldTime: 
                        Time = (int(Time)+1)
                    elif Time == (int(oldTime)-1):
                        Time = (int(Time)+2)
                        
                    Time = str(Time)
                    '''
                    Underneath: 
                    Get all the variables from raw fil
                    '''
                    NetCDF_filename = (filename[:(filename.index('-')+1)]+Time+'.nc')
                    
                    mrualphax = FileData.mrualphax[0]
                                                        
                    sampleinterval = FileData.PingData[increment].sampleinterval[0]
    
                    gaintx = FileData.PingData[increment].gaintx     #Vect
                    
                    datatype = FileData.PingData[increment].datatype[0]        
    
                    ncomplexpersample = FileData.PingData[increment].ncomplexpersample[0]   
    
                    frequency = FileData.PingData[increment].frequency[0]                  
    
                    transmitpower = FileData.PingData[increment].transmitpower[0]       
    
                    pulslength = FileData.PingData[increment].pulslength[0]     
    
                    bandwidth = FileData.PingData[increment].bandwidth[0] 
    
                    soundvelocity = FileData.PingData[increment].soundvelocity[0]
    
                    absorptioncoefficient = FileData.PingData[increment].absorptioncoefficient[0]
    
                    heave = FileData.PingData[increment].heave[0]
    
                    roll = FileData.PingData[increment].roll[0]
    
                    pitch = FileData.PingData[increment].pitch[0]
    
                    temperature = FileData.PingData[increment].temperature[0]
    
                    heading = FileData.PingData[increment].heading[0]
    
                    transmitmode = FileData.PingData[increment].transmitmode[0]
    
                    pulseform = FileData.PingData[increment].pulseform[0]
    
                    dirx = FileData.PingData[increment].dirx
    
                    diry = FileData.PingData[increment].diry
    
                    dirz = FileData.PingData[increment].dirz
    
                    for i_v in range(len(dirz)): 
                        dirz[i_v] = FileData.dirz[0]
    
                    gainrx = FileData.PingData[increment].gainrx
    
                    sacorrection = FileData.PingData[increment].sacorrection
    
                    equivalentbeamangle = FileData.PingData[increment].equivalentbeamangle
    
                    beammode = FileData.PingData[increment].beammode[0]
    
                    BeamAmplitudeData =  FileData.PingData[increment].BeamAmplitudeData

                    Longitude = Lon[IdxNMEAtime][0]

                    Latitude = Lat[IdxNMEAtime][0]

            
                    try: 
                        [R, L] = BeamAmplitudeData.shape
                        
                        if L ==64: 
                            '''
                            Create the netcdf file
                            In future this should include several of pings in 
                            one file
                            '''
                            f=Dataset(NetCDF_filename,'w',format='NETCDF4')
                            
                            f.createDimension('we',1)
                            
                            f.createDimension('be',numberOfBeams)   #Make this more general
                            
                            f.createDimension('re',len(BeamAmplitudeData))
                            
                            
                            
                            '''
                            Underneath: 
                                make all the veriables in netcdf format
                            '''
                            innput = f.createVariable('sampleinterval','d',('we','we'))
                            innput[:] = sampleinterval
                            
                            
                            innput = f.createVariable('mrualphax','d',('we','we'))
                            innput[:] = mrualphax
                            
                            innput = f.createVariable('gaintx','d',('be','we'))
                            innput[:] = gaintx
                            
                            innput = f.createVariable('datatype','d',('we','we'))
                            innput[:] = datatype
                            
                            innput = f.createVariable('ncomplexpersample','d',('we','we'))
                            innput[:] = ncomplexpersample
                            
                            innput = f.createVariable('frequency','d',('we','we'))
                            innput[:] = frequency
                            
                            innput = f.createVariable('transmitpower','d',('we','we'))
                            innput[:] = transmitpower
                            
                            innput = f.createVariable('pulslength','d',('we','we'))
                            innput[:] = pulslength
                            
                            innput = f.createVariable('bandwidth','d',('we','we'))
                            innput[:] = bandwidth
                            
                            innput = f.createVariable('soundvelocity','d',('we','we'))
                            innput[:] = soundvelocity
        
                            innput = f.createVariable('absorptioncoefficient','d',('we','we'))
                            innput[:] = absorptioncoefficient
        
                            innput = f.createVariable('heave','d',('we','we'))
                            innput[:] = heave
        
                            innput = f.createVariable('roll','d',('we','we'))
                            innput[:] = roll
        
                            innput = f.createVariable('pitch','d',('we','we'))
                            innput[:] = pitch
        
                            innput = f.createVariable('temperature','d',('we','we'))
                            innput[:] = temperature
        
                            innput = f.createVariable('heading','d',('we','we'))
                            innput[:] = heading
        
                            innput = f.createVariable('transmitmode','d',('we','we'))
                            innput[:] = transmitmode
        
                            innput = f.createVariable('pulseform','d',('we','we'))
                            innput[:] = pulseform
        
                            innput = f.createVariable('dirx','d',('be','we'))
                            innput[:] = dirx
        
                            innput = f.createVariable('diry','d',('be','we'))
                            innput[:] = diry
        
                            innput = f.createVariable('dirz','d',('be','we'))
                            innput[:] = dirz
        
                            innput = f.createVariable('gainrx','d',('be','we'))
                            innput[:] = gainrx
        
                            innput = f.createVariable('sacorrection','d',('be','we'))
                            innput[:] = sacorrection
        
                            innput = f.createVariable('beammode','d',('we','we'))
                            innput[:] = beammode
        
                            innput = f.createVariable('equivalentbeamangle','d',('be','we'))
                            innput[:] = equivalentbeamangle
        
                            innput = f.createVariable('BeamAmplitudeData','f8',('re','be'))
                            innput[:,:] = BeamAmplitudeData
        
                            innput = f.createVariable('Longitude','d',('we','we'))
                            innput[:] = Longitude
        
                            innput = f.createVariable('Latitude','d',('we','we'))
                            innput[:] = Latitude
        
                            innput = f.createVariable('Time','d',('we','we'))
                            innput[:] = Time
        
        
                            '''This is perhaps redundant. Will be deleted 
                            when properly tested
                            '''
                            oldTime = Time
        
                            f.close()            
                
                    except ValueError: 
                        print('Bad Ping')
                        
            f_id = open(DirectoryToNCProg+'/ConvertedFiles.txt','w')
            f_id.write(filename)
            f_id.close()
            

            
