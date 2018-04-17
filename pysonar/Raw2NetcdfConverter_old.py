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


#import som modual
#
#
#
#
#def Raw2NetcdfConverter(directoryToRaw,vessel_name,platform_type,
#                        platform_code,maxPingInFile,MaxNumberOfFilesInNC,
#                        WhereNMEA,WhereNMEAheading,directory): 
#    
#    
#    import glob,os, datetime, pynmea2,os.path, pytz
#    from netCDF4 import Dataset
#    import numpy as np
#    from shutil import copyfile
#    from RawConverter import ReadRawData
#    
#    
#    #Get current direction
#    #This is used to get back to the template 
#    cwd=os.getcwd()
#    
#    
#    
#    #Check if netcdf folder exist
#    #If not create the directory
#    if (not os.path.isdir(directory)): 
#        os.makedirs(directory)
#        
#        
#    
#    #Version controll of the conversion, spesific for this script
#    ConversionSoftwareName = 'IMR_RAW_2_NETCDF4'
#    ConversionSoftwareVersion = '1.0.0'
#    
#    
#
#    #Some bokkeeping
#    NumberOfFilesInNC = 0
#    Frequency_index = []
#
#    
#
#    #Change the directory to the raw files
#    os.chdir(directoryToRaw)
#    
#    
#    #Go through each .raw file in folder
#    for filename in sorted(glob.glob("*.raw")):
#        
#        
#        #Display msg to user
#        print('    -Converting '+filename+' to netcdf',end='\r')
#        
#        
#        #Change the directory back to .raw folder. 
#        #Later it will be send to the netcdf folder
#        os.chdir(directoryToRaw)
#        
#        
#        
#        #open .raw file
#        fid = open(filename,mode='rb')
#        
#       
#        #Read the .raw file according to SIMRAD spesifications
#        FileData = ReadRawData.ReadRawData(fid,12)
#        
#        #Get the transduser rotation
#        #This is a correction of the true direction of the beams
#        transducer_rotation =  FileData.dirz
#                    
#        
#        #Set the directory to where the netcdf should be stored
#        os.chdir(directory)
#        
#        
#        #If first file in .nc file
#        if NumberOfFilesInNC == 0:
#        
#            
#            #Set the name of the .nc file to the same as the first .raw file                        
#            ncfilename=filename[:-4]+'.nc'
#        
#    
#            #Copy template to current folder, and rename it
#            copyfile(cwd+'\\src\\netcdf4.nc',ncfilename)
#            
#            
#            #add index of current .raw files in .nc
#            NumberOfFilesInNC = NumberOfFilesInNC+1
#            
#            
#            #Some bookkeeping                    
#            vertical_ping = 0
#            horizontal_ping = 0
#            annotation_index = 0 
#            old_time = 0                       
#            LoopIndex1 = 0
#            LoopIndex2 = 0
#            LoopIndex3= 0
#            Frequency_index = []
#                             
#            
#            #open new .nc file
#            f =Dataset(ncfilename,'a') 
#            
#            
#            #Create a new file
#            #fileID=Dataset('test.nc','w') 
#            #fileID.setncattr('Conventions','CF-1.6')
#            #fileID.setncattr('file_created',str(pytz.timezone('Europe/Oslo').localize(datetime.datetime.utcnow())))
#            #fileID.setncattr('file_format_authority','ICES')
#            #fileID.setncattr('file_format_name','SONAR-netCDF4')
#            #fileID.setncattr('file_format_version','1.0')
#            #fileID.setncattr('file_licensing','No licensing applied')
#            #fileID.setncattr('file_rights','Unrestricted rights')
#            
#            
#        else:
#            #Increase file index
#            NumberOfFilesInNC = NumberOfFilesInNC+1
#            
#            
#            #if last file to be added to .nc
#            if NumberOfFilesInNC >= MaxNumberOfFilesInNC: 
#                NumberOfFilesInNC=0
#                Frequency_index = []
#
#
#        #Msg to user about progress
##        print('     - Write NMEA')
#        
##        try: 
##            fileID.groups['Provenance']
##        except KeyError :
##            grp = fileID.createGroup('Provenance')
##            grp.setncattr('conversion_software_name',ConversionSoftwareName)
##            grp.setncattr('conversion_software_version',ConversionSoftwareVersion)
##            grp.setncattr('conversion_time',str(pytz.timezone('Europe/Oslo').localize(datetime.datetime.utcnow())))
#        
#        
#            
#            
#       
#        #Set the time of creation to the .nc file
#        filetime=pytz.timezone('Europe/Oslo').localize(datetime.datetime.utcnow())
#        
#        
#        #Write creation of file
#        f.file_created = str(filetime) 
#        
#        
#        #Add .raw file name to list in Provenance
#        #Add this so it is more easier to see which files that is included
#        addvar = f.groups['Provenance'].variables['source_filenames']
#        addvar[NumberOfFilesInNC-1] = filename
#
#     
#
#        #Write Provenance information about the new file
#        f.groups['Provenance'].conversion_software_name= ConversionSoftwareName 
#        f.groups['Provenance'].conversion_software_version= ConversionSoftwareVersion
#        f.groups['Provenance'].conversion_time= str(filetime)
#        
#        
#        #Add information of the sonar
#        f.groups['Sonar'].sonar_model = str(FileData.soundername)
#        f.groups['Sonar'].sonar_software_version = str(FileData.version)
#        f.groups['Sonar'].sonar_serial_number = 'NaN'
#        
#
#        #Add information of the platform
#        f.groups['Platform'].platform_code_ICES = platform_code 
#        f.groups['Platform'].platform_name = vessel_name
#        f.groups['Platform'].platform_type = platform_type
#        
#        
#
#                
##        print('    Write Platform')
#                
#                
#                
#        #Run through each NMEA information 
#        for LoopIndex0 in range(len(FileData.NMEA_time)): 
#            
#            
#            #Load  NMEA time and info
#            Time = FileData.NMEA_time[LoopIndex0]
#            info = FileData.NMEA_info[LoopIndex0]
#
#
#            #A bug fix
#            if Time[0]>0: 
#                if LoopIndex0==0: 
#                    
#                    #Add some information of the annotation, may be deleted later                    
#                    addvar = f.groups['Annotation'].variables['time']
#                    addvar[annotation_index] = Time[0]
#
#                    addvar = f.groups['Annotation'].variables['annotation_category']
#                    addvar[annotation_index] = 'Update'
#                    
#                    addvar = f.groups['Annotation'].variables['annotation_text']
#                    addvar[annotation_index] = 'Converted '+filename+' to' +ncfilename
#
#                    annotation_index = annotation_index+1
#                    
#                    
#                    
#                #Open Platform group
#                
#                addvardir = f.groups['Platform'].groups['NMEA']
#                
#
#                #Add info to platform group 
#                addvar = addvardir.variables['time']                            
#                addvar[LoopIndex3]=  Time[0]
#                addvar = addvardir.variables['NMEA_datagram']                            
#                addvar[LoopIndex3]=  info
#                LoopIndex3=LoopIndex3+1
#                
#
#                #Get the preffered telegram, 'H
#                if info.split(',')[0][-3:]==WhereNMEA[-3:]: 
#                    
#                    msg = pynmea2.parse(info,check= False)
#                    
#                    
#                    #Check if nmea telegram is corrupted
#                    
#                    
#                    try:
#                        t= msg.latitude
#                        
#                        #Parse the msg
#                        
#                        if (Time[0]-old_time)>0: 
#                        
#                            #Add NMEA information to platform
#                            addvar = f.groups['Platform'].variables['latitude']                            
#                            addvar[LoopIndex1]=  msg.latitude
#                            
#                            addvar = f.groups['Platform'].variables['longitude']                            
#                            addvar[LoopIndex1]=  msg.longitude
#                            
#                            addvar = f.groups['Platform'].variables['time1']                            
#                            addvar[LoopIndex1]=  Time[0]
#    
#    
#                            old_time = Time[0]
#                            
#                            
#                            LoopIndex1 = LoopIndex1+1 
#                        else: 
#                            print('         *NMEA was replicated',end='/r')
#                            
#                    except AttributeError:
#                        print('         *Bad GPS',end='/r')
#                        
#                        
#                        
#                        
#                #Get the information of vessel heading from the preferred nmea telegram                     
#                elif info.split(',')[0][-3:]==WhereNMEAheading[-3:]:
#                    
#                    'ok'
#                    #Parse information to get heading information
#                    msg = pynmea2.parse(info)
#                    
#                    
#                    #Write information
#                    addvar = f.groups['Platform'].variables['heading']                            
#                    addvar[LoopIndex2]=  msg.heading
#                    addvar = f.groups['Platform'].variables['time2']                            
#                    addvar[LoopIndex2]=  Time[0]
#                    LoopIndex2=LoopIndex2+1
#                    
#                    
#        
#        #Set nan to unused values
##        f.groups['Platform'].variables['distance'][0] = np.float32(0.0)
#        f.groups['Platform'].variables['speed_ground'][0] = np.nan
##        f.groups['Platform'].variables['speed_relative'][0] = np.float32(0.0)
##        f.groups['Platform'].variables['time4'][0] = np.int64(0)
#                    
#
#        
#        
#        
#        #Loop through every ping in raw fil
#        for LoopIndex0 in range(0,len(FileData.PingData)): 
#            
#            
#            #If there is information inside the ping continue
#            if not FileData.NMEA[LoopIndex0] ==[]:
#                
#                
#                #Get the time of the ping from nmea
#                Time = FileData.NMEA[LoopIndex0]
#                
#
#                #Get relevant information from raw data
#                sampleinterval = FileData.PingData[LoopIndex0].sampleinterval[0]
#
#                frequency = FileData.PingData[LoopIndex0].frequency[0]                  
#
#                transmitpower = FileData.PingData[LoopIndex0].transmitpower[0]       
#
#                pulslength = FileData.PingData[LoopIndex0].pulslength[0]     
#
#                bandwidth = FileData.PingData[LoopIndex0].bandwidth[0] 
#
#                heave = FileData.PingData[LoopIndex0].heave[0]
#
#                roll = FileData.PingData[LoopIndex0].roll[0]
#
#                pitch = FileData.PingData[LoopIndex0].pitch[0]
#        
#                transducer_gain = FileData.PingData[LoopIndex0].gainrx+FileData.PingData[LoopIndex0].gaintx
#                
#                transmitmode = FileData.PingData[LoopIndex0].transmitmode[0]
#
#                dirx = FileData.PingData[LoopIndex0].dirx
#
#                diry = FileData.PingData[LoopIndex0].diry
#
#                equivalentbeamangle = FileData.PingData[LoopIndex0].equivalentbeamangle
#
#                BeamAmplitudeData =  FileData.PingData[LoopIndex0].BeamAmplitudeData
#                
#                BeamAmplitudeData_imaginary =  FileData.PingData[LoopIndex0].BeamAmplitudeData_imaginary
#
#                absorption = FileData.PingData[LoopIndex0].absorptioncoefficient[0]
#
#                soundvelocity = FileData.PingData[LoopIndex0].soundvelocity[0]
#
#                BWalong = FileData.PingData[LoopIndex0].beamwidthalongshiprx
#                
#                BWatwh=FileData.PingData[LoopIndex0].beamwidthathwartshiprx
#
#                BWalongTX = FileData.PingData[LoopIndex0].beamwidthhorizontaltx[0]
#
#                BWatwhTX=FileData.PingData[LoopIndex0].beamwidthverticaltx[0]
#
#                NoiseFilter=FileData.PingData[LoopIndex0].noisefilter[0]
#
#               
#                #Convert beam direction into unit vectors
#                dirx_u = np.cos((dirx)*np.pi/180)*np.cos((diry+transducer_rotation[0])*np.pi/180)
#                
#                diry_u = np.cos((dirx)*np.pi/180)*np.sin((diry+transducer_rotation[0])*np.pi/180)
#                
#                dirz_u = np.sin((dirx)*np.pi/180)
#                
#                
#
#                #Get all unique frequencies and store it
#                if not frequency in Frequency_index: 
#
#                    Frequency_index = np.hstack((Frequency_index,frequency))
#                    
#                    
#                    #Write Environment information
#                    f.groups['Environment'].variables['frequency'][len(Frequency_index)-1]=frequency
#                    f.groups['Environment'].variables['absorption_indicative'][len(Frequency_index)-1]=absorption
#                    f.groups['Environment'].variables['sound_speed_indicative'][len(Frequency_index)-1]=soundvelocity
#
#
#                #Write platform data to group
#                addvar = f.groups['Platform'].variables['vertical_offset']                            
#                addvar[LoopIndex0]=  heave
#                
#                addvar = f.groups['Platform'].variables['roll']                            
#                addvar[LoopIndex0]=  roll
#                
#                addvar = f.groups['Platform'].variables['pitch']                            
#                addvar[LoopIndex0]=  pitch
#                
#                addvar = f.groups['Platform'].variables['time3']                            
#                addvar[LoopIndex0]=  Time[0]
#
#
#                #Find if the beam group is horizontal or vertical 
#                #This is a temporarly fix, other solution may be applied
#                if diry[0] ==2.8125:
#                    addvardir = f.groups['Sonar'].groups['Beam_group2']
#                    ping = horizontal_ping
#                
#                else: 
#                    addvardir = f.groups['Sonar'].groups['Beam_group1']
#                    ping = vertical_ping
#             
#                
#                #Bug fix to see if there is data avaliable in ping
#                if len(BeamAmplitudeData.shape)==2: 
#                    
#                    
#                    #Go though each beam
#                    for beamnum in range(0,64): 
#                        
#                        #Write beam data
#                        addvar = addvardir.variables['backscatter_r']
#                        addvar[ping,beamnum] = np.array(BeamAmplitudeData[:,beamnum])
#                        
#                        addvar = addvardir.variables['backscatter_i']
#                        addvar[ping,beamnum] = np.array(BeamAmplitudeData_imaginary[:,beamnum])
#                        
#                        
#                    
#                    #Write other sonar data
#                    addvar = addvardir.variables['equivalent_beam_angle']
#                    addvar[ping,:] = 10**(equivalentbeamangle/10)
#
#                    addvar = addvardir.variables['transducer_gain']
#                    addvar[ping,:] = transducer_gain
#
#                    addvar = addvardir.variables['beamwidth_receive_major']
#                    addvar[ping,:] = BWalong
#
#                    addvar = addvardir.variables['beamwidth_receive_minor']
#                    addvar[ping,:] = BWatwh
#
#                    addvar = addvardir.variables['beamwidth_transmit_major']
#                    addvar[ping,:] = BWalongTX
#
#                    addvar = addvardir.variables['beamwidth_transmit_minor']
#                    addvar[ping,:] = BWatwhTX
#                    
#                    addvar = addvardir.variables['transmit_bandwidth']
#                    addvar[ping] = bandwidth
#                    
#                    addvar = addvardir.variables['sample_interval']
#                    addvar[ping] = sampleinterval
#                    
#                    addvar = addvardir.variables['transmit_frequency_start']
#                    addvar[ping] = frequency
#                    
#                    addvar = addvardir.variables['transmit_frequency_stop']
#                    addvar[ping] = frequency
#                    
#                    addvar = addvardir.variables['transmit_duration_equivalent']
#                    addvar[ping] = pulslength
#                    
#                    addvar = addvardir.variables['transmit_duration_nominal']
#                    addvar[ping] = pulslength
#                    
#                    addvar = addvardir.variables['transmit_power']
#                    addvar[ping] =transmitpower
#                    
#                    addvar = addvardir.variables['sample_time_offset']
#                    addvar[ping] =0.004
#                    
#                    addvar = addvardir.variables['beam_type']
#                    addvar[ping] =0
#                    
#                    addvar = addvardir.variables['transmit_type']
#                    addvar[ping] =transmitmode
#                    
#                    addvar = addvardir.variables['ping_time']
#                    addvar[ping] = Time[0]
#                    
#                    addvar = addvardir.variables['beam_direction_vector']
#                    addvar[ping,:,0] = dirx_u
#                    addvar[ping,:,1] = diry_u
#                    addvar[ping,:,2] = dirz_u
#                    
#                    addvar=addvardir.variables['beam_stabilisation']
#                    addvar[ping] = 1
#                    
#                    addvar=addvardir.variables['non_quantitative_processing']
#                    addvar[ping] = NoiseFilter
#                    
#                    
#                    #For bookeeping
#                    #Use the same beam identificator as above
#                    if diry[0] ==2.8125:
#                        horizontal_ping = horizontal_ping + 1
#                        
#                    else: 
#                        vertical_ping = vertical_ping + 1
#                    
#                        
#                    
#        if horizontal_ping >=maxPingInFile:
#            NumberOfFilesInNC=0
#            Frequency_index = []
#        
#
#        #if last file to be added to .nc
#        if NumberOfFilesInNC >= MaxNumberOfFilesInNC: 
#            NumberOfFilesInNC=0
#        
#        if NumberOfFilesInNC == 0: 
#            f.close()
#            
#            
#


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
            print('Start converting')
            
            
        if LastFile == True: 
            os.chdir(directory)
            
            #Change to raw folder and display text to user
            print('Loading '+filename+'        ')
            
            
            #open .raw file
            fid = open(filename,mode='rb')
            
            
            #Read the .raw file
            FileData = ReadRawData.ReadRawData(fid,headerlength)
            
            print('File '+filename+' is loaded' , end='\r')
            
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
                            
#                    print(np.min(np.abs(NMEATime-int(str(Time)[8:]))))
                            
#                        time2 =  fulldate + datetime.timedelta(milliseconds=int(FileData.NMEA_time[0]/10000))
#                        time2 = str(time).replace("-","").replace(" ","").replace(":","").replace(".","")
#                        NMEATime = np.hstack((NMEATime,time2))
                    
#                    NMEATime = [str(time = fulldate + datetime.timedelta(milliseconds=int(FileData.NMEA_time[i]/10000))).replace("-","").replace(" ","").replace(":","").replace(".","") for i in range(len(FileData.NMEA_time[:]))]
                    
            
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
            
            print('File '+filename+' is converted' , end='\r')

            
