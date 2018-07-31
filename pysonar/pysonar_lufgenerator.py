# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 15:03:21 2018

@author: sindrev
"""

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
      
      


def LUFICEStemplate(report_file): 
    k=1

    
    
    
    

def LUF20template(directory2Data): 
    import os
    import scipy.io as sc
    from xml.etree import ElementTree as ET
    
    
    
    #Loop through each work file
    first = True
    
    
    #Loop through all files
    for report in os.listdir(directory2Data.dir_result+'/Vertical/'): 
        
        
        #Load report file
        report_file = sc.loadmat(directory2Data.dir_result+'/Vertical/'+report)
        
        
        #If this is the first file
        if first == True:
            first = False
        
            #Start the xml document
            root = ET.Element("echosounder_dataset")
            ET.SubElement(root,'report_time').text = report_file['Report_time'][0]
            ET.SubElement(root,'lsss_version').text = report_file['DataProcessingMethod']['SoftwareName'][0][0][0]+' V'+report_file['DataProcessingMethod']['SoftwareVersion'][0][0][0]
            ET.SubElement(root,'nation').text =  report_file['Cruice']['Nation'][0][0][0]
            ET.SubElement(root,'platform').text =  report_file['Cruice']['Platform'][0][0][0]
       
        
            distancelist = ET.SubElement(root,'distancelist')
    
            
            
        #Start on the distance list and add informatinon of each log distance
        distance = ET.SubElement(distancelist,'distance')
        distance.set('log_start',report_file['Cruice']['Log'][0][0]['Distance'][0][0][0])
        distance.set('start_time',report_file['Cruice']['Log'][0][0]['TimeStart'][0][0][0])
        ET.SubElement(distance,'integrator_dist').text =  str(report_file['Cruice']['Log'][0][0]['Sample'][0][0]['PingAxisInterval'][0][0][0][0])
        ET.SubElement(distance,'pel_ch_thickness').text =  str(report_file['Cruice']['Log'][0][0]['Sample'][0][0]['ChannelThickness'][0][0][0][0])
        ET.SubElement(distance,'include_estimate').text =  '1'
        ET.SubElement(distance,'lat_start').text =  str(report_file['Cruice']['Log'][0][0]['Latitude_start'][0][0][0])
        ET.SubElement(distance,'lon_start').text =  str(report_file['Cruice']['Log'][0][0]['Longitude_start'][0][0][0])
        ET.SubElement(distance,'lat_stop').text = str(report_file['Cruice']['Log'][0][0]['Latitude_stop'][0][0][0])
        ET.SubElement(distance,'lon_stop').text =str(report_file['Cruice']['Log'][0][0]['Longitude_stop'][0][0][0])
        ET.SubElement(distance,'stop_time').text =  str(report_file['Cruice']['Log'][0][0]['TimeStop'][0][0][0])
        
        
        
        #Start on the frequency list and add all the frequency stuff
        frequency = ET.SubElement(distance,'frequency')
        frequency.set('frequency',str(int(report_file['Instrument']['Frequency'][0][0][0][0])))
        frequency.set('tranceiver',report_file['Instrument']['TransducerSerial'][0][0][0])
        
        ET.SubElement(frequency,'quality').text =  '2'
        ET.SubElement(frequency,'bubble_corr').text =  '0'
        ET.SubElement(frequency,'threshold').text =str(report_file['Cruice']['Log'][0][0]['Sample'][0][0]['SvThreshold'][0][0])
        ET.SubElement(frequency,'num_pel_ch').text =  str(len(report_file['Cruice']['Log'][0][0]['Sample'][0][0]['DataValue'][0][0]))
        ET.SubElement(frequency,'upper_interpret_depth').text =  '0'
        ET.SubElement(frequency,'upper_integrator_depth').text =  '0'
        if (len(report_file['Cruice']['Log'][0][0]['Sample'][0][0]['DataValue'][0][0])*report_file['Cruice']['Log'][0][0]['Sample'][0][0]['ChannelThickness'][0][0][0][0])<abs(report_file['Cruice']['Log'][0][0]['Bottom_depth'][0][0][0][0]):
            ET.SubElement(frequency,'lower_interpret_depth').text = str(len(report_file['Cruice']['Log'][0][0]['Sample'][0][0]['DataValue'][0][0])*report_file['Cruice']['Log'][0][0]['Sample'][0][0]['ChannelThickness'][0][0][0][0])
            ET.SubElement(frequency,'lower_integrator_depth').text =  str(len(report_file['Cruice']['Log'][0][0]['Sample'][0][0]['DataValue'][0][0])*report_file['Cruice']['Log'][0][0]['Sample'][0][0]['ChannelThickness'][0][0][0][0])
            
        else: 
            ET.SubElement(frequency,'lower_interpret_depth').text = str(abs(report_file['Cruice']['Log'][0][0]['Bottom_depth'][0][0][0][0]))
            ET.SubElement(frequency,'lower_integrator_depth').text =  str(abs(report_file['Cruice']['Log'][0][0]['Bottom_depth'][0][0][0][0]))
        
        
    
            
            
            
        #Start on the acoustic data stuff
        #The acocat must be generalised
        chtype = ET.SubElement(frequency,'ch_type')
        chtype.set('type','P')
        sa_by_acocat = ET.SubElement(chtype,'ch_type')
        sa_by_acocat.set('acocat','12')
        
        
        #Print all NASC values but skipp 0-s
        sa = report_file['Cruice']['Log'][0][0]['Sample'][0][0]['DataValue'][0][0]
        for i in range(len(sa)): 
            if sa[i][0]>0: 
                sa_value = ET.SubElement(sa_by_acocat,'sa')
                sa_value.set('ch',str(i+i))
                sa_value.text = str(sa[i][0])
        
    
    #Fix the xml structure
    indent(root)
    
    
    #Write the report file
    tree = ET.ElementTree(root)
    tree.write(directory2Data.dir_result+'/ListUserFile20_SU90_vertical.xml')
                
    

    

    
    

def reader(directory2Data, LUF20 = False, LUFICES = False):
    
    
    if LUF20 == True: 
        LUF20template(directory2Data)
    if LUFICES == True:         
        LUFICEStemplate(directory2Data)

