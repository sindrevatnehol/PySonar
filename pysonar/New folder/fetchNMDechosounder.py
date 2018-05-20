# -*- coding: utf-8 -*-
#import urllib
import urllib.request as req
from xml.dom import minidom   
#import sys
#sys.modules[__name__].__dict__.clear()
#
#
def geturl(string): #Forskningsfart%C3%B8y/2015/G%20O%20Sars-LMEL/2015106 

    with req.urlopen(string) as response:
       html = response.read()#.decode("utf-8") #.replace('b\'','')
   
       
       
    doc = minidom.parseString(html.decode("utf-8"))#("Output.xml")
    return doc

test = 'http://tomcat7-test.imr.no:8080/apis/nmdapi/cruise/v2/Leiefart%C3%B8y/2016/2016851/'
print('start scan')
string = 'http://tomcat7-staging.imr.no:8080/apis/nmdapi/reference/v2/dataset/cruiseseries?version=2.0'
string = 'http://tomcat7-staging.imr.no:8080/apis/nmdapi/reference/v2/dataset/cruiseseries?version=2.0'
doc = geturl(string)

    

#for row in doc.getElementsByTagName("code"): 
#    if int(row.firstChild.nodeValue) == 18:#20: 
#        for samples in doc.getElementsByTagName("samples"): 
#            for sample in samples.getElementsByTagName("sample"): 
#                for cruises in sample.getElementsByTagName("cruises"):
#                    for cruise in cruises.getElementsByTagName("cruise"):
#                        cruise_numer = cruise.getElementsByTagName("cruisenr")
#                        vessel = cruise.getElementsByTagName("shipName")
#                        
#                        print(cruise_numer[0].firstChild.nodeValue,vessel[0].firstChild.nodeValue)
#        for ii in i.getElementsByTagName('samples'): 
#            print(ii)



#for i in doc.getElementsByTagName("element"): 
#    missiontype = i.firstChild.nodeValue
#    missiontype = missiontype.replace('ø','%C3%B8')
#    
#    doc2 = geturl(string+missiontype)
#    
#        
#    for i2 in doc2.getElementsByTagName("element"): 
#        year = i2.firstChild.nodeValue
#        year = year.replace(' ','')
#        
#        if int(year)>2010: 
#            doc3 = geturl(string+missiontype+'/'+str(year))
#            
#            for i3 in doc3.getElementsByTagName("element"): 
#                if i3.getAttribute('name') == 'path': 
#                    vesselPath = i3.firstChild.nodeValue
#                    doc4 = geturl(string+missiontype+'/'+str(year)+'/'+vesselPath.replace(' ','%20').replace('ø','%C3%B8').replace('å','%C3%A5').replace('æ','%C3%A6').replace('Æ','%C3%86').replace('Ø','%C3%98'))
#                    
#                        
#                    for i4 in doc4.getElementsByTagName("element"): 
#                        delivery = str(i4.firstChild.nodeValue)
#                        print(string+missiontype+'/'+str(year)+'/'+vesselPath.replace(' ','%20').replace('ø','%C3%B8').replace('å','%C3%A5').replace('æ','%C3%A6').replace('Æ','%C3%86').replace('Ø','%C3%98')+'/'+delivery+'/dataset?version=2.0')
#                        doc5 = geturl(string+missiontype+'/'+str(year)+'/'+vesselPath.replace(' ','%20').replace('ø','%C3%B8').replace('å','%C3%A5').replace('æ','%C3%A6').replace('Æ','%C3%86').replace('Ø','%C3%98')+'/'+delivery+'/dataset?version=2.0')
            
            
        
#        missiontype = missiontype.replace('ø','%C3%B8')


#from rpy2.robjects.packages import importr
#import rpy2
#
#
##Load utils
#utils = importr('utils')
#
#
#
##install Rstox
##R_pack = ['data.table','ggplot2','pbapply','rgdal','rgeos','rJava','sp','XML',
##"ftp://ftp.imr.no/StoX/Download/Rstox/Rstox_1.9.tar.gz"]
##
##[utils.install_packages(i) for i  in R_pack]
#
##utils.install_packages("CRAN")
##utils.install_packages('http://cran.us.r-project.org/bin/windows/contrib/3.4/data.table_1.10.4-3.zip')
##
#utils.install_packages("data.table")
##data_table = importr('data.table')
##print(data_table)
#
#utils.install_packages("ftp://ftp.imr.no/StoX/Download/Rstox/Rstox_1.9.tar.gz")
#
##import RStox
#importr('Rstox')
#
