# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 09:36:41 2018

@author: sindrev
"""


from graphviz import Digraph, Graph
import os

os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'



g = Digraph('G', filename='cluster')
g.format = 'tiff'





with g.subgraph(name='cluster_0') as c:
    c.attr(style='filled')
    c.attr(color='lightgrey')
    c.node_attr.update(style='filled', color='white')
    
    c.node('.raw data',url = 'http://bt.no')
    c.node('.dat data')
    c.node('read .raw',color = 'green')
    c.node('read .dat',color = 'red')
    
    c.node('convert to .nc',color = 'green')
    c.node('.nc',color = 'green')
    c.node('.nc (single ping)')
    
    
    c.edges([('read .raw','convert to .nc')])
    c.edges([('convert to .nc','.nc')])
    c.edges([('convert to .nc','.nc (single ping)')])
    c.edges([('.raw data','read .raw')])
    c.edges([('.dat data','read .dat')])
    c.edges([('read .dat','convert to .nc')])
    
    c.attr(label='RawToNetcdfConverter.py')

    


with g.subgraph(name='cluster_1') as c:
    c.attr(style='filled')
    c.attr(color='lightgrey')
    c.node_attr.update(style='filled', color='white')
    
    c.node('ACOUSTIC_DATA/SX/ORIGINAL_RAWDATA',color = 'green')
    c.node('ACOUSTIC_DATA/SX/RAWDATA',color = 'green')
    c.node('ACOUSTIC_DATA/SX/temp')
    
    
    c.attr(label='Structuring data')

    
    

with g.subgraph(name='cluster_2') as c:
    c.attr(style='filled')
    c.attr(color='lightgrey')
    c.node_attr.update(style='filled', color='white')
    
    c.node('Make Search Matrix',color = 'yellow')
    c.node('Make Work Files',color = 'yellow')
    c.node('Echo integrate -> nc',color = 'yellow')
    c.node('nc -> LUF20',color = 'green')
    
    
    c.attr(label='Data Processing')

    
    
    
    
    
    
    


g.edges([('.raw data','ACOUSTIC_DATA/SX/ORIGINAL_RAWDATA')])
g.edges([('.nc','ACOUSTIC_DATA/SX/RAWDATA')])
g.edges([('.nc (single ping)','ACOUSTIC_DATA/SX/temp')])

g.edges([('ACOUSTIC_DATA/SX/temp','Make Search Matrix')])
g.edges([('ACOUSTIC_DATA/SX/RAWDATA','Make Search Matrix')])
g.edges([('Make Search Matrix','Make Work Files')])
g.edges([('Make Work Files','Echo integrate -> nc')])
g.edges([('Make Work Files','.nc')])
g.edges([('Echo integrate -> nc','nc -> LUF20')])


    

g.render('cluster', view=True,cleanup=True)  