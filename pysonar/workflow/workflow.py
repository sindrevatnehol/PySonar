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






with g.subgraph(name='cluster_00') as c:
    c.attr(style='filled')
    c.attr(color='lightgrey')
    c.node_attr.update(style='filled', color='white')
    
    c.node('Initialize')
    c.node('Load timeseries info \n from NMDAPI',color = 'red')
    c.node('Make folder structure',color = 'green')
    c.node('Copy .raw/.dat data',color = 'green')
    c.node('Copy LUF20',color = 'red')
    
    
    
    
    c.attr(label='Initialization')

    
    
    

with g.subgraph(name='cluster_0') as c:
    c.attr(style='filled')
    c.attr(color='lightgrey')
    c.node_attr.update(style='filled', color='white')
    
    c.node('read .raw',color = 'green')
    c.node('read .dat',color = 'red')
    
    c.node('convert to .nc',color = 'green')
    c.node('.nc',color = 'green')
    
    
    
    c.attr(label='RawToNetcdfConverter.py')

    


with g.subgraph(name='cluster_1') as c:
    c.attr(style='filled')
    c.attr(color='lightgrey')
    c.node_attr.update(style='filled', color='white')
    
    c.node('.../ORIGINAL_RAWDATA',color = 'green')
    c.node('.../RAWDATA',color = 'green')
    c.node('.../psonar/ekluf20',color = 'green')
    c.node('.../psonar/result',color = 'green')
    
    c.attr(label='Folder structure')

    
with g.subgraph(name='cluster_2') as c:
    c.attr(style='filled')
    c.attr(color='lightgrey')
    c.node_attr.update(style='filled', color='white')
    
    c.node('Get list of files on each transect',color = 'green')
    c.node('pysonar.makesearch()',color = 'yellow')
    c.node('serchmatrix.mat',color = 'yellow')
    
    c.attr(label='Make Search Matrix')

    
    
    
    
    
    
    
    
    
g.edges([('Initialize','Load timeseries info \n from NMDAPI')])
g.edges([('Load timeseries info \n from NMDAPI','Make folder structure')])

g.edges([('.../psonar/ekluf20','Get list of files on each transect')])
g.edges([('Get list of files on each transect','pysonar.makesearch()')])
g.edges([('.../RAWDATA','pysonar.makesearch()')])
g.edges([('pysonar.makesearch()','serchmatrix.mat')])
    
    
    
g.edges([('read .raw','convert to .nc')])
g.edges([('convert to .nc','.nc')])
g.edges([('.../ORIGINAL_RAWDATA','read .raw')])
g.edges([('.../ORIGINAL_RAWDATA','read .dat')])
g.edges([('read .dat','convert to .nc')])
    
    

g.edges([('Make folder structure','Copy .raw/.dat data')])
g.edges([('Make folder structure','Copy LUF20')])
g.edges([('Copy LUF20','.../psonar/ekluf20')])

g.edges([('Copy .raw/.dat data','.../ORIGINAL_RAWDATA')])
g.edges([('.nc','.../RAWDATA')])

g.edges([('serchmatrix.mat','.../psonar/result')])
g.edges([('Make Work Files','Echo integrate -> nc')])
g.edges([('Make Work Files','.nc')])


g.edges([('.../psonar/result','Make Work Files')])
g.edges([('.../RAWDATA','Make Work Files')])


g.edges([('Echo integrate -> nc','nc -> LUF20')])
g.edges([('.../RAWDATA','nc -> LUF20')])


    

g.render('cluster', view=True,cleanup=True)  