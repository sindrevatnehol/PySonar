import numpy as np
import scipy.io as sc

    
def processinnput(Work,x__school_traj,X_traj,y__school_traj,Y_traj,case,corr): 
    
    if len(case['BananaTool'])==1:
        BananaTool = case['BananaTool'][0]
    else: 
        BananaTool = case['BananaTool']


    for ii in range(len(x__school_traj)):   
        R = (x__school_traj[ii]-X_traj)**2+(y__school_traj[ii]-Y_traj)**2
        if ii == 0: 
            r_test = R
        else: 
            r_test = np.minimum(r_test,R)
    
    x_i,y_i,z_i = np.where(np.sqrt(r_test)<=(int(BananaTool[0])+corr))
    Work[x_i,y_i,z_i]=1



    return Work
        
    
    
    
    
def MakeIndex(case,Work,x__school_traj_port,y__school_traj_port,X_traj,Y_traj,
              x__school_traj_stb,y__school_traj_stb,filename): 
    '''Make index file '''
    
    
    #Correction function when the distance of beams increases
    
    corr =  case['MakeRange']*np.sin(2*np.pi/64)

    print('    Start (1/2)')
    Work = processinnput(Work,x__school_traj_port,X_traj,y__school_traj_port,Y_traj,case,corr)

    
    print('    Start (2/2)')
    Work = processinnput(Work,x__school_traj_stb,X_traj,y__school_traj_stb,Y_traj,case,corr)

    #Save the work file for future use
    sc.savemat(filename,mdict={'Work':Work})
    
    

    