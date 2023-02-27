import os
from os.path import join
import sys
import matplotlib.pyplot as plt
import json
import neuron
import time as tm
import numpy as np

import LFPy

if __name__ == '__main__':
    Vshifts = [-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50]
    cm_factor = 1.0 
    iamps = [0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.03,0.031,0.032,0.033,0.034,0.035,0.036,0.037,0.038,0.039,0.04,0.041,0.042,0.043,0.044,0.045,0.046,0.047,0.048,0.049,0.05,0.051,0.052,0.053,0.054,0.055,0.056,0.057,0.058,0.059,0.06]
    idur   = 1000 # 5000 # 
    idelay = 100
    tstop  = idur+idelay+10. #2.
    dtexp  = -7
    v_init = -86 #-78.16 # -70 #-60 # 
    t_before_rec = -600.
    Ni  = len(iamps)
    
    somasize = 10
    
    model_folder = '' # We are in the model folder
    
    varyE = 0
    varymech = 'None' # 'Epas' # 'EK' # 'ENa' # 
    namestring = ''
    if varymech=='ENa':
        varyE = 63 #[40,50,60,70]
        namestring = namestring + 'ENa'+str(varyE)
    elif varymech=='EK':
        varyE = -97#-107
        namestring = namestring + 'EK'+str(varyE)
    elif varymech=='Epas': 
        varyE = -20 # Vary by shifts
        namestring = namestring + 'Epasshift'+str(varyE)
    
    vary_naf    = False
    vary_kaf    = False
    vary_canin  = True # False
    vary_gpas   = False #
    
    gnaf    = 1.0
    gkaf    = 1.0
    gcanins = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.25,1.5,2.0,2.5,3.0,5.0,10.0]
    gpas    = 1.0
    Ng = len(gcanins)
    
    
    print('iamps:',iamps)
    for Vshift in Vshifts:
        for i in range(Ni):
            iamp = iamps[i]
            folder = 'Vshift_%i/Results/Soma%i/current_idur%i_iamp'%(Vshift,somasize,idur)+str(iamp)+'/'
            folder_total = folder
            ###folder_total = 'C:/Users/Kine/Documents/Projects_PhD/P4_PNN_channels/Allen/'+folder
            if os.path.exists(folder_total)==False:
                os.mkdir(folder_total)
            print('Step ', i+1, ' out of ', Ni)
            for j in range(Ng):
                gcanin = gcanins[j]
                
                namestring = ''
                if vary_naf==True:
                    namestring = namestring + '_gnaf'+str(gnaf)+'p'    
                if vary_kaf==True:
                    namestring = namestring + '_gkaf'+str(gkaf)+'p'    
                if vary_canin==True:
                    namestring = namestring + '_gcanin'+str(gcanin)+'p'
                if vary_gpas==True: 
                    namestring = namestring + '_gpas'+str(gpas)+'p'    
                namestring = namestring +'_'
                # Delete results:
                filename = folder+namestring+'somaonly_cm'+str(cm_factor)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'_V.txt'
                try:
                    os.remove(filename)    
                except:
                    kine = 'me'
    
