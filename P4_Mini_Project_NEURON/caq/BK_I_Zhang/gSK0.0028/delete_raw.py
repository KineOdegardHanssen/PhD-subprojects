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
    cm_factors = [1.0]#[0.5,1,1.25,1.5]
    iamps =  [0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.03,0.031,0.032,0.033,0.034,0.035,0.036,0.037,0.038,0.039,0.04,0.041,0.042,0.043,0.044,0.045,0.046,0.047,0.048,0.049,0.05,0.051,0.052,0.053,0.054,0.055,0.056,0.057,0.058,0.059,0.06]
    idur   = 1000 # 5000 # 
    idelay = 100
    tstop  = idur+idelay+10. #2.
    dtexp  = -7
    v_init = -86 #-78.16 # -70 #-60 # 
    t_before_rec = -600.
    Ncm = len(cm_factors)
    Ni  = len(iamps)
    
    somasize = 10
    
    model_folder = '' # We are in the model folder
    
    
    vary_naf    = False
    vary_kaf    = False
    vary_Ca_HVA = True # False
    vary_bk     = True # False
    vary_gpas   = False #
    
    gnaf    = 1.0
    gkaf    = 1.0
    gcahvas = [0.1,0.5,1.0,2,5.0,10]
    gbks    = [0.1,0.5,1.0,2,5.0,10]
    gpas    = 1.0
    
    for gcahva in gcahvas:
        for gbk in gbks:
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
        
            if vary_naf==True:
                namestring = namestring + '_gnaf'+str(gnaf)+'p'
            if vary_kaf==True:
                namestring = namestring + '_gkaf'+str(gkaf)+'p'
            if vary_bk==True:
                namestring = namestring + '_gbk'+str(gbk)+'p'
            if vary_Ca_HVA==True:
                namestring = namestring + '_gCaHVA'+str(gcahva)+'p'
            if vary_gpas==True: 
                namestring = namestring + '_gpas'+str(gpas)+'p'
            namestring = namestring +'_'
            
            print('iamps:',iamps)
            for i in range(Ni):        
                iamp = iamps[i]
                folder = 'Results/Soma%i/current_idur%i_iamp'%(somasize,idur)+str(iamp)+'/'
                folder_total = folder
                ###folder_total = 'C:/Users/Kine/Documents/Projects_PhD/P4_PNN_channels/Allen/'+folder
                if os.path.exists(folder_total)==False:
                    os.mkdir(folder_total)
                print('Step ', i+1, ' out of ', Ni)
                for j in range(Ncm):
                    cm_factor = cm_factors[j]
                    # Save results:
                    filename = folder+namestring+'somaonly_cm'+str(cm_factor)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'_V.txt'
                    try:
                        os.remove(filename)    
                    except:
                        kine = 'me'
