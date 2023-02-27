import os
from os.path import join
import sys
import matplotlib.pyplot as plt
import json
import neuron
import time as tm
import numpy as np

import LFPy

def return_allen_cell_model(model_folder, cm_factor, somasize, idur,iamp,idelay,tstop,dtexp,v_init,t_before_rec):

    params = json.load(open(join(model_folder, "fit_parameters.json"), 'r'))

    Ra = params["passive"][0]["ra"]
    
    e_pas = params["passive"][0]["e_pas"]
    celsius = params["conditions"][0]["celsius"]
    cm = params["passive"][0]["cm"]
    reversal_potentials = params["conditions"][0]["erev"]
    v_init = params["conditions"][0]["v_init"]
    active_mechs = params["genome"]
    Ra = params["passive"][0]["ra"] 
    neuron.h.celsius = celsius
    # print(Ra, celsius, v_init)
    # print(reversal_potentials)
    # print(active_mechs)
    cm = cm[0]["cm"]
    
    print('cm:',cm)
    print('e_pas:',e_pas)

    mydt = 2.**dtexp
    # Define cell parameters
    cell_params = {
        'morphology': "somaonly.hoc",
        'v_init': v_init,    # initial membrane potential # Was -86 mV
        'passive': False,   # turn on NEURONs passive mechanism for all sections
        'nsegs_method': 'fixed_length', # spatial discretization method
        'max_nsegs_length': 20.,
        'dt': mydt,      # simulation time step size
        'tstart': t_before_rec,    # start time of simulation, recorders start at t=0
        'tstop': tstop,     # stop simulation at 100 ms.
    }
    
    cell = LFPy.Cell(**cell_params)
    
    
    for sec in neuron.h.allsec():
        sec.insert("pas")
        sectype = sec.name().split("[")[0]
        for sec_dict in active_mechs:
            if sec_dict["section"] == sectype:
                if sectype=="soma": # Works
                    sec.cm   = cm*cm_factor
                    sec.L    = somasize
                    sec.diam = somasize
                if not sec_dict["mechanism"] == "":
                    sec.insert(sec_dict["mechanism"])
                exec("sec.{} = {}".format(sec_dict["name"], sec_dict["value"]))
    
    stim_idx = 0
    stim_params = {
                'idx': stim_idx,
                'record_current': True,
                'pptype': 'IClamp',
                'amp': iamp,
                'dur': idur,
                'delay': idelay,
            }
    
    stimulus = LFPy.StimIntElectrode(cell, **stim_params)
    #syn.set_spike_times(np.array([1]))
    cell.simulate(rec_vmem=True, rec_imem=True)
    t, v = cell.tvec.copy(), cell.vmem[0].copy()
    
    if cm_factor==1 and iamp==0.02:
        for sec in neuron.h.allsec():
            neuron.h.psection()
    
    cell.__del__()
    return t, v


if __name__ == '__main__':
    cm_factors = [1.0]#[0.5,1,1.25,1.5]
    iamps = [0.01]# [0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.03,0.031,0.032,0.033,0.034,0.035,0.036,0.037,0.038,0.039,0.04,0.041,0.042,0.043,0.044,0.045,0.046,0.047,0.048,0.049,0.05,0.051,0.052,0.053,0.054,0.055,0.056,0.057,0.058,0.059,0.06]#[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2]#[0.4,0.5,0.6,0.7,0.8,0.9,1.0]#[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2]#,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8] #[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2]#[0.4,0.5,0.6,0.7,0.8,0.9,1.0]#[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2]#,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8]
    idur   = 100 #1000
    idelay = 10
    tstop  = idur+idelay+20. #2.
    dtexp  = -7
    v_init = -86 #-78.16 # -70 #-60 # 
    t_before_rec = -600.
    Ncm = len(cm_factors)
    Ni  = len(iamps)
    
    somasize = 10
    
    model_folder = '' # We are in the model folder
    
    print('iamps:',iamps)
    for i in range(Ni):
        iamp = iamps[i]
        folder = 'Results/Soma%i/current_idur%i_iamp'%(somasize,idur)+str(iamp)+'/'
        folder_total = folder
        if os.path.exists(folder_total)==False:
            os.mkdir(folder_total)
        print('Step ', i+1, ' out of ', Ni)
        for j in range(Ncm):
            cm_factor = cm_factors[j]
            t, v = return_allen_cell_model(model_folder,cm_factor,somasize,idur,iamp,idelay,tstop,dtexp,v_init,t_before_rec)
            
            # Save results:
            filename = folder+'somaonly_cm'+str(cm_factor)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'_V.txt'
            plotname   = folder+'somaonly_cm'+str(cm_factor)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'_V.png'
            
            Nt = len(t)
            file = open(filename,'w')
            for k in range(Nt):
                file.write('%.16e %.16e\n' % (t[k],v[k]))
            file.close()   
            
            '''
            plt.figure(figsize=(6,5))
            plt.plot(t,v)
            plt.xlabel('Time (ms)')
            plt.ylabel('Membrane potential (mV)')
            plt.title('Voltage vs time, %.2f*Cm' % cm_factor)
            plt.savefig(plotname)
            plt.show()'''
