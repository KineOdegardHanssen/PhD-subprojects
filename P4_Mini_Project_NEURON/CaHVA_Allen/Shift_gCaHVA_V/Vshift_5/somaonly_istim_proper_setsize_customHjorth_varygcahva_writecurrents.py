import os
from os.path import join
import sys
import matplotlib.pyplot as plt
import json
import neuron
import time as tm
import numpy as np

import LFPy

def return_allen_cell_model(model_folder, cm_factor, somasize, idur,iamp,idelay,tstop,dtexp,gnaf,gkaf,gcahva,gpas,vary_naf,vary_kaf,vary_Ca_HVA,vary_gpas,varyE,varymech,v_init,t_before_rec):

    params = json.load(open(join(model_folder, "fit_parameters.json"), 'r'))

    Ra = params["passive"][0]["ra"]
    
    e_pas = params["passive"][0]["e_pas"]
    celsius = params["conditions"][0]["celsius"]
    cm = params["passive"][0]["cm"]
    reversal_potentials = params["conditions"][0]["erev"]
    #v_init = params["conditions"][0]["v_init"]
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
        'v_init': v_init,    # initial membrane potential # -86.8 was the default before
        'passive': False,   # turn on NEURONs passive mechanism for all sections
        'nsegs_method': 'fixed_length', # spatial discretization method
        'max_nsegs_length': 20.,
        'dt': mydt,      # simulation time step size
        'tstart': t_before_rec,    # start time of simulation, recorders start at t=0 #
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
                    print('sec_dict[mechanism]:',sec_dict["mechanism"])
                    sec.insert(sec_dict["mechanism"])
                exec("sec.{} = {}".format(sec_dict["name"], sec_dict["value"]))
    
    for sec in neuron.h.allsec():
        sectype = sec.name().split("[")[0]
        # First: Change mechanisms 
        # Mechanisms that can be anywhere
        if varymech=='Epas':
            sec.e_pas += varyE
        if sectype=='soma':
            if varymech=='ENa':
                sec.ena = varyE
            elif varymech=='EK':
                sec.ek = varyE
            if vary_naf==True:
                sec.gbar_naf *= gnaf
            if vary_kaf==True:
                sec.gbar_kaf *= gkaf
            if vary_Ca_HVA==True:
                sec.gbar_Ca_HVA *= gcahva
            if vary_gpas==True: 
                sec.g_pas *= gpas
    
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
    cell.simulate(rec_variables=["ica_Ca_HVA","ina_naf","ik_kaf","cai"],rec_vmem=True,  rec_imem=True, rec_icap=True, rec_ipas=True)
    print('cell.stimulate() on')
    t, v = cell.tvec.copy(), cell.vmem[0].copy()
    
    # Operate with current names
    I_Ca_HVA = cell.rec_variables["ica_Ca_HVA"][0]
    I_naf = cell.rec_variables["ina_naf"][0]
    I_kaf = cell.rec_variables["ik_kaf"][0]
        
    print('sum(I_Ca_HVA):',sum(I_Ca_HVA))
    print('sum(I_naf):',sum(I_naf))
    print('sum(I_SK):',sum(I_kaf))
    
    if cm_factor==1 and iamp==0.02:
        for sec in neuron.h.allsec():
            neuron.h.psection()
    
    cell.__del__()
    return t, v, I_naf, I_kaf, I_Ca_HVA, cell.icap[0], cell.imem[0], cell.ipas[0]


if __name__ == '__main__':
    vshift = 0
    cm_factor = 1.0 
    iamps = [0.002,0.005,0.01,0.015,0.02,0.04]
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
    vary_Ca_HVA = True # False
    vary_gpas   = False #
    
    gnaf    = 1.0
    gkaf    = 1.0
    gcahvas = [0.1,0.5,1.0,1.5,2.0]
    gpas    = 1.0
    Ng = len(gcahvas)
    
    
    print('iamps:',iamps)
    for i in range(Ni):
        iamp = iamps[i]
        folder = 'Results/Soma%i/current_idur%i_iamp'%(somasize,idur)+str(iamp)+'/'
        folder_total = folder
        ###folder_total = 'C:/Users/Kine/Documents/Projects_PhD/P4_PNN_channels/Allen/'+folder
        if os.path.exists(folder_total)==False:
            os.mkdir(folder_total)
        print('Step ', i+1, ' out of ', Ni)
        for j in range(Ng):
            gcahva = gcahvas[j]
            
            namestring = ''
            if vary_naf==True:
                namestring = namestring + '_gnaf'+str(gnaf)+'p'
            if vary_kaf==True:
                namestring = namestring + '_gkaf'+str(gkaf)+'p'
            if vary_Ca_HVA==True:
                namestring = namestring + '_gCaHVA'+str(gcahva)+'p'
            if vary_gpas==True: 
                namestring = namestring + '_gpas'+str(gpas)+'p'
            namestring = namestring +'_'
            t, v, I_naf, I_kaf, I_Ca_HVA, icap, imem, ipas = return_allen_cell_model(model_folder,cm_factor,somasize,idur,iamp,idelay,tstop,dtexp,gnaf,gkaf,gcahva,gpas,vary_naf,vary_kaf,vary_Ca_HVA,vary_gpas,varyE,varymech,v_init,t_before_rec)
            
            # Save results:
            filename = folder+namestring+'somaonly_cm'+str(cm_factor)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'_V.txt'
            currname = folder+namestring+'somaonly_cm'+str(cm_factor)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'_currents.txt'
            plotname = folder+namestring+'somaonly_cm'+str(cm_factor)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'V.png'
            currplotname = folder+namestring+'somaonly_cm'+str(cm_factor)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'V.png'
            
            Nt = len(t)
            file = open(filename,'w')
            currfile = open(currname,'w')
            for k in range(Nt):
                file.write('%.16e %.16e\n' % (t[k],v[k]))
                currfile.write('%.16e %.16e %.16e %.16e %.16e %.16e %.16e\n' % (t[k],I_naf[k],I_kaf[k],I_Ca_HVA[k],icap[k],imem[k],ipas[k]))
            file.close()   
            
            plt.figure(figsize=(6,5))
            plt.plot(t,v)
            plt.xlabel('Time (ms)')
            plt.ylabel('Membrane potential (mV)')
            plt.title('Voltage vs time, %.2f*Cm' % cm_factor)
            plt.savefig(plotname)
            
            plt.figure(figsize=(6,5))
            plt.plot(t,imem, label=r'$I_{\mathregular{mem}}$')
            plt.plot(t,icap, label=r'$I_{\mathregular{cap}}$')
            plt.plot(t,ipas, label=r'$I_{\mathregular{pas}}$')
            plt.plot(t,I_naf, label=r'$I_{\mathregular{naf}}$')
            plt.plot(t,I_kaf, label=r'$I_{\mathregular{kaf}}$')
            plt.plot(t,I_Ca_HVA, label=r'$I_{\mathregular{CaHVA}}$')
            plt.xlabel('Time (ms)')
            plt.ylabel('Current (nA)')
            plt.title(r'Currents vs time, $V_{\mathregular{shift}}$=%i' % vshift)
            plt.savefig(currplotname)
            #plt.show()
            
            print('v[-1]:',v[-1])
    