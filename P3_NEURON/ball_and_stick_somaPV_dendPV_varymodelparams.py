import os
from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import json
import neuron
import LFPy

def run_sim(varymech, varyE, varyg, cm=1.0, idur=1.0, iamp=1.0, idelay=1.0, v_init=-86.5, tstart=0., tstop=50., Ra=150, testmodel=496497595):
    
    ### Copy-pasted, more or less, from the test_allan_cell_models-script
    testmodelname = 'neur_%i' % testmodel
    model_folder = 'cell_models/%s/' % testmodelname
    params = json.load(open(join(model_folder, "fit_parameters.json"), 'r'))
    
    e_pas = params["passive"][0]["e_pas"]
    celsius = params["conditions"][0]["celsius"]
    reversal_potentials = params["conditions"][0]["erev"]
    #v_init = params["conditions"][0]["v_init"] # This might be wrong, will set it by hand above
    active_mechs = params["genome"]
    Ra = params["passive"][0]["ra"]  # This was not here before
    neuron.h.celsius = celsius
    
    dt = 2.**-6 #-4
    cell_params = {          # various cell parameters,
                'morphology': "ballandstick.hoc",
                'v_init' : v_init,    # initial crossmembrane potential
                'cm': 1.0,
                'Ra': Ra,
                'passive': False,
                'nsegs_method' : "lambda_f",
                'lambda_f': 200.,
                'dt' : dt,   # [ms] dt's should be in powers of 2
                'tstart' : tstart,
                'tstop' : tstop,
            }


    cell = LFPy.Cell(**cell_params)
    
    for sec in neuron.h.allsec():
        sec.insert("pas")
        sec.e_pas = e_pas
        sec.Ra = Ra
        sectype = sec.name().split("[")[0]
        
        for sec_dict in active_mechs:
            if sec_dict["section"] == sectype:
                if sectype=="soma": # Works
                    exec("sec.cm = {}".format(cm))
                if sectype=="axon": # Works
                    exec("sec.cm = {}".format(cm_axon))
                if not sec_dict["mechanism"] == "":
                    sec.insert(sec_dict["mechanism"])
                exec("sec.{} = {}".format(sec_dict["name"], sec_dict["value"]))

        for sec_dict in reversal_potentials:
            if sec_dict["section"] == sectype:
                # print(sectype, sec_dict)
                for key in sec_dict.keys():
                    if not key == "section":
                        exec("sec.{} = {}".format(key, sec_dict[key]))

    stim_idx = 0
    stim_params = {
            'idx': stim_idx,
            'record_current': True,
            'pptype': 'IClamp',
            'amp': iamp,
            'dur': idur,
            'delay': idelay,
        }    
    
    # Vary properties:
    for sec in neuron.h.soma:
        if varymech=='NaV':
            sec.ena = varyE
            sec.gbar_NaV = varyg
        elif varymech=='pas':
            sec.e_pas = varyE
            sec.g_pas = varyg
        elif varymech=='Kd':
            sec.ek = varyE
            sec.gbar_Kd = varyg
        elif varymech=='Kv2like':
            sec.ek = varyE
            sec.gbar_Kv2like = varyg
        elif varymech=='Kv3_1':
            sec.ek = varyE
            sec.gbar_Kv3_1 = varyg
        elif varymech=='SK':
            sec.ek = varyE
            sec.gbar_SK = varyg
        elif varymech=='K_T':
            sec.ek = varyE
            sec.gbar_K_T = varyg
        elif varymech=='Im_v2':
            sec.ek = varyE
            sec.gbar_Im_v2 = varyg



    stimulus = LFPy.StimIntElectrode(cell, **stim_params)
    cell.simulate(rec_vmem=True, rec_imem=True)
    t, v = cell.tvec.copy(), cell.vmem[0].copy()
    
    # Try to get information on compartments:
    # print out section information: # Works even though I do everything through LFPy
    if cm==1:
        for sec in neuron.h.allsec():
            neuron.h.psection()
            for seg in sec:
                print('seg.cm:',seg.cm, '; seg.x:', seg.x, 'seg.diam:', seg.diam, 'seg.g_pas:',seg.g_pas, 'seg.e_pas:',seg.e_pas)
 
        idxs = cell.get_idx(section="dend")
        print('idxs:',idxs)
    
    cell.__del__()
    return t, v


if __name__ == '__main__':
    # varymech: Difficult to make this make sense for other ions than Na (and leak)...
    # 478513407: gbar_Kv2like=0.0854285
    varymech = 'NaV' #'SK'#'Im_v2'#'Kv2like'#'Kd'#'pas' #
    varyE = -53 # Default...407Kd-107 Default478513407pas: -83.6528;  Default478513407Na: 53
    varyg = 0.09 # Default...407Kd2.94396e-010 # Default478513407pas: 0.000362109; Default478513407Na: 0.0409177
    if varymech=='NaV':
        folderstring = 'VaryNa/' 
    elif varymech=='pas':
        folderstring = 'VaryPas/'
    elif varymech=='Kd':
        folderstring = 'VaryKd/'
    elif varymech=='Kv2like':
        folderstring = 'VaryKv2like/'
    elif varymech=='Kv3_1':
        folderstring = 'VaryKv3_1/'
    elif varymech=='SK':
        folderstring = 'VarySK/'
    elif varymech=='K_T':
        folderstring = 'VaryK_T/'
    elif varymech=='Im_v2':
        folderstring = 'VaryIm_v2/'
    #
    testmodel = 478513407# 488462965#478513437# ### Weirdo:480633479#  ###AA:496497595
    smallCm = True
    cm = 1.0
    iamps = [0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18]
    idur = 1000
    idelay  = 100.0
    afteri  = 100.0
    tstart  = -500.
    tstop   = idur+afteri+idelay
    v_init  = -86.5
    Ra      = 150
    outfolder_base = 'Results/%i/IStim/' % testmodel + folderstring
    for iamp in iamps:
        print('iamp:',iamp)
        t, v = run_sim(varymech,varyE,varyg,cm, idur, iamp, idelay, v_init, tstart, tstop, Ra, testmodel)
        currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
        outfolder = outfolder_base+currentfolder
        outfilename = outfolder+'basps_cm'+str(cm)+'_E'+str(varyE)+'_g'+str(varyg)+'_Ra%i_vinit' %Ra+str(v_init)+'_V.txt' 
        plotname = outfolder+'basps_cm'+str(cm)+'_E'+str(varyE)+'_g'+str(varyg)+'_Ra%i_vinit' %Ra+str(v_init)+'_V.png' 
        outfile = open(outfilename,'w')
        # Write to 'Results/IStim/'

        Nt   = len(t)
        for i in range(Nt):
            outfile.write('%.16f %.16f\n' % (t[i],v[i]))
        outfile.close()
        plt.figure(figsize=(6,5))
        plt.plot(t, v)
        plt.xlabel('Time [ms]')
        plt.ylabel('Voltage [mV]')
        plt.title(r'Voltage vs time (current input)')
        plt.savefig(plotname)