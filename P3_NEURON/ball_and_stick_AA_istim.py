import os
from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import json
import neuron
import LFPy

def run_sim(cm=1.0, idur=1.0, iamp=1.0, idelay=1.0, v_init=-86.5, tstart=0., tstop=50., Ra=150, testmodel=496497595):
    
    ### Copy-pasted, more or less, from the test_allan_cell_models-script
    testmodelname = 'neur_%i' % testmodel
    model_folder = 'cell_models/%s/' % testmodelname
    params = json.load(open(join(model_folder, "fit_parameters.json"), 'r'))
    
    celsius = params["conditions"][0]["celsius"]
    reversal_potentials = params["conditions"][0]["erev"]
    #v_init = params["conditions"][0]["v_init"] # This might be wrong, will set it by hand above
    active_mechs = params["genome"]
    Ra = params["passive"][0]["ra"]  # This was not here before
    neuron.h.celsius = celsius
    
    dt = 2.**-4
    cell_params = {          # various cell parameters,
                'morphology': "ballandstick.hoc",
                'v_init' : v_init,    # initial crossmembrane potential
                'cm': cm,
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
        sectype = sec.name().split("[")[0]
        for sec_dict in active_mechs:
            if sec_dict["section"] == sectype:
                if sectype=="soma": # Works
                    sec.cm = cm
                if sectype=="dend": # Works
                    sec.cm = cm
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

    stimulus = LFPy.StimIntElectrode(cell, **stim_params)
    cell.simulate(rec_vmem=True, rec_imem=True)
    t, v = cell.tvec.copy(), cell.vmem[0].copy()
    
    # Try to get information on compartments:
    # print out section information: # Works even though I do everything through LFPy
    if cm==1:
        for sec in neuron.h.allsec():
            neuron.h.psection()
        idxs = cell.get_idx(section="dend")
        print('idxs:',idxs)
    
    cell.__del__()
    return t, v


if __name__ == '__main__':
    testmodel = 496497595
    cms = [0.5,1,2,3,4,4.5,5,5.5,6,7,8,9,10]
    iamp = -0.5
    idur = 100
    idelay  = 100.0
    afteri  = 100.0
    tstart  = -500.
    tstop   = idur+afteri+idelay
    v_init  = -86.5
    Ra      = 150
    outfolder = 'Results/%i/IStim/' % testmodel
    currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    outfolder = outfolder+currentfolder
    plotname = outfolder+'basaa_varycm_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit'%Ra+str(v_init)+'_V.png' 
    for cm in cms:
        print('cm:',cm)
        t, v = run_sim(cm, idur, iamp, idelay, v_init, tstart, tstop, Ra, testmodel)
        outfilename = outfolder+'basaa_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_V.txt' 
        outfile = open(outfilename,'w')
        # Write to 'Results/IStim/'

        Nt   = len(t)
        for i in range(Nt):
            outfile.write('%.16f %.16f\n' % (t[i],v[i]))
        outfile.close()    
        plt.plot(t, v, label="c$_m$= {:1.2f} µF/cm²".format(cm))
    plt.xlabel('Time [ms]')
    plt.ylabel('Voltage [mV]')
    plt.title(r'Voltage vs time for different cell capacitances $c_m$ (current input)')
    plt.legend()
    plt.savefig(plotname)