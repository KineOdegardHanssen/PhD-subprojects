import os
from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import json
import neuron
import LFPy
import math

def run_sim(dendlen=1000, cm=1.0, idur=1.0, iamp=1.0, idelay=1.0, v_init=-86.5, tstart=0., tstop=50., Ra=150, testmodel=496497595,nsegments=65):
    
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
                    if not sec_dict["mechanism"] == "": # and sectype=="soma" # Merge with the one above.
                        sec.insert(sec_dict["mechanism"])
                    exec("sec.{} = {}".format(sec_dict["name"], sec_dict["value"]))
                if sectype=="dend":
                    sec.L = dendlen
                    sec.nseg = nsegments

        for sec_dict in reversal_potentials:
            if sec_dict["section"] == sectype:
                # print(sectype, sec_dict)
                for key in sec_dict.keys():
                    if not key == "section":
                        if sectype=="soma":
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
    testmodel = 478513437#478513407#488462965# ### Weirdo:480633479#  ###AA:496497595
    smallCm = True
    dendlens = [10,100,500,1000,2000,5000,10000]
    if smallCm==True:
        cms = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]#[]
    else:
        cms = [0.5,1,2,3,4,4.5,5,5.5,6,7,8,9,10]
    iamp = 0.6
    idur = 1000
    idelay  = 100.0
    afteri  = 100.0
    tstart  = -500.
    tstop   = idur+afteri+idelay
    v_init  = -86.5
    Ra      = 150
    outfolder = 'Results/%i/IStim/Vary_dendrite_length/' % testmodel
    currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    outfolder = outfolder+currentfolder
    
    #lambda = 180 #np.sqrt(1000*2e-4/(4*150))
    maxlen = 15. # Improving the accuracy somewhat. 8.3333..% of lambda. About the accuracy of the original model
     
    for dendlen in dendlens:
        plt.figure(figsize=(6,5))
        outfolder_dl = outfolder + 'dendlen%i/' % dendlen
        plotname = outfolder_dl+'basps_varycm_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit'%Ra+str(v_init)+'_dendlen%i_V.png'% dendlen
        nsegments = int(math.ceil(dendlen/maxlen))
        for cm in cms:
            t, v = run_sim(dendlen,cm, idur, iamp, idelay, v_init, tstart, tstop, Ra, testmodel,nsegments)
            outfilename = outfolder_dl+'basps_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_dendlen%i_V.txt' % dendlen
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