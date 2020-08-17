import os
from os.path import join
import sys
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import json
import neuron
import numpy as np

import LFPy

# Defaulting to original values:
# DO NOT TOUCH THESE!
# SET THEM BELOW INSTEAD!
cm_soma = 1.14805
cm_dend = 9.98231
cm_axon = 3.00603

# Changing values of membrane capacitance:
cm_soma = 0.01
cm_dend = 0.01
cm_axon = 0.01

# Change current:
idur = 100 # ms
iamp = 0.5 # nA

#
tstop_i = idur+20.

ca_on = False # Switch to not generate unneccessary plots

def return_allen_cell_model(model_folder):

    params = json.load(open(join(model_folder, "fit_parameters.json"), 'r'))

    celsius = params["conditions"][0]["celsius"]
    reversal_potentials = params["conditions"][0]["erev"]
    v_init = params["conditions"][0]["v_init"]
    active_mechs = params["genome"]
    neuron.h.celsius = celsius
    # print(Ra, celsius, v_init)
    # print(reversal_potentials)
    # print(active_mechs)


    # Define cell parameters
    cell_parameters = {
        'morphology' : join(model_folder, 'reconstruction.swc'),
        'v_init' : v_init,    # initial membrane potential
        'passive' : False,   # turn on NEURONs passive mechanism for all sections
        'nsegs_method' : 'lambda_f', # spatial discretization method
        'lambda_f' : 200.,           # frequency where length constants are computed
        'dt' : 2.**-5,      # simulation time step size
        'tstart' : -100.,      # start time of simulation, recorders start at t=0
        'tstop' : tstop_i,     # stop simulation at 100 ms.
        # 'custom_code': ['remove_axon.hoc']
    }

    cell = LFPy.Cell(**cell_parameters)
    cell.set_rotation(z=np.pi/1.25)

    for sec in neuron.h.allsec():
        sec.insert("pas")
        sectype = sec.name().split("[")[0]
        for sec_dict in active_mechs:
            if sec_dict["section"] == sectype:
                #print(sectype, sec_dict)
                if sectype=="soma": # Works
                    #print('sectype==soma')
                    sec.cm = cm_soma
                if sectype=="dend": # Works
                    #print('sectype==dend')
                    sec.cm = cm_dend
                if sectype=="axon": # Works
                    #print('sectype==axon')
                    sec.cm = cm_axon
                if not sec_dict["mechanism"] == "":
                    #print(sectype, sec_dict)
                    sec.insert(sec_dict["mechanism"])
                exec("sec.{} = {}".format(sec_dict["name"], sec_dict["value"]))

        for sec_dict in reversal_potentials:
            if sec_dict["section"] == sectype:
                # print(sectype, sec_dict)
                for key in sec_dict.keys():
                    if not key == "section":
                        exec("sec.{} = {}".format(key, sec_dict[key]))
    return cell

cell_models_folder = "cell_models" #"cell_models_Arkhipov_et_al" #"all_cell_models"

all_models = [f for f in os.listdir(cell_models_folder) if f.startswith("neur")
              and os.path.isdir(join(cell_models_folder, f))]
mod_folder = "Allen_test" #"modfiles" #"Allen_test/modfiles" #"all_mods"
print("mod_folder:", mod_folder)
if "win64" in sys.platform:
    print('Detected sys.platform as win64')
    warn("no autompile of NMODL (.mod) files on Windows. " 
         + "Run mknrndll from NEURON bash in the folder cells and rerun example script")
    if not pth in neuron.nrn_dll_loaded:
        neuron.h.nrn_load_dll(mod_folder+"/nrnmech.dll")
    neuron.nrn_dll_loaded.append(mod_folder)

############## Only holds for cell-id 496497595 #############
N_plotdendrites = 5
dendrite_id     = ['A','D','E','K','M']
dendrite_idxs   = [[0,4,6,8,11],[0,21,23,29,35],[0,36,41,46],[0,84,87],[0,99]]
#############################################################

for model_idx in range(len(all_models)): # Change this?
    model_name = all_models[model_idx]
    #print(model_idx, model_name)
    model_folder = join("cell_models", model_name)

    cell = return_allen_cell_model(model_folder)

    pointprocess = {
            'idx': 0,
            'record_current': True,
            'pptype': 'IClamp',
            'amp': iamp,
            'dur': idur,
            'delay': 1,
        }
    stimulus = LFPy.StimIntElectrode(cell,**pointprocess)
    cell.simulate(rec_vmem=True, rec_variables=['cai'])
    
    for i in range(N_plotdendrites): # NB NB NB!: Put current info in the file name #NAMES NOT CHANGED! Put (custom) dendrite name on them too
        idxs = dendrite_idxs[i]
        this_dendrite_id = dendrite_id[i]

        fig = plt.figure(figsize=[12, 8])
        [plt.plot(cell.tvec, cell.vmem[idx, :], label='%i' % idx) for idx in idxs]
        plt.xlabel('Time (ms)')
        plt.ylabel('Potential (mV)')
        plt.title('Membrane potential along dendrite %s' % this_dendrite_id)
        plt.legend(loc='upper right')
        fig.savefig(join("figures", "%s" % model_name, "dend", "dend%s" % this_dendrite_id, '{}_{}_idur{}_iamp{}_cmsoma{}_cmdend{}_cmaxon{}_dend{}.png'.format(model_idx, model_name,idur,iamp,cm_soma,cm_dend,cm_axon,this_dendrite_id)))
        
        # NB! No [Ca2+]-recording (no cai) in dendrites for cell with cell-id 496497595
        if ca_on==True:
            fig = plt.figure(figsize=[12, 8])
            [plt.plot(cell.tvec, cell.rec_variables['cai'][idx, :], label='Dendrite seg. %i' % idx) for idx in idxs]
            plt.xlabel('Time (ms)')
            plt.ylabel(r'Ca$^{2+}$-concentration (mM)')
            plt.title('Ca$^{2+}$-concentration along dendrite %s' % this_dendrite_id)
            plt.legend(loc='lower right')
            fig.savefig(join("figures", "%s" % model_name, "dend", "dend%s" % this_dendrite_id, '{}_{}_idur{}_iamp{}_cmsoma{}_cmdend{}_cmaxon{}_dend{}_Ca.png'.format(model_idx, model_name,idur,iamp,cm_soma,cm_dend,cm_axon,this_dendrite_id)))
    
    sys.exit()
