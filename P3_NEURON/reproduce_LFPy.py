import os
from os.path import join
import sys
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import json
import neuron
import h5py
import numpy as np

import LFPy


def return_allen_cell_model(model_file, morph_file):
    print('Calling json')
    params = json.load(open(model_file, 'r'))
    print('json called')

    Ra = params["passive"][0]["ra"]

    e_pas = params["passive"][0]["e_pas"]
    celsius = params["conditions"][0]["celsius"]
    cms = params["passive"][0]["cm"]
    reversal_potentials = params["conditions"][0]["erev"]
    v_init = params["conditions"][0]["v_init"]
    active_mechs = params["genome"]
    neuron.h.celsius = celsius
    # print(Ra, celsius, v_init)
    # print(reversal_potentials)
    # print(active_mechs)


    # Define cell parameters
    cell_parameters = {
        'morphology' : morph_file,
        'v_init' : -80,    # initial membrane potential
        'passive' : False,   # turn on NEURONs passive mechanism for all sections
        # 'nsegs_method' : 'fixed_length', # spatial discretization method
        'max_nsegs_length': 20.,
        'lambda_f' : 200.,           # frequency where length constants are computed
        'dt' : 0.1,      # simulation time step size
        'tstart' : 0,      # start time of simulation, recorders start at t=0
        'tstop' : 1200.,     # stop simulation at 2000 ms.
        'custom_code': [join("..", 'remove_axon.hoc')]
    }

    print('Setting cell') # My line
    cell = LFPy.Cell(**cell_parameters)
    print('Cell set')
    # cell.set_rotation(z=np.pi/1.25)

    for sec in neuron.h.allsec():
        sec.insert("pas")
        sec.e_pas = e_pas
        sec.Ra = Ra
        sectype = sec.name().split("[")[0]

        for cm_dict in cms:
            if cm_dict["section"] == sectype:
                exec("sec.cm = {}".format(cm_dict["cm"]))

        for sec_dict in active_mechs:
            if sec_dict["section"] == sectype:
                # print(sectype, sec_dict)
                if not sec_dict["mechanism"] == "":

                    if not sec.has_membrane(sec_dict["mechanism"]):
                        print('sec_dict["mechanism"]:', sec_dict["mechanism"]) # My line
                        sec.insert(sec_dict["mechanism"])
                        # print("Inserted ", sec_dict["mechanism"])
                exec("sec.{} = {}".format(sec_dict["name"], sec_dict["value"]))

        for sec_dict in reversal_potentials:
            if sec_dict["section"] == sectype:
                # print(sectype, sec_dict)
                for key in sec_dict.keys():
                    if not key == "section":
                        exec("sec.{} = {}".format(key, sec_dict[key]))

    # for sec in neuron.h.allsec():
    #     if hasattr(sec, "eca"):
    #         print(sec.cao)
    #         sec.cao = 2
            # print(sec.name(), sec.eca)
    return cell


cell_models_folder = join("sim_ch01", "components", "biophysical_neuron_models")


model_id = "491623973"#"496497595"#"472363762"
model_file = join(cell_models_folder, "%s_fit.json" % model_id)
morph_file = join("sim_ch01", "components", "morphologies", 'Pvalb_490387590_m.swc')
#morph_file = join("sim_ch01", "components", "morphologies", 'Pvalb_487667205_m.swc')
#morph_file = join("sim_ch01", "components", "morphologies", 'Scnn1a_473845048_m.swc')

model_idx = 0

mod_folder = join("sim_ch01", "components", "mechanisms")

if "win64" in sys.platform: # Not sure any of this works...
    print('Detected sys.platform as win64')
    warn("no autompile of NMODL (.mod) files on Windows. " 
         + "Run mknrndll from NEURON bash in the folder cells and rerun example script")
    if not pth in neuron.nrn_dll_loaded:
        neuron.h.nrn_load_dll(mod_folder+"/nrnmech.dll")
    neuron.nrn_dll_loaded.append(mod_folder)
''' # This works?
if "win32" in sys.platform:
    print('Detected sys.platform as win32')
    print("no autompile of NMODL (.mod) files on Windows. " 
         + "Run mknrndll from NEURON bash in the folder cells and rerun example script")
    if not pth in neuron.nrn_dll_loaded:
        neuron.h.nrn_load_dll(mod_folder+"/nrnmech.dll")
    neuron.nrn_dll_loaded.append(mod_folder)
'''
#print('Loading NEURON mechanisms')
#neuron.load_mechanisms(mod_folder)
print('NEURON mechanisms loaded')


cell = return_allen_cell_model(model_file, morph_file)


current_stim_params = {
    'idx': 0,
    'pptype': 'IClamp',
    'amp': 0.61,
    'dur': 1000.0,
    'delay': 100.0,

}

stim = LFPy.StimIntElectrode(cell, **current_stim_params)

cell.simulate(rec_imem=True, rec_vmem=True)


plt.close("all")
#fig = plt.figure(figsize=[6, 6])
fig = plt.figure(figsize=[10, 6])
fig.subplots_adjust(wspace=0.5)

ax_v = fig.add_subplot(111, ylabel="soma membrane potential (mV)", xlabel="time (ms)", ylim=[-100, 45])

allen_v_file = join("sim_ch01", "output", "v_report.h5")
v_allen = h5py.File(allen_v_file, 'r')["report"]["mcortex"]["data"]

l1, = ax_v.plot(cell.tvec, cell.somav)
l2, = ax_v.plot(cell.tvec[:-1], v_allen, c='gray', ls=':')
fig.legend([l1, l2], ["LFPy", "BMTK"])

fig.savefig(join('reproduce_allen_{}.png'.format(model_id)))

# Count spikes:
Nspikes = 0
vprev = 0
for i in range(len(cell.somav)):
    v = cell.somav[i]
    if v>-40 and vprev<-40:
            Nspikes+=1
    vprev = v
print(Nspikes)
print(max(cell.somav))
