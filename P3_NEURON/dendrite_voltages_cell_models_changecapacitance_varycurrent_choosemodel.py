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

## Choose model(s) to test:
#testmodel = 515175347 #Axonabundance (Hundreds of axon sections, branched dends, at least one part of the neuron is not connected to the rest) # AXON CONNECTED TO DENDRITE!!!
#testmodel = 501282059 #Originalplayer (Not everything is connected to the rest. 156 axon sec., 35 dend. sec.) # AXON CONNECTED TO DENDRITE!!!
testmodel = 496497595 #Developmentcell (One short axon, branched dendrites, everything seems to be connected, no errors)
#testmodel = 497232392
#testmodel = 496538958
testmodelname = 'neur_%i' % testmodel
all_models    = [testmodelname]

v_init = -86.5#-90#

# Change this! Or come up with some automated solution.
dendnr = 6
idxs = [37,38,39,40,41,42,43,44,45,46]
Nidx = len(idxs)

# Defaulting to original values:
# DO NOT TOUCH THESE!
# SET THEM BELOW INSTEAD!
cm_soma = 1.14805
cm_dend = 9.98231
cm_axon = 3.00603

# Changing values of membrane capacitance:
#cm_soma = 2.0
#cm_dend = 15.0
#cm_axon = 0.01

# Change current:
idur = 1 #100    #1000 # ms #100 #
iamp = 1.0 #0.271  #0.41 # nA #0.5 #

idelay = 1  #100 # ms # 1
afteri = 10#20 #100 # ms # 20

#
tstop_i = idur+afteri+idelay

def return_allen_cell_model(model_folder):

    params = json.load(open(join(model_folder, "fit_parameters.json"), 'r'))

    celsius = params["conditions"][0]["celsius"]
    reversal_potentials = params["conditions"][0]["erev"]
    #v_init = params["conditions"][0]["v_init"] # This might be wrong, will set it by hand above
    active_mechs = params["genome"]
    Ra = params["passive"][0]["ra"]  # This was not here before
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
        'tstop' : tstop_i,     # stop simulation at idur+20 ms.
        'Ra':Ra,
        # 'custom_code': ['remove_axon.hoc']
    }

    cell = LFPy.Cell(**cell_parameters)
    cell.set_rotation(z=np.pi/1.25)

    for sec in neuron.h.allsec():
        sec.insert("pas")
        sectype = sec.name().split("[")[0]
        for sec_dict in active_mechs:
            if sec_dict["section"] == sectype:
                if sectype=="soma": # Works
                    sec.cm = cm_soma
                if sectype=="dend": # Works
                    sec.cm = cm_dend
                if sectype=="axon": # Works
                    sec.cm = cm_axon
                if not sec_dict["mechanism"] == "":
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

#all_models = [f for f in os.listdir(cell_models_folder) if f.startswith("neur")
#              and os.path.isdir(join(cell_models_folder, f))]
mod_folder = "Allen_test_changecapacitance/cell_models/"+testmodelname+"/modfiles"
if "win64" in sys.platform:
    print('Detected sys.platform as win64')
    warn("no autompile of NMODL (.mod) files on Windows. " 
         + "Run mknrndll from NEURON bash in the folder cells and rerun example script")
    if not pth in neuron.nrn_dll_loaded:
        neuron.h.nrn_load_dll(mod_folder+"/nrnmech.dll")
    neuron.nrn_dll_loaded.append(mod_folder)

for model_idx in range(len(all_models)):
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
            'delay': idelay,
        }
    stimulus = LFPy.StimIntElectrode(cell,**pointprocess)
    cell.simulate(rec_vmem=True, rec_variables=['cai'])
    
    # print out section information: # Works even though I do everything through LFPy
    for sec in neuron.h.allsec():
        neuron.h.psection()
    
    time = cell.tvec # time
    Nt   = len(time)
    
    folder = "figures/%i/current_idur%i_iamp" % (testmodel,idur) + str(iamp)+"/dendritepropagation/"
    for j in range(Nidx):
        outfilename = folder+"idur%i_iamp" % idur + str(iamp)+"_cmsoma" + str(cm_soma) + "_cmdend" + str(cm_dend) + "_cmaxon"+ str(cm_axon) + "_vinit"+str(v_init)+"_addedRa_V_dendsec%i_dend%i.txt" % (j,dendnr)
        outfile = open(outfilename,'w')
        vmem = cell.vmem[idxs[j],:]   
        
        for i in range(Nt):
            outfile.write('%.16f %.16f\n' % (time[i],vmem[i]))
        outfile.close()  
    plotname = folder+"idur%i_iamp" % idur + str(iamp)+"_cmsoma" + str(cm_soma) + "_cmdend" + str(cm_dend) + "_cmaxon"+ str(cm_axon) + "_vinit"+str(v_init)+"_addedRa_V_dend%i.png" % dendnr 
    fig = plt.figure(figsize=[12, 8])
    [plt.plot(cell.tvec, cell.vmem[idx, :], label='%i' % idx) for idx in idxs]
    plt.xlabel('Time (ms)')
    plt.ylabel('Potential (mV)')
    plt.title('Membrane potential along dendrite %i' % dendnr)
    plt.legend(loc='upper right')
    plt.savefig(plotname)  

    sys.exit()
