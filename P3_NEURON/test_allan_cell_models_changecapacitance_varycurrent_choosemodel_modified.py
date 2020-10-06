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

# Change this! Or come up with some automated solution.
idxs = [0,6,25,39,58]

# Defaulting to original values:
# DO NOT TOUCH THESE!
# SET THEM BELOW INSTEAD!
cm_soma = 1.14805
cm_dend = 9.98231
cm_axon = 3.00603

# Changing values of membrane capacitance:
#cm_soma = 0.01
#cm_dend = 0.01
#cm_axon = 10.0

# Change current:
idur = 1000 # ms
iamp = 0.41 # nA

idelay = 100 # ms # 1
afteri = 100 # ms # 20

#
tstop_i = idur+afteri+idelay

def return_allen_cell_model(model_folder):

    params = json.load(open(join(model_folder, "fit_parameters.json"), 'r'))

    celsius = params["conditions"][0]["celsius"]
    reversal_potentials = params["conditions"][0]["erev"]
    v_init = params["conditions"][0]["v_init"]
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
    
    #idxs = cell.get_idx(section="axon")
    
    plt.close("all")
    fig = plt.figure(figsize=[12, 8])
    fig.subplots_adjust(wspace=0.5)
    ax1 = fig.add_subplot(231, aspect=1, xlabel="x", ylabel="z")
    ax2 = fig.add_subplot(234, aspect=1, xlabel="x", ylabel="y")
    ax3 = fig.add_subplot(132, xlabel="ms", ylabel="mV")
    ax4 = fig.add_subplot(133, xlabel="ms", ylabel="nA")

    [ax1.plot([cell.xstart[idx], cell.xend[idx]],
              [cell.zstart[idx], cell.zend[idx]], c='k') for idx in range(cell.totnsegs)]

    [ax2.plot([cell.xstart[idx], cell.xend[idx]],
              [cell.ystart[idx], cell.yend[idx]], c='k') for idx in range(cell.totnsegs)]

    ax3.plot(cell.tvec, cell.somav)

    ax4.plot(cell.tvec, stimulus.i)

    fig.savefig(join("figures", "%i" % testmodel,"current_idur%i_iamp" % idur + str(iamp),  'idur{}_iamp{}_{}_{}_cmsoma{}_cmdend{}_cmaxon{}_addedRa.png'.format(idur, iamp, model_idx, model_name,cm_soma,cm_dend,cm_axon)))
    
    # print out section information: # Works even though I do everything through LFPy
    #for sec in neuron.h.allsec():
    #    neuron.h.psection()

    fig = plt.figure(figsize=[12, 8])
    [plt.plot(cell.tvec, cell.vmem[idx, :], label='%i' % idx) for idx in idxs]
    plt.xlabel('Time (ms)')
    plt.ylabel('Potential (mV)')
    plt.title('Membrane potential along axon')
    plt.legend(loc='upper right')
    fig.savefig(join("figures", "%i" % testmodel,"current_idur%i_iamp" % idur + str(iamp), "axon", 'idur{}_iamp{}_{}_{}_cmsoma{}_cmdend{}_cmaxon{}_axon_addedRa.png'.format(idur, iamp, model_idx, model_name,cm_soma,cm_dend,cm_axon)))
    
    fig = plt.figure(figsize=[12, 8])
    #plt.plot(cell.tvec, cell.rec_variables['cai'][0, :], label='Soma') # Soma is high
    [plt.plot(cell.tvec, cell.rec_variables['cai'][idx, :], label='Axon seg. %i' % idx) for idx in idxs] # Worth a shot # Is this even the axon? Should it not be the soma?
    plt.xlabel('Time (ms)')
    plt.ylabel(r'Ca$^{2+}$-concentration (mM)')
    plt.title('Ca$^{2+}$-concentration along axon')
    plt.legend(loc='lower right')
    fig.savefig(join("figures", "%i" % testmodel, "current_idur%i_iamp" % idur + str(iamp), "axon", 'idur{}_iamp{}_{}_{}_cmsoma{}_cmdend{}_cmaxon{}_axon_Ca_addedRa.png'.format(idur, iamp, model_idx, model_name,cm_soma,cm_dend,cm_axon)))
    
    fig = plt.figure(figsize=[12, 8])
    plt.plot(cell.tvec, cell.rec_variables['cai'][0, :])
    plt.xlabel('Time (ms)')
    plt.ylabel(r'Ca$^{2+}$-concentration (mM)')
    plt.title('Ca$^{2+}$-concentration in soma')
    fig.savefig(join("figures", "%i" % testmodel,"current_idur%i_iamp" % idur + str(iamp), 'idur{}_iamp{}_{}_{}_cmsoma{}_cmdend{}_cmaxon{}_Ca_addedRa.png'.format(idur, iamp, model_idx, model_name,cm_soma,cm_dend,cm_axon)))
    
    outfilename = "figures/%i/current_idur%i_iamp" % (testmodel,idur) + str(iamp)+"/idur%i_iamp" % idur + str(iamp)+"_cmsoma" + str(cm_soma) + "_cmdend" + str(cm_dend) + "_cmaxon"+ str(cm_axon) + "_addedRa.txt"
    outfile = open(outfilename,'w')
    vmem_soma = cell.vmem[0,:]#[idx, :] # Have time array too...
    Ca_soma   = cell.rec_variables['cai'][0, :] # Don't know if I need this...
    print('vmem_soma:',vmem_soma)
    
    time = cell.tvec # time
    Nt   = len(time)
    
    for i in range(Nt):
        outfile.write('%.16f %.16f %.16f\n' % (time[i],vmem_soma[i],Ca_soma[i]))
    outfile.close()    
    
    sys.exit()
