import os
from os.path import join
import sys
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import json
import neuron
import time as tm
import numpy as np

import LFPy

t0 = tm.clock()
## Choose model(s) to test:
#testmodel = 515175347 #Axonabundance (Hundreds of axon sections, branched dends, at least one part of the neuron is not connected to the rest) # AXON CONNECTED TO DENDRITE!!!
#testmodel = 501282059 #Originalplayer (Not everything is connected to the rest. 156 axon sec., 35 dend. sec.) # AXON CONNECTED TO DENDRITE!!!
#testmodel = 496497595 #Developmentcell (One short axon, branched dendrites, everything seems to be connected, no errors)
#testmodel = 497232392
#testmodel = 496538958
#testmodel = 497230463
#testmodel = 497233075
#testmodel = 497233271 # Better Cm's, but I worry this won't run properly ######
#testmodel = 488462965 # PERISOMATIC model of Developmentcell
#testmodel = 478513407 # Perisomatic version of 497233271, but only yields one peak
#testmodel = 480633479 # Perisomatic version of 497230463, only one peak ####
testmodel = 478513437 # Perisomatic version of 497233075
#testmodel = 489931686 # NEW perisomatic model
testmodelname = 'neur_%i' % testmodel
all_models    = [testmodelname]

if testmodel==480633479:
    v_init = -96.8#-83.7#-90#-86.5# # Have a test here too
elif testmodel==496497595:
    v_init = -86.5
elif testmodel==488462965:
    v_init = -86.5 # Maybe I should have changed this...
elif testmodel==497230463:
    v_init = -90
elif testmodel==497233075:
    v_init = -90
elif testmodel==478513437:
    v_init = -86.8
elif testmodel==478513407:
    v_init = -83.7
elif testmodel==497233271:
    v_init = -90
elif testmodel==489931686:
    v_init = -95.7
#v_init = -93.3

# Change this! Or come up with some automated solution.
idxs = [0,6,25,39,58]

# Defaulting to original values:
# DO NOT TOUCH THESE!
# SET THEM BELOW INSTEAD!
if testmodel==496497595:
    cm_soma = 1.14805
    cm_dend = 9.98231
    cm_axon = 3.00603
elif testmodel==497233271:
    cm_soma = 0.783229
    cm_dend = 1.94512
    cm_axon = 8.25387
elif testmodel==497230463:
    cm_soma = 1.23729
    cm_dend = 2.57923
    cm_axon = 5.02697
elif testmodel==497233075:
    cm_soma = 1.64168
    cm_dend = 2.83035
    cm_axon = 9.98442
elif testmodel==488462965:
    cm_soma = 3.31732779736 # Strange values, right?
    cm_dend = 3.31732779736
    cm_axon = 3.31732779736
elif testmodel==478513407:
    cm_soma = 1.0
    cm_dend = 1.0
    cm_axon = 1.0
elif testmodel==480633479:
    cm_soma = 0.704866 # 0.704866118957
    cm_dend = 0.704866 # 0.704866118957
    cm_axon = 0.704866 # 0.704866118957
elif testmodel==478513437:
    cm_soma = 2.34539964752
    cm_dend = 2.34539964752
    cm_axon = 2.34539964752
elif testmodel==489931686:
    cm_soma = 1.66244903951
    cm_dend = 1.66244903951
    cm_axon = 1.66244903951
    
# Changing values of membrane capacitance:
cm_changecm = 2.0

# Change current:
idur = 1000 # 200 #2 # 100 # ms # 1 #
iamp = 0 #0.41 #-0.5 # 1.0 # 0.2023 # -0.01 #  0.264 # 0.5 # nA 

idelay = 100 # 1 #     ms #
afteri = 100 # 20 # 1 #     ms # 

#
tstop_i = idur+afteri+idelay

def return_allen_cell_model(model_folder):

    params = json.load(open(join(model_folder, "fit_parameters.json"), 'r'))

    Ra = params["passive"][0]["ra"]

    e_pas = params["passive"][0]["e_pas"]
    celsius = params["conditions"][0]["celsius"]
    #cms = params["passive"][0]["cm"]
    reversal_potentials = params["conditions"][0]["erev"]
    v_init = params["conditions"][0]["v_init"]
    active_mechs = params["genome"]
    neuron.h.celsius = celsius
    # print(Ra, celsius, v_init)
    # print(reversal_potentials)
    # print(active_mechs)


    # Define cell parameters
    cell_parameters = {
        'morphology': join(model_folder, 'reconstruction.swc'),
        'v_init': -93.3,    # initial membrane potential
        'passive': False,   # turn on NEURONs passive mechanism for all sections
        'nsegs_method': 'fixed_length', # spatial discretization method
        'max_nsegs_length': 20.,
        #'lambda_f' : 200.,           # frequency where length constants are computed
        'dt': 2.**-7,      # simulation time step size
        'tstart': -100.,      # start time of simulation, recorders start at t=0
        'tstop': tstop_i,     # stop simulation at 100 ms.
        'custom_code': ['remove_axon.hoc']
    }

    cell = LFPy.Cell(**cell_parameters)
    cell.set_rotation(z=np.pi/1.25)
    
    pnncutoff = 0
    for sec in neuron.h.allsec(): # Extracting soma length
        sectype = sec.name().split("[")[0]
        if sectype=="soma":
            somasec = sec
            pnncutoff = 3.5*sec.L # What Sverre said

    for sec in neuron.h.allsec():
        sec.insert("pas")
        sec.e_pas = e_pas
        sec.Ra = Ra
        sectype = sec.name().split("[")[0]
        
        for sec_dict in active_mechs:
            if sec_dict["section"] == sectype:
                if sectype=="soma": # Works
                    exec("sec.cm = {}".format(cm_changecm))
                if sectype=="dend": # Works
                    ### TEST FOR DISTANCE # Weird to test from middle of soma? # Hard to do it any other way...
                    section = sec
                    dist = neuron.h.distance(somasec(0.5),section(1))
                    if dist<=pnncutoff:
                        exec("sec.cm = {}".format(cm_changecm))
                    else:
                        exec("sec.cm = {}".format(cm_dend))
                        for seg in sec:
                            if neuron.h.distance(somasec(0.5),seg)<=pnncutoff: # Have no idea if any of this will work...
                                exec("seg.cm = {}".format(cm_changecm)) # 
                                
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

    fig.savefig(join("figures", "%i" % testmodel,"current_idur%i_iamp" % idur + str(iamp),  'idur{}_iamp{}_{}_changecm{}_cmdend{}_cmaxon{}_vinit{}_addedRa.png'.format(idur, iamp, testmodel,cm_changecm,cm_dend,cm_axon,v_init)))
    
    # print out section information: # Works even though I do everything through LFPy
    for sec in neuron.h.allsec():
        neuron.h.psection()

    fig = plt.figure(figsize=[12, 8])
    plt.plot(cell.tvec, cell.vmem[0, :])
    plt.xlabel('Time (ms)')
    plt.ylabel('Potential (mV)')
    plt.title('Membrane potential in soma')
    plt.legend(loc='upper right')
    fig.savefig(join("figures", "%i" % testmodel,"current_idur%i_iamp" % idur + str(iamp),  'idur{}_iamp{}_{}_changecm{}_cmdend{}_cmaxon{}_vinit{}_addedRa_big.png'.format(idur, iamp, testmodel,cm_changecm,cm_dend,cm_axon,v_init)))
    
    '''
    fig = plt.figure(figsize=[12, 8])
    [plt.plot(cell.tvec, cell.vmem[idx, :], label='%i' % idx) for idx in idxs]
    plt.xlabel('Time (ms)')
    plt.ylabel('Potential (mV)')
    plt.title('Membrane potential along axon')
    plt.legend(loc='upper right')
    fig.savefig(join("figures", "%i" % testmodel,"current_idur%i_iamp" % idur + str(iamp), "axon", 'idur{}_iamp{}_{}_changecm{}_cmdend{}_cmaxon{}_vinit{}_ax_wRa.png'.format(idur, iamp,testmodel,cm_changecm,cm_dend,cm_axon,v_init)))
    
    fig = plt.figure(figsize=[12, 8])
    #plt.plot(cell.tvec, cell.rec_variables['cai'][0, :], label='Soma') # Soma is high
    [plt.plot(cell.tvec, cell.rec_variables['cai'][idx, :], label='Axon seg. %i' % idx) for idx in idxs] # Worth a shot # Is this even the axon? Should it not be the soma?
    plt.xlabel('Time (ms)')
    plt.ylabel(r'Ca$^{2+}$-concentration (mM)')
    plt.title('Ca$^{2+}$-concentration along axon')
    plt.legend(loc='lower right')
    fig.savefig(join("figures", "%i" % testmodel, "current_idur%i_iamp" % idur + str(iamp), "axon", "Ca", 'idur{}_iamp{}_{}_chcm{}_cmde{}_cmax{}_vinit{}_axCaRa.png'.format(idur, iamp,testmodel,cm_changecm,cm_dend,cm_axon,v_init)))
    
    fig = plt.figure(figsize=[12, 8])
    plt.plot(cell.tvec, cell.rec_variables['cai'][0, :])
    plt.xlabel('Time (ms)')
    plt.ylabel(r'Ca$^{2+}$-concentration (mM)')
    plt.title('Ca$^{2+}$-concentration in soma')
    fig.savefig(join("figures", "%i" % testmodel,"current_idur%i_iamp" % idur + str(iamp), 'Ca', 'idur{}_iamp{}_{}_changecm{}_cmdend{}_cmaxon{}_vinit{}_Ca_wRa.png'.format(idur, iamp, testmodel,cm_changecm,cm_dend,cm_axon,v_init)))
    '''
    
    outfilename = "figures/%i/current_idur%i_iamp" % (testmodel,idur) + str(iamp)+"/idur%i_iamp" % idur + str(iamp)+"_changecm" + str(cm_changecm) + "_cmdend" + str(cm_dend) + "_cmaxon"+ str(cm_axon) + "_vinit"+str(v_init)+"_addedRa.txt"
    outfile = open(outfilename,'w')
    vmem_soma = cell.vmem[0,:]#[idx, :] # Have time array too...
    Ca_soma   = cell.rec_variables['cai'][0, :] # Don't know if I need this...
    #print('vmem_soma[0]:',vmem_soma[0])
    
    time = cell.tvec # time
    Nt   = len(time)
    
    for i in range(Nt):
        outfile.write('%.16f %.16f %.16f\n' % (time[i],vmem_soma[i],Ca_soma[i]))
    outfile.close()    
    
    vmax = max(vmem_soma) 
    vmin = min(vmem_soma) 
    deltav = vmax-vmin
    vthr  = vmax-0.15*deltav # If there is a peak above this value, we count it
    vprev = vthr-40 # A peak never kicks in at initiation, even if I change vthr
    Npeaks = 0
    for i in range (1,len(vmem_soma)-1):  
        #print(vmem_soma[i])
        if vmem_soma[i-1]<vmem_soma[i] and vmem_soma[i+1]<vmem_soma[i] and vmem_soma[i]>vthr:
            Npeaks+=1
    print(Npeaks, ' peaks for model ', testmodel, ', current ', iamp)
    print('vmax:', vmax)

    outfilename = "figures/%i/current_idur%i_iamp" % (testmodel,idur) + str(iamp)+"/idur%i_iamp" % idur + str(iamp)+"_changecm" + str(cm_changecm) + "_cmdend" + str(cm_dend) + "_cmaxon"+ str(cm_axon) + "_vinit"+str(v_init)+"_addedRa_vmax.txt"
    outfile = open(outfilename,'w') 
    outfile.write('%.5f' % vmax)
    outfile.close()   
    
    outfilename = "figures/%i/current_idur%i_iamp" % (testmodel,idur) + str(iamp)+"/idur%i_iamp" % idur + str(iamp)+"_changecm" + str(cm_changecm) + "_cmdend" + str(cm_dend) + "_cmaxon"+ str(cm_axon) + "_vinit"+str(v_init)+"_addedRa_Npeaks.txt"
    outfile = open(outfilename,'w') 
    outfile.write('%i' % Npeaks)
    outfile.close()   
    
    #print('vmem_soma[0]:',vmem_soma[0])
    
    print('iamp:',iamp)
    t1 = tm.clock()
    print('Run time:', t1-t0)
    print('changecm:',cm_changecm, '; cm_dend:', cm_dend, '; cm_axon:', cm_axon)

    pnncutoff = 0
    for sec in neuron.h.allsec(): # Extracting soma length
        sectype = sec.name().split("[")[0]
        if sectype=="soma":
            somasec = sec
            pnncutoff = 3.5*sec.L # What Sverre said
    
    for sec in neuron.h.allsec(): # Testing
        sectype = sec.name().split("[")[0]
        if sectype=="soma":
            somasec = sec
        if sectype=="dend": # Works
            section = sec
            dist = neuron.h.distance(somasec(0.5),section(1))
            if dist<=pnncutoff:
                print(sec.name(), '; dist:',dist, '; Cm:',sec.cm)
            else:
                segnr = 0
                for seg in sec:
                    segnr+=1
                    dist2 = neuron.h.distance(somasec(0.5),seg)
                    print(sec.name(), '; distsegm:',dist2, '; segnr:',segnr, ' of ',sec.nseg,'; Cm:',seg.cm)
                print('---------------------')
                            
    sys.exit()
