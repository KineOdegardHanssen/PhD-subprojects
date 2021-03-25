import os
from os.path import join
import sys
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import json
import neuron
import numpy as np
import math

import LFPy

plotit = True # Bool for plotting
## Choose model(s) to test:
#testmodel = 515175347 #Axonabundance (Hundreds of axon sections, branched dends, at least one part of the neuron is not connected to the rest) # AXON CONNECTED TO DENDRITE!!!
#testmodel = 501282059 #Originalplayer (Not everything is connected to the rest. 156 axon sec., 35 dend. sec.) # AXON CONNECTED TO DENDRITE!!!
testmodel = 496497595 #Developmentcell (One short axon, branched dendrites, everything seems to be connected, no errors)
#testmodel = 497232392
#testmodel = 496538958
testmodel = 488462965 # PERISOMATIC model of Developmentcell
#testmodel = 480633479

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

# Defaulting to original values:
# DO NOT TOUCH THESE!
# SET THEM BELOW INSTEAD!
if testmodel==496497595:
    cm_soma = 1.14805
    cm_dend = 9.98231
    cm_axon = 3.00603
    Ndendsecs = 18 # 24? Plot and neuron.psection()-output different.
elif testmodel==488462965:
    cm_soma = 3.31732779736 # Strange values, right?
    cm_dend = 3.31732779736
    cm_axon = 3.31732779736
elif testmodel==480633479:
    cm_soma = 0.704866 # 0.704866118957
    cm_dend = 0.704866 # 0.704866118957
    cm_axon = 0.704866 # 0.704866118957


# Changing values of membrane capacitance:
cm_soma = 0.1
#cm_dend = cm_soma # 15.0
#cm_axon = cm_soma # 0.01

# Change current:
idur = 2 # 1 #100    #1000 # ms #100 #
iamp = 1.0 #0.271  #0.41 # nA #0.5 #

idelay = 1  #100 # ms # 1
afteri = 10#20 #100 # ms # 20

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

    Ndendsecs=0
    for sec in neuron.h.allsec():
        sec.insert("pas")
        sec.e_pas = e_pas
        sec.Ra = Ra
        sectype = sec.name().split("[")[0]
        if sectype=="dend":
            Ndendsecs+=1
        
        for sec_dict in active_mechs:
            if sec_dict["section"] == sectype:
                if sectype=="soma": # Works
                    exec("sec.cm = {}".format(cm_soma))
                if sectype=="dend": # Works
                    exec("sec.cm = {}".format(cm_dend))
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
    return cell, Ndendsecs

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

    cell, Ndendsecs = return_allen_cell_model(model_folder)
    print('Ndendsecs:',Ndendsecs)

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
    lengths=[]
    for sec in neuron.h.allsec():
        neuron.h.psection()
        sectype = sec.name().split("[")[0]
        if sectype=="dend":
            lengths.append(sec.L)
    
    time = cell.tvec # time
    Nt   = len(time)
    
    #idxs = cell.get_idx(section="dend[6]")
    #Nidx = len(idxs)
    #k = 2
    #kstr = str(k)
    #idxsk = cell.get_idx(section="dend[%s]"%kstr) # Need this awkward notation.
    #idxs2 = cell.get_idx(section="dend[2]")
    
    print('Ndendsecs:',Ndendsecs)
    folder = "figures/%i/current_idur%i_iamp" % (testmodel,idur) + str(iamp)+"/dendritepropagation/"
    metafolder = 'figures/%i/' % testmodel # this is all I need.
    metafilename = metafolder+'%i_dendrites_section_segments.txt' % testmodel
    metafile = open(metafilename,'w')
    for k in range(Ndendsecs):
        kstr = str(k)
        idxs = cell.get_idx(section="dend[%s]"%kstr) # Will this work?
        Nidx = len(idxs)
        metafile.write('%i %i %.16f\n' % (k,Nidx,lengths[k]))
        print('dendrite', k, ', Nidx:',Nidx)
        for j in range(Nidx):
            outfilename = folder+"idur%i_iamp" % idur + str(iamp)+"_cms" + str(cm_soma) + "_cmd" + str(cm_dend) + "_cma"+ str(cm_axon) + "_vinit"+str(v_init)+"_wRa_V_dsec%i_d%i.txt" % (j,k)
            outfile = open(outfilename,'w')
            vmem = cell.vmem[idxs[j],:]   
            
            for i in range(Nt):
                outfile.write('%.16f %.16f\n' % (time[i],vmem[i]))
            outfile.close()  
        plotname = folder+"idur%i_iamp" % idur + str(iamp)+"_cms" + str(cm_soma) + "_cmd" + str(cm_dend) + "_cma"+ str(cm_axon) + "_vin"+str(v_init)+"_wRa_V_d%i.png" % k
        plotname_few = folder+"idur%i_iamp" % idur + str(iamp)+"_cms" + str(cm_soma) + "_cmd" + str(cm_dend) + "_cma"+ str(cm_axon) + "_vin"+str(v_init)+"_wRa_V_d%i_few.png" % k
        plotname_withpos = folder+"idur%i_iamp" % idur + str(iamp)+"_cms" + str(cm_soma) + "_cmd" + str(cm_dend) + "_cma"+ str(cm_axon) + "_vin"+str(v_init)+"_wRa_V_d%i_wp.png" % k
        
        print('idxs:',idxs)
        if len(idxs)==0:
            break
        if plotit==True:
            fig = plt.figure(figsize=[12, 8])
            [plt.plot(cell.tvec, cell.vmem[idx, :], label='%i' % idx) for idx in idxs]
            plt.xlabel('Time (ms)')
            plt.ylabel('Potential (mV)')
            plt.title('Membrane potential along dendrite %i' % k)
            plt.legend(loc='upper right')
            plt.savefig(plotname)  # Too many plots?
            
            fig = plt.figure(figsize=[12, 8])
            plt.plot(cell.tvec, cell.vmem[0, :], label='Soma')
            plt.plot(cell.tvec, cell.vmem[idxs[0], :], label='%i' % idxs[0])
            plt.plot(cell.tvec, cell.vmem[idxs[int(math.floor(Nidx/2))], :], label='%i' % idxs[int(math.floor(Nidx/2))])
            plt.plot(cell.tvec, cell.vmem[idxs[Nidx-1], :], label='%i' % idxs[Nidx-1])
            plt.xlabel('Time (ms)')
            plt.ylabel('Potential (mV)')
            plt.title('Membrane potential along dendrite %i' % k)
            plt.legend(loc='upper right')
            plt.savefig(plotname_few)  
            
            fig = plt.figure(figsize=[19, 9])
            fig.subplots_adjust(wspace=0.5)
            ax_morph = fig.add_subplot(121, aspect=1, xlabel="x", ylabel="y", title="Morphology and section names")
            ax_v = fig.add_subplot(122, ylabel="membrane potential (mV)", xlabel="time (ms)", title="Membrane potential vs time")
            
            [ax_morph.plot([cell.xstart[idx], cell.xend[idx]],
                           [cell.ystart[idx], cell.yend[idx]], c='k', lw=cell.diam[idx])
            for idx in range(cell.totnsegs)]
            
            ax_v.plot(cell.tvec, cell.somav)
            
            path_idx_clrs = lambda num: plt.cm.viridis(num / len(idxs))
            
            for num, idx in enumerate(idxs):
                ax_v.plot(cell.tvec, cell.vmem[idx], c=path_idx_clrs(num))
                ax_morph.plot(cell.xmid[idx], cell.ymid[idx], 'o', c=path_idx_clrs(num), ms=12)
            
            for sec in neuron.h.allsec():
                sec_name = sec.name()
                sec_idxs = cell.get_idx(sec_name)
                plot_idx = sec_idxs[-1]
                ax_morph.text(cell.xend[plot_idx], cell.yend[plot_idx], sec.name(), color='r')
                
            fig.savefig(plotname_withpos)
    metafile.close()
    print('testmodel:',testmodel)
    print('idxs, dend6:', cell.get_idx(section="dend[6]"))
    print('Ndendsecs:',Ndendsecs)

    sys.exit()
