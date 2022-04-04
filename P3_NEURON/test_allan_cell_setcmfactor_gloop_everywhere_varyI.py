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
testmodel = 478513407 # Perisomatic version of 497233271, but only yields one peak
#testmodel = 480633479 # Perisomatic version of 497230463, only one peak ####
#testmodel = 478513437 # Perisomatic version of 497233075
#testmodel = 489931686 # NEW perisomatic model
#testmodel = 485694403
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
elif testmodel==485694403:
    v_init = -88.8

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
elif testmodel==485694403:
    cm_soma = 0.763348896
    cm_dend = 0.763348896
    cm_axon = 0.763348896
    
# Changing values of membrane capacitance:
cm_changecmf = 1.0

varymech = 'None' #'Na' # 'K' # 'pas'
varyE_bool = False # True # Easiest to do this. I think.
namestringfirst = ''
if varymech=='Na':
    varyE = 40 #[40,50,53,60,70] # 53 not run
    namestringfirst = namestringfirst + 'ENa'+str(varyE)
elif varymech=='K':
    varyE = -100 
    namestring = namestringfirst + 'EK'+str(varyE)
elif varymech=='pas': # THIS CAN VARY A LOT, RIGHT? 
    varyE = -83 #[-70,-73,-77,-80,-90]
    namestringfirst = namestringfirst + 'EK'+str(varyE)
elif varymech=='None':
    varyE = 'None'
varygbool = True # False # 
varyIh       = False # True # 
vary_NaV     = False # True # 
vary_Kd      = False # True # 
vary_Kv2like = True # False # True # 
vary_Kv3_1   = False # True # 
vary_K_T     = False # True # 
vary_Im_v2   = False # True # 
vary_SK      = False # True #   
vary_Ca_HVA  = False # True # 
vary_Ca_LVA  = False # True # 
vary_gpas    = False # True # 

varyglist = [6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9]#[21.0,22.0,23.0,24.0]#[22.0,25.0,27.0,30.0,35.0,40.0,45.0]#[11.1,11.2,11.3,11.4,11.5,11.6,11.7,11.8,11.9]#[4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9]#[0.00001,0.0001]#[25.0,50.0,75.0,100.0]#[15.0,20.0,50.0,75.0,100.0]#[25.0,30.0,40.0,50.0,100.0]#[0.001,0.01,0.05]#[11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,50.0]#[20,50,100]#[4.0,5.0,6.0,7.0,8.0,9.0,10.0]#[1.0] #[0.1,0.3,0.5,0.7,0.9,1.1,1.5,2.0,3.0] # 'None' # 
 # Placeholder, figure this out

# "gbar_Ih"	"gbar_NaV"	"gbar_Kd"	"gbar_Kv2like"		"gbar_Kv3_1"	"gbar_K_T"	"gbar_Im_v2"	"gbar_SK"	"gbar_Ca_HVA"		"gbar_Ca_LVA"	
# "g_pas" (set individually in soma, dend, axon)
# Change current:
idur = 2000 # 10 # 200 #2 100 #  # ms # 1 #
iamp = 0.4

idelay = 100  #     ms #
afteri = 100  # 1 #     ms # 

tstop_i = idur+afteri+idelay

if varyIh==True:
    namestringfirst = namestringfirst + '_gIh'
if vary_NaV==True:
    namestringfirst = namestringfirst + '_gNaV'
if vary_Kd==True:
    namestringfirst = namestringfirst + '_gKd'
if vary_Kv2like==True:
    namestringfirst = namestringfirst + '_gKv2like'
if vary_Kv3_1==True:
    namestringfirst = namestringfirst + '_gKv31'
if vary_K_T==True:
    namestringfirst = namestringfirst + '_gKT'
if vary_Im_v2==True:
    namestringfirst = namestringfirst + '_gImv2'
if vary_SK==True:
    namestringfirst = namestringfirst + '_gSK'
if vary_Ca_HVA==True:
    namestringfirst = namestringfirst + '_gCaHVA'
if vary_Ca_LVA==True:
    namestringfirst = namestringfirst + '_gCaLVA'
if vary_gpas==True: 
    namestringfirst = namestringfirst + '_gpas'

def return_allen_cell_model(model_folder,gfac):

    params = json.load(open(join(model_folder, "fit_parameters.json"), 'r'))

    Ra = params["passive"][0]["ra"]

    e_pas = params["passive"][0]["e_pas"]
    celsius = params["conditions"][0]["celsius"]
    cm_base = params["passive"][0]["cm"][0]["cm"]
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
        'v_init': v_init,    # initial membrane potential
        'passive': False,   # turn on NEURONs passive mechanism for all sections
        'nsegs_method': 'fixed_length', # spatial discretization method
        'max_nsegs_length': 20.,
        #'lambda_f' : 200.,           # frequency where length constants are computed
        'dt': 2.**-7,      # simulation time step size
        'tstart': -600.,      # start time of simulation, recorders start at t=0
        'tstop': tstop_i,     # stop simulation at 100 ms.
        'custom_code': ['remove_axon.hoc']
    }

    cell = LFPy.Cell(**cell_parameters)
    cell.set_rotation(z=np.pi/1.25)

    for sec in neuron.h.allsec():
        sec.insert("pas")
        sec.e_pas = e_pas
        sec.Ra = Ra
        sectype = sec.name().split("[")[0]
        
        cm_new = cm_base*cm_changecmf
        for sec_dict in active_mechs:
            if sec_dict["section"] == sectype:
                exec("sec.cm = {}".format(cm_new))
                if not sec_dict["mechanism"] == "":
                    sec.insert(sec_dict["mechanism"])
                exec("sec.{} = {}".format(sec_dict["name"], sec_dict["value"]))

        for sec_dict in reversal_potentials:
            if sec_dict["section"] == sectype:
                print(sectype, sec_dict)
                for key in sec_dict.keys():
                    if not key == "section":
                        exec("sec.{} = {}".format(key, sec_dict[key]))
        #for seg in sec: # Successfully changing cm everywhere
        #    print('sectype',sectype,'; seg.cm:', seg.cm)
        
     
    for sec in neuron.h.allsec():
        sectype = sec.name().split("[")[0]
        warning = 'Appropriate amount of gs set (zero or one)'
        gnumber = 0
        # First: Change mechanisms 
        # Mechanisms that can be anywhere
        if varymech=='pas':
            if varyE!='None':
                sec.e_pas = varyE
        if vary_gpas==True: ### NB! THIS CHANGES DEPENDING ON WHERE ON THE NEURON WE ARE!
            sec.g_pas *= gfac
            gnumber+=1
        if sectype=='soma':
            if varymech=='Na':
                if varyE!='None':
                    sec.ena = varyE
            elif varymech=='K':
                if varyE!='None':
                    sec.ek = varyE
            if varygbool==True: ## REMEMBER TO IMPLEMENT!: ISSUE A WARNING IF MORE g's HAVE BEEN CHANGED
                if varyIh==True:
                    sec.gbar_Ih *= gfac 
                    gnumber+=1
                if vary_NaV==True:
                    sec.gbar_NaV *= gfac
                    gnumber+=1
                if vary_Kd==True:
                    sec.gbar_Kd *= gfac 
                    gnumber+=1
                if vary_Kv2like==True:
                    sec.gbar_Kv2like *= gfac
                    gnumber+=1
                if vary_Kv3_1==True:
                    sec.gbar_Kv3_1 *= gfac
                    gnumber+=1
                if vary_K_T==True:
                    sec.gbar_K_T *= gfac 
                    gnumber+=1
                if vary_Im_v2==True:
                    sec.gbar_Im_v2 *= gfac 
                    gnumber+=1
                if vary_SK==True:
                    sec.gbar_SK *= gfac
                    gnumber+=1
                if vary_Ca_HVA==True:
                    sec.gbar_Ca_HVA *= gfac 
                    gnumber+=1
                if vary_Ca_LVA==True:
                    sec.gbar_Ca_LVA *= gfac 
                    gnumber+=1
                if gnumber>1:
                    warning = 'WARNING! %i gs set, max. 1 appropriate!' % gnumber
    
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

for g in varyglist:
    namestring = namestringfirst + str(g)+'p'
    for model_idx in range(len(all_models)):
        model_name = all_models[model_idx]
        #print(model_idx, model_name)
        model_folder = join("cell_models", model_name)
    
        cell = return_allen_cell_model(model_folder,g)
    
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

        fig.savefig(join("figures", "%i" % testmodel,"current_idur%i_iamp" % idur + str(iamp),  '{}_changecmf{}_everywhere_vinit{}_addedRa.png'.format(namestring,cm_changecmf,v_init)))
    
        # print out section information: # Works even though I do everything through LFPy
        for sec in neuron.h.allsec():
            neuron.h.psection()

        fig = plt.figure(figsize=[12, 8])
        plt.plot(cell.tvec, cell.vmem[0, :])
        plt.xlabel('Time (ms)')
        plt.ylabel('Potential (mV)')
        plt.title('Membrane potential in soma')
        plt.legend(loc='upper right')
        fig.savefig(join("figures", "%i" % testmodel,"current_idur%i_iamp" % idur + str(iamp),  '{}_changecmf{}_everywhere_vinit{}_addedRa_big.png'.format(namestring,cm_changecmf,v_init)))
    
        '''
        fig = plt.figure(figsize=[12, 8])
        [plt.plot(cell.tvec, cell.vmem[idx, :], label='%i' % idx) for idx in idxs]
        plt.xlabel('Time (ms)')
        plt.ylabel('Potential (mV)')
        plt.title('Membrane potential along axon')
        plt.legend(loc='upper right')
        fig.savefig(join("figures", "%i" % testmodel,"current_idur%i_iamp" % idur + str(iamp), "axon", '{}_changecmf{}_everywhere_vinit{}_ax_wRa.png'.format(namestring,cm_changecmf,v_init)))
    
        fig = plt.figure(figsize=[12, 8])
        #plt.plot(cell.tvec, cell.rec_variables['cai'][0, :], label='Soma') # Soma is high
        [plt.plot(cell.tvec, cell.rec_variables['cai'][idx, :], label='Axon seg. %i' % idx) for idx in idxs] # Worth a shot # Is this even the axon? Should it not be the soma?
        plt.xlabel('Time (ms)')
        plt.ylabel(r'Ca$^{2+}$-concentration (mM)')
        plt.title('Ca$^{2+}$-concentration along axon')
        plt.legend(loc='lower right')
        fig.savefig(join("figures", "%i" % testmodel, "current_idur%i_iamp" % idur + str(iamp), "axon", "Ca", '{}_chcmf{}_everywhere_vinit{}_axCaRa.png'.format(namestring,cm_changecmf,v_init)))
        
        fig = plt.figure(figsize=[12, 8])
        plt.plot(cell.tvec, cell.rec_variables['cai'][0, :])
        plt.xlabel('Time (ms)')
        plt.ylabel(r'Ca$^{2+}$-concentration (mM)')
        plt.title('Ca$^{2+}$-concentration in soma')
        fig.savefig(join("figures", "%i" % testmodel,"current_idur%i_iamp" % idur + str(iamp), 'Ca', '{}_changecmf{}_everywhere_vinit{}_Ca_wRa.png'.format(namestring,cm_changecmf,v_init)))
        '''
        
        outfilename = "figures/%i/current_idur%i_iamp" % (testmodel,idur) + str(iamp)+"/"+namestring+"_changecmf" + str(cm_changecmf) + "_everywhere_vinit"+str(v_init)+"_addedRa.txt"
        outfile = open(outfilename,'w')
        vmem_soma = cell.vmem[0,:]#[idx, :] # Have time array too...
        Ca_soma   = cell.rec_variables['cai'][0, :] # Don't know if I need this...
        #print('vmem_soma[0]:',vmem_soma[0])
        
        time = cell.tvec # time
        Nt   = len(time)
        
        for i in range(Nt):
            outfile.write('%.16f %.16f\n' % (time[i],vmem_soma[i]))
        outfile.close()    
        
        vmax = max(vmem_soma) 
        vmin = min(vmem_soma) 
        deltav = vmax-vmin
        vthr  = -40 #vmax-0.15*deltav # If there is a peak above this value, we count it
        vprev = vthr-40 # A peak never kicks in at initiation, even if I change vthr
        Npeaks = 0
        for i in range (1,len(vmem_soma)-1):  
            #print(vmem_soma[i])
            if vmem_soma[i-1]<vmem_soma[i] and vmem_soma[i+1]<vmem_soma[i] and vmem_soma[i]>vthr:
                Npeaks+=1
        print(Npeaks, ' peaks for model ', testmodel, ',' , namestring)
        print('iamp:',iamp)
        print(namestring)
        print('vmax:', vmax)

        outfilename = "figures/%i/current_idur%i_iamp" % (testmodel,idur) + str(iamp)+"/"+namestring+"_changecmf" + str(cm_changecmf) + "_everywhere_vinit"+str(v_init)+"_addedRa_vmax.txt"
        outfile = open(outfilename,'w') 
        outfile.write('%.5f' % vmax)
        outfile.close()   
        
        outfilename = "figures/%i/current_idur%i_iamp" % (testmodel,idur) + str(iamp)+"/"+namestring+"_changecmf" + str(cm_changecmf) + "_everywhere_vinit"+str(v_init)+"_addedRa_Npeaks.txt"
        outfile = open(outfilename,'w') 
        outfile.write('%i' % Npeaks)
        outfile.close()   
    
        #print('vmem_soma[0]:',vmem_soma[0])
    
        print(namestring)
        t1 = tm.clock()
        print('Run time:', t1-t0)
        print('changecmf:',cm_changecmf)
    cell.__del__()                            
sys.exit()
