import uncertainpy as un
import os
from os.path import join
import sys
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import json
import neuron
import chaospy as cp
import numpy as np
import efel
    
import LFPy
    

if __name__ == '__main__': # EVERYTHING inside this?

    def runit(cm_soma,cm_dend,cm_axon,testmodel):
        testmodelname = 'neur_%i' % testmodel
        all_models    = [testmodelname]
        def return_allen_cell_model(model_folder,testmodel):
            
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
            
            params = json.load(open(join(model_folder, "fit_parameters.json"), 'r'))
            
            celsius = params["conditions"][0]["celsius"]
            reversal_potentials = params["conditions"][0]["erev"]
            #v_init = params["conditions"][0]["v_init"] # This might be wrong, will set it by hand above
            active_mechs = params["genome"]
            Ra = params["passive"][0]["ra"]  # This was not here before
            neuron.h.celsius = celsius
        
            # Define cell parameters
            cell_parameters = {
                'morphology' : join(model_folder, 'reconstruction.swc'),
                'v_init' : v_init,    # initial membrane potential
                'passive' : False,   # turn on NEURONs passive mechanism for all sections
                'nsegs_method' : 'lambda_f', # spatial discretization method
                'lambda_f' : 200.,           # frequency where length constants are computed
                'dt' : 2.**-5,      # simulation time step size
                'tstart' : -100.,      # start time of simulation, recorders start at t=0
                'tstop' : 1200.,     # stop simulation at this time.
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
        mod_folder = "Allen_test_changecapacitance/cell_models/"+testmodelname+"/modfiles"
        if "win64" in sys.platform:
            print('Detected sys.platform as win64')
            warn("no autompile of NMODL (.mod) files on Windows. " 
                 + "Run mknrndll from NEURON bash in the folder cells and rerun example script")
            if not pth in neuron.nrn_dll_loaded:
                neuron.h.nrn_load_dll(mod_folder+"/nrnmech.dll")
            neuron.nrn_dll_loaded.append(mod_folder)
    
        model_idx = 0
        model_name = all_models[model_idx]
        model_folder = join("cell_models", model_name)
    
        cell = return_allen_cell_model(model_folder,testmodel)
        
        pointprocess = {
                'idx': 0,
                'record_current': True,
                'pptype': 'IClamp',
                'amp': -0.5,
                'dur': 1000.,
                'delay': 100.,
            }
        stimulus = LFPy.StimIntElectrode(cell,**pointprocess)
        cell.simulate(rec_vmem=True, rec_variables=['cai'])
        
        
        time = cell.tvec # time
        Nt   = len(time)
        voltage = cell.vmem[0,:]
        
        #sys.exit()
        
        # Trying to find Cm # Ugh, this can lead to a lot of crashes...
        
        # Now we will construct the datastructure that will be passed to eFEL
    
        # A 'trace' is a dictionary
        trace1 = {}
        
        # Set the 'T' (=time) key of the trace
        trace1['T'] = time
        
        # Set the 'V' (=voltage) key of the trace
        trace1['V'] = voltage
        
        Nt = len(time)
        current2 = []
        current = np.zeros(Nt) # Or set this through LFPy?
        for i in range(Nt):
            if time[i]>100. and time[i]<1100.:
                 current[i] = -0.5
                 current2.append(-0.5)
        current2 = np.array(current2)
        
        # Set the 'I' (=current) key of the trace # Hope this helps # Or do I need the stimulus current?
        trace1['I'] = current2 #!!!!
        print('current:',current)
        print('np.sum(current):',np.sum(current))
        
        # Fix this:
        # Stim Current sample index for detecting current of stimulation step
        trace1['stimulus_current'] = current2 #-0.5 #current
    
        # Set the 'stim_start' (time at which a stimulus starts, in ms)
        # key of the trace
        # Warning: this need to be a list (with one element)
        trace1['stim_start'] = [100.]
        
        # Set the 'stim_end' (time at which a stimulus end) key of the trace
        # Warning: this need to be a list (with one element)
        trace1['stim_end'] = [1100.]
        
        # Multiple traces can be passed to the eFEL at the same time, so the
        # argument should be a list
        traces = [trace1]
        
        v_deflect = (efel.getFeatureValues(traces,['voltage_deflection']))
        trace1['voltage_deflection']= v_deflect[0]['voltage_deflection']
        
        input_res = (efel.getFeatureValues(traces,['ohmic_input_resistance']))
        input_res = (input_res[0]['ohmic_input_resistance'])#*1000
        
        # Time constant, see Efel docs.
        time_constant = efel.getFeatureValues(traces,['time_constant'])
        time_constant = time_constant[0]['time_constant']
        
        # Capacitance in Farads.
        capacitance_out = (time_constant*10**-3)/(input_res*10**6)
        capacitance_out = capacitance_out*10**12 # in pF
        
        # Do I need this?
        # Stim Current sample index for detecting current of stimulation step
        #stim = trace1['I'][1100.]
        return_time = 1100. # 'cause it needs time, for some reason
        
        return return_time, capacitance_out
    
    
    testmodel = 496497595
    # Defaulting to original values:
    # DO NOT TOUCH THESE!
    # SET THEM BELOW INSTEAD!    # Do I need these now? Will I run Uncertainpy for all?
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
        cm_soma = 0.704866118957
        cm_dend = 0.704866118957
        cm_axon = 0.704866118957
    elif testmodel==478513437:
        cm_soma = 2.34539964752
        cm_dend = 2.34539964752
        cm_axon = 2.34539964752
    
    model = un.Model(run=runit)
    cm_soma = cp.Uniform(0.5,5.0)
    cm_dend = cp.Uniform(0.5,5.0)
    cm_axon = cp.Uniform(0.5,5.0)
    
    parameters = {"cm_soma":cm_soma,"cm_dend":cm_dend,"cm_axon":cm_axon,"testmodel": testmodel}    
    
    UQ = un.UncertaintyQuantification(model=model, parameters=parameters)
    data = UQ.quantify(seed=10) # To easier be able to reproduce the result # Does this fix ALL the randomness?
    
    print('data:', data) # No clue if this will work
    
    # How do I plot/visualize? # Can test running first