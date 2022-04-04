import numpy as np
import matplotlib.pyplot as plt

import neuron
from neuron import h
import LFPy
################################ HH #############################
# Default HH values:
ena = 50
ek = -77
el_hh = -54.3
gnabar_hh = 0.12
gkbar_hh = 0.036
gl_hh = 0.0003

gnafactor = 5.0

### Change HH values here: ####
#ena = 49.3
#ek = -85
#el_hh = -70
gnabar_hh *= gnafactor
#gkbar_hh = 0.020
#gl_hh = 0

######################### Other params ##########################
# iamps: nA
iamps = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5]
idur = 1000
idelay = 10
tstop = idur+idelay+10. #2.
v_init = -70 #-86.5         # Was -70

cm = 1.5
Ra = 100. #150.
somasize = 10  ### Change this in h("""...""") too!!!

for iamp in iamps:
    # Need to manually change somasize in the following section?
    h("""
    create soma
    objref all
    
    soma.L = 10
    soma.diam = 10  
    
    all = new SectionList()
    soma all.append()
    
    forall {nseg = 1}
    
    Ra = 100.
    cm = 1.5
    Rm = 30000
    
    forall {
        insert hh
        }
    """) # Had insert pas, that didn't give any spiking # Somehow the cm that is set here doesn't work
    
    # Vary HH properties:
    for sec in h.soma:
        sec.ena       = ena
        sec.ek        = ek
        sec.el_hh     = el_hh
        sec.gnabar_hh = gnabar_hh
        sec.gkbar_hh  = gkbar_hh
        sec.gl_hh     = gl_hh
    
    dt = 2**-6 #-3
    print("dt: ", dt)
    cell_params = {          # various cell parameters,
                'morphology': h.all, # look this up
                'delete_sections': False,
                'v_init' : v_init,    # initial crossmembrane potential
                'cm': cm,
                'Ra': Ra,
                'passive' : False,   # switch on passive mechs
                'nsegs_method' : None,
                'dt' : dt,   # [ms] dt's should be in powers of 2 for both,
                'tstart' : -300.,    # start time of simulation, recorders start at t=0
                'tstop' : tstop,   # stop simulation at 200 ms. These can be overridden
            }
    
    
    cell = LFPy.Cell(**cell_params)
    
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
    #syn.set_spike_times(np.array([1]))
    cell.simulate(rec_vmem=True, rec_imem=True)
    
    for sec in neuron.h.allsec():
        neuron.h.psection()
    
    #ena, ek, el_hh, gnabar_hh, gkbar_hh, gl_hh
    folder = 'Results/IStim/Soma%i/current_idur%i_iamp'%(somasize,idur)+str(iamp)+'/'
    hhstring = '_gnafactor'+str(gnafactor)
    filename = folder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+hhstring+'_Ra%i_vinit' %Ra+str(v_init)+'_V.txt' 
    plotname = folder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+hhstring+'_Ra%i_vinit' %Ra+str(v_init)+'_V.png' 
    
    plt.figure(figsize=(6,5))
    plt.plot(cell.tvec, cell.vmem[0,:])
    
    plt.savefig(plotname, dpi=300)

    V = cell.vmem[0,:]
    
    file = open(filename,'w')
    for i in range(len(V)):
        file.write('%.16e %.16e\n' % (cell.tvec[i],V[i]))
    file.close()
    
    cell.__del__()      