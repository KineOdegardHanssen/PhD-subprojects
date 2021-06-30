import numpy as np
import matplotlib.pyplot as plt

import neuron
from neuron import h
import LFPy


iamp = 0.0 #10.0
idur = 1000
idelay = 10
tstop = idur+idelay+10. #2.
v_init = -70 #-86.5         # Was -70

cm = 1.0
Ra = 100. #150.
somasize = 15

# Need to manually change somasize in the following section?
h("""
create soma
objref all

soma.L = 15
soma.diam = 15  

all = new SectionList()
soma all.append()

forall {nseg = 1}

Ra = 100.
cm = 1.0
Rm = 30000

forall {
    insert hh
    insert pas
    }
""")

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
            'tstart' : -100.,    # start time of simulation, recorders start at t=0
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

folder = 'Results/IStim/Soma%i/Withpas/current_idur%i_iamp'%(somasize,idur)+str(iamp)+'/'
filename = folder+'somaonly_withpas_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_V.txt' 
plotname = folder+'somaonly_withpas_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_V.png' 

plt.plot(cell.tvec, cell.vmem[0,:])

plt.savefig(plotname, dpi=300)

V = cell.vmem[0,:]

file = open(filename,'w')
for i in range(len(V)):
    file.write('%.16e %.16e\n' % (cell.tvec[i],V[i]))
file.close()