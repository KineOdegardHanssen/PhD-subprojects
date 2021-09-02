import numpy as np
import matplotlib.pyplot as plt

import neuron
from neuron import h
import LFPy

dtexp = -8

iamp = -0.1
idur = 100
idelay = 10
tstop = idur+idelay+10.

v_init = -65
cm = 3.0
Ra = 100.
somasize = 10
gpas = 0.0003
epas = -65

h("""
create soma[1]
objref all

soma.L = 10
soma.diam = 10 

all = new SectionList()
soma[0] all.append()

forall {nseg = 1}

Ra = 100.
cm = 0.5
Rm = 1./0.0003

forall {
    insert pas
    }
""")

dt = 2**dtexp #6 #3
print("dt: ", dt)
cell_params = {          # various cell parameters,
            'morphology': h.all,
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

# Set length
for sec in h.allsec():
    sec.insert("pas")
    sectype = sec.name().split("[")[0]
    sec.cm = cm
    sec.L = somasize
    sec.diam = somasize
    for seg in sec:
        seg.g_pas = gpas #0.0003: Standard value (HH)
        seg.e_pas = epas

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

folder = 'Results/IStim/Soma%i/current_idur%i_iamp'%(somasize,idur)+str(iamp)+'/dtexp%i/' % dtexp
filename = folder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_gpas'+str(gpas)+'_epas'+str(epas)+'_V.txt' 
plotname = folder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_gpas'+str(gpas)+'_epas'+str(epas)+'_V.png' 

plt.plot(cell.tvec, cell.vmem[0,:])

plt.savefig(plotname, dpi=300)

V = cell.vmem[0,:]

file = open(filename,'w')
for i in range(len(V)):
    file.write('%.16e %.16e\n' % (cell.tvec[i],V[i]))
file.close()

print('neuron.h.soma[0].L:',neuron.h.soma[0].L)
