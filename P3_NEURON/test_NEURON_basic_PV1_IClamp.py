import matplotlib.pyplot as plt
from neuron import h, gui

h.load_file('import3d.hoc')

### No AllenSDK

### Found this on https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3257
### Courtesy of R. A. McDougal
def instantiate_swc(filename):
    """load an swc file and instantiate it"""
    
    ### I have done this already: ###########################################
    ## load the NEURON library (just in case h is defined otherwise elsewhere)
    #from neuron import h, gui
    
    ## a helper library, included with NEURON
    #h.load_file('import3d.hoc')
    #########################################################################
    
    # load the data. Use Import3d_SWC_read for swc, Import3d_Neurolucida3 for
    # Neurolucida V3, Import3d_MorphML for MorphML (level 1 of NeuroML), or
    # Import3d_Eutectic_read for Eutectic. (There is also an 
    # Import3d_Neurolucida_read for old Neurolucida files, but I've never seen one
    # in practice; try Import3d_Neurolucida3 first.)
    cell = h.Import3d_SWC_read()
    cell.input(filename)

    # easiest to instantiate by passing the loaded morphology to the Import3d_GUI
    # tool; with a second argument of 0, it won't display the GUI, but it will allow
    # use of the GUI's features
    i3d = h.Import3d_GUI(cell, 0)
    i3d.instantiate(None)


# Hmmm, syntax.
basepath  = 'C:/Users/Kine/Documents/Projects_PhD/P3_PNN_Capacitance/Cellemodeller/' #'C:\Users\Kine\Documents\Projects_PhD\P3_PNN_Capacitance\Cellemodeller\'
##modelpath = basepath + 'Arkhipov_etal_2018/SI_1/example_1/cell_models/'
##PV1_path  = modelpath + '330080937/'
##PV2_path  = modelpath + '318331342/'
PV1_morph = 'Pvalb-IRES-Cre_Ai14_IVSCC_-176847.04.02.01_470522102_m.swc'

print('PV1_morph:',PV1_morph)

## Tried standard tools, they did not work: ##
#h.xopen(PV1_morph)
#neuron.h.xopen(PV1_morph)
##############################################

print('Importing morphology file:')

## Found this online:
instantiate_swc(PV1_morph)

h.psection()

print('Have run h.psection()')
##################################################################
# Set biophysical parameters
##################################################################
###Changing capacitance (Need to figure out how to do this just for the soma)###

for sec in neuron.h.allsec():
    sec.Ra = 100       # Axial resistivity in Ohm*m
    sec.cm = 1         # membrane capacitance in uF/cm2

### Do I need to do something like this?: ###################
# insert 'passive' membrane mechanism, adjust parameters.
# None: without a leak mechanism, the neuron will be a
# perfect integrator
#soma.insert('pas')          
#soma(0.5).pas.g = 0.0002    # membrane conducance in S/cm2
#soma(0.5).pas.e = -65.       # leak reversal potential in mV
#############################################################

################################################################################
# Model instrumentation
################################################################################
## Could use IClamp instead...
# Attach current clamp to the neuron
iclamp = neuron.h.IClamp(0.5, sec=h.soma[0])
iclamp.delay = 100. # current delay period in ms
iclamp.dur = 1.   # duration of stimulus current in ms
iclamp.amp = 0.2    # amplitude of current in nA

##################################################################
# Set up recording of variables
##################################################################
# NEURON variables can be recorded using Vector objects. Here, we
# set up recordings of time, voltage and stimulus current with the
# record attributes.
t = neuron.h.Vector()   
v = neuron.h.Vector()
i = neuron.h.Vector()
# recordable variables must be preceded by '_ref_':
t.record(neuron.h._ref_t)   
v.record(h.soma[0](0.5)._ref_v)
i.record(iclamp._ref_i)


# print out section information again
for sec in neuron.h.allsec():
    neuron.h.psection()

##################################################################
# Simulation control
##################################################################
neuron.h.dt = 0.1          # simulation time resolution
tstop = 500.        # simulation duration
v_init = -65        # membrane voltage(s) at t = 0

print('Have set parameters')

### From Gaute's course:
def initialize():
    '''
    initializing function, setting the membrane voltages to v_init and
    resetting all state variables
    '''
    h.finitialize(v_init)
    h.fcurrent()

def integrate():
    '''
    run the simulation up until the simulation duration
    '''
    while h.t < tstop:
        h.fadvance()

# run simulation
initialize()
integrate()

################################################################################
# Plot simulated output
################################################################################
fig, axes = plt.subplots(2)
fig.suptitle('stimulus current and point-neuron response')
axes[0].plot(t, i, 'r', lw=2)
axes[0].set_ylabel('current (nA)')

axes[1].plot(t, v, 'r', lw=2)
axes[1].set_ylabel('voltage (mV)')
axes[1].set_xlabel('time (ms)')

# tight layout
for ax in axes: ax.axis(ax.axis('tight'))

fig.savefig('Test_results/test_IClamp_currdur%i_ampl%.2f_allC%.1f_modchanged.pdf' % (iclamp.dur,iclamp.amp,sec.cm))
plt.show()


#h.run()
