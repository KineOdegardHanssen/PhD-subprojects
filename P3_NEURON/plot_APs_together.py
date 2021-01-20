import matplotlib.pyplot as plt
import numpy as np

##### Adjustable parameters/specifications #########

# Choose capacitance values to vary. cm1='Default' yields value in original model
cm1 = 'Default'
cm2 = 15.0#0.01

# Bools for varying Cmsoma or Cmdend:
varysoma = False # True:soma,  False:dendrite

# File/simulation selection:
testmodel = 496497595
idur = 1 # ms
iamp = 1.0 # nA
v_init = -86.5 # mV

####################################################

# Defaulting to original values:
# DO NOT TOUCH THESE!
cm_soma = 1.14805
cm_dend = 9.98231
cm_axon = 3.00603

cms = []
if varysoma==True:
    if cm1!='Default':
        cm_soma = cm1
    cms.append(cm_soma)
else:
    if cm1!='Default':
        cm_dend = cm1
    cms.append(cm_dend)
     
    
times1 = []
times2 = []
V1     = []
V2     = []

folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/'
figname = '/model%i_idur%i_iamp' %(testmodel,idur)+str(iamp)
if varysoma==True:
    folder = folder + 'Varycm_soma/'
else:
    folder = folder + 'Varycm_dend/'

### First file ###
filename = folder+'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) + '_vinit'+str(v_init)+'_addedRa.txt'

infile = open(filename,'r')
lines   = infile.readlines()
for line in lines:
    words = line.split()
    if len(words)!=0: # words[0]:time, words[1]:V, words[2]:[Ca^2+]
        times1.append(float(words[0]))
        V1.append(float(words[1]))
infile.close()


if varysoma==True:
    cm_soma = cm2
    cms.append(cm_soma)
else:
    cm_dend = cm2
    cms.append(cm_dend)

### Second file ###
filename = folder+'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) + '_vinit'+str(v_init)+'_addedRa.txt'

infile = open(filename,'r')
lines   = infile.readlines()
for line in lines:
    words = line.split()
    if len(words)!=0: # words[0]:time, words[1]:V, words[2]:[Ca^2+]
        times2.append(float(words[0]))
        V2.append(float(words[1]))
infile.close()

# Plot names:
folder = folder+'Compare/'
if varysoma==True:
    figname = folder+figname + 'somacm%.5f_%.5f' % (cms[0],cms[1])
    figname_2 = figname +'_cut.png'
    figname   = figname +'.png'
else:
    figname = folder+figname + 'dendcm%.5f_%.5f.png' % (cms[0],cms[1])
    figname_2 = figname +'_cut.png'
    figname   = figname +'.png'

# Plotting:

plt.figure(figsize=(6,5))
plt.plot(times1, V1, label=r'$C_m$=%.5f $\mu$F/cm$^2$' % cms[0])
plt.plot(times2, V2, label=r'$C_m$=%.5f $\mu$F/cm$^2$' % cms[1])
plt.xlabel('Time [ms]')
plt.ylabel('Voltage [mV]')
if varysoma==True:
    plt.title('Action potentials for different capacitances of the soma')
else:
    plt.title('Action potentials for different capacitances of the basal dendrites')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(figname)

Nind = len(times1) # Hacky solution, maybe improve later
N2   = int(round(Nind/2)) 
plt.figure(figsize=(6,5))
plt.plot(times1[:N2], V1[:N2], label=r'$C_m$=%.5f $\mu$F/cm$^2$' % cms[0])
plt.plot(times2[:N2], V2[:N2], label=r'$C_m$=%.5f $\mu$F/cm$^2$' % cms[1])
plt.xlabel('Time [ms]')
plt.ylabel('Voltage [mV]')
if varysoma==True:
    plt.title('Action potentials for different capacitances of the soma')
else:
    plt.title('Action potentials for different capacitances of the basal dendrites')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(figname_2)