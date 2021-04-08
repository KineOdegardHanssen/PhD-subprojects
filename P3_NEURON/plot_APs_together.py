import matplotlib.pyplot as plt
import numpy as np

##### Adjustable parameters/specifications #########

# Choose capacitance values to vary. cm1='Default' yields value in original model
cm1 = 0.1 #'Default'
cm2 = 0.2 #0.01 #5.0 #15.0#0.01

# Bools for varying Cmsoma or Cmdend:
varysoma = True # True:soma,  False:dendrite
varyall  = False # Hacky solution for varying Cm everywhere. This overrides varysoma

# File/simulation selection:
testmodel = 480633479#478513437#478513407#488462965 #496497595
if testmodel==496497595:
    idur = 1 # ms
else: #elif testmodel==488462965:
    idur = 2 # ms
idur = 1000
iamp = 0.41 #1.0 # nA
v_init = -86.5 # mV

####################################################

# Defaulting to original values:
# DO NOT TOUCH THESE!
if testmodel==496497595:
    cm_soma = 1.14805
    cm_dend = 9.98231
    cm_axon = 3.00603
elif testmodel==488462965:
    cm_soma = 3.31732779736
    cm_dend = 3.31732779736
    cm_axon = 3.31732779736
    v_init = -86.5 # mV
elif testmodel==480633479:
    cm_soma = 0.704866 # 0.704866118957
    cm_dend = 0.704866 # 0.704866118957
    cm_axon = 0.704866 # 0.704866118957
    v_init = -96.8
elif testmodel==478513407: # Should have I=0.17
    cm_soma = 1.0
    cm_dend = 1.0
    cm_axon = 1.0
    v_init = -83.7
elif testmodel==478513437:
    cm_soma = 2.34539964752
    cm_dend = 2.34539964752
    cm_axon = 2.34539964752
    v_init = -86.8

cms = []
if varysoma==True and varyall==False:
    if cm1!='Default':
        cm_soma = cm1
    cms.append(cm_soma)
elif varysoma==False and varyall==False:
    if cm1!='Default':
        cm_dend = cm1
    cms.append(cm_dend)
else:
    cms = []
    if cm1!='Default':
        cm_soma = cm1
        cm_dend = cm1
        cm_axon = cm1
    cms.append(cm_dend)
     
    
times1 = []
times2 = []
V1     = []
V2     = []

folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/'
figname = '/model%i_idur%i_iamp' %(testmodel,idur)+str(iamp)
if varysoma==True and varyall==False:
    folder = folder + 'Varycm_soma/'
elif varysoma==False and varyall==False:
    folder = folder + 'Varycm_dend/'
else:
    folder = folder + 'Varycm_all/'

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


if varysoma==True and varyall==False:
    cm_soma = cm2
    cms.append(cm_soma)
elif varysoma==False and varyall==False:
    cm_dend = cm2
    cms.append(cm_dend)
else:
    cm_soma = cm2
    cm_dend = cm2
    cm_axon = cm2
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
if varysoma==True and varyall==False:
    figname = folder+figname + '_somacm%.5f_%.5f' % (cms[0],cms[1])
    figname_2 = figname +'_cut.png'
    figname   = figname +'.png'
elif varysoma==False and varyall==False:
    figname = folder+figname + '_dendcm%.5f_%.5f' % (cms[0],cms[1])
    figname_2 = figname +'_cut.png'
    figname   = figname +'.png'
else:
    figname = folder+figname + '_allcm%.5f_%.5f' % (cms[0],cms[1])
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
if varyall==True:
    plt.title('Action potentials for different capacitances of the neuron')
if idur>10:
    plt.legend(loc='lower center')
else:
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
if varyall==True:
    plt.title('Action potentials for different capacitances of the neuron')
if idur>10:
    plt.legend(loc='lower center')
else:
    plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(figname_2)
plt.show()