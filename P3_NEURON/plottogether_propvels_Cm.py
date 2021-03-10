import matplotlib.pyplot as plt
import numpy as np
import math
import sys

# Varycm: All, soma or dendrite
#varywhichcm = 'a' # All # Only for perisomatic so far
#varywhichcm = 's' # Soma
varywhichcm = 'd' # Dendrite

##### Adjustable parameters/specifications #########
# File/simulation selection:
testmodel = 496497595 # 488462965 # 
iduraa = 2
if testmodel==496497595:
    idur = iduraa # ms
elif testmodel==488462965:
    idur = 2 # ms
iamp = 1.0 # nA
v_init = -86.5 # mV


####################################################

# Defaulting to original values:
# DO NOT TOUCH THESE!
if testmodel==496497595:
    cm_soma = 1.14805
    cm_dend = 9.98231
    cm_axon = 3.00603
    if varywhichcm=='a':
        cms = [0.01,0.1,0.5,1.0,2.0,3.0,5.0,10.0]
    elif varywhichcm=='s':
        cms = [0.01,0.1,0.5,1.0,1.14805,2.0,3.0,5.0,10.0]
    elif varywhichcm=='d':
        cms = [0.01,0.1,0.5,1.0,2.0,3.0,5.0,9.98231,10.0]
elif testmodel==488462965:
    cm_soma = 3.31732779736
    cm_dend = 3.31732779736
    cm_axon = 3.31732779736
    cms = [0.01,0.1,0.5,1.0,2.0,3.0,5.0,10.0] # Might change. We'll see

Ncm = len(cms)
propvels = np.zeros(Ncm)
propvels_rms = np.zeros(Ncm)
folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/dendritepropagation/' # Yeah, this works for all.


if varywhichcm=='a':
    textsnippet = 'whole neuron'
    outfilename = folder+'idur%i_iamp' % idur+str(iamp)+'_varycmall_vi'+str(v_init)+ '_wRa_propvels.txt'
    plotname = folder+'idur%i_iamp' % idur+str(iamp)+'_varycmall_vi'+str(v_init)+ '_wRa_propvels.png'
elif varywhichcm=='s':
    textsnippet = 'soma'
    outfilename = folder+'idur%i_iamp' % idur+str(iamp)+'_varycmsoma_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_propvels.txt'
    plotname = folder+'idur%i_iamp' % idur+str(iamp)+'_varycmsoma_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_propvels.png'
elif varywhichcm=='d':
    textsnippet = 'dendrite'
    outfilename = folder+'idur%i_iamp' % idur+str(iamp)+'_varycmdend' + str(cm_soma) +'_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_propvels.txt'
    plotname = folder+'idur%i_iamp' % idur+str(iamp)+'_varycmdend' + str(cm_soma) +'_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_propvels.png'

outfile = open(outfilename,'w')

for i in range(Ncm):
    if varywhichcm=='a':
        cm_soma = cms[i]
        cm_dend = cms[i]
        cm_axon = cms[i]
    elif varywhichcm=='s':
        cm_soma = cms[i]
    elif varywhichcm=='d':
        cm_dend = cms[i]
    infilename = folder+'idur%i_iamp' % idur+str(iamp)+'_cms' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_avgpropvel.txt'
    infile = open(infilename,'r')
    line = infile.readline() # We only need the first line
    propvels[i] = float(line.split()[0])
    propvels_rms[i] = float(line.split()[1])
    outfile.write('%.2f %.16f %.16f\n' % (cms[i],propvels[i],propvels_rms[i]))
    infile.close()
outfile.close()

plt.figure(figsize=(6,5))
plt.errorbar(cms, propvels, yerr=propvels_rms,capsize=2)
plt.xlabel(r'$C_m$ of %s [$\mu$F/cm$^2$]' % textsnippet)
plt.ylabel(r'Propagation velocity of signal [m/s]')
plt.title(r'Propagation velocity of signal vs $C_m$ of %s' % textsnippet)
plt.tight_layout()
plt.savefig(plotname)