import matplotlib.pyplot as plt
import numpy as np


# Varycm: All, soma or dendrite
varywhichcm = 'a' # All # Only for perisomatic so far
#varywhichcm = 's' # Soma
#varywhichcm = 'd' # Dendrite

##### Adjustable parameters/specifications #########
iduraa = 2
# File/simulation selection:
testmodel = 496497595 #488462965 # 
if testmodel==496497595:
    idur = iduraa # ms
elif testmodel==488462965:
    idur = 2 # ms
iamp = 1.0 # nA
v_init = -86.5 # mV
dendnr = 6

if testmodel==496497595:
    cm_soma = 1.14805
    cm_dend = 9.98231
    cm_axon = 3.00603
    if varywhichcm=='a':
        varycm = [0.01,0.1,0.5,1.0,2.0,3.0]
    elif varywhichcm=='s':
        varycm = [0.01,0.1,0.5,1.0,1.14805,2.0,3.0]
    else: # Dendrite
        varycm  = [0.1,0.5,1.0,2.0,3.0,5.0,9.98231,15.0] # 0.01 yields inf
elif testmodel==488462965:
    cm_soma = 3.31732779736
    cm_dend = 3.31732779736
    cm_axon = 3.31732779736
    if varywhichcm!='a':
        varycm  = [0.01,0.1,0.5,1.0,2.0,3.31732779736,5.0]
    else:
        varycm  = [0.1,0.5,1.0,2.0,3.31732779736,5.0] # Cm=5.0 maybe a bit high for dend
N = len(varycm)

if varywhichcm=='a':
    subfolder = 'Varycm_all/'
    textsnippet = 'varycmall'
    plottitletext = 'whole neuron'
elif varywhichcm=='s':
    subfolder = 'Varycm_soma/'
    textsnippet = 'varycmsoma'
    plottitletext = 'soma'
else:
    subfolder = 'Varycm_dend/'
    textsnippet = 'varycmdend'
    plottitletext = 'dendrites'


folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/dendritepropagation/'
outfolder = folder+subfolder
outfilename = outfolder +'idur%i_iamp' % idur+str(iamp)+ '_'+str(testmodel)+ '_'+ textsnippet+ '_vinit'+str(v_init)+'_addedRa_propvel.txt'
plotname = outfolder +'idur%i_iamp' % idur+str(iamp)+ '_'+str(testmodel)+ '_'+ textsnippet+ '_vinit'+str(v_init)+'_addedRa_propvel.png'
outfile = open(outfilename,'w')
propvels = np.zeros(N)

for i in range(N):
    # Update cm
    if varywhichcm=='a':
        cm_soma = varycm[i]
        cm_dend = varycm[i]
        cm_axon = varycm[i]
    elif varywhichcm=='s':
        cm_soma = varycm[i]
    else:
        cm_dend = varycm[i]
    # Update file name
    if testmodel==496497595:
        infilename = folder+'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) +'_vinit'+str(v_init)+ '_addedRa_V_dend%i_propvel.txt' % dendnr
    else:
        infilename = folder+'idur%i_iamp' % idur+str(iamp)+'_cms' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_V_dend%i_propvel.txt' % dendnr # Needed to cut down on the name length
    infile = open(infilename,'r')
    line = infile.readline()
    propvel = float(line.split()[0])
    propvels[i] = propvel
    outfile.write('%.11f %.16f\n' % (varycm[i],propvel))
    infile.close()
outfile.close()

# Plot results
plt.figure(figsize=(6,5))
plt.plot(varycm,propvels,'-o')
plt.xlabel(r'$C_{m}$ of %s [$\mu$ F/cm$^2$]' % plottitletext)
plt.ylabel(r'Signal velocity along dendrite')
plt.title(r'Signal velocity vs %s $C_{m}$' % plottitletext)
plt.tight_layout()
plt.savefig(plotname)