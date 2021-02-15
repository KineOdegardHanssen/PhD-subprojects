import matplotlib.pyplot as plt
import numpy as np

# Varycm: All, soma or dendrite
varywhichcm = 'a' # All
#varywhichcm = 's' # Soma
#varywhichcm = 'd' # Dendrite


# File/simulation selection:
testmodel = 488462965 #496497595
idur = 1000 # ms
iamp = -0.5 # nA
idelay = 100
idur   = 1000
v_init = -86.5 # mV

if varywhichcm=='a':
    subfolder = 'Varycm_all/'
    textsnippet = 'varycmall'
    plottitletext = 'all compartments'
elif varywhichcm=='s':
    subfolder = 'Varycm_soma/'
    textsnippet = 'varycmsoma'
    plottitletext = 'soma'
else:
    subfolder = 'Varycm_dend/'
    textsnippet = 'varycmdend'
    plottitletext = 'dendrites'

varycm = []
if testmodel==496497595:
    A = 479.039 
    cm_soma = 1.14805
    cm_dend = 9.98231
    cm_axon = 3.00603
    if varywhichcm=='s':
        varycm = [0.01,1.14805,3]
    elif varywhichcm=='d':
        varycm = [1.0,9.98231,15.0]
elif testmodel==488462965:
    A = 479.039 
    cm_soma = 3.31732779736 
    cm_dend = 3.31732779736 
    cm_axon = 3.31732779736
    if varywhichcm!='a':
        varycm  = [0.01,0.1,0.5,1.0,2.0,3.31732779736,5.0,10.0] # This crashes for all.
    else:
        varycm  = [0.5,1.0,2.0,3.31732779736,5.0,10.0]
N = len(varycm)
cms = np.zeros(N)

folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/'
outfolder = folder +subfolder
outfilename = outfolder +'idur%i_iamp' % idur+str(iamp)+ '_'+str(testmodel)+ '_'+ textsnippet+ '_vinit'+str(v_init)+'_addedRa_Cm.txt'
plotname    = outfolder +'idur%i_iamp' % idur+str(iamp)+ '_'+str(testmodel)+ '_'+ textsnippet+ '_vinit'+str(v_init)+'_addedRa_Cm.png'
outfile = open(outfilename,'w')
for i in range(N):
    # Pick somas
    if varywhichcm=='a':
        cm_soma = varycm[i]
        cm_dend = varycm[i]
        cm_axon = varycm[i]
    elif varywhichcm=='s':
        cm_soma = varycm[i]
    else:
        cm_dend = varycm[i]
    # Input file name
    filename = folder +'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) + '_vinit'+str(v_init)+'_addedRa_Cmtotal.txt'
    infile = open(filename,'r')
    line = infile.readline()
    cmtot = float(line.split()[0])
    infile.close()
    cm = cmtot/A
    cms[i] = cm
    outfile.write('%.11f %.16f\n' % (varycm[i],cm))
outfile.close()

# Plot results
plt.figure(figsize=(6,5))
plt.plot(varycm,cms,'-o')
plt.xlabel(r'Parameter $C_{m}$ [$\mu$ F/cm$^2$]')
plt.ylabel(r'Output $C_{m}$ [$\mu$ F/cm$^2$]')
plt.title(r'Input $C_{m}$ vs measured $C_{m}$, %s' % plottitletext)
plt.tight_layout()
plt.savefig(plotname)
