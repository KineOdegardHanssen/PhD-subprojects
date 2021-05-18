import numpy as np
import matplotlib.pyplot as plt

# Varycm: All, soma or dendrite
#varywhichcm = 'a' # All # Only for perisomatic so far
varywhichcm = 's' # Soma
#varywhichcm = 'd' # Dendrite

testmodels = [478513437,478513407,488462965]
v_inits    = [-86.8,-83.7,-86.5]
Nmodels    = len(testmodels)
idur   = 2 # ms # Different idur?
iduraa = 2
idelay = 20
iamp   = 1.0 # nA

textsnippetp = 'Cmsomaprox'
plottitletext = 'soma and prox. dend.'

outfolder = 'Comparemodels/'
cellmodelstring = ''
for i in range(Nmodels):
    cellmodelstring = cellmodelstring +str(testmodels[i])+'_' # Will get a trailing underscore, but whatever
cellmodelstring = cellmodelstring + 'ball-and-stick'
outfolder = outfolder + cellmodelstring +'/'

plotname  = outfolder+'cellmodels'+cellmodelstring+'current_idur%i_iamp'% idur+str(iamp) +'_propvels_vs_'+textsnippetp+'.png'

cms = []
propvels = []
propvels_rms = []
# Input: Cm feature
for j in range(Nmodels):
    cm_temp = []
    propvels_temp = []
    propvels_rms_temp = []
    v_init = v_inits[j]
    testmodel = testmodels[j]
    if testmodel==496497595: # Will this cause a problem? Run all active for idur=2 too?
        idur = iduraa # ms
    elif testmodel==488462965:
        idur = 2 # ms
    folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp'% (testmodel,idur)+str(iamp)+'/dendritepropagation/'+subfolder
    infilename = folder +'idur%i_iamp' % idur+str(iamp)+ '_cmspr_'+str(testmodel)+ '_'+'_vinit'+str(v_init)+'_addedRa_somaprox_avgpropvel_ends.txt'
    # Open files
    infile = open(infilename,'r')
    # Read lines
    lines = infile.readlines()
    # Extract values
    for i in range(len(lines)):
        # Npeaks
        words = lines[i].split()
        cm_this = float(words[0])
        propvels_this = float(words[1])
        propvels_rms_this = float(words[2])
        # Append to array
        cm_temp.append(cm_this)
        propvels_temp.append(propvels_this)
        propvels_rms_temp.append(propvels_rms_this)
    cms.append(cm_temp)
    propvels.append(propvels_temp)
    propvels_rms.append(propvels_rms_temp)
    infile.close()

plt.figure(figsize=(6,5))
for i in range(Nmodels):
    plt.errorbar(cms[i], propvels[i], yerr=propvels_rms[i], capsize=2, label='Model %i' % testmodels[i])
plt.xlabel(r'$C_{m}$ of %s [$\mu$ F/cm$^2$]' % plottitletext)
plt.ylabel(r'Signal velocity [m/s]')
plt.title(r'Signal velocity vs capacitance of %s' % plottitletext)
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(plotname)
plt.show()
