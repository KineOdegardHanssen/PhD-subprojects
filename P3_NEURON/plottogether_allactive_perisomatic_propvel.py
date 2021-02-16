import numpy as np
import matplotlib.pyplot as plt

# Varycm: All, soma or dendrite
#varywhichcm = 'a' # All # Only for perisomatic so far
#varywhichcm = 's' # Soma
varywhichcm = 'd' # Dendrite

testmodels = [488462965,496497595]
Nmodels    = len(testmodels)
idur   = 2 # ms # Different idur?
iduraa = 2
idelay = 20
iamp   = 1.0 # nA
v_init = -86.5 # mV

if varywhichcm=='a': # All is probably not that relevant
    subfolder = 'Varycm_all/'
    textsnippet = 'varycmall'
    textsnippetp = 'Cmall'
    plottitletext = 'whole neuron'
elif varywhichcm=='s':
    subfolder = 'Varycm_soma/'
    textsnippet = 'varycmsoma'
    textsnippetp = 'Cmsoma'
    plottitletext = 'soma'
else:
    subfolder = 'Varycm_dend/'
    textsnippet = 'varycmdend'
    textsnippetp = 'Cmdend'
    plottitletext = 'dendrites'

outfolder = 'Allen_test_changecapacitance/figures/Comparemodels/'
cellmodelstring = ''
for i in range(Nmodels):
    cellmodelstring = cellmodelstring +str(testmodels[i])+'_' # Will get a trailing underscore, but whatever
outfolder = outfolder + cellmodelstring +'/'+subfolder
plotname  = outfolder+'cellmodels'+cellmodelstring+'current_idur%i_iamp'% idur+str(iamp) +'_propvels_vs_'+textsnippetp+'.png'

cms = []
propvels = []
# Input: Cm feature
for j in range(Nmodels):
    cm_temp = []
    propvels_temp = []
    testmodel = testmodels[j]
    if testmodel==496497595: # Will this cause a problem? Run all active for idur=2 too?
        idur = iduraa # ms
    elif testmodel==488462965:
        idur = 2 # ms
    folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp'% (testmodel,idur)+str(iamp)+'/dendritepropagation/'+subfolder
    infilename = folder +'idur%i_iamp' % idur+str(iamp)+ '_'+str(testmodel)+ '_'+ textsnippet+ '_vinit'+str(v_init)+'_addedRa_propvel.txt'
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
        # Append to array
        cm_temp.append(cm_this)
        propvels_temp.append(propvels_this)
    cms.append(cm_temp)
    propvels.append(propvels_temp)
    infile.close()

plt.figure(figsize=(6,5))
for i in range(Nmodels):
    plt.plot(cms[i], propvels[i], '-o', label='Model %i' % testmodels[i])
plt.xlabel(r'$C_{m}$ of %s [$\mu$ F/cm$^2$]' % plottitletext)
plt.ylabel(r'Signal velocity')
plt.title(r'Signal velocity vs capacitance of %s' % plottitletext)
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(plotname)
