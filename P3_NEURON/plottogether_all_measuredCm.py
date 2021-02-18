import numpy as np
import matplotlib.pyplot as plt

# Varycm: All, soma or dendrite
#varywhichcm = 'a' # All # Only for perisomatic so far
varywhichcm = 's' # Soma
#varywhichcm = 'd' # Dendrite

testmodels = [488462965,496497595]
Nmodels    = len(testmodels)
idur   = 1000 # ms
idelay = 100
iamp   = -0.5 # nA
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

outfolder = 'Comparemodels/'
cellmodelstring = ''
for i in range(Nmodels):
    cellmodelstring = cellmodelstring +str(testmodels[i])+'_' # Will get a trailing underscore, but whatever
cellmodelstring = cellmodelstring + 'ball-and-stick'
outfolder = outfolder + cellmodelstring +'/'


plotname = outfolder +'idur%i_iamp' % idur+str(iamp)+ '_'+cellmodelstring+ textsnippet+ '_vinit'+str(v_init)+'_addedRa_Cm.png'

cms_in  = []
cms_out = []
# Input: Cm_in Cm_out
for j in range(Nmodels):
    cm_in_temp = []
    cm_out_temp = []
    testmodel = testmodels[j]
    folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp'% (testmodel,idur)+str(iamp)+'/'+subfolder
    infilename = folder +'idur%i_iamp' % idur+str(iamp)+ '_'+str(testmodel)+ '_'+ textsnippet+ '_vinit'+str(v_init)+'_addedRa_Cm.txt'
    # Open files
    infile = open(infilename,'r')
    # Read lines
    lines = infile.readlines()
    # Extract values
    for i in range(len(lines)):
        # Npeaks
        words = lines[i].split()
        cm_in_this = float(words[0])
        cm_out_this = float(words[1])
        # Append to array
        cm_in_temp.append(cm_in_this)
        cm_out_temp.append(cm_out_this)
    cms_in.append(cm_in_temp)
    cms_out.append(cm_out_temp)
    infile.close()

# Finding results from ball-and-stick model:
## Extract from ball-and-stick
cms_in_bas  = []
cms_out_bas = []

# Set names
Ra = 150
v_init = -65 
folder = 'Ball_and_stick/Results/IStim/current_idur%.1f_iamp'%idur+str(iamp)+'/'
infilename = folder +'bas_varycm_idur%.1f_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_Cm.txt' 
infile = open(infilename,'r')

lines = infile.readlines()
Nlines = len(lines)

# Could have made arrays here... well, well

for i in range(Nlines):
    # Npeaks
    words = lines[i].split()
    if len(words)!=0:
        cms_in_bas.append(float(words[0]))
        cms_out_bas.append(float(words[1]))
infile.close()

plt.figure(figsize=(6,5))
for i in range(Nmodels):
    plt.plot(cms_in[i], cms_out[i], '-o', label='Model %i' % testmodels[i])
plt.plot(cms_in_bas, cms_out_bas, '-o', label='Ball-and-stick')
plt.xlabel(r'Parameter $C_{m}$ [$\mu$ F/cm$^2$]')
plt.ylabel(r'Output $C_{m}$ [$\mu$ F/cm$^2$]')
plt.title(r'Input $C_{m}$ vs measured $C_{m}$, %s' % plottitletext)
plt.legend(loc='center right')
plt.tight_layout()
plt.savefig(plotname)
