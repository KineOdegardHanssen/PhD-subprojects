import matplotlib.pyplot as plt
import numpy as np

## EDIT!!!!

# Varycm: All, soma or dendrite
#varywhichcm = 'a' # All # Only for perisomatic so far
varywhichcm = 's' # Soma
#varywhichcm = 'd' # Dendrite

##### Adjustable parameters/specifications #########
smallCm = True
iduraa = 2
# File/simulation selection:
testmodels = [478513437,478513407,488462965]
v_inits    = [-86.8,-83.7,-86.5]
idur = 2 # ms
iamp = 1.0 # nA
v_init = -86.5 # mV
Nmodels = len(testmodels)
Nt      = Nmodels
Ncm     = 21 # Stupid to hard-code this, but I ran into some trouble.
vmax_totavgs = np.zeros(Ncm)
vmax_totrmss = np.zeros(Ncm)

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

cms_storage = []
vmax_avgs_storage = []
vmax_rmss_storage = []
plt.figure(figsize=(6,5))
outfolder = 'Comparemodels/'
cellmodelstring = ''
for i in range(Nmodels):
    print('i:',i)
    cellmodelstring = cellmodelstring +str(testmodels[i])+'_' # Will get a trailing underscore, but whatever
outfolder = outfolder + cellmodelstring +'/'
plotname  = outfolder +'idur%i_iamp' % idur+str(iamp)+ '_allmodels_'+ textsnippet+ '_vinit'+str(v_init)+'_addedRa_avgVmax_smallCm_ends.png'
plotname_avgrms  = outfolder +'idur%i_iamp' % idur+str(iamp)+ '_allmodels_'+ textsnippet+ '_vinit'+str(v_init)+'_addedRa_avgVmax_avgrms_smallCm_ends.png'
for i in range(Nt):
    testmodel = testmodels[i]
    v_init    = v_inits[i]
    folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/dendritepropagation/' # Yeah, this works for all.
    infolder = folder+subfolder
    infilename_avg = outfolder+'idur%i_iamp' % idur+str(iamp)+'_cms' + '_'+str(testmodel)+ '_'+     textsnippet+ '_vinit'+str(v_init)+'_addedRa_avgVmax_ends.txt'
    if smallCm==True:
        infilename_avg = infolder+'idur%i_iamp' % idur+str(iamp)+'_cms' + '_'+str(testmodel)+ '_'+ textsnippet+ '_vinit'+str(v_init)+'_addedRa_avgVmax_smallCm_ends.txt'        
    
    cms       = []
    vmax_avgs = []
    vmax_rmss = []
    
    infile = open(infilename_avg,'r')
    lines = infile.readlines()
    j = 0
    for line in lines:
        words = line.split()
        cms.append(float(words[0]))
        vmax_avg_this = float(words[1])
        vmax_rms_this = float(words[2])
        vmax_avgs.append(vmax_avg_this)
        vmax_rmss.append(vmax_rms_this)
        vmax_totavgs[j] += vmax_avg_this
        vmax_totrmss[j] += vmax_rms_this**2
        j+=1
    infile.close()
    
    plt.errorbar(cms,vmax_avgs,yerr=vmax_rmss,capsize=2,label='%i' % testmodel)
    vmax_avgs_storage.append(vmax_avgs)
    vmax_rmss_storage.append(vmax_rmss)
plt.xlabel(r'$C_{m}$ of %s [$\mu$ F/cm$^2$]' % plottitletext)
plt.ylabel(r'Maximum membrane potential $V_{max}$ in dendrite ends [mV]')
plt.title(r'$V_{max}$ vs %s $C_{m}$' % plottitletext)
plt.tight_layout()
plt.savefig(plotname)
print('Done with read-in and first plot')

print('vmax_avgs_storage:',vmax_avgs_storage)

vmax_totavgs /= Nt
vmax_totrmss  = np.sqrt(vmax_totrmss/(Nt-1))

# Plot results
plt.figure(figsize=(6,5))
plt.errorbar(cms,vmax_totavgs,yerr=vmax_totrmss,capsize=2)
plt.xlabel(r'$C_{m}$ of %s [$\mu$ F/cm$^2$]' % plottitletext)
plt.ylabel(r'Maximum membrane potential $V_{max}$ in dendrite ends [mV]')
plt.title(r'$V_{max}$ vs %s $C_{m}$' % plottitletext)
plt.tight_layout()
plt.savefig(plotname_avgrms)
