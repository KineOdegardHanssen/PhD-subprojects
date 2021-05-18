import matplotlib.pyplot as plt
import numpy as np

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

plottitletext = 'soma and prox. dend.'

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

plotname  = outfolder +'idur%i_iamp' % idur+str(iamp)+ '_allmodels_vinit'+str(v_init)+'_addedRa_somaprox_avgVmax_smallCm_ends.png'
plotname_avgrms  = outfolder +'idur%i_iamp' % idur+str(iamp)+ '_allmodels_vinit'+str(v_init)+'_addedRa_somaprox_avgVmax_avgrms_smallCm_ends.png'
for i in range(Nt):
    testmodel = testmodels[i]
    v_init    = v_inits[i]
    folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/dendritepropagation/' # Yeah, this works for all.
    infolder = folder+subfolder
    infilename_avg = outfolder+'idur%i_iamp' % idur+str(iamp)+'_cmspr' + '_'+str(testmodel)+ '_somaprox_vinit'+str(v_init)+'_addedRa_avgVmax_ends.txt'   
    
    cms       = []
    vmax_avgs = []
    vmax_rmss = []
    
    infile = open(infilename_avg,'r')
    lines = infile.readlines()
    for line in lines:
        words = line.split()
        cms.append(float(words[0]))
        vmax_avgs.append(float(words[1]))
        vmax_rmss.append(float(words[2]))
    infile.close()
    
    plt.errorbar(cms,vmax_avgs,yerr=vmax_rmss,capsize=2,label='%i' % testmodel)
    vmax_avgs_storage.append(vmax_avgs)
    vmax_rmss_storage.append(vmax_rmss)
plt.xlabel(r'$C_{m}$ of %s [$\mu$ F/cm$^2$]' % plottitletext)
plt.ylabel(r'Maximum membrane potential $V_{max}$ in dendrite ends [mV]')
plt.title(r'$V_{max}$ in dendrite ends vs %s $C_{m}$' % plottitletext)
plt.tight_layout()
plt.savefig(plotname)
print('Done with read-in and first plot')

print('vmax_avgs_storage:',vmax_avgs_storage)

Ncm  = len(cms) # Assuming same length
vmax_totavgs = np.zeros(Ncm)
vmax_totrmss = np.zeros(Ncm)
print('Before loop: vmax_totavgs:',vmax_totavgs,'; vmax_totrmss:',vmax_totrmss,)
for i in range(Nt):
    vmax_totavgs_thismodel = vmax_avgs_storage[i]
    vmax_totrmss_thismodel = vmax_rmss_storage[i]
    print('vmax_totavgs_thismodel:',vmax_totavgs_thismodel)
    print('vmax_totrmss_thismodel:',vmax_totrmss_thismodel)
    for j in range(Ncm):
        print('i:',i,'; j:',j)
        print('vmax_totavgs_thismodel[j]:',vmax_totavgs_thismodel[j])
        print('vmax_totrmss_thismodel[j]:',vmax_totrmss_thismodel[j])
        vmax_totavgs[j] += vmax_totavgs_thismodel[j]
        vmax_totrmss[j] += vmax_totrmss_thismodel[j]**2
        print('In loop: vmax_totavgs:',vmax_totavgs,'; vmax_totrmss:',vmax_totrmss,)
vmax_totavgs /= Nt
vmax_totrmss = np.sqrt(vmax_totrmss/(Nt-1))

    
# Plot results
plt.figure(figsize=(6,5))
plt.errorbar(cms,vmax_totavgs,yerr=vmax_totrmss,capsize=2)
plt.xlabel(r'$C_{m}$ of %s [$\mu$ F/cm$^2$]' % plottitletext)
plt.ylabel(r'Maximum membrane potential $V_{max}$ in dendrite ends [mV]')
plt.title(r'Avg. $V_{max}$ in dendrite ends vs %s $C_{m}$' % plottitletext)
plt.tight_layout()
plt.savefig(plotname_avgrms)