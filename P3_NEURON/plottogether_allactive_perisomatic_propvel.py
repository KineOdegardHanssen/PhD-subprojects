import numpy as np
import matplotlib.pyplot as plt

# Varycm: All, soma or dendrite
#varywhichcm = 'a' # All # Only for perisomatic so far
#varywhichcm = 's' # Soma
varywhichcm = 'd' # Dendrite

testmodels = [488462965,496497595]
Nmodels    = len(testmodels)
idur   = 1000 # ms # Different idur?
idelay = 100
iamp   = 0.41 # nA
v_init = -86.5 # mV

if varywhichcm=='a': # All is probably not that relevant
    subfolder = 'Varycm_all/'
    textsnippet = 'Cmall'
    plottitletext = 'whole neuron'
elif varywhichcm=='s':
    subfolder = 'Varycm_soma/'
    textsnippet = 'Cmsoma'
    plottitletext = 'soma'
else:
    subfolder = 'Varycm_dend/'
    textsnippet = 'Cmdend'
    plottitletext = 'dendrites'

outfolder = 'Allen_test_changecapacitance/figures/Comparemodels/'
cellmodelstring = ''
for i in range(Nmodels):
    cellmodelstring = cellmodelstring +str(testmodels[i])+'_' # Will get a trailing underscore, but whatever
outfolder = outfolder + cellmodelstring +'/'+subfolder

cms_in = []
cms_out = []
# Input: Cm feature
for j in range(Nmodels):
    cm_in_temp = []
    cm_out_temp = []
    testmodel = testmodels[j]
    folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp'% (testmodel,idur)+str(iamp)+'/'+subfolder
    infilename_Nspikes = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_Nspikes_vs_'+textsnippet+'.txt'
    infilename_APampl  = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APampl_vs_'+textsnippet+'.txt'
    infilename_APdhw   = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APdurhalfwidth_vs_'+textsnippet+'.txt'
    if testmodel==496497595: # Will this cause a problem? Run all active for idur=2 too?
        idur = 1 # ms
    elif testmodel==488462965:
        idur = 2 # ms
    folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp'% (testmodel,idur)+str(iamp)+'/dendritepropagation/'
    infilename = folder +'idur%i_iamp' % idur+str(iamp)+ '_'+str(testmodel)+ '_'+ textsnippet+ '_vinit'+str(v_init)+'_addedRa_propvel.txt'
    # Open files
    infile_Nspikes = open(infilename_Nspikes,'r')
    infile_APampl  = open(infilename_APampl,'r')
    infile_APdhw   = open(infilename_APdhw,'r')
    # Read lines
    lines_Nspikes = infile_Nspikes.readlines()
    lines_APampl = infile_APampl.readlines()
    lines_APdhw = infile_APdhw.readlines()
    # Extract values
    for i in range(len(lines_Nspikes)):
        # Npeaks
        words = lines_Nspikes[i].split()
        cm_this = float(words[0])
        Npeaks_this = float(words[1])
        # AP Amplitude
        words = lines_APampl[i].split()
        APampl_this = float(words[1])
        APampl_rms_this = float(words[2])
        # AP duration at half width
        words = lines_APdhw[i].split()
        APdhw_this = float(words[1])
        APdhw_rms_this = float(words[2])
        # Append to array
        cm_temp.append(cm_this)
        Npeaks_temp.append(Npeaks_this)
        APampl_temp.append(APampl_this)
        APdhw_temp.append(APdhw_this)
        APampl_rms_temp.append(APampl_rms_this)
        APdhw_rms_temp.append(APdhw_rms_this)
    cms.append(cm_temp)
    Npeaks.append(Npeaks_temp)
    APampl.append(APampl_temp)
    APdhw.append(APdhw_temp)
    APampl_rms.append(APampl_rms_temp)
    APdhw_rms.append(APdhw_rms_temp)
    infile_Nspikes.close()
    infile_APampl.close()
    infile_APdhw.close()

plt.figure(figsize=(6,5))
for i in range(Nmodels):
    plt.plot(cms[i], Npeaks[i], '-o', label='Model %i' % testmodels[i])
plt.xlabel(r'$C_{m}$ of %s [$\mu$ F/cm$^2$]' % plottitletext)
plt.ylabel(r'Number of spikes')
plt.title(r'Number of spikes vs capacitance of %s' % plottitletext)
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(plotname_Nspikes)
