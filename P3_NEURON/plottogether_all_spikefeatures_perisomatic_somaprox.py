import numpy as np
import matplotlib.pyplot as plt


# Varycm: All, soma or dendrite
varywhichcm = 'a' # All

testmodels = [478513437,478513407,488462965]
iamps      = [0.41,0.41,0.41]
Nmodels    = len(testmodels)
idur   = 1000 # ms
idelay = 100
iamp   = 0.41 # nA
v_init = -86.5 # mV

plottitletext = 'soma and proximal dendrites'

outfolder = 'Comparemodels/'
cellmodelstring = ''
for i in range(Nmodels):
    cellmodelstring = cellmodelstring +str(testmodels[i])+'_' # Will get a trailing underscore, but whatever
cellmodelstring = cellmodelstring + 'ball-and-stick'
outfolder = outfolder + cellmodelstring +'/'

# No use for text outputfiles, I think. Maybe?
plotname_Nspikes    = outfolder+'cellmodels'+cellmodelstring+'current_idur%i_Nspikes_vs_'% idur+'_somaprox.png'
plotname_APampl     = outfolder+'cellmodels'+cellmodelstring+'current_idur%i_APampl_vs_'% idur+'_somaprox.png'
plotname_APdhw      = outfolder+'cellmodels'+cellmodelstring+'current_idur%i_APdurhalfwidth_vs_'% idur+'_somaprox.png' 

cms = []
Npeaks = []
APampl = []
APdhw = []
APampl_rms = []
APdhw_rms  = []
# Input: Cm feature (rms)
for j in range(Nmodels):
    iamp = iamps[j]
    cm_temp = []
    Npeaks_temp = []
    APampl_temp = []
    APdhw_temp  = []
    APampl_rms_temp = []
    APdhw_rms_temp  = []
    testmodel = testmodels[j]
    folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp'% (testmodel,idur)+str(iamp)+'/'
    infilename_Nspikes = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_Nspikes_vs_Cmsomaprox.txt'
    infilename_APampl  = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APampl_vs_Cmsomaprox.txt'
    infilename_APdhw   = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APdurhalfwidth_vs_Cmsomaprox.txt'
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
    plt.plot(cms[i], Npeaks[i], '-o', label='Model %i, I=%.2f nA' % (testmodels[i],iamps[i]))
plt.xlabel(r'$C_{m}$ of %s [$\mu$F/cm$^2$]' % plottitletext)
plt.ylabel(r'Number of spikes')
plt.title(r'Number of spikes vs capacitance of %s' % plottitletext)
plt.legend(loc='center right')
plt.tight_layout()
plt.savefig(plotname_Nspikes)

plt.figure(figsize=(6,5))
for i in range(Nmodels):
    plt.errorbar(cms[i], APampl[i], yerr=APampl_rms[i], capsize=2, label='Model %i, I=%.2f nA' % (testmodels[i],iamps[i]))
plt.xlabel(r'$C_{m}$ of %s [$\mu$F/cm$^2$]' % plottitletext)
plt.ylabel(r'Amplitude [mV]')
plt.title(r'Amplitude vs capacitance of %s' % plottitletext)
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(plotname_APampl)

plt.figure(figsize=(6,5))
for i in range(Nmodels):
    plt.errorbar(cms[i], APdhw[i], yerr=APdhw_rms[i], capsize=2, label='Model %i, I=%.2f nA' % (testmodels[i],iamps[i]))
plt.xlabel(r'$C_{m}$ of %s [$\mu$F/cm$^2$]' % plottitletext)
plt.ylabel(r'AP width at half amplitude [ms]')
plt.title(r'AP width at half amplitude vs capacitance of %s' % plottitletext)
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(plotname_APdhw)
plt.show()
