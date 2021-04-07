import numpy as np
import matplotlib.pyplot as plt


# Varycm: All, soma or dendrite
varywhichcm = 'a' # All

testmodels = [478513437,478513407,480633479,488462965]
iamps      = [0.41,0.41,0.41,0.41]
Nmodels    = len(testmodels)
idur   = 1000 # ms
idelay = 100
iamp   = 0.41 # nA

subfolder = 'Varycm_soma/'
textsnippet = 'Cmsoma'
plottitletext = 'soma'

outfolder = 'Cell_average/'
cellmodelstring = ''
for i in range(Nmodels):
    cellmodelstring = cellmodelstring +str(testmodels[i])+'_' # Will get a trailing underscore, but whatever
outfolder = outfolder + cellmodelstring +'/'

# No use for text outputfiles, I think. Maybe?
plotname_Nspikes    = outfolder+'cellmodels'+cellmodelstring+'current_idur%i_Nspikes_vs_'% idur+textsnippet+'_smallCm_avgrms.png'
plotname_APampl     = outfolder+'cellmodels'+cellmodelstring+'current_idur%i_APampl_vs_'% idur+textsnippet+'_smallCm_avgrms.png'
plotname_APdhw      = outfolder+'cellmodels'+cellmodelstring+'current_idur%i_APdurhalfwidth_vs_'% idur+textsnippet+'_smallCm_avgrms.png' 

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
    folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp'% (testmodel,idur)+str(iamp)+'/'+subfolder
    infilename_Nspikes = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_Nspikes_vs_'+textsnippet+'_smallCm.txt'
    infilename_APampl  = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APampl_vs_'+textsnippet+'_smallCm.txt'
    infilename_APdhw   = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APdurhalfwidth_vs_'+textsnippet+'_smallCm.txt'
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


cms = cms[0] # Run for the same values
Ncms = len(cms)
Npeaks_all = np.zeros((Ncms,Nmodels))
Npeaks_avg = np.zeros(Ncms)
APampl_avg = np.zeros(Ncms)
APdhw_avg  = np.zeros(Ncms)
Npeaks_rms = np.zeros(Ncms)
APampl_all_rms = np.zeros(Ncms)
APdhw_all_rms  = np.zeros(Ncms)

for i in range(Nmodels):
    Npeaks_these     = Npeaks[i]
    APampl_these     = APampl[i]
    APdhw_these      = APdhw[i]
    APampl_these_rms = APampl_rms[i]
    APdhw_these_rms  = APdhw_rms[i]
    for j in range(Ncms):
        Npeaks_this = Npeaks_these[j]
        APampl_this = APampl_these[j]
        APdhw_this  = APdhw_these[j]
        Npeaks_avg[j] += Npeaks_this
        APampl_avg[j] += APampl_this*Npeaks_this # Weight it
        APdhw_avg[j]  += APdhw_this*Npeaks_this  # Weight it
        APampl_all_rms[j] += (APampl_these_rms[j]*Npeaks_this)**2
        APdhw_all_rms[j]  += (APdhw_these_rms[j]*Npeaks_this)**2
        Npeaks_all[j,i] = Npeaks_this


for i in range(Ncms):
    APampl_avg[i] /= Npeaks_avg[i] # Divide by number of peaks
    APdhw_avg[i]  /= Npeaks_avg[i]
    APampl_all_rms[i] /= Npeaks_avg[i]**2
    APdhw_all_rms[i]  /= Npeaks_avg[i]**2
    APampl_all_rms[i]  = np.sqrt(APampl_all_rms[i])
    APdhw_all_rms[i]   = np.sqrt(APdhw_all_rms[i])
    Npeaks_avg[i] /= Nmodels       # Get average number of peaks
    Npeaks_temp = Npeaks_all[i,:]
    for j in range(Nmodels):
        Npeaks_rms[i] += (Npeaks_avg[i]-Npeaks_temp[j])**2
    Npeaks_rms[i] = np.sqrt(Npeaks_rms[i]/(Nmodels-1))

plt.figure(figsize=(6,5))
plt.errorbar(cms, Npeaks_avg, yerr=Npeaks_rms, capsize=2)
plt.xlabel(r'$C_{m}$ of %s [$\mu$F/cm$^2$]' % plottitletext)
plt.ylabel(r'Number of spikes')
plt.title(r'Number of spikes vs capacitance of %s' % plottitletext)
plt.tight_layout()
plt.savefig(plotname_Nspikes)

plt.figure(figsize=(6,5))
plt.errorbar(cms, APampl_avg, yerr=APampl_all_rms, capsize=2)
plt.xlabel(r'$C_{m}$ of %s [$\mu$F/cm$^2$]' % plottitletext)
plt.ylabel(r'Amplitude [mV]')
plt.title(r'Amplitude vs capacitance of %s' % plottitletext)
plt.tight_layout()
plt.savefig(plotname_APampl)

plt.figure(figsize=(6,5))
plt.errorbar(cms, APdhw_avg, yerr=APdhw_all_rms, capsize=2)
plt.xlabel(r'$C_{m}$ of %s [$\mu$F/cm$^2$]' % plottitletext)
plt.ylabel(r'AP width at half amplitude [ms]')
plt.title(r'AP width at half amplitude vs capacitance of %s' % plottitletext)
plt.tight_layout()
plt.savefig(plotname_APdhw)
plt.show()
