import numpy as np
import matplotlib.pyplot as plt
import math

######################################### MULTIPLE MODELS #############################################
######################################### IMPORTING DATA ##############################################
# Varycm: All, soma or dendrite
varywhichcm = 'a' # All

testmodels = [478513437,478513407,480633479,488462965]
iamps      = [0.41,0.41,0.41,0.41]#[0.41,0.17,0.41,0.41] # Should be 0.17 for 478513407
Nmodels    = len(testmodels)
idur   = 1000 # ms
idelay = 100
iamp   = 0.41 # nA
v_init = -86.5 # mV

subfolder = 'Varycm_soma/'
textsnippet = 'Cmsoma'
plottitletext = 'soma'

outfolder = 'Comparemodels/'
cellmodelstring = ''
for i in range(Nmodels):
    cellmodelstring = cellmodelstring +str(testmodels[i])+'_' # Will get a trailing underscore, but whatever
cellmodelstring = cellmodelstring + 'ball-and-stick'
outfolder = outfolder + cellmodelstring +'/'

# No use for text outputfiles, I think. Maybe?
plotname_Nspikes    = outfolder+'cellmodels'+cellmodelstring+'current_idur%i_Nspikes_vs_'% idur+textsnippet+'.png'
plotname_APampl     = outfolder+'cellmodels'+cellmodelstring+'current_idur%i_APampl_vs_'% idur+textsnippet+'.png'
plotname_APdhw      = outfolder+'cellmodels'+cellmodelstring+'current_idur%i_APdurhalfwidth_vs_'% idur+textsnippet+'.png' 

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
    infilename_Nspikes = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_Nspikes_vs_'+textsnippet+'.txt'
    infilename_APampl  = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APampl_vs_'+textsnippet+'.txt'
    infilename_APdhw   = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APdurhalfwidth_vs_'+textsnippet+'.txt'
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

## Extract from ball-and-stick
cms_bas = []
Npeaks_bas = []
APampl_bas = []
APdhw_bas = []
APampl_rms_bas = []
APdhw_rms_bas  = []

# Set names
folder = 'Ball_and_stick/Results/IStim/current_idur%.1f_iamp'%idur+str(iamp)+'/'
infilename_Nspikes = folder+'bas_idur%.1f_iamp'% (idur)+str(iamp) +'_Nspikes_vs_Cm.txt'
infilename_APampl  = folder+'bas_idur%.1f_iamp'% (idur)+str(iamp) +'_APampl_vs_Cm.txt'
infilename_APdhw   = folder+'bas_idur%.1f_iamp'% (idur)+str(iamp) +'_APdurhalfwidth_vs_Cm.txt'

infile_Nspikes = open(infilename_Nspikes,'r')
infile_APampl = open(infilename_APampl,'r')
infile_APdhw = open(infilename_APdhw,'r')

lines_Nspikes = infile_Nspikes.readlines()
lines_APampl = infile_APampl.readlines()
lines_APdhw = infile_APdhw.readlines()
Nlines = len(lines_Nspikes)

# Could have made arrays here... well, well

for i in range(Nlines):
    # Npeaks
    words = lines_Nspikes[i].split()
    if len(words)!=0:
        cms_bas.append(float(words[0]))
        Npeaks_bas.append(float(words[1]))
        # AP Amplitude
        words = lines_APampl[i].split()
        APampl_bas.append(float(words[1]))
        APampl_rms_bas.append(float(words[2]))
        # AP duration at half width
        words = lines_APdhw[i].split()
        APdhw_bas.append(float(words[1]))
        APdhw_rms_bas.append(float(words[2]))
infile_Nspikes.close()
infile_APampl.close()
infile_APdhw.close()

#################################### SKELETON FOR SUBPLOTTING #########################################
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10,4)) # ???
###fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4)) # ???
ax1.set_title('Number of spikes vs capacitance of %s' % plottitletext)#, fontsize=28, y=1.03)
for i in range(Nmodels):
    ax1.plot(cms[i], Npeaks[i], '-o', label='Model %i, I=%.2f nA' % (testmodels[i],iamps[i]))
ax1.plot(cms_bas, Npeaks_bas, '-o', label='Ball-and-stick')
ax1.set(xlabel=r'$C_{m}$ of %s [$\mu$F/cm$^2$]' % plottitletext, ylabel='Number of spikes')#, fontsize=20)
ax1.legend(loc='center right')

ax2.set_title(r'Amplitude vs capacitance of %s' % plottitletext)
for i in range(Nmodels):
    ax2.errorbar(cms[i], APampl[i], yerr=APampl_rms[i], capsize=2, label='Model %i, I=%.2f nA' % (testmodels[i],iamps[i]))
ax2.errorbar(cms_bas, APampl_bas, yerr=APampl_rms_bas, capsize=2, label='Ball-and-stick')
ax2.set(xlabel=r'$C_{m}$ of %s [$\mu$F/cm$^2$]' % plottitletext, ylabel='Amplitude [mV]')#, fontsize=20)

ax3.set_title('AP width at half amplitude vs capacitance of %s' % plottitletext)
for i in range(Nmodels):
    ax3.errorbar(cms[i], APdhw[i], yerr=APdhw_rms[i], capsize=2, label='Model %i, I=%.2f nA' % (testmodels[i],iamps[i]))
ax3.errorbar(cms_bas, APdhw_bas, yerr=APdhw_rms_bas, capsize=2, label='Ball-and-stick')
ax3.set(xlabel=r'$C_{m}$ of %s [$\mu$F/cm$^2$]' % plottitletext, ylabel='AP width at half amplitude [ms]')

######################################### AVERAGE AND RMS #############################################
######################################### IMPORTING DATA ##############################################
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
Nampl  = []
Ndur   = []
APampl = []
APdhw = []
APampl_rms = []
APdhw_rms  = []
# Input: Cm feature (rms)
for j in range(Nmodels):
    iamp = iamps[j]
    cm_temp = []
    Npeaks_temp = []
    Nampl_temp  = []
    Ndur_temp = []
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
    # Extract values ############################
    for i in range(len(lines_Nspikes)):
        # Npeaks
        words = lines_Nspikes[i].split()
        cm_this = float(words[0])
        Npeaks_this = float(words[1])
        # AP Amplitude
        words = lines_APampl[i].split()
        APampl_this = float(words[1])
        APampl_rms_this = float(words[2])
        Nampl_this      = float(words[3])
        # AP duration at half width
        words = lines_APdhw[i].split()
        APdhw_this = float(words[1])
        APdhw_rms_this = float(words[2])
        Ndur_this = float(words[3])
        # Append to array
        cm_temp.append(cm_this)
        Npeaks_temp.append(Npeaks_this)
        Nampl_temp.append(Nampl_this)
        Ndur_temp.append(Ndur_this)
        APampl_temp.append(APampl_this)
        APampl_rms_temp.append(APampl_rms_this)
        APdhw_temp.append(APdhw_this)
        APdhw_rms_temp.append(APdhw_rms_this)
    cms.append(cm_temp)
    Npeaks.append(Npeaks_temp)
    Nampl.append(Nampl_temp)
    Ndur.append(Ndur_temp)
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
Nampl_avg  = np.zeros(Ncms)
Ndur_avg   = np.zeros(Ncms)
APampl_avg = np.zeros(Ncms)
APdhw_avg  = np.zeros(Ncms)
Npeaks_rms = np.zeros(Ncms)
APampl_all_rms = np.zeros(Ncms)
APdhw_all_rms  = np.zeros(Ncms)

for i in range(Nmodels):
    Npeaks_these     = Npeaks[i]
    Nampl_these      = Nampl[i]
    Ndur_these       = Ndur[i]
    APampl_these     = APampl[i]
    APdhw_these      = APdhw[i]
    APampl_these_rms = APampl_rms[i]
    APdhw_these_rms  = APdhw_rms[i]
    for j in range(Ncms):
        Npeaks_this = Npeaks_these[j]
        Nampl_this  = Nampl_these[j]
        Ndur_this   = Ndur_these[j]
        APampl_this = APampl_these[j]
        APdhw_this  = APdhw_these[j]
        Npeaks_avg[j] += Npeaks_this
        Nampl_avg[j]  += Nampl_this
        Ndur_avg[j]   += Ndur_this
        APampl_avg[j] += APampl_this*Nampl_this # Weight it
        APdhw_avg[j]  += APdhw_this*Ndur_this   # Weight it
        APampl_all_rms[j] += (APampl_these_rms[j]*Nampl_this)**2
        APdhw_all_rms[j]  += (APdhw_these_rms[j]*Ndur_this)**2
        Npeaks_all[j,i] = Npeaks_this


for i in range(Ncms):
    APampl_avg[i] /= Nampl_avg[i] # Divide by number of peaks
    APdhw_avg[i]  /= Ndur_avg[i]
    APampl_all_rms[i] /= Nampl_avg[i]**2
    APdhw_all_rms[i]  /= Ndur_avg[i]**2
    APampl_all_rms[i]  = np.sqrt(APampl_all_rms[i])
    APdhw_all_rms[i]   = np.sqrt(APdhw_all_rms[i])
    Npeaks_avg[i] /= Nmodels       # Get average number of peaks
    Nampl_avg[i]  /= Nmodels       # Get average number of amplitudes ## Don't think I need these two
    Ndur_avg[i]   /= Nmodels       # Get average number of durations  ## Don't think I need these two
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

ax4.set_title('Number of spikes vs capacitance of %s' % plottitletext)
ax4.errorbar(cms, Npeaks_avg, yerr=Npeaks_rms, capsize=2)
ax4.set(xlabel=r'$C_{m}$ of %s [$\mu$F/cm$^2$]',ylabel='Number of spikes')





#ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#fig.tight_layout()
#plt.legend(loc='lower right')
plt.show()
#fig.savefig(plotname_twoinone)

plt.show()