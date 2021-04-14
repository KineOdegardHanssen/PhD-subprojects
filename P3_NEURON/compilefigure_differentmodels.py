import numpy as np
import matplotlib.pyplot as plt
import math

def avg_and_rms(x):
    N = len(x)
    avgx = np.mean(x)
    print('avgx:',avgx)
    rmsx = 0
    for i in range(N):
        rmsx += (avgx-x[i])**2
    rmsx = np.sqrt(rmsx/(N-1))
    return avgx,rmsx

######################################### MULTIPLE MODELS #############################################
######################################### IMPORTING DATA ##############################################
# Varycm: All, soma or dendrite
varywhichcm = 'a' # All

testmodels = [478513437,478513407,488462965] #[478513437,478513407,480633479,488462965]
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

## Extract from ball-and-stick
cms_bas = []
Npeaks_bas = []
APampl_bas = []
APdhw_bas = []
APampl_rms_bas = []
APdhw_rms_bas  = []

# Set names
folder = 'Ball_and_stick/Results/IStim/current_idur%.1f_iamp'%idur+str(iamp)+'/'
infilename_Nspikes = folder+'bas_idur%.1f_iamp'% (idur)+str(iamp) +'_Nspikes_vs_Cm_smallCm.txt'
infilename_APampl  = folder+'bas_idur%.1f_iamp'% (idur)+str(iamp) +'_APampl_vs_Cm_smallCm.txt'
infilename_APdhw   = folder+'bas_idur%.1f_iamp'% (idur)+str(iamp) +'_APdurhalfwidth_vs_Cm_smallCm.txt'

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

#################################### SUBPLOTTING, PART 1 ###############################################
plottitletext = 'soma'
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(17,15)) # ???
ax1.set_title('Number of spikes vs capacitance of %s for different models' % plottitletext)#, fontsize=28, y=1.03)
for i in range(Nmodels):
    ax1.plot(cms[i], Npeaks[i], '-o', label='Model %i, I=%.2f nA' % (testmodels[i],iamps[i]))
#ax1.plot(cms_bas, Npeaks_bas, '-o', label='Ball-and-stick')
ax1.set(ylabel='Number of spikes')#, fontsize=20)
ax1.legend(loc='lower left')

ax2.set_title(r'Amplitude vs capacitance of %s for different models' % plottitletext)
for i in range(Nmodels):
    ax2.errorbar(cms[i], APampl[i], yerr=APampl_rms[i], capsize=2, label='Model %i, I=%.2f nA' % (testmodels[i],iamps[i]))
ax2.set(ylabel='Amplitude [mV]')#, fontsize=20)

ax3.set_title('AP width at half amplitude vs capacitance of %s for different models' % plottitletext)
for i in range(Nmodels):
    ax3.errorbar(cms[i], APdhw[i], yerr=APdhw_rms[i], capsize=2, label='Model %i, I=%.2f nA' % (testmodels[i],iamps[i]))
ax3.set(xlabel=r'$C_{m}$ of %s [$\mu$F/cm$^2$]' % plottitletext, ylabel='AP width at half amplitude [ms]')



######################################### DENDPROP ####################################################
######################################## ALL MODELS ###################################################
# Varycm: All, soma or dendrite
#varywhichcm = 'a' # All # Only for perisomatic so far
varywhichcm = 's' # Soma
#varywhichcm = 'd' # Dendrite

testmodels = [478513437,478513407,488462965]#[478513437,478513407,480633479,488462965]
v_inits    = [-86.8,-83.7,-86.5]#[-86.8,-83.7,-96.8,-86.5]
Nmodels    = len(testmodels)
idur   = 2 # ms # Different idur?
iduraa = 2
idelay = 20
iamp   = 1.0 # nA

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
    infilename = folder +'idur%i_iamp' % idur+str(iamp)+ '_cms_'+str(testmodel)+ '_'+ textsnippet+ '_vinit'+str(v_init)+'_addedRa_avgpropvel_smallCm.txt'
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

## Read in values from ball-and-stick model
cms_bas  = []
propvels_bas = []

# Set names
Ra = 150
v_init = -65 
folder = 'Ball_and_stick/Results/IStim/current_idur%.1f_iamp'%idur+str(iamp)+'/dendritepropagation/'
infilename = folder+'bas_varycm_idur%.1f_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_propvel_vs_cm_smallCm.txt'
infile = open(infilename,'r')

lines = infile.readlines()
Nlines = len(lines)

# Could have made arrays here... well, well

for i in range(Nlines):
    # Npeaks
    words = lines[i].split()
    if len(words)!=0:
        cms_bas.append(float(words[0]))
        propvels_bas.append(float(words[1]))
infile.close()


ax4.set_title(r'Signal velocity vs capacitance of %s for different models' % plottitletext)
for i in range(Nmodels):
    ax4.errorbar(cms[i], propvels[i], yerr=propvels_rms[i], capsize=2, label='Model %i' % testmodels[i])
ax4.set(xlabel=r'$C_{m}$ of %s [$\mu$F/cm$^2$]' % plottitletext, ylabel=r'Signal velocity [m/s]')
plt.legend(loc='upper right')

fig.tight_layout()
plt.show()
#fig.savefig(plotname_twoinone)