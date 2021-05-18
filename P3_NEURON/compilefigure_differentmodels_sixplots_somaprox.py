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
testmodels = [478513437,478513407,488462965]
iamps      = [0.41,0.41,0.41,0.41]
Nmodels    = len(testmodels)
idur   = 1000 # ms
idelay = 100
iamp   = 0.41 # nA
v_init = -86.5 # mV

textsnippet = 'Cmsomaprox'
plottitletext = 'soma and prox. dend.'

outfolder = 'Comparemodels/'
cellmodelstring = ''
for i in range(Nmodels):
    cellmodelstring = cellmodelstring +str(testmodels[i])+'_' # Will get a trailing underscore, but whatever
cellmodelstring = cellmodelstring + 'ball-and-stick'
outfolder = outfolder + cellmodelstring +'/'

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

#################################### SUBPLOTTING, PART 1 ###############################################
plottitletext = 'soma'
fig, ((ax1, ax2), (ax3, ax4),(ax5,ax6)) = plt.subplots(3, 2, figsize=(15,15)) # ???
ax1.set_title('Number of spikes vs capacitance of %s for different models' % plottitletext)#, fontsize=28, y=1.03)
for i in range(Nmodels):
    ax1.plot(cms[i], Npeaks[i], '-o', label='Model %i, I=%.2f nA' % (testmodels[i],iamps[i]))
ax1.set(ylabel='Number of spikes')#, fontsize=20)
ax1.legend(loc='lower left')

ax2.set_title(r'Amplitude vs capacitance of %s for different models' % plottitletext)
for i in range(Nmodels):
    ax2.errorbar(cms[i], APampl[i], yerr=APampl_rms[i], capsize=2, label='Model %i, I=%.2f nA' % (testmodels[i],iamps[i]))
ax2.set(ylabel='Amplitude [mV]')#, fontsize=20)

ax3.set_title('AP width at half amplitude vs capacitance of %s for different models' % plottitletext)
for i in range(Nmodels):
    ax3.errorbar(cms[i], APdhw[i], yerr=APdhw_rms[i], capsize=2, label='Model %i, I=%.2f nA' % (testmodels[i],iamps[i]))
ax3.set(ylabel='AP width at half amplitude [ms]')



######################################### DENDPROP ####################################################
######################################## ALL MODELS ###################################################
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

textsnippet = 'Cmsomaprox'
textsnippetp = 'Cmsomaprox'
plottitletext = 'soma and dend. prop.'

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
    folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp'% (testmodel,idur)+str(iamp)+'/dendritepropagation/'
    infilename = folder +'idur%i_iamp' % idur+str(iamp)+'_cmspr' + '_'+str(testmodel)+ '_'+'_vinit'+str(v_init)+'_addedRa_somaprox_avgpropvel_ends.png'
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

ax4.set_title(r'Signal velocity vs capacitance of %s for different models' % plottitletext)
for i in range(Nmodels):
    ax4.errorbar(cms[i], propvels[i], yerr=propvels_rms[i], capsize=2, label='Model %i' % testmodels[i])
ax4.set(ylabel=r'Signal velocity [m/s]')

#################### Vmax, dend ends #################################
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
outfolder = 'Comparemodels/'
cellmodelstring = ''
for i in range(Nmodels):
    print('i:',i)
    cellmodelstring = cellmodelstring +str(testmodels[i])+'_' # Will get a trailing underscore, but whatever
outfolder = outfolder + cellmodelstring +'/'

for i in range(Nt):
    testmodel = testmodels[i]
    v_init    = v_inits[i]
    folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/dendritepropagation/' # Yeah, this works for all.
    infolder = folder
    infilename_avg = outfolder+'idur%i_iamp' % idur+str(iamp)+'_cms' + '_'+str(testmodel)+ '_vinit'+str(v_init)+'_addedRa_somaprox_avgVmax_ends.txt'
    
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
    
    ax5.errorbar(cms,vmax_avgs,yerr=vmax_rmss,capsize=2,label='%i' % testmodel)
    vmax_avgs_storage.append(vmax_avgs)
    vmax_rmss_storage.append(vmax_rmss)
ax5.set(xlabel=r'$C_{m}$ of %s [$\mu$F/cm$^2$]' % plottitletext,ylabel=r'Maximum membrane potential $V_{max}$ in dendrite ends [mV]')
ax5.set_title(r'$V_{max}$ in dendrite ends vs %s $C_{m}$' % plottitletext)

###################### SPIKING THRESHOLD ##############################################################

## Redo!

# Cm: Soma
cm_478513437 = np.array([0.01,0.1,0.5,1.0,2.0])
spthr_478513437 = np.array([0.1995,0.1997,0.2002,0.2008,0.2019])

cm_488462965 = np.array([0.01,0.1,0.5,1.0,2.0])
spthr_488462965 = np.array([0.1406,0.1407,0.1410,0.1413,0.1420])

cm_478513407 = np.array([0.01,0.1,0.5,1.0,2.0])
spthr_478513407 = np.array([0.1509,0.1510,0.1511,0.1514,0.1519])

ax6.plot(cm_478513437,spthr_478513437, '-o')
ax6.plot(cm_478513407,spthr_478513407, '-o')
ax6.plot(cm_488462965,spthr_488462965, '-o')
ax6.set(xlabel='Soma capacitance, $C_m$ [$\mu$F/cm$^2$]',ylabel='Spiking threshold [nA]')
ax6.set_title('Spiking threshold vs soma capacitance')


fig.tight_layout()
plt.show()
#fig.savefig(plotname_twoinone)