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
fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4, 2, figsize=(20,16)) # ???
ax1.set_title('Number of spikes vs capacitance of %s for different models' % plottitletext)#, fontsize=28, y=1.03)
for i in range(Nmodels):
    ax1.plot(cms[i], Npeaks[i], '-o', label='Model %i, I=%.2f nA' % (testmodels[i],iamps[i]))
#ax1.plot(cms_bas, Npeaks_bas, '-o', label='Ball-and-stick')
ax1.set(ylabel='Number of spikes')#, fontsize=20)
ax1.legend(loc='lower left')

ax3.set_title(r'Amplitude vs capacitance of %s for different models' % plottitletext)
for i in range(Nmodels):
    ax3.errorbar(cms[i], APampl[i], yerr=APampl_rms[i], capsize=2, label='Model %i, I=%.2f nA' % (testmodels[i],iamps[i]))
#ax3.errorbar(cms_bas, APampl_bas, yerr=APampl_rms_bas, capsize=2, label='Ball-and-stick')
ax3.set(ylabel='Amplitude [mV]')#, fontsize=20)

ax5.set_title('AP width at half amplitude vs capacitance of %s for different models' % plottitletext)
for i in range(Nmodels):
    ax5.errorbar(cms[i], APdhw[i], yerr=APdhw_rms[i], capsize=2, label='Model %i, I=%.2f nA' % (testmodels[i],iamps[i]))
#ax5.errorbar(cms_bas, APdhw_bas, yerr=APdhw_rms_bas, capsize=2, label='Ball-and-stick')
ax5.set(ylabel='AP width at half amplitude [ms]')

######################################### AVERAGE AND RMS #############################################
######################################### IMPORTING DATA ##############################################
# Varycm: All, soma or dendrite
varywhichcm = 'a' # All

testmodels = [478513437,478513407,488462965]#[478513437,478513407,480633479,488462965]
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



#################################### SUBPLOTTING, PART 2 ###############################################

plottitletext = 'soma'

ax2.set_title('Average number of spikes vs capacitance of %s' % plottitletext)
ax2.errorbar(cms, Npeaks_avg, yerr=Npeaks_rms, capsize=2)
ax2.set(ylabel='Number of spikes')


ax4.set_title(r'Average amplitude vs capacitance of %s' % plottitletext)
ax4.errorbar(cms, APampl_avg, yerr=APampl_all_rms, capsize=2)
ax4.set(ylabel='Amplitude [mV]')


ax6.set_title(r'Average AP width at half amplitude vs capacitance of %s' % plottitletext)
ax6.errorbar(cms, APdhw_avg, yerr=APdhw_all_rms, capsize=2)
ax6.set(ylabel='AP width at half amplitude [ms]')



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


ax7.set_title(r'Signal velocity vs capacitance of %s for different models' % plottitletext)
for i in range(Nmodels):
    ax7.errorbar(cms[i], propvels[i], yerr=propvels_rms[i], capsize=2, label='Model %i' % testmodels[i])
#ax7.plot(cms_bas, propvels_bas, '-o', label='Ball-and-stick')
ax7.set(xlabel=r'$C_{m}$ of %s [$\mu$ F/cm$^2$]' % plottitletext, ylabel=r'Signal velocity [m/s]')
plt.legend(loc='upper right')

######################################### AVG AND RMS ##############################################

# Varycm: All, soma or dendrite
#varywhichcm = 'a' # All # Only for perisomatic so far
varywhichcm = 's' # Soma
#varywhichcm = 'd' # Dendrite

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

##### Adjustable parameters/specifications #########
iduraa = 2
# File/simulation selection:
testmodels = [478513437,478513407,488462965] #[478513437,478513407,480633479,488462965] 
Nmodels = len(testmodels)
idur = 2 # ms
iamp = 1.0 # nA

varycm = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
N = len(varycm)
propvels = np.zeros(N)
propvel_rmss = np.zeros(N)

outfolder = 'Cell_average/'
cellmodelstring = ''
for i in range(Nmodels):
    cellmodelstring = cellmodelstring +str(testmodels[i])+'_' # Will get a trailing underscore, but whatever
outfolder = outfolder + cellmodelstring +'/'

outfilename_avg = outfolder+'idur%i_iamp' % idur+str(iamp)+ '_'+ textsnippet+ '_addedRa_avgrmsallpropvel.txt'
plotname = outfolder +'idur%i_iamp' % idur+str(iamp)+ '_'+ textsnippet+ '_addedRa_avgrmsallpropvel.png'
outfile = open(outfilename_avg,'w')

for i in range(N):
    propvels_all = []
    for k in range(Nmodels):
        testmodel = testmodels[k]
        if testmodel==496497595:
            cm_soma = 1.14805
            cm_dend = 9.98231
            cm_axon = 3.00603
            v_init = -86.5 # mV
        elif testmodel==488462965:
            cm_soma = 3.31732779736
            cm_dend = 3.31732779736
            cm_axon = 3.31732779736
            v_init = -86.5 # mV
        elif testmodel==480633479:
            cm_soma = 0.704866 # 0.704866118957
            cm_dend = 0.704866 # 0.704866118957
            cm_axon = 0.704866 # 0.704866118957
            v_init = -96.8
        elif testmodel==478513407:
            cm_soma = 1.0
            cm_dend = 1.0
            cm_axon = 1.0
            v_init = -83.7
        elif testmodel==478513437:
            cm_soma = 2.34539964752
            cm_dend = 2.34539964752
            cm_axon = 2.34539964752
            v_init = -86.8
        # Update cm
        if varywhichcm=='a':
            cm_soma = varycm[i]
            cm_dend = varycm[i]
            cm_axon = varycm[i]
        elif varywhichcm=='s':
            cm_soma = varycm[i]
        else:
            cm_dend = varycm[i]
        
        folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/dendritepropagation/'
        # Update file name
        infilename = folder+'idur%i_iamp' % idur+str(iamp)+'_cms'+str(cm_soma) + '_cmd'+str(cm_dend)+'_cma'+str(cm_axon)+'_vi'+str(v_init)+ '_wRa_allpropvels.txt' # Needed to cut down on the name length
        infile = open(infilename,'r')
        lines = infile.readlines()
        for line in lines:
            propvel = float(line.split()[0])
            propvels_all.append(propvel)
        infile.close()
    propvels[i],propvel_rmss[i]=avg_and_rms(propvels_all)
    outfile.write('%.2f %.16f %.16f\n' % (varycm[i],propvels[i],propvel_rmss[i]))
outfile.close()

# Plot results
#ax8.set_title(r'Average signal velocity vs %s $C_{m}$' % plottitletext)
ax8.set_title(r'Average signal velocity vs capacitance of %s' % plottitletext)
ax8.errorbar(varycm,propvels,yerr=propvel_rmss,capsize=2)
ax8.set(xlabel=r'$C_{m}$ of %s [$\mu$ F/cm$^2$]' % plottitletext,ylabel='Signal velocity along dendrites [m/s]')

fig.tight_layout()
plt.show()
#fig.savefig(plotname_twoinone)