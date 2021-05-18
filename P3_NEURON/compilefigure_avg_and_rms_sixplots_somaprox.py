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

######################################### AVERAGE AND RMS #############################################
######################################### IMPORTING DATA ##############################################
# Varycm: All, soma or dendrite
varywhichcm = 'a' # All

testmodels = [478513437,478513407,488462965]
iamps      = [0.41,0.41,0.41]
Nmodels    = len(testmodels)
idur   = 1000 # ms
idelay = 100
iamp   = 0.41 # nA

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
    folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp'% (testmodel,idur)+str(iamp)+'/'
    infilename_Nspikes = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_Nspikes_vs_Cmsomaprox.txt'
    infilename_APampl  = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APampl_vs_Cmsomaprox.txt'
    infilename_APdhw   = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APdurhalfwidth_vs__Cmsomaprox.txt'
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



#################################### SUBPLOTTING, PART 1 ###############################################
plottitletext = 'soma'
fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(15,15))

ax1.set_title('Average number of spikes vs capacitance of %s' % plottitletext)
ax1.errorbar(cms, Npeaks_avg, yerr=Npeaks_rms, capsize=2)
ax1.set(ylabel='Number of spikes')


ax2.set_title(r'Average amplitude vs capacitance of %s' % plottitletext)
ax2.errorbar(cms, APampl_avg, yerr=APampl_all_rms, capsize=2)
ax2.set(ylabel='Amplitude [mV]')


ax3.set_title(r'Average AP width at half amplitude vs capacitance of %s' % plottitletext)
ax3.errorbar(cms, APdhw_avg, yerr=APdhw_all_rms, capsize=2)
ax3.set(ylabel='AP width at half amplitude [ms]')

######################################### DENDPROP ####################################################
######################################### AVG AND RMS ##############################################

textsnippet = 'Cmsomaprox'
textsnippetp = 'Cmsomaprox'
plottitletext = 'soma and dend. prop.'

##### Adjustable parameters/specifications #########
iduraa = 2
# File/simulation selection:
testmodels = [478513437,478513407,488462965]
Nmodels = len(testmodels)
idur = 2 # ms
iamp = 1.0 # nA

varycm = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0]
N = len(varycm)
propvels = np.zeros(N)
propvel_rmss = np.zeros(N)

outfolder = 'Cell_average/'
cellmodelstring = ''
for i in range(Nmodels):
    cellmodelstring = cellmodelstring +str(testmodels[i])+'_' # Will get a trailing underscore, but whatever
outfolder = outfolder + cellmodelstring +'/'

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
        cm_soma = varycm[i]
        
        folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/dendritepropagation/'
        # Update file name
        infilename = folder+'idur%i_iamp' % idur+str(iamp)+'_cms'+str(cm_soma) + '_cmd'+str(cm_dend)+'_cma'+str(cm_axon)+'_vi'+str(v_init)+ '_wRa_somaprox_allpropvels_ends.txt' # Needed to cut down on the name length
        infile = open(infilename,'r')
        lines = infile.readlines()
        for line in lines:
            propvel = float(line.split()[0])
            propvels_all.append(propvel)
        infile.close()
    propvels[i],propvel_rmss[i]=avg_and_rms(propvels_all)

# Plot results
ax4.set_title(r'Average signal velocity vs capacitance of %s' % plottitletext)
ax4.errorbar(varycm,propvels,yerr=propvel_rmss,capsize=2)
ax4.set(ylabel='Signal velocity along dendrites [m/s]')

########################### Vmax, dend ends ###################################################
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

textsnippet = 'Cmsomaprox'
plottitletext = 'soma and dend. prop.'

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
    
    vmax_avgs_storage.append(vmax_avgs)
    vmax_rmss_storage.append(vmax_rmss)
print('Done with read-in and first plot')

print('vmax_avgs_storage:',vmax_avgs_storage)

Ncm  = len(cms) # Assuming same length
vmax_totavgs = np.zeros(Ncm)
vmax_totrmss = np.zeros(Ncm)
print('Before loop: vmax_totavgs:',vmax_totavgs,'; vmax_totrmss:',vmax_totrmss)
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
ax5.errorbar(cms,vmax_totavgs,yerr=vmax_totrmss,capsize=2)
ax5.set(xlabel=r'$C_{m}$ of %s [$\mu$F/cm$^2$]' % plottitletext,ylabel=r'Maximum membrane potential $V_{max}$ in dendrite ends [mV]')
ax5.set_title(r'Avg. $V_{max}$ in dendrite ends vs %s$C_{m}$' % plottitletext)

################################ SPIKING THRESHOLD ###################################################
# Could have done this better, thb.

# Cm: Soma
cm = np.array([0.01,0.1,0.5,1.0,2.0])
spthr_478513437 = np.array([0.1995,0.1997,0.2002,0.2008,0.2019])
spthr_488462965 = np.array([0.1406,0.1407,0.1410,0.1413,0.1420])
spthr_478513407 = np.array([0.1509,0.1510,0.1511,0.1514,0.1519])

Ncm = len(cm)

spthr_avg = (spthr_478513437+spthr_488462965+spthr_478513407)/3.
spthr_rms = np.zeros(Ncm)
for i in range(Ncm):
    spthr_rms[i] = (spthr_avg[i]-spthr_478513437[i])**2+(spthr_avg[i]-spthr_488462965[i])**2+(spthr_avg[i]-spthr_478513407[i])**2
spthr_rms = np.sqrt(spthr_rms/(Ncm-1))


ax6.errorbar(cm,spthr_avg,yerr=spthr_rms,capsize=2)
ax6.set(xlabel='Soma capacitance, $C_m$ [$\mu$F/cm$^2$]',ylabel='Spiking threshold [nA]')
ax6.set_title('Spiking threshold vs soma capacitance')

fig.tight_layout()
plt.show()