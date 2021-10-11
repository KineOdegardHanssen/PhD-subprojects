import numpy as np
import matplotlib.pyplot as plt

testmodels = [485694403,478513437,489931686,488462965,478513407]
idur       = 1000 #100 # ms
idelay     = 100
iamp       = 0.4 # nA
v_init     = -86.5 # mV
Ra         = 150
somasize   = 10 # 15 # 
dendlen    = 1000
denddiam   = 2
nsegments  = 200 
Nmodels    = len(testmodels)

varymech = 'Kd'#'NaV' #'SK'#'Im_v2'#'Kv2like'#'pas' #
varyE_bool = True
varyE = -107 #[-90,-80]#[-130,-107,-100,-70]#[40,50,53,65,70]#[-60,-53,-30,-10,10,30,60] # Default...407Kd-107 Default478513407pas: -83.6528;  Default478513407Na: 53
varyg = 'None' # Default...407Kd2.94396e-010 # Default478513407pas: 0.000362109; Default478513407Na: 0.0409177

varylist = [] # Should be redundant
plotstring = '_vary'
if varyE_bool==True:
    varylist = varyE
    plotstring = plotstring + 'E'
else:
    varylist = varyg
    plotstring = plotstring + 'g'

if varymech=='NaV':
    folderstring = 'VaryNa/' 
    plotstring   = plotstring + '_NaV'
elif varymech=='pas':
    folderstring = 'VaryPas/'
    plotstring   = plotstring + '_Pas'
elif varymech=='Kd':
    folderstring = 'VaryKd/'
    plotstring   = plotstring + '_Kd'
elif varymech=='Kv2like':
    folderstring = 'VaryKv2like/'
    plotstring   = plotstring + '_Kv2like'
elif varymech=='Kv3_1':
    folderstring = 'VaryKv3_1/'
    plotstring   = plotstring + '_Kv3_1'
elif varymech=='SK':
    folderstring = 'VarySK/'
    plotstring   = plotstring + '_SK'
elif varymech=='K_T':
    folderstring = 'VaryK_T/'
    plotstring   = plotstring + '_K_T'
elif varymech=='Im_v2':
    folderstring = 'VaryIm_v2/'
    plotstring   = plotstring + '_Im_v2'

changestring =''
if varyE_bool==True:
    changestring = changestring+'_E'+str(varyE)+'_gdf'
else:
    changestring = changestring+'_Edf_g'+str(varyg)

cm_master = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0]
    
NCms = len(cm_master)

Nspikes_all          = []
avg_AP_ampl_all      = []
rms_AP_ampl_all      = []
avg_AP_halfwidth_all = []
rms_AP_halfwidth_all = []
Cm_Nspikes_all       = []
Cm_AP_ampl_all       = []
Cm_AP_halfwidth_all  = []
    
for testmodel in testmodels:
    Nspikes          = []
    Cm_Nspikes       = []
    Cm_AP_ampl       = []
    Cm_AP_halfwidth  = []
    avg_AP_ampl      = []
    rms_AP_ampl      = []
    avg_AP_halfwidth = []
    rms_AP_halfwidth = []
    
    # Set names
    infolder = 'Results/%i/IStim/Soma%i/dendlen%i/denddiam'% (testmodel,somasize,dendlen)+str(denddiam)+'/'+folderstring
    currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    infolder = infolder+currentfolder
    infilename_Nspikes = infolder+'baspv_current_idur%i_iamp'% (idur)+str(iamp) +'_cmfs_f_vs_Cmall.txt'
    infilename_APampl  = infolder+'baspv_current_idur%i_iamp'% (idur)+str(iamp) +'_cmfs_APampl_vs_Cmsprx.txt'
    infilename_APdhw   = infolder+'baspv_current_idur%i_iamp'% (idur)+str(iamp) +'_cmfs_APdurhalfwidth_vs_Cmsprx.txt'
    # Read files
    infile_Nspikes = open(infilename_Nspikes,'r')
    infile_APampl  = open(infilename_APampl,'r')
    infile_APdhw   = open(infilename_APdhw,'r')
    
    lines_Nspikes = infile_Nspikes.readlines()
    lines_APampl  = infile_APampl.readlines()
    lines_APdhw   = infile_APdhw.readlines()
    Nlines = len(lines_Nspikes) # All data types have the same length
    
    for i in range(Nlines):
        words_Nspikes = lines_Nspikes[i].split()
        words_APampl  = lines_APampl[i].split()
        words_APdhw   = lines_APdhw[i].split()
        if len(words_Nspikes)>0:
            Cm_Nspikes.append(float(words_Nspikes[0]))
            Nspikes.append(float(words_Nspikes[1]))
        if len(words_APampl)>0: # One test should be enough, but better safe than sorry
            Cm_AP_ampl.append(float(words_APampl[0]))
            avg_AP_ampl.append(float(words_APampl[1]))
            rms_AP_ampl.append(float(words_APampl[2]))
        if len(words_APdhw)>0:
            Cm_AP_halfwidth.append(float(words_APdhw[0]))
            avg_AP_halfwidth.append(float(words_APdhw[1]))
            rms_AP_halfwidth.append(float(words_APdhw[2]))
    
    infile_Nspikes.close()
    infile_APampl.close()
    infile_APdhw.close()
    
    Nspikes_all.append(Nspikes)
    avg_AP_ampl_all.append(avg_AP_ampl)
    rms_AP_ampl_all.append(rms_AP_ampl)
    avg_AP_halfwidth_all.append(avg_AP_halfwidth)
    rms_AP_halfwidth_all.append(rms_AP_halfwidth)
    Cm_Nspikes_all.append(Cm_Nspikes)
    Cm_AP_ampl_all.append(Cm_AP_ampl)
    Cm_AP_halfwidth_all.append(Cm_AP_halfwidth)

## Should I take avg. and rms. of all models? Might be difficult since I allow for different Cms.
# Check against cm_master
Nspikes     = np.zeros(NCms)
Nspikes_rms = np.zeros(NCms)
avg_AP_ampl = np.zeros(NCms)
rms_AP_ampl = np.zeros(NCms)
avg_AP_halfwidth = np.zeros(NCms)
rms_AP_halfwidth = np.zeros(NCms)


Nspikes_bucket   = np.zeros((NCms,Nmodels))


Nspikes_collections = []
avg_AP_ampl_collections = []
rms_AP_ampl_collections = []
avg_AP_halfwidth_collections = []
rms_AP_halfwidth_collections = []

fcounter    = np.zeros(NCms)
amplcounter = np.zeros(NCms)
hwcounter   = np.zeros(NCms)

for i in range(Nmodels):
    Nspikes_this = Nspikes_all[i]
    avg_AP_ampl_this = avg_AP_ampl_all[i]
    rms_AP_ampl_this = rms_AP_ampl_all[i]
    avg_AP_halfwidth_this = avg_AP_halfwidth_all[i]
    rms_AP_halfwidth_this = rms_AP_halfwidth_all[i]
    Cm_Nspikes_this       = Cm_Nspikes_all[i]
    Cm_AP_ampl_this       = Cm_AP_ampl_all[i]
    Cm_AP_halfwidth_this  = Cm_AP_halfwidth_all[i]
    Nspikes_collections = []
    avg_AP_ampl_collections = []
    rms_AP_ampl_collections = []
    avg_AP_halfwidth_collections = []
    rms_AP_halfwidth_collections = []
    for j in range(NCms):
        for k in range(len(Cm_Nspikes_this)):
            if Cm_Nspikes_this[k]==cm_master[j]: # Does this make any sense?
                Nspikes[j] += Nspikes_this[k]
                Nspikes_bucket[j,i] = Nspikes_this[k]
                fcounter[j]+=1 
                continue
    for j in range(NCms):
        for k in range(len(Cm_AP_ampl_this)):
            if Cm_AP_ampl_this[k]==cm_master[j]: # Does this make any sense?
                avg_AP_ampl[j] += avg_AP_ampl_this[k]
                rms_AP_ampl[j] += rms_AP_ampl_this[k]**2 # Square now, square root at the end
                amplcounter[j]+=1 
                continue
    for j in range(NCms):
        for k in range(len(Cm_AP_halfwidth_this)):
            if Cm_AP_halfwidth_this[k]==cm_master[j]: # Does this make any sense?
                avg_AP_halfwidth[j] += avg_AP_halfwidth_this[k]
                rms_AP_halfwidth[j] += rms_AP_halfwidth_this[k]**2 # Square now, square root at the end
                hwcounter[j]+=1 
                continue

# Average and finish rms
for j in range(NCms):
    if fcounter[j]!=0:
        Nspikes[j] /= fcounter[j]
        for i in range(int(fcounter[j])):
            Nspikes_rms[j]+=(Nspikes[j]-Nspikes_bucket[j,i])**2
        Nspikes_rms[j] = np.sqrt(Nspikes_rms[j]/(fcounter[j]-1))
for j in range(NCms):
    if amplcounter[j]!=0:
        avg_AP_ampl[j] /=amplcounter[j]
        rms_AP_ampl[j]  = np.sqrt(rms_AP_ampl[j]/(amplcounter[j]-1))
for j in range(NCms):
    if hwcounter[j]!=0:
        avg_AP_halfwidth[j] /= hwcounter[j]
        rms_AP_halfwidth[j]  = np.sqrt(rms_AP_halfwidth[j]/(hwcounter[j]-1))

# Plotting

plotfolder = 'Results/Comparemodels/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+currentfolder
plotname_all = plotfolder+'features_all.png'
plotname_avgrms = plotfolder+'features_avgrms.png'



## avg and rms:
plottitletext = 'soma and proximal dendrites'
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(17,15))

ax1.set_title(r'$f$ vs $C_m$ of %s' % plottitletext)
ax1.errorbar(cm_master, Nspikes_rms,yerr=Nspikes_rms,capsize=2)
ax1.set(ylabel='$f$')

ax2.set_title(r'Amplitude vs $C_m$ of %s' % plottitletext)
ax2.errorbar(cm_master, avg_AP_ampl, yerr=rms_AP_ampl, capsize=2)
ax2.set(ylabel='Amplitude [mV]')

ax3.set_title(r'Average AP width at half amplitude vs capacitance of %s' % plottitletext)
ax3.errorbar(cm_master, avg_AP_halfwidth, yerr=rms_AP_halfwidth, capsize=2)
ax3.set(xlabel=r'$C_{m}$ of %s [$\mu$ F/cm$^2$]' % plottitletext,ylabel='AP width at half amplitude [ms]')
plt.savefig(plotname_avgrms)

## All:
plottitletext = 'soma and proximal dendrites'
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(17,15))

ax1.set_title(r'$f$ vs $C_m$ of %s' % plottitletext)
for i in range(Nmodels):
    ax1.plot(Cm_Nspikes_all[i], Nspikes_all[i],label='%i' % testmodels[i])
ax1.set(ylabel='$f$')
ax1.legend(loc='upper right')

ax2.set_title(r'Amplitude vs $C_m$ of %s' % plottitletext)
for i in range(Nmodels):
    ax2.errorbar(Cm_AP_ampl_all[i], avg_AP_ampl_all[i], yerr=rms_AP_ampl_all[i], capsize=2,label='%i' % testmodels[i])
ax2.set(ylabel='Amplitude [mV]')
#ax2.legend(loc='upper right')

ax3.set_title(r'Average AP width at half amplitude vs capacitance of %s' % plottitletext)
for i in range(Nmodels):
    ax3.errorbar(Cm_AP_halfwidth_all[i], avg_AP_halfwidth_all[i], yerr=rms_AP_halfwidth_all[i], capsize=2,label='%i' % testmodels[i])
ax3.set(xlabel=r'$C_{m}$ of %s [$\mu$ F/cm$^2$]' % plottitletext,ylabel='AP width at half amplitude [ms]')
plt.savefig(plotname_all)
#ax3.legend(loc='upper right')

plt.show()
    