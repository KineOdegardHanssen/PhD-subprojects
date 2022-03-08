import numpy as np
import matplotlib.pyplot as plt

# change the default font family
plt.rcParams.update({'font.family':'Arial'})

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

mylinewidth = 2

idur       = 1000 #100 # ms
idelay     = 100
iamp       = 0.3 # nA
v_init     = -86.5 # mV
Ra         = 150
somasize   = 10 # 15 # 
dendlen    = 1000
denddiam   = 1
nsegments  = 200 
spikedurat = -40

varymech = 'Na'#'Kd'#'NaV' #'SK'#'Im_v2'#'Kv2like'#'pas' #
varyE_bool = True
varyE = 50 #[-90,-80]#[-130,-107,-100,-70]#[40,50,53,65,70]#[-60,-53,-30,-10,10,30,60] # Default...407Kd-107 Default478513407pas: -83.6528;  Default478513407Na: 53
varyg = 'None' # Default...407Kd2.94396e-010 # Default478513407pas: 0.000362109; Default478513407Na: 0.0409177

varylist = [] # Should be redundant
plotstring = '_vary'
if varyE_bool==True:
    varylist = varyE
    plotstring = plotstring + 'E'+str(varyE)
else:
    varylist = varyg
    plotstring = plotstring + 'g'+str(varyg)

if varymech=='Na':
    folderstring = 'VaryNa/' 
    plotstring   = plotstring + '_Na'
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

Nspikes          = []
avg_AP_halfwidth = []
rms_AP_halfwidth = []
Cm_Nspikes       = []
Cm_AP_halfwidth  = []
    
# Set names
infolder = 'Ball-and-stick models/BAS_somaHH_dendpassive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+folderstring
currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
vrestfolder = infolder
infolder = infolder+currentfolder
infilename_Nspikes = infolder+'basHHdpas_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmfs_Nspikes_vs_Cmall.txt'
infilename_APdhw   = infolder+'basHHdpas_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmfs_sdurat%s_vs_Cmall.txt'% str(spikedurat) 
# Read files
infile_Nspikes = open(infilename_Nspikes,'r')
infile_APdhw   = open(infilename_APdhw,'r')

lines_Nspikes = infile_Nspikes.readlines()
lines_APdhw   = infile_APdhw.readlines()
Nlines_Nspikes = len(lines_Nspikes)
Nlines_APdhw   = len(lines_APdhw)

print('infilename_APdhw:',infilename_APdhw)
print('Nlines_Nspikes:',Nlines_Nspikes)
print('lines_Nspikes:',lines_Nspikes)

for i in range(Nlines_Nspikes):
    words_Nspikes = lines_Nspikes[i].split()
    if len(words_Nspikes)>0:
        Cm_Nspikes.append(float(words_Nspikes[0]))
        Nspikes.append(float(words_Nspikes[1]))
        print('Appending to Nspikes:',float(words_Nspikes[1]))

for i in range(Nlines_APdhw):
    words_APdhw   = lines_APdhw[i].split()
    if len(words_APdhw)>0:
        Cm_AP_halfwidth.append(float(words_APdhw[0]))
        avg_AP_halfwidth.append(float(words_APdhw[1]))
        rms_AP_halfwidth.append(float(words_APdhw[2]))

infile_Nspikes.close()
infile_APdhw.close()


print('Nspikes, all:',Nspikes)

# Plotting

plotfolder = 'Ball-and-stick models/BAS_somaHH_dendpassive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+currentfolder+'/Comparemodels'
plotname_all = plotfolder+'features_cmcombi_BASHHdpas_all_iamp'+str(iamp)+'.png'

## avg and rms:
#fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(17,15),dpi=300)
fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(15,15),dpi=300)

plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('legend', fontsize=16)

ax1.set_title(r'Frequency $f$ vs $C_m$',fontsize=16)
ax1.plot(Cm_Nspikes, Nspikes,label=r'$C_{m,\mathregular{all}}$', linewidth=mylinewidth)
ax1.set_ylabel('$f$ [Hz]',fontsize=14)

ax3.set_title(r'Spike duration at %s mV vs $C_m$' % str(spikedurat),fontsize=16)
ax3.errorbar(Cm_AP_halfwidth, avg_AP_halfwidth, yerr=rms_AP_halfwidth, capsize=2, linewidth=mylinewidth)
ax3.set_ylabel('Spike duration [ms]',fontsize=14)


############# SOMAPROX ##################################
Nspikes          = []
avg_AP_halfwidth = []
rms_AP_halfwidth = []
Cm_Nspikes       = []
Cm_AP_halfwidth  = []
    
# Set names
infolder = 'Ball-and-stick models/BAS_somaHH_dendpassive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+folderstring
currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
vrestfolder   = infolder
infolder      = infolder+currentfolder
infilename_Nspikes = infolder+'basHHdpas_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmfs_Nspikes_vs_Cmsprx.txt'
infilename_APdhw   = infolder+'basHHdpas_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmfs_sdurat%s_vs_Cmsprx.txt'% str(spikedurat) 
# Read files
infile_Nspikes = open(infilename_Nspikes,'r')
infile_APdhw   = open(infilename_APdhw,'r')

lines_Nspikes = infile_Nspikes.readlines()
lines_APdhw   = infile_APdhw.readlines()
Nlines_Nspikes = len(lines_Nspikes)
Nlines_APdhw   = len(lines_APdhw)

print('infilename_APdhw:',infilename_APdhw)
print('Nlines_Nspikes:',Nlines_Nspikes)
print('lines_Nspikes:',lines_Nspikes)

for i in range(Nlines_Nspikes):
    words_Nspikes = lines_Nspikes[i].split()
    if len(words_Nspikes)>0:
        Cm_Nspikes.append(float(words_Nspikes[0]))
        Nspikes.append(float(words_Nspikes[1]))
        print('Appending to Nspikes:',float(words_Nspikes[1]))


for i in range(Nlines_APdhw):
    words_APdhw   = lines_APdhw[i].split()
    if len(words_APdhw)>0:
        Cm_AP_halfwidth.append(float(words_APdhw[0]))
        avg_AP_halfwidth.append(float(words_APdhw[1]))
        rms_AP_halfwidth.append(float(words_APdhw[2]))

infile_Nspikes.close()
infile_APdhw.close()

print('infilename_APdhw:',infilename_APdhw)
print('Nspikes, somaprox:',Nspikes)

print('cm_master:',cm_master) 
print('Nspikes:',Nspikes)
ax1.plot(Cm_Nspikes, Nspikes,label=r'$C_{m,\mathregular{somaprox.}}$', linewidth=mylinewidth)
ax3.errorbar(Cm_AP_halfwidth, avg_AP_halfwidth, yerr=rms_AP_halfwidth, capsize=2, linewidth=mylinewidth)

print('Vrest, BAS sprx:',avg_Vrest)
#################### Soma only, Hodgkin-Huxley #########################################
Nspikes          = []
avg_AP_halfwidth = []
rms_AP_halfwidth = []
Cm_Nspikes       = []
Cm_AP_halfwidth  = []

# Default HH values:
ena = 50
ek = -77
el_hh = -54.3
gnabar_hh = 0.12
gkbar_hh = 0.036
gl_hh = 0.0003

infolder_shh = 'Somaonly/Results/IStim/Soma%i/' % somasize
vrestfolder  = infolder_shh 
currentfolder = 'current_idur%i_iamp'%idur+str(iamp)+'/'
infolder_shh  = infolder_shh+currentfolder
hhstring = '_ena'+str(ena)+'_ek'+str(ek)+'_el'+str(el_hh)+'_gnabar'+str(gnabar_hh)+'_gkbar'+str(gkbar_hh)+'_gl'+str(gl_hh)
infilename_Nspikes = infolder_shh+'somaonlyHH_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmfs_Nspikes_vs_Cm.txt'
infilename_APdhw   = infolder_shh+'somaonlyHH_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmfs_sdurat%s_vs_Cm.txt' % str(spikedurat)
infile_Nspikes = open(infilename_Nspikes,'r')
infile_APdhw   = open(infilename_APdhw,'r')

lines_Nspikes = infile_Nspikes.readlines()
lines_APdhw   = infile_APdhw.readlines()
Nlines = len(lines_Nspikes) # All data types have the same length

for i in range(Nlines):
    words_Nspikes = lines_Nspikes[i].split()
    words_APampl  = lines_APampl[i].split()
    words_APmins  = lines_APmins[i].split()
    words_APdhw   = lines_APdhw[i].split()
    words_ISI     = lines_ISI[i].split()
    if len(words_Nspikes)>0:
        Cm_Nspikes.append(float(words_Nspikes[0]))
        Nspikes.append(float(words_Nspikes[1]))
    if len(words_APdhw)>0:
        Cm_AP_halfwidth.append(float(words_APdhw[0]))
        avg_AP_halfwidth.append(float(words_APdhw[1]))
        rms_AP_halfwidth.append(float(words_APdhw[2]))

infile_Nspikes.close()
infile_APdhw.close()

print('Vrest, onecomp:',avg_Vrest)
################## ADDING THRESHOLDS AT THE END: #######################################
Ncm_somaHH        = 6
infolder_thr      = 'Ball-and-stick models/BAS_somaHH_dendpassive/Results/IStim/Soma%i/dendlen%i/denddiam%i/' % (somasize,dendlen,denddiam)+folderstring
infolder_thr_soma = 'Somaonly/Results/IStim/Soma%i/' % somasize
infilename_all    = infolder_thr + 'thresholds_bashhdpas_everywhere_varyfactor'+plotstring+'.txt' 
infilename_sprx   = infolder_thr + 'thresholds_bashhdpas_somaprox_varyfactor'+plotstring+'.txt' 
infilename_soma   = infolder_thr_soma + 'thresholds_somaonlyHH_varycmfactor_%icms_idur%i' % (Ncm_somaHH,idur)+hhstring+'.txt'

infile_all  = open(infilename_all,'r')
infile_sprx = open(infilename_sprx,'r')
infile_soma = open(infilename_soma,'r')
lines_all   = infile_all.readlines()
lines_sprx  = infile_sprx.readlines()
lines_soma  = infile_soma.readlines()

cm_all    = []
cm_sprx   = []
cm_soma   = []
thrs_all  = []
thrs_sprx = []
thrs_soma = []

for line in lines_all:
    words = line.split()
    if len(words)>0:
        cm_all.append(float(words[0]))
        thrs_all.append(float(words[1]))

for line in lines_sprx:
    words = line.split()
    if len(words)>0:
        cm_sprx.append(float(words[0]))
        thrs_sprx.append(float(words[1]))

for line in lines_soma:
    words = line.split()
    if len(words)>0:
        cm_soma.append(float(words[0]))
        thrs_soma.append(float(words[1]))

cm_all    = np.array(cm_all)
cm_sprx   = np.array(cm_sprx)
cm_soma   = np.array(cm_soma)
thrs_all  = np.array(thrs_all)
thrs_sprx = np.array(thrs_sprx)
thrs_soma = np.array(thrs_soma)


# Plotting

plotname_all = 'Comparemodels/BAS_vs_somaonly_HH/features_cmcombi_BASHHdpas_somaHH_all_iamp'+str(iamp)+'.png'

fig.suptitle('I = %s nA' % str(iamp),fontsize=18)

ax1.set_title(r'A',loc='left',fontsize=16)
ax2.set_title(r'B',loc='left',fontsize=16)
ax3.set_title(r'C',loc='left',fontsize=16)
ax4.set_title(r'D',loc='left',fontsize=16)
ax5.set_title(r'E',loc='left',fontsize=16)
ax6.set_title(r'F',loc='left',fontsize=16)

ax1.plot(Cm_Nspikes, Nspikes,label=r'One comp.', linewidth=mylinewidth)
ax2.errorbar(Cm_AP_ampl, avg_AP_ampl, yerr=rms_AP_ampl, capsize=2, color='#2ca02c', label=r'$V_{\mathregular{max}}$', linewidth=mylinewidth)
ax2.errorbar(Cm_AP_mins, avg_AP_mins, yerr=rms_AP_mins, capsize=2, color='#2ca02c', label=r'$V_{\mathregular{min}}$', linewidth=mylinewidth)
ax3.errorbar(Cm_AP_halfwidth, avg_AP_halfwidth, yerr=rms_AP_halfwidth, capsize=2, linewidth=mylinewidth)
ax5.errorbar(Cm_Vrest, avg_Vrest, yerr=rms_Vrest, capsize=2, linewidth=mylinewidth)
ax5.axis([-0.3, 3.3,-66,-64])
ax6.plot(cm_all,thrs_all,label=r'$C_{\mathregular{all}}$', linewidth=mylinewidth)
ax6.plot(cm_sprx,thrs_sprx,label=r'$C_{\mathregular{somaprox}}$', linewidth=mylinewidth)
ax6.plot(cm_soma,thrs_soma,label=r'One comp.', linewidth=mylinewidth)
ax1.legend(loc='upper right')

ax6.set_title(r'Threshold current vs $C_m$',fontsize=16)
ax6.set_xlabel(r'$C_{m}$ [$\mu$ F/cm$^2$]',fontsize=14)
ax6.set_ylabel(r'Threshold current [nA]',fontsize=14)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])

#### PLOTTING ALLEN #####

testmodels = [478513437,488462965,478513407] #485694403,489931686,
idur       = 2000 #100 # ms
idelay     = 100
iamp       = 0.2 # nA
v_init     = -86.5 # mV
Ra         = 150
somasize   = 10 # 15 # 
dendlen    = 1000
denddiam   = 1
nsegments  = 200 
Nmodels    = len(testmodels)
spikedurat = -40

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
avg_AP_halfwidth_all = []
rms_AP_halfwidth_all = []
Cm_Nspikes_all       = []
Cm_AP_halfwidth_all  = []
    
for testmodel in testmodels:
    Nspikes          = []
    Cm_Nspikes       = []
    Cm_AP_halfwidth  = []
    avg_AP_halfwidth = []
    rms_AP_halfwidth = []
    
    # Set names
    infolder      = 'figures/%i/' % (testmodel)
    vrestfolder   = infolder 
    currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    infolder = infolder+currentfolder
    infilename_Nspikes = infolder+'%i_idur%i_iamp'% (testmodel,idur)+str(iamp)+'_manual_cmfs_Nspikes_vs_Cmall.txt'
    infilename_APdhw   = infolder+'%i_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_manual_cmfs_sdurat%s_vs_Cmall.txt' % str(spikedurat)
    # Read files
    infile_Nspikes = open(infilename_Nspikes,'r')
    infile_APdhw   = open(infilename_APdhw,'r')
    
    lines_Nspikes = infile_Nspikes.readlines()
    lines_APdhw   = infile_APdhw.readlines()
    Nlines = len(lines_Nspikes) # All data types have the same length
    
    for i in range(Nlines):
        words_Nspikes = lines_Nspikes[i].split()
        words_APdhw   = lines_APdhw[i].split()
        if len(words_Nspikes)>0:
            Cm_Nspikes.append(float(words_Nspikes[0]))
            Nspikes.append(float(words_Nspikes[1]))
        if len(words_APdhw)>0:
            Cm_AP_halfwidth.append(float(words_APdhw[0]))
            avg_AP_halfwidth.append(float(words_APdhw[1]))
            rms_AP_halfwidth.append(float(words_APdhw[2]))
    
    infile_Nspikes.close()
    infile_APdhw.close()
    
    Nspikes_all.append(Nspikes)
    avg_AP_halfwidth_all.append(avg_AP_halfwidth)
    rms_AP_halfwidth_all.append(rms_AP_halfwidth)
    Cm_Nspikes_all.append(Cm_Nspikes)
    Cm_AP_halfwidth_all.append(Cm_AP_halfwidth)

## Should I take avg. and rms. of all models? Might be difficult since I allow for different Cms.
# Check against cm_master
Nspikes     = np.zeros(NCms)
Nspikes_rms = np.zeros(NCms)
avg_AP_halfwidth = np.zeros(NCms)
rms_AP_halfwidth = np.zeros(NCms)

Nspikes_bucket   = np.zeros((NCms,Nmodels))

Nspikes_collections = []
avg_AP_halfwidth_collections = []
rms_AP_halfwidth_collections = []

fcounter    = np.zeros(NCms)
isicounter  = np.zeros(NCms)
amplcounter = np.zeros(NCms)
minscounter = np.zeros(NCms)
hwcounter   = np.zeros(NCms)
vrcounter   = np.zeros(NCms)

for i in range(Nmodels):
    Nspikes_this = Nspikes_all[i]
    avg_AP_halfwidth_this = avg_AP_halfwidth_all[i]
    rms_AP_halfwidth_this = rms_AP_halfwidth_all[i]
    Cm_Nspikes_this       = Cm_Nspikes_all[i]
    Cm_AP_halfwidth_this  = Cm_AP_halfwidth_all[i]
    Nspikes_collections = []
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
    if hwcounter[j]!=0:
        avg_AP_halfwidth[j] /= hwcounter[j]
        rms_AP_halfwidth[j]  = np.sqrt(rms_AP_halfwidth[j]/(hwcounter[j]-1))

# Plotting

plotfolder = 'figures/Comparemodels/'+currentfolder
plotname_avgrms = plotfolder+'features_combi_avgrms.png'

#Nspikes_bucket   = np.zeros((NCms,Nmodels))
cm_short = []
Nspikes_short = []
for k in range(Nmodels):
    cm_short_this = []
    Nspikes_short_this = []
    for i in range(NCms): # Is this gonna work?
        if Nspikes_bucket[i,k]!=0:
            cm_short_this.append(cm_master[i])
            Nspikes_short_this.append(Nspikes_bucket[i,k])
        else:
            break
    cm_short.append(cm_short_this)
    Nspikes_short.append(Nspikes_short_this)
    

## avg and rms:
#fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(17,15))
#fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(15,15))

ax1.set_title(r'A',loc='left',fontsize=16)
ax2.set_title(r'B',loc='left',fontsize=16)
ax3.set_title(r'C',loc='left',fontsize=16)
ax4.set_title(r'D',loc='left',fontsize=16)
ax5.set_title(r'E',loc='left',fontsize=16)
ax6.set_title(r'F',loc='left',fontsize=16)

ax1.set_title(r'Frequency $f$ vs $C_m$',fontsize=16)
ax1.errorbar(cm_master, Nspikes,yerr=Nspikes_rms,label=r'$C_{m,\mathregular{all}}$',capsize=2, linewidth=2)
ax1.set_ylabel('$f$ [Hz]',fontsize=14)

ax3.set_title(r'Spike duration at %s mV vs $C_m$' % str(spikedurat),fontsize=16)
ax3.errorbar(cm_master, avg_AP_halfwidth, yerr=rms_AP_halfwidth, capsize=2, linewidth=2)
ax3.set_ylabel('Spike duration at %s mV [ms]' % str(spikedurat),fontsize=14)


## Saving arrays for next figure:
all_cm_short             = cm_short
all_Nspikes_short        = Nspikes_short
all_Cm_AP_halfwidth_all  = Cm_AP_halfwidth_all
all_avg_AP_halfwidth_all = avg_AP_halfwidth_all
all_rms_AP_halfwidth_all = rms_AP_halfwidth_all

############################################### SOMAPROX #######################################################
Nspikes_all          = []
avg_AP_halfwidth_all = []
rms_AP_halfwidth_all = []
Cm_Nspikes_all       = []
Cm_AP_halfwidth_all  = []
    
for testmodel in testmodels:
    Nspikes          = []
    Cm_Nspikes       = []
    Cm_AP_ampl       = []
    Cm_AP_mins       = []
    Cm_AP_halfwidth  = []
    avg_AP_halfwidth = []
    rms_AP_halfwidth = []
    
    # Set names
    infolder      = 'figures/%i/' % (testmodel)
    vrestfolder   = infolder
    currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    infolder = infolder+currentfolder
    infilename_Nspikes = infolder+'%i_idur%i_iamp'% (testmodel,idur)+str(iamp)+'_manual_cmfs_Nspikes_vs_Cmsprx.txt'
    infilename_APdhw   = infolder+'%i_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_manual_cmfs_sdurat%s_vs_Cmsprx.txt' % str(spikedurat)
    # Read files
    infile_Nspikes = open(infilename_Nspikes,'r')
    infile_APdhw   = open(infilename_APdhw,'r')
    
    lines_Nspikes = infile_Nspikes.readlines()
    lines_APdhw   = infile_APdhw.readlines()
    Nlines = len(lines_Nspikes) # All data types have the same length
    
    for i in range(Nlines):
        words_Nspikes = lines_Nspikes[i].split()
        words_APdhw   = lines_APdhw[i].split()
        if len(words_Nspikes)>0:
            Cm_Nspikes.append(float(words_Nspikes[0]))
            Nspikes.append(float(words_Nspikes[1]))
        if len(words_APdhw)>0:
            Cm_AP_halfwidth.append(float(words_APdhw[0]))
            avg_AP_halfwidth.append(float(words_APdhw[1]))
            rms_AP_halfwidth.append(float(words_APdhw[2]))
    
    infile_Nspikes.close()
    infile_APdhw.close()
    
    Nspikes_all.append(Nspikes)
    avg_AP_halfwidth_all.append(avg_AP_halfwidth)
    rms_AP_halfwidth_all.append(rms_AP_halfwidth)
    Cm_Nspikes_all.append(Cm_Nspikes)
    Cm_AP_halfwidth_all.append(Cm_AP_halfwidth)

## Should I take avg. and rms. of all models? Might be difficult since I allow for different Cms.
# Check against cm_master
Nspikes     = np.zeros(NCms)
Nspikes_rms = np.zeros(NCms)
avg_AP_halfwidth = np.zeros(NCms)
rms_AP_halfwidth = np.zeros(NCms)

Nspikes_bucket   = np.zeros((NCms,Nmodels))

Nspikes_collections = []
avg_AP_halfwidth_collections = []
rms_AP_halfwidth_collections = []

fcounter    = np.zeros(NCms)
isicounter  = np.zeros(NCms)
amplcounter = np.zeros(NCms)
minscounter = np.zeros(NCms)
hwcounter   = np.zeros(NCms)
vrcounter   = np.zeros(NCms)

for i in range(Nmodels):
    Nspikes_this = Nspikes_all[i]
    avg_AP_halfwidth_this = avg_AP_halfwidth_all[i]
    rms_AP_halfwidth_this = rms_AP_halfwidth_all[i]
    Cm_Nspikes_this       = Cm_Nspikes_all[i]
    Cm_AP_halfwidth_this  = Cm_AP_halfwidth_all[i]
    Nspikes_collections = []
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
    if hwcounter[j]!=0:
        avg_AP_halfwidth[j] /= hwcounter[j]
        rms_AP_halfwidth[j]  = np.sqrt(rms_AP_halfwidth[j]/(hwcounter[j]-1))

# Plotting

plotfolder = 'figures/Comparemodels/'+currentfolder
#plotname_all = plotfolder+'features_vs_cm_Allen_all_iamp'+str(iamp)+'.png'
#plotname_all_II = plotfolder+'features_vs_cm_Allen_all_iamp'+str(iamp)+'_noerrorbars.png'
#plotname_avgrms = plotfolder+'features_vs_cm_Allen_all_iamp'+str(iamp)+'_avgrms.png'

#Nspikes_bucket   = np.zeros((NCms,Nmodels))
cm_short = []
Nspikes_short = []
for k in range(Nmodels):
    cm_short_this = []
    Nspikes_short_this = []
    for i in range(NCms): # Is this gonna work?
        if Nspikes_bucket[i,k]!=0:
            cm_short_this.append(cm_master[i])
            Nspikes_short_this.append(Nspikes_bucket[i,k])
        else:
            break
    cm_short.append(cm_short_this)
    Nspikes_short.append(Nspikes_short_this)

## avg and rms:
ax1.errorbar(cm_master, Nspikes,yerr=Nspikes_rms,label=r'$C_{m, \mathregular{somaprox.}}$',ls='--',capsize=2)
ax1.legend(loc='upper right')
ax3.errorbar(cm_master, avg_AP_halfwidth, yerr=rms_AP_halfwidth,ls='--', capsize=2, linewidth=2)

fig.tight_layout(rect=[0, 0.03, 1, 0.95])
#plt.savefig(plotname_avgrms) ## ?? Which one?

##### Last plot: #####################################################################################

# This will have threshold data in it too # Doesn't make much sense to average the thresholds, I guess... #
# Or maybe averaging the thresholds do make sense... But then I'd... No.
#### Thresholds: ####

all_cm_all    = []
all_cm_sprx   = []
all_thrs_all  = []
all_thrs_sprx = []
for model in testmodels:
    infilename_all  = 'figures/threshold_Allen_everywhere_models%s.txt' %str(model) 
    infilename_sprx = 'figures/threshold_Allen_somaprox_models%s.txt' %str(model) 
    
    infile_all  = open(infilename_all,'r')
    infile_sprx = open(infilename_sprx,'r')
    lines_all   = infile_all.readlines()
    lines_sprx  = infile_sprx.readlines()
    
    cm_all    = []
    cm_sprx   = []
    thrs_all  = []
    thrs_sprx = []
    
    for line in lines_all:
        words = line.split()
        if len(words)>0:
            cm_all.append(float(words[0]))
            thrs_all.append(float(words[1]))

    for line in lines_sprx:
        words = line.split()
        if len(words)>0:
            cm_sprx.append(float(words[0]))
            thrs_sprx.append(float(words[1]))
        
    cm_all    = np.array(cm_all)
    cm_sprx   = np.array(cm_sprx)
    thrs_all  = np.array(thrs_all)
    thrs_sprx = np.array(thrs_sprx)
    
    # Appending for plotting
    all_cm_all.append(cm_all)
    all_cm_sprx.append(cm_sprx)
    all_thrs_all.append(thrs_all)
    all_thrs_sprx.append(thrs_sprx)

ax1.set_title(r'A',loc='left',fontsize=18)
ax2.set_title(r'B',loc='left',fontsize=18)
ax3.set_title(r'C',loc='left',fontsize=18)
ax4.set_title(r'D',loc='left',fontsize=18)
ax5.set_title(r'E',loc='left',fontsize=18)
ax6.set_title(r'F',loc='left',fontsize=18)

ax1.set_title(r'Frequency $f$ vs $C_m$',fontsize=16)
for i in range(Nmodels):
    colors = colorarray[i]
    ax1.plot(all_cm_short[i], all_Nspikes_short[i],label='%i: $C_{m,\mathregular{all}}$' % testmodels[i],color=colors[0], linewidth=2)
    ax1.plot(cm_short[i], Nspikes_short[i],ls='--',label='%i: $C_{m,\mathregular{somaprox.}}$' % testmodels[i],color=colors[1], linewidth=2)
ax1.set_ylabel('$f$ [Hz]',fontsize=14)
ax1.legend(loc='upper right')

ax3.set_title(r'Spike duration at %s mV vs $C_m$' % str(spikedurat),fontsize=16)
for i in range(Nmodels):
    colors = colorarray[i]
    ax3.errorbar(all_Cm_AP_halfwidth_all[i], all_avg_AP_halfwidth_all[i], yerr=all_rms_AP_halfwidth_all[i], capsize=2,label='%i,$C_{m,\mathregular{all}}$' % testmodels[i], color=colors[0], linewidth=2)
    ax3.errorbar(Cm_AP_halfwidth_all[i], avg_AP_halfwidth_all[i], yerr=rms_AP_halfwidth_all[i], capsize=2,ls='--',label='%i,$C_{m,\mathregular{somaprox.}}$' % testmodels[i], color=colors[1], linewidth=2)
#ax3.set(xlabel=r'$C_{m}$ [$\mu$ F/cm$^2$]',ylabel='Spike duration at %s mV [ms]' % str(spikedurat))
ax3.set_ylabel('Spike duration [ms]',fontsize=14)

ax6.set_title(r'Threshold current vs $C_m$',fontsize=16)
for i in range(Nmodels):
    colors = colorarray[i]
    ax6.plot(all_cm_all[i], all_thrs_all[i], label='%i,$C_{m,all}$' % testmodels[i],color=colors[0], linewidth=2)
    ax6.plot(all_cm_sprx[i], all_thrs_sprx[i], label='%i,$C_{m,somaprox.}$' % testmodels[i], color=colors[1],ls='--', linewidth=2)
ax6.set_xlabel(r'$C_{m}$ [$\mu$ F/cm$^2$]',fontsize=14)
ax6.set_ylabel('Threshold current [nA]',fontsize=14)

fig.tight_layout(rect=[0, 0.03, 1, 0.95])
#plt.savefig(plotname_all) # ???


#### NO ERROR BARS ###
ax1.set_title(r'A',loc='left',fontsize=16)
ax2.set_title(r'B',loc='left',fontsize=16)
ax3.set_title(r'C',loc='left',fontsize=16)
ax4.set_title(r'D',loc='left',fontsize=16)
ax5.set_title(r'E',loc='left',fontsize=16)
ax6.set_title(r'F',loc='left',fontsize=16)
fig.suptitle('I = %s nA' % str(iamp),fontsize=18)

ax1.set_title(r'Frequency $f$ vs $C_m$',fontsize=16)
for i in range(Nmodels):
    colors = colorarray[i]
    ax1.plot(all_cm_short[i], all_Nspikes_short[i],label='%i: $C_{m,\mathregular{all}}$' % testmodels[i],color=colors[0], linewidth=2)
    ax1.plot(cm_short[i], Nspikes_short[i],ls='--',label='%i: $C_{m,\mathregular{somaprox.}}$' % testmodels[i],color=colors[1], linewidth=2)
ax1.set_ylabel('$f$ [Hz]',fontsize=14)
ax1.legend(loc='upper right')

ax3.set_title(r'Spike duration at %s mV vs $C_m$' % str(spikedurat),fontsize=16)
for i in range(Nmodels):
    colors = colorarray[i]
    ax3.plot(all_Cm_AP_halfwidth_all[i], all_avg_AP_halfwidth_all[i],label='%i,$C_{m,all}$' % testmodels[i], color=colors[0], linewidth=2)
    ax3.plot(Cm_AP_halfwidth_all[i], avg_AP_halfwidth_all[i],ls='--',label='%i,$C_{m,somaprox.}$' % testmodels[i], color=colors[1], linewidth=2)
#ax3.set(xlabel=r'$C_{m}$ [$\mu$ F/cm$^2$]',ylabel='Spike duration at %s mV [ms]' % str(spikedurat))
ax3.set_ylabel('Spike duration [ms]',fontsize=14)

ax6.set_title(r'Threshold current vs $C_m$',fontsize=16)
for i in range(Nmodels):
    colors = colorarray[i]
    ax6.plot(all_cm_all[i], all_thrs_all[i], label='%i,$C_{m,all}$' % testmodels[i],color=colors[0], linewidth=2)
    ax6.plot(all_cm_sprx[i], all_thrs_sprx[i], ls='--', label='%i,$C_{m,somaprox.}$' % testmodels[i], color=colors[1], linewidth=2)
ax6.set_xlabel(r'$C_{m}$ [$\mu$ F/cm$^2$]',fontsize=14)
ax6.set_ylabel('Threshold current [nA]',fontsize=14)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.savefig(plotname_all)

print('plotname_all:',plotname_all)
plt.show()
    