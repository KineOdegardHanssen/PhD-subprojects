import numpy as np
import matplotlib.pyplot as plt

# change the default font family
plt.rcParams.update({'font.family':'Arial'})

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

def reldiffs(fother,f1,iother,i1):
    iout   = []
    iout2  = []
    rdout  = []
    rdout2 = []
    N1     = len(f1)
    Nother = len(fother)
    for i in range(Nother):
        ithis = iother[i]
        for j in range(N1):
            ibasic = i1[j]
            if ithis==ibasic:
                falt   = fother[i]
                fbasic = f1[j]
                if fbasic!=0 and falt!=0:
                    iout.append(ithis)
                    rdout.append((fbasic-falt)/float(fbasic))
                    rdout2.append((fbasic-falt)/float(falt))
                break
    iout   = np.array(iout)
    rdout  = np.array(rdout)
    return iout, rdout, rdout2

mylinewidth = 2

idur       = 1000 # ms
idelay     = 100
v_init     = -86.5 # mV
Ra         = 150
somasize   = 10 # 15 # 
dendlen    = 1000
denddiam   = 1
nsegments  = 200 
spikedurat = -40

varymech = 'Na' # Not what we want to test here
varyE_bool = True
varyE = 50 
varyg = 'None' 

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

i_master = [0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5]
cms = [0.1,0.5,0.75,1.0,1.25,1.5,1.99]
    
NI  = len(i_master)
Ncm = len(cms) # Need this later

i_master_everywhere         = []
Nspikes_everywhere          = []
avg_AP_ampl_everywhere      = []
rms_AP_ampl_everywhere      = []
avg_AP_mins_everywhere      = []
rms_AP_mins_everywhere      = []
avg_AP_halfwidth_everywhere = []
rms_AP_halfwidth_everywhere = []
avg_ISI_everywhere          = []
rms_ISI_everywhere          = []
I_Nspikes_everywhere        = []
I_AP_ampl_everywhere        = []
I_AP_mins_everywhere        = []
I_AP_halfwidth_everywhere   = []
I_ISI_everywhere            = []

i_master_sprx         = []
Nspikes_sprx          = []
avg_AP_ampl_sprx      = []
rms_AP_ampl_sprx      = []
avg_AP_mins_sprx      = []
rms_AP_mins_sprx      = []
avg_AP_halfwidth_sprx = []
rms_AP_halfwidth_sprx = []
avg_ISI_sprx          = []
rms_ISI_sprx          = []
I_Nspikes_sprx        = []
I_AP_ampl_sprx        = []
I_AP_mins_sprx        = []
I_AP_halfwidth_sprx   = []
I_ISI_sprx            = []

i_master_onecomp         = []
Nspikes_onecomp          = []
avg_AP_ampl_onecomp      = []
rms_AP_ampl_onecomp      = []
avg_AP_mins_onecomp      = []
rms_AP_mins_onecomp      = []
avg_AP_halfwidth_onecomp = []
rms_AP_halfwidth_onecomp = []
avg_ISI_onecomp          = []
rms_ISI_onecomp          = []
I_Nspikes_onecomp        = []
I_AP_ampl_onecomp        = []
I_AP_mins_onecomp        = []
I_AP_halfwidth_onecomp   = []
I_ISI_onecomp            = []


for cm in cms:
    Nspikes   = []
    I_Nspikes = []
    
    # Set names
    infolder = 'Ball-and-stick models/BAS_somaHH_dendpassive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+folderstring
    infilename_Nspikes = infolder+'basHHdpas_cmall'+str(cm)+'_idur%i_varyiamp'% (idur) +'_manual_Nspikes_vs_I.txt'
    # Read files
    infile_Nspikes = open(infilename_Nspikes,'r')

    lines_Nspikes = infile_Nspikes.readlines()
    Nlines_Nspikes = len(lines_Nspikes)

    for i in range(Nlines_Nspikes):
        words_Nspikes = lines_Nspikes[i].split()
        if len(words_Nspikes)>0:
            I_Nspikes.append(float(words_Nspikes[0]))
            Nspikes.append(float(words_Nspikes[1]))

    infile_Nspikes.close()
    
    Nspikes_everywhere.append(Nspikes)
    I_Nspikes_everywhere.append(I_Nspikes)

    ############# SOMAPROX ##################################
    Nspikes   = []
    I_Nspikes = []
        
    # Set names
    infolder = 'Ball-and-stick models/BAS_somaHH_dendpassive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+folderstring
    infilename_Nspikes = infolder+'basHHdpas_csprx'+str(cm)+'_idur%i_varyiamp'% (idur) +'_manual_Nspikes_vs_I.txt'
    # Read files
    infile_Nspikes = open(infilename_Nspikes,'r')
    lines_Nspikes  = infile_Nspikes.readlines()
    Nlines_Nspikes = len(lines_Nspikes)

    for i in range(Nlines_Nspikes):
        words_Nspikes = lines_Nspikes[i].split()
        if len(words_Nspikes)>0:
            I_Nspikes.append(float(words_Nspikes[0]))
            Nspikes.append(float(words_Nspikes[1]))

    infile_Nspikes.close()
    
    Nspikes_sprx.append(Nspikes)
    I_Nspikes_sprx.append(I_Nspikes)

    #print('i_master:',i_master) 
    #print('Nspikes:',Nspikes)

    #################### Soma only, Hodgkin-Huxley #########################################
    Nspikes   = []
    I_Nspikes = []

    # Default HH values:
    ena = 50
    ek = -77
    el_hh = -54.3
    gnabar_hh = 0.12
    gkbar_hh = 0.036
    gl_hh = 0.0003

    infolder_shh = 'Somaonly/Results/IStim/Soma%i/' % somasize
    hhstring = '_ena'+str(ena)+'_ek'+str(ek)+'_el'+str(el_hh)+'_gnabar'+str(gnabar_hh)+'_gkbar'+str(gkbar_hh)+'_gl'+str(gl_hh)
    infilename_Nspikes = infolder_shh+'somaonlyHH_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I.txt'
    # Read files
    infile_Nspikes = open(infilename_Nspikes,'r')
    lines_Nspikes  = infile_Nspikes.readlines()
    Nlines_Nspikes = len(lines_Nspikes)

    for i in range(Nlines_Nspikes):
        words_Nspikes = lines_Nspikes[i].split()
        if len(words_Nspikes)>0:
            I_Nspikes.append(float(words_Nspikes[0]))
            Nspikes.append(float(words_Nspikes[1]))

    infile_Nspikes.close()
    
    Nspikes_onecomp.append(Nspikes)
    I_Nspikes_onecomp.append(I_Nspikes)


# Plotting
# These are out of use for the time being... Need to add more colors.
color_bashhall  = ['#1f77b4','#d62728']
color_bashhsprx = ['#ff7f0e','#9467bd'] # Skip this?
color_somahh    = ['#2ca02c','#8c564b']

plotfolder = 'Comparemodels/All/'
plotname_all = plotfolder+'fI_allmodels.png'
plotname_rdf = plotfolder+'fI_rdiff_allmodels.png'
plotname_rdf_maxdiff  = plotfolder+'fI_rdiff_maxdiff_allmodels.png'
plotname_rdf_lastdiff = plotfolder+'fI_rdiff_lastdiff_allmodels.png'


fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(15,15),dpi=300)
#fig, ((ax1, ax2), (ax3, ax4), (ax5)) = plt.subplots(3, 2, figsize=(15,15),dpi=300)

ax1.set_title(r'A',loc='left',fontsize=18)
ax2.set_title(r'B',loc='left',fontsize=18)
ax3.set_title(r'C',loc='left',fontsize=18)
ax4.set_title(r'D',loc='left',fontsize=18)
ax5.set_title(r'E',loc='left',fontsize=18)


plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('legend', fontsize=16)

ax1.set_title('One compartment HH neuron',fontsize=18)
for i in range(Ncm):
    ax1.plot(I_Nspikes_onecomp[i], Nspikes_onecomp[i],label=r'One comp., $C_m$ = %.2f' % cms[i], linewidth=mylinewidth)
ax1.set_ylabel('$f$ [Hz]',fontsize=16)

ax2.set_title('Ball-and-stick HH neuron',fontsize=18)
for i in range(Ncm):
    ax2.plot(I_Nspikes_everywhere[i], Nspikes_everywhere[i],label=r'$C_{m,all}$ = %.2f' % cms[i], linewidth=mylinewidth)
    #ax2.plot(I_Nspikes_sprx[i], Nspikes_sprx[i],label=r'$C_{m,somaprox.}$ = %.2f' % cms[i], linewidth=mylinewidth)

## One panel for each Allen model too.

# Allen stuff:

testmodels = [478513437,488462965,478513407]
idur       = 2000 #100 # ms
idelay     = 100
v_init     = -86.5 # mV
Ra         = 150
somasize   = 10 # 15 # 
dendlen    = 1000
denddiam   = 1
nsegments  = 200 
Nmodels    = len(testmodels)
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

i_master = [0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5]
    
NI  = len(i_master)

i_master_everywhere_all   = []
Nspikes_everywhere_all    = []
I_Nspikes_everywhere_all  = []

i_master_sprx_all         = []
Nspikes_sprx_all          = []
I_Nspikes_sprx_all        = []

for testmodel in testmodels:
    i_master_everywhere         = []
    Nspikes_everywhere_allen    = []
    I_Nspikes_everywhere_allen  = []
    
    i_master_sprx               = []
    Nspikes_sprx_allen          = []
    I_Nspikes_sprx_allen        = []

    infolder      = 'Allen_test_changecapacitance/figures/%i/' % (testmodel)
    vrestfolder   = infolder 
    for cm in cms:
        Nspikes   = []
        I_Nspikes = []
        
        # Set names
        infilename_Nspikes = infolder+'%i_cmfall'%testmodel+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
        # Read files
        infile_Nspikes = open(infilename_Nspikes,'r')
        lines_Nspikes = infile_Nspikes.readlines()
        Nlines_Nspikes = len(lines_Nspikes)
    
        for i in range(Nlines_Nspikes):
            words_Nspikes = lines_Nspikes[i].split()
            if len(words_Nspikes)>0:
                I_Nspikes.append(float(words_Nspikes[0]))
                Nspikes.append(float(words_Nspikes[1]))
    
        infile_Nspikes.close()

        Nspikes_everywhere_allen.append(Nspikes)
        I_Nspikes_everywhere_allen.append(I_Nspikes)
    
        ############# SOMAPROX ##################################
        Nspikes         = []
        I_Nspikes       = []
        
        # Set names
        infilename_Nspikes = infolder+'%i_cmsprx'%testmodel+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
        # Read files
        infile_Nspikes = open(infilename_Nspikes,'r')
    
        lines_Nspikes = infile_Nspikes.readlines()
        Nlines_Nspikes = len(lines_Nspikes)
    
        for i in range(Nlines_Nspikes):
            words_Nspikes = lines_Nspikes[i].split()
            if len(words_Nspikes)>0:
                I_Nspikes.append(float(words_Nspikes[0]))
                Nspikes.append(float(words_Nspikes[1]))
    
        infile_Nspikes.close()
        
        Nspikes_sprx_allen.append(Nspikes)
        I_Nspikes_sprx_allen.append(I_Nspikes)
    Nspikes_everywhere_all.append(Nspikes_everywhere_allen)
    I_Nspikes_everywhere_all.append(I_Nspikes_everywhere_allen)
    
    Nspikes_sprx_all.append(Nspikes_sprx_allen)
    I_Nspikes_sprx_all.append(I_Nspikes_sprx_allen)

# Plotting
    
#color_pv_all   = ['#1f77b4','#d62728','#2ca02c']
#color_pv_sprx  = ['#ff7f0e','#9467bd','#8c564b']

color_pv_all   = ['#1f77b4','#d62728','#2ca02c']
#color_pv_sprx  = ['darkblue','#9467bd','darkgreen']
color_pv_sprx  = ['rebeccapurple','#9467bd','olive']

## avg and rms:

ax3.set_title('Allen model %s' % testmodels[0],fontsize=18)
ax4.set_title('Allen model %s' % testmodels[1],fontsize=18)
ax5.set_title('Allen model %s' % testmodels[2],fontsize=18)
I_Nspikes_everywhere_this0 = I_Nspikes_everywhere_all[0]
I_Nspikes_everywhere_this1 = I_Nspikes_everywhere_all[1]
I_Nspikes_everywhere_this2 = I_Nspikes_everywhere_all[2]
I_Nspikes_sprx_this0       = I_Nspikes_sprx_all[0]
I_Nspikes_sprx_this1       = I_Nspikes_sprx_all[1]
I_Nspikes_sprx_this2       = I_Nspikes_sprx_all[2]
Nspikes_everywhere_this0  = Nspikes_everywhere_all[0]
Nspikes_everywhere_this1  = Nspikes_everywhere_all[1]
Nspikes_everywhere_this2  = Nspikes_everywhere_all[2]
Nspikes_sprx_this0        = Nspikes_sprx_all[0]
Nspikes_sprx_this1        = Nspikes_sprx_all[1]
Nspikes_sprx_this2        = Nspikes_sprx_all[2]
for i in range(Ncm):
    ax3.plot(I_Nspikes_everywhere_this0[i], Nspikes_everywhere_this0[i],
label=r'$C_{m,\mathregular{all}}$=%.2f' % cms[i], linewidth=mylinewidth)
    ax4.plot(I_Nspikes_everywhere_this1[i], Nspikes_everywhere_this1[i],
label=r'$C_{m,\mathregular{all}}$=%.2f' % cms[i], linewidth=mylinewidth)
    ax5.plot(I_Nspikes_everywhere_this2[i], Nspikes_everywhere_this2[i],label=r'$C_{m,\mathregular{all}}$=%.2f' % cms[i], linewidth=mylinewidth)
    ### Somaprox:
    #ax3.plot(I_Nspikes_sprx_this0[i], Nspikes_sprx_this0[i],'--',label=r'$C_{m,\mathregular{somaprox.}}$=%.2f' % cms[i], linewidth=mylinewidth)
    #ax4.plot(I_Nspikes_sprx_this1[i], Nspikes_sprx_this1[i],'--',label=r'$C_{m,\mathregular{somaprox.}}$=%.2f' % cms[i], linewidth=mylinewidth)
    #ax5.plot(I_Nspikes_sprx_this2[i], Nspikes_sprx_this2[i],'--',label=r'$C_{m,\mathregular{somaprox.}}$=%.2f' % cms[i], linewidth=mylinewidth)

ax1.legend(loc='lower right')
ax2.legend(loc='lower right')
ax3.legend(loc='lower right')
ax4.legend(loc='upper left')
ax5.legend(loc='lower right')
ax3.set_xlabel(r'$I$ (nA)',fontsize=16)
ax3.set_ylabel(r'$f$ (Hz)',fontsize=16)
ax4.set_xlabel(r'$I$ (nA)',fontsize=16)
ax4.set_ylabel(r'$f$ (Hz)',fontsize=16)
ax5.set_xlabel(r'$I$ (nA)',fontsize=16)
ax5.set_ylabel(r'$f$ (Hz)',fontsize=16)


fig.tight_layout()
plt.savefig(plotname_all)

############ Differences #############
# Store the max values. Plot the relative difference
fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(15,15),dpi=300)
currat0      = []
currat1      = []
currat2      = [] 
maxdiff0     = []
maxdiff0_v1  = []
maxdiff0_v2  = []
maxdiff1     = []
maxdiff1_v1  = []
maxdiff1_v2  = []
maxdiff2     = [] 
maxdiff2_v1  = [] 
maxdiff2_v2  = [] 
maxdiff_tot  = []
lastdiff1_m0 = []
lastdiff1_m1 = []
lastdiff1_m2 = []
for i in range(Ncm): # Should I do this for somaprox too?
    if i!=3:
        #print('---------------------------')
        im0, rd1_m0, rd2_m0 = reldiffs(Nspikes_everywhere_this0[i],Nspikes_everywhere_this0[3],I_Nspikes_everywhere_this0[i],I_Nspikes_everywhere_this0[3]) 
        im1, rd1_m1, rd2_m1 = reldiffs(Nspikes_everywhere_this1[i],Nspikes_everywhere_this1[3],I_Nspikes_everywhere_this1[i],I_Nspikes_everywhere_this1[3]) 
        im2, rd1_m2, rd2_m2 = reldiffs(Nspikes_everywhere_this2[i],Nspikes_everywhere_this2[3],I_Nspikes_everywhere_this2[i],I_Nspikes_everywhere_this2[3])
        if len(rd1_m0)>0:
            maxrd1_m0 = max(rd1_m0) # 1
            minrd1_m0 = min(rd1_m0)
            if maxrd1_m0<abs(minrd1_m0): # 1
                maxrd1_m0 = minrd1_m0
            maxdiff0_v1.append(maxrd1_m0)
            currat0.append(im0[np.argmax(abs(rd1_m0))])
        if len(rd1_m1)>0:
            maxrd1_m1 = max(rd1_m1) # 2
            minrd1_m1 = min(rd1_m1)
            if maxrd1_m1<abs(minrd1_m1): # 2
                maxrd1_m1 = minrd1_m1
            maxdiff1_v1.append(maxrd1_m1)
            currat1.append(im1[np.argmax(abs(rd1_m1))])
        if len(rd1_m2)>0:
            maxrd1_m2 = max(rd1_m2) # 3
            minrd1_m2 = min(rd1_m2)
            if maxrd1_m2<abs(minrd1_m2): # 3
                maxrd1_m2 = minrd1_m2
            maxdiff2_v1.append(maxrd1_m2)
            currat2.append(im2[np.argmax(abs(rd1_m2))])
        if len(rd2_m0)>0:
            maxrd2_m0 = max(rd2_m0) # 4
            minrd2_m0 = min(rd2_m0)
            if maxrd2_m0<abs(minrd2_m0): # 4
                maxrd2_m0 = minrd2_m0
            maxdiff0_v2.append(maxrd2_m0)
        if len(rd2_m1)>0:
            maxrd2_m1 = max(rd2_m1) # 5
            minrd2_m1 = min(rd2_m1)
            if maxrd2_m1<abs(minrd2_m1): # 5
                maxrd2_m1 = minrd2_m1
            maxdiff1_v2.append(maxrd2_m1)
        if len(rd2_m2)>0:
            maxrd2_m2 = max(rd2_m2) # 6
            minrd2_m2 = min(rd2_m2)
            if maxrd2_m2<abs(minrd2_m2): # 6
                maxrd2_m2 = minrd2_m2
            maxdiff2_v2.append(maxrd2_m2)
        maxes = []
        if len(rd1_m0)>0 and len(rd2_m0)>0:
            maxtot_m0 = maxrd1_m0
            if abs(maxrd1_m0)<abs(maxrd2_m0):
                maxtot_m0 = maxrd2_m0
            maxdiff0.append(maxtot_m0)
            maxes.append(maxtot_m0)
        if len(rd1_m1)>0 and len(rd2_m1)>0:
            maxtot_m1 = maxrd1_m1
            if abs(maxrd1_m1)<abs(maxrd2_m1):
                maxtot_m1 = maxrd2_m1
            maxdiff1.append(maxtot_m1)
            maxes.append(maxtot_m1)
        if len(rd1_m2)>0 and len(rd2_m2)>0:
            maxtot_m2 = maxrd1_m2
            if abs(maxrd1_m2)<abs(maxrd2_m2):
                maxtot_m2 = maxrd2_m2
            maxdiff2.append(maxtot_m2)
            maxes.append(maxtot_m2)
        maxes = np.array(maxes)
        maxtot = max(maxes)
        #print('rd1_m0:',rd1_m0)
        #print('np.argmax(abs(rd1_m0)):',np.argmax(abs(rd1_m0)))
        #print('im0:',im0)
        #print('im0[np.argmax(abs(rd1_m0)):',im0[np.argmax(abs(rd1_m0))])
        maxdiff_tot.append(maxtot)
        #print('---------------------------')
        # Plotting:
        if len(rd1_m0)>0:
            ax3.plot(im0, rd1_m0, label='$C_m$=%s' % str(cms[i]))
            lastdiff1_m0.append(rd1_m0[-1])
        if len(rd1_m1)>0:
            ax4.plot(im1, rd1_m1, label='$C_m$=%s' % str(cms[i]))
            lastdiff1_m1.append(rd1_m1[-1])
        if len(rd1_m2)>0:
            ax5.plot(im2, rd1_m2, label='$C_m$=%s' % str(cms[i]))
            lastdiff1_m2.append(rd1_m2[-1])
ax1.legend(loc='lower right')
ax2.legend(loc='lower right')
ax3.legend(loc='lower right')
ax4.legend(loc='lower right')
ax5.legend(loc='lower right')
ax1.set_title(r'A',loc='left',fontsize=18)
ax2.set_title(r'B',loc='left',fontsize=18)
ax3.set_title(r'C',loc='left',fontsize=18)
ax4.set_title(r'D',loc='left',fontsize=18)
ax5.set_title(r'E',loc='left',fontsize=18)
ax3.set_title(r'Allen model %i' % testmodels[0],fontsize=18)
ax4.set_title(r'Allen model %i' % testmodels[1],fontsize=18)
ax5.set_title(r'Allen model %i' % testmodels[2],fontsize=18)
ax3.set_ylabel(r'Relative difference in $f$',fontsize=16)
ax4.set_xlabel(r'$I$ (nA)',fontsize=16)
ax4.set_ylabel(r'Relative difference in $f$',fontsize=16)
ax5.set_ylabel(r'Relative difference in $f$',fontsize=16)
ax5.set_xlabel(r'$I$ (nA)',fontsize=16)
# Convert to arrays?
currat_bas    = []
maxdiff_bas   = []
maxdiff1_bas  = []
maxdiff2_bas  = []
lastdiff1_bas = []
for i in range(Ncm): # Should I do this for somaprox too?
    if i!=3:
        im_bas, rd1_bas, rd2_bas = reldiffs(Nspikes_everywhere[i],Nspikes_everywhere[3],I_Nspikes_everywhere[i],I_Nspikes_everywhere[3]) 
        if len(rd1_bas)>0:
            maxrd1 = max(rd1_bas)
            minrd1 = min(rd1_bas)
            if maxrd1<abs(minrd1):
                maxrd1 = minrd1
            maxdiff1_bas.append(maxrd1)
        if len(rd2_bas)>0:
            maxrd2 = max(rd2_bas)
            minrd2 = min(rd2_bas)
            if maxrd2<abs(minrd2):
                maxrd2 = minrd2
            maxdiff2_bas.append(maxrd2)
        if len(rd1_bas)>0 and len(rd2_bas)>0:
            if abs(maxrd1)<abs(maxrd2):
                maxrd = maxrd2
            else:
                maxrd = maxrd1
            maxdiff_bas.append(maxrd)
        if len(rd1_bas)>0 or len(rd2_bas)>0:
            currat_bas.append(im_bas[np.argmax(abs(rd1_bas))])
        if len(rd1_bas)>0:
            ax2.plot(im_bas, rd1_bas,label='$C_m$=%s' % str(cms[i]))
            lastdiff1_bas.append(rd1_bas[-1])
ax2.set_title(r'Ball-and-stick HH',fontsize=18)

currat_oc    = []
maxdiff_oc   = []
maxdiff1_oc  = []
maxdiff2_oc  = []
lastdiff1_oc = []
for i in range(Ncm):
    if i!=3:
        im_oc, rd1_oc, rd2_oc = reldiffs(Nspikes_onecomp[i],Nspikes_onecomp[3],I_Nspikes_onecomp[i],I_Nspikes_onecomp[3]) 
        if len(rd1_oc)>0:
            maxrd1 = max(rd1_oc)
            minrd1 = min(rd1_oc)
            if maxrd1<abs(minrd1):
                maxrd1 = minrd1
            maxdiff1_oc.append(maxrd1)
        if len(rd2_oc)>0:
            maxrd2 = max(rd2_oc)
            minrd2 = min(rd2_oc)
            if maxrd2<abs(minrd2):
                maxrd2 = minrd2
            maxdiff2_oc.append(maxrd2)
        if len(rd1_oc)>0 and len(rd2_oc):
            if abs(maxrd1)<abs(maxrd2):
                maxrd = maxrd2
            else:
                maxrd = maxrd1
            maxdiff_oc.append(maxrd)
        if len(rd1_oc)>0 or len(rd2_oc):
            currat_oc.append(im_oc[np.argmax(abs(rd1_oc))])
        if len(rd1_oc)>0:
            ax1.plot(im_oc, rd1_oc,label='$C_m$=%s' % str(cms[i]))
            lastdiff1_oc.append(rd1_oc[-1])
ax1.set_title(r'One compartment HH',fontsize=18)
ax1.set_ylabel(r'Relative difference in $f$',fontsize=16)

fig.tight_layout()
plt.savefig(plotname_rdf)

print('------- Last difference --------')
print('cms:',cms)
print('lastdiff_oc:',lastdiff1_oc)
print('lastdiff1_bas:',lastdiff1_bas)
print('lastdiff1_m0:',lastdiff1_m0)
print('lastdiff1_m1:',lastdiff1_m1)
print('lastdiff1_m2:',lastdiff1_m2)

print('------- Max difference --------')
print('cms:',cms)
print('maxdiff1_oc:',maxdiff1_oc)
print('maxdiff1_bas:',maxdiff1_bas)
print('maxdiff0_v1:',maxdiff0_v1)
print('maxdiff1_v1:',maxdiff1_v1)
print('maxdiff2_v1:',maxdiff2_v1)


fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(15,15),dpi=300)
ax1.plot(cms[:-1],maxdiff1_oc)
ax2.plot(cms[:-2],maxdiff1_bas)
ax3.plot(cms[:-1],maxdiff0_v1)
ax4.plot(cms[:-2],maxdiff1_v1)
ax5.plot(cms[:-1],maxdiff2_v1)
ax1.set_title(r'A',loc='left',fontsize=18)
ax2.set_title(r'B',loc='left',fontsize=18)
ax3.set_title(r'C',loc='left',fontsize=18)
ax4.set_title(r'D',loc='left',fontsize=18)
ax5.set_title(r'E',loc='left',fontsize=18)
ax1.set_title(r'One compartment HH',fontsize=18)
ax1.set_title(r'Ball-and-stick HH',fontsize=18)
ax3.set_title(r'Allen model %i' % testmodels[0],fontsize=18)
ax4.set_title(r'Allen model %i' % testmodels[1],fontsize=18)
ax5.set_title(r'Allen model %i' % testmodels[2],fontsize=18)
ax1.set_ylabel(r'Max. relative difference in $f$',fontsize=16)
ax3.set_ylabel(r'Max. relative difference in $f$',fontsize=16)
ax4.set_xlabel(r'$C_m$ ($\mu$F/cm$^2$)',fontsize=16)
ax4.set_ylabel(r'Max. relative difference in $f$',fontsize=16)
ax5.set_ylabel(r'Max. relative difference in $f$',fontsize=16)
ax5.set_xlabel(r'$C_m$ ($\mu$F/cm$^2$)',fontsize=16)

plt.savefig(plotname_rdf_maxdiff)

fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(15,15),dpi=300)
ax1.plot(cms[:-1],lastdiff1_oc)
ax2.plot(cms[:-2],lastdiff1_bas)
ax3.plot(cms[:-1],lastdiff1_m0)
ax4.plot(cms[:-2],lastdiff1_m1)
ax5.plot(cms[:-1],lastdiff1_m2)
ax1.set_title(r'A',loc='left',fontsize=18)
ax2.set_title(r'B',loc='left',fontsize=18)
ax3.set_title(r'C',loc='left',fontsize=18)
ax4.set_title(r'D',loc='left',fontsize=18)
ax5.set_title(r'E',loc='left',fontsize=18)
ax1.set_title(r'One compartment HH',fontsize=18)
ax1.set_title(r'Ball-and-stick HH',fontsize=18)
ax3.set_title(r'Allen model %i' % testmodels[0],fontsize=18)
ax4.set_title(r'Allen model %i' % testmodels[1],fontsize=18)
ax5.set_title(r'Allen model %i' % testmodels[2],fontsize=18)
ax1.set_ylabel(r'Relative difference in $f$ at highest stim.',fontsize=16)
ax3.set_ylabel(r'Relative difference in $f$ at highest stim.',fontsize=16)
ax4.set_xlabel(r'$C_m$ ($\mu$F/cm$^2$)',fontsize=16)
ax4.set_ylabel(r'Relative difference in $f$ at highest stim.',fontsize=16)
ax5.set_ylabel(r'Relative difference in $f$ at highest stim.',fontsize=16)
ax5.set_xlabel(r'$C_m$ ($\mu$F/cm$^2$)',fontsize=16)

plt.savefig(plotname_rdf_lastdiff)

plt.show()
