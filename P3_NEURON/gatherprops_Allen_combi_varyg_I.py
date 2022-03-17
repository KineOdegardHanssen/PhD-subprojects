import numpy as np
import matplotlib.pyplot as plt

# change the default font family
plt.rcParams.update({'font.family':'Arial'})

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend', fontsize=13)

mylinewidth = 2

testmodels = [488462965]#[478513437,488462965,478513407]
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

varyg_bool = True # Easiest to do this. I think.
namestringfirst = ''
varylist = [] # Should be redundant
namestringfirst = ''
plotstring   = ''
varygbool    = True
varyIh       = False
vary_NaV     = False
vary_Kd      = False
vary_Kv2like = True
vary_Kv3_1   = False
vary_K_T     = False
vary_Im_v2   = False
vary_SK      = False
vary_Ca_HVA  = False
vary_Ca_LVA  = False
vary_gpas    = False 

if varyIh==True:
    namestringfirst = namestringfirst + '_gIh'
    plotstring = plotstring + '_gIh'
    plotlabel = r'$\bar{g}_{\mathregular{Ih}}$'
if vary_NaV==True:
    namestringfirst = namestringfirst + '_gNaV'
    plotstring = plotstring + '_gNaV'
    plotlabel = r'$\bar{g}_{\mathregular{NaV}}$'
if vary_Kd==True:
    namestringfirst = namestringfirst + '_gKd'
    plotstring = plotstring + '_gKd'
    plotlabel = r'$\bar{g}_{\mathregular{Kd}}$'
if vary_Kv2like==True:
    namestringfirst = namestringfirst + '_gKv2like'
    plotstring = plotstring + '_gKv2like'
    plotlabel = r'$\bar{g}_{\mathregular{Kv2like}}$'
if vary_Kv3_1==True:
    namestringfirst = namestringfirst + '_gKv31'
    plotstring = plotstring + '_gKv31'
    plotlabel = r'$\bar{g}_{\mathregular{Kv3 1}}$'
if vary_K_T==True:
    namestringfirst = namestringfirst + '_gKT'
    plotstring = plotstring + '_gKT'
    plotlabel = r'$\bar{g}_{\mathregular{KT}}$'
if vary_Im_v2==True:
    namestringfirst = namestringfirst + '_gImv2'
    plotstring = plotstring + '_gImv2'
    plotlabel = r'$\bar{g}_{\mathregular{Im v2}}$'
if vary_SK==True:
    namestringfirst = namestringfirst + '_gSK'
    plotstring = plotstring + '_gSK'
    plotlabel = r'$\bar{g}_{\mathregular{SK}}$'
if vary_Ca_HVA==True:
    namestringfirst = namestringfirst + '_gCaHVA'
    plotstring = plotstring + '_gCaHVA'
    plotlabel = r'$\bar{g}_{\mathregular{Ca HVA}}$'
if vary_Ca_LVA==True:
    namestringfirst = namestringfirst + '_gCaLVA'
    plotstring = plotstring + '_gCaLVA'
    plotlabel = r'$\bar{g}_{\mathregular{Ca LVA}}$'
if vary_gpas==True: 
    namestringfirst = namestringfirst + '_gpas'
    plotstring = plotstring + '_gpas'
    plotlabel = r'$\bar{g}_{\mathregular{pas}}$'

varyg = [0.0004,0.0008,0.0016,0.0032,0.0064]




changestring =''
changestring = changestring+'_g'+str(varyg)+'_gdf'

i_master = [0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5]
cm = 1.0 # Keep this a scalar for now
    
NI  = len(i_master)
Ng = len(varyg)

i_master_everywhere_all         = []
Nspikes_everywhere_all          = []
avg_AP_ampl_everywhere_all      = []
rms_AP_ampl_everywhere_all      = []
avg_AP_mins_everywhere_all      = []
rms_AP_mins_everywhere_all      = []
avg_AP_halfwidth_everywhere_all = []
rms_AP_halfwidth_everywhere_all = []
avg_ISI_everywhere_all          = []
rms_ISI_everywhere_all          = []
I_Nspikes_everywhere_all        = []
I_AP_ampl_everywhere_all        = []
I_AP_mins_everywhere_all        = []
I_AP_halfwidth_everywhere_all   = []
I_ISI_everywhere_all            = []

i_master_sprx_all         = []
Nspikes_sprx_all          = []
avg_AP_ampl_sprx_all      = []
rms_AP_ampl_sprx_all      = []
avg_AP_mins_sprx_all      = []
rms_AP_mins_sprx_all      = []
avg_AP_halfwidth_sprx_all = []
rms_AP_halfwidth_sprx_all = []
avg_ISI_sprx_all          = []
rms_ISI_sprx_all          = []
I_Nspikes_sprx_all        = []
I_AP_ampl_sprx_all        = []
I_AP_mins_sprx_all        = []
I_AP_halfwidth_sprx_all   = []
I_ISI_sprx_all            = []

for testmodel in testmodels:
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

    infolder      = 'figures/%i/' % (testmodel)
    vrestfolder   = infolder 
    for g in varyg:
        Nspikes          = []
        avg_AP_ampl      = []
        rms_AP_ampl      = []
        avg_AP_mins      = []
        rms_AP_mins      = []
        avg_AP_halfwidth = []
        rms_AP_halfwidth = []
        avg_ISI          = []
        rms_ISI          = []
        I_Nspikes       = []
        I_AP_ampl       = []
        I_AP_mins       = []
        I_AP_halfwidth  = []
        I_ISI           = []
        
        namestring = namestringfirst+str(g)
        
        # Set names
        infilename_Nspikes = infolder+'%s_%i_cmfall'%(namestring,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
        infilename_APampl = infolder+'%s_%i_cmfall'%(namestring,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Vmax_vs_I.txt'
        infilename_APmins = infolder+'%s_%i_cmfall'%(namestring,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Vmin_vs_I.txt'
        infilename_APdhw  = infolder+'%s_%i_cmfall'%(namestring,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_sdurat%s_vs_I.txt' % str(spikedurat)
        infilename_ISI   = infolder+'%s_%i_cmfall'%(namestring,testmodel)+'_idur%i_varyiamp'% idur+'_manual_ISI_vs_I.txt'
        # Read files
        infile_Nspikes = open(infilename_Nspikes,'r')
        infile_APampl  = open(infilename_APampl,'r')
        infile_APmins  = open(infilename_APmins,'r')
        infile_APdhw   = open(infilename_APdhw,'r')
        infile_ISI     = open(infilename_ISI,'r')
    
        lines_Nspikes = infile_Nspikes.readlines()
        lines_APampl  = infile_APampl.readlines()
        lines_APmins  = infile_APmins.readlines()
        lines_APdhw   = infile_APdhw.readlines()
        lines_ISI     = infile_ISI.readlines()
        Nlines_Nspikes = len(lines_Nspikes)
        Nlines_APampl  = len(lines_APampl)
        Nlines_APmins  = len(lines_APmins)
        Nlines_APdhw   = len(lines_APdhw)
        Nlines_ISI     = len(lines_ISI)
    
        for i in range(Nlines_Nspikes):
            words_Nspikes = lines_Nspikes[i].split()
            if len(words_Nspikes)>0:
                I_Nspikes.append(float(words_Nspikes[0]))
                Nspikes.append(float(words_Nspikes[1]))
    
        for i in range(Nlines_APampl):
            words_APampl  = lines_APampl[i].split()
            if len(words_APampl)>0: # One test should be enough, but better safe than sorry
                I_AP_ampl.append(float(words_APampl[0]))
                avg_AP_ampl.append(float(words_APampl[1]))
                rms_AP_ampl.append(float(words_APampl[2]))

        for i in range(Nlines_APmins):
            words_APmins  = lines_APmins[i].split()
            if len(words_APmins)>0:
                I_AP_mins.append(float(words_APmins[0]))
                avg_AP_mins.append(float(words_APmins[1]))
                rms_AP_mins.append(float(words_APmins[2]))
    
        for i in range(Nlines_APdhw):
            words_APdhw   = lines_APdhw[i].split()
            if len(words_APdhw)>0:
                I_AP_halfwidth.append(float(words_APdhw[0]))
                avg_AP_halfwidth.append(float(words_APdhw[1]))
                rms_AP_halfwidth.append(float(words_APdhw[2]))
    
        for i in range(Nlines_ISI):
            words_ISI     = lines_ISI[i].split()
            if len(words_ISI)>0:
                I_ISI.append(float(words_ISI[0]))
                avg_ISI.append(float(words_ISI[1]))
                rms_ISI.append(float(words_ISI[2]))

        infile_Nspikes.close()
        infile_APampl.close()
        infile_APmins.close()
        infile_APdhw.close()
        infile_ISI.close()
        
        Nspikes_everywhere.append(Nspikes)
        avg_AP_ampl_everywhere.append(avg_AP_ampl)
        rms_AP_ampl_everywhere.append(rms_AP_ampl)
        avg_AP_mins_everywhere.append(avg_AP_mins)
        rms_AP_mins_everywhere.append(rms_AP_mins)
        avg_AP_halfwidth_everywhere.append(avg_AP_halfwidth)
        rms_AP_halfwidth_everywhere.append(rms_AP_halfwidth)
        avg_ISI_everywhere.append(avg_ISI)
        rms_ISI_everywhere.append(rms_ISI)
        I_Nspikes_everywhere.append(I_Nspikes)
        I_AP_ampl_everywhere.append(I_AP_ampl)
        I_AP_mins_everywhere.append(I_AP_mins)
        I_AP_halfwidth_everywhere.append(I_AP_halfwidth)
        I_ISI_everywhere.append(I_ISI)
    
    
        ############# SOMAPROX ##################################
        Nspikes          = []
        avg_AP_ampl      = []
        rms_AP_ampl      = []
        avg_AP_mins      = []
        rms_AP_mins      = []
        avg_AP_halfwidth = []
        rms_AP_halfwidth = []
        avg_ISI          = []
        rms_ISI          = []
        I_Nspikes       = []
        I_AP_ampl       = []
        I_AP_mins       = []
        I_AP_halfwidth  = []
        I_ISI           = []
        
        # Set names
        infilename_Nspikes = infolder+'%s_%i_cmfsprx'%(namestring,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
        infilename_APampl = infolder+'%s_%i_cmfsprx'%(namestring,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Vmax_vs_I.txt'
        infilename_APmins = infolder+'%s_%i_cmfsprx'%(namestring,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Vmin_vs_I.txt'
        infilename_APdhw  = infolder+'%s_%i_cmfsprx'%(namestring,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_sdurat%s_vs_I.txt' % str(spikedurat)
        infilename_ISI   = infolder+'%s_%i_cmfsprx'%(namestring,testmodel)+'_idur%i_varyiamp'% idur+'_manual_ISI_vs_I.txt'
        # Read files
        infile_Nspikes = open(infilename_Nspikes,'r')
        infile_APampl  = open(infilename_APampl,'r')
        infile_APmins  = open(infilename_APmins,'r')
        infile_APdhw   = open(infilename_APdhw,'r')
        infile_ISI     = open(infilename_ISI,'r')
    
        lines_Nspikes = infile_Nspikes.readlines()
        lines_APampl  = infile_APampl.readlines()
        lines_APmins  = infile_APmins.readlines()
        lines_APdhw   = infile_APdhw.readlines()
        lines_ISI     = infile_ISI.readlines()
        Nlines_Nspikes = len(lines_Nspikes)
        Nlines_APampl  = len(lines_APampl)
        Nlines_APmins  = len(lines_APmins)
        Nlines_APdhw   = len(lines_APdhw)
        Nlines_ISI     = len(lines_ISI)
    
        for i in range(Nlines_Nspikes):
            words_Nspikes = lines_Nspikes[i].split()
            if len(words_Nspikes)>0:
                I_Nspikes.append(float(words_Nspikes[0]))
                Nspikes.append(float(words_Nspikes[1]))
    
        for i in range(Nlines_APampl):
            words_APampl  = lines_APampl[i].split()
            if len(words_APampl)>0: # One test should be enough, but better safe than sorry
                I_AP_ampl.append(float(words_APampl[0]))
                avg_AP_ampl.append(float(words_APampl[1]))
                rms_AP_ampl.append(float(words_APampl[2]))

        for i in range(Nlines_APmins):
            words_APmins  = lines_APmins[i].split()
            if len(words_APmins)>0:
                I_AP_mins.append(float(words_APmins[0]))
                avg_AP_mins.append(float(words_APmins[1]))
                rms_AP_mins.append(float(words_APmins[2]))
    
        for i in range(Nlines_APdhw):
            words_APdhw   = lines_APdhw[i].split()
            if len(words_APdhw)>0:
                I_AP_halfwidth.append(float(words_APdhw[0]))
                avg_AP_halfwidth.append(float(words_APdhw[1]))
                rms_AP_halfwidth.append(float(words_APdhw[2]))
    
        for i in range(Nlines_ISI):
            words_ISI     = lines_ISI[i].split()
            if len(words_ISI)>0:
                I_ISI.append(float(words_ISI[0]))
                avg_ISI.append(float(words_ISI[1]))
                rms_ISI.append(float(words_ISI[2]))

        infile_Nspikes.close()
        infile_APampl.close()
        infile_APmins.close()
        infile_APdhw.close()
        infile_ISI.close()
        
        Nspikes_sprx.append(Nspikes)
        avg_AP_ampl_sprx.append(avg_AP_ampl)
        rms_AP_ampl_sprx.append(rms_AP_ampl)
        avg_AP_mins_sprx.append(avg_AP_mins)
        rms_AP_mins_sprx.append(rms_AP_mins)
        avg_AP_halfwidth_sprx.append(avg_AP_halfwidth)
        rms_AP_halfwidth_sprx.append(rms_AP_halfwidth)
        avg_ISI_sprx.append(avg_ISI)
        rms_ISI_sprx.append(rms_ISI)
        I_Nspikes_sprx.append(I_Nspikes)
        I_AP_ampl_sprx.append(I_AP_ampl)
        I_AP_mins_sprx.append(I_AP_mins)
        I_AP_halfwidth_sprx.append(I_AP_halfwidth)
        I_ISI_sprx.append(I_ISI)
    Nspikes_everywhere_all.append(Nspikes_everywhere)
    avg_AP_ampl_everywhere_all.append(avg_AP_ampl_everywhere)
    rms_AP_ampl_everywhere_all.append(rms_AP_ampl_everywhere)
    avg_AP_mins_everywhere_all.append(avg_AP_mins_everywhere)
    rms_AP_mins_everywhere_all.append(rms_AP_mins_everywhere)
    avg_AP_halfwidth_everywhere_all.append(avg_AP_halfwidth_everywhere)
    rms_AP_halfwidth_everywhere_all.append(rms_AP_halfwidth_everywhere)
    avg_ISI_everywhere_all.append(avg_ISI_everywhere)
    rms_ISI_everywhere_all.append(rms_ISI_everywhere)
    I_Nspikes_everywhere_all.append(I_Nspikes_everywhere)
    I_AP_ampl_everywhere_all.append(I_AP_ampl_everywhere)
    I_AP_mins_everywhere_all.append(I_AP_mins_everywhere)
    I_AP_halfwidth_everywhere_all.append(I_AP_halfwidth_everywhere)
    I_ISI_everywhere_all.append(I_ISI_everywhere)
    
    Nspikes_sprx_all.append(Nspikes_sprx)
    avg_AP_ampl_sprx_all.append(avg_AP_ampl_sprx)
    rms_AP_ampl_sprx_all.append(rms_AP_ampl_sprx)
    avg_AP_mins_sprx_all.append(avg_AP_mins_sprx)
    rms_AP_mins_sprx_all.append(rms_AP_mins_sprx)
    avg_AP_halfwidth_sprx_all.append(avg_AP_halfwidth_sprx)
    rms_AP_halfwidth_sprx_all.append(rms_AP_halfwidth_sprx)
    avg_ISI_sprx_all.append(avg_ISI_sprx)
    rms_ISI_sprx_all.append(rms_ISI_sprx)
    I_Nspikes_sprx_all.append(I_Nspikes_sprx)
    I_AP_ampl_sprx_all.append(I_AP_ampl_sprx)
    I_AP_mins_sprx_all.append(I_AP_mins_sprx)
    I_AP_halfwidth_sprx_all.append(I_AP_halfwidth_sprx)
    I_ISI_sprx_all.append(I_ISI_sprx)

# Plotting
    
color_pv_all   = [['#1f77b4','darkblue','#7f7f7f','b','tab:blue'],['#d62728','#ff7f0e','#d62728','m','tab:pink'],['#2ca02c','#8c564b','#9467bd','g','tab:brown']]
color_pv_sprx  = [['rebeccapurple','#e377c2','#8c564b','c','tab:purple'],['#9467bd','#bcbd22','#17becf','xkcd:sky blue','tab:orange'],['olive','darkgreen','mediumseagreen','tab:green','tab:olive']]

plotfolder = 'figures/Comparemodels/'
plotname = plotfolder+'features_vs_I_cmgraphs_combi_AllenPV.png'

## avg and rms:
#fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(17,15))
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15,15))
fig.suptitle('Vary %s' % plotlabel,fontsize=20)

ax1.set_title(r'Frequency $f$ vs $I$',fontsize=16)
for i in range(Nmodels):
    color_these_all = color_pv_all[i]
    color_these_spx = color_pv_sprx[i]
    I_Nspikes_everywhere_this = I_Nspikes_everywhere_all[i]
    I_Nspikes_sprx_this       = I_Nspikes_sprx_all[i]
    Nspikes_everywhere_this   = Nspikes_everywhere_all[i]
    Nspikes_sprx_this         = Nspikes_sprx_all[i]
    for j in range(Ng):
    ### Everywhere:
        ax1.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],
color=color_these_all[j],label=r'%i all: %s=%.2e' % (testmodels[i],plotlabel,varyg[j]), linewidth=mylinewidth)
        ### Somaprox:
        ax1.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',
color=color_these_spx[j],label=r'%i somaprox: %s=%.2e' % (testmodels[i],plotlabel,varyg[j]), linewidth=mylinewidth)
    ''' # May rewrite.
    f1 = 0
    f2_all = 0
    f2_sprx = 0
    for j in range(len(I_Nspikes_everywhere_this1)):
        Ithis = I_Nspikes_everywhere_this1
        if Ithis[j]==0.44:
            f1 = Nspikes_everywhere_this1[j]
    for j in range(len(I_Nspikes_everywhere_this2)):
        Ithis = I_Nspikes_everywhere_this2
        if Ithis[j]==0.44:
            f2_all = Nspikes_everywhere_this2[j]
    for j in range(len(I_Nspikes_sprx_this2)):
        Ithis = I_Nspikes_sprx_this2
        if Ithis[j]==0.44:
            f2_sprx = Nspikes_sprx_this2[j]
    if f1!=0:
        if f2_all!=0:
            print('Model:',testmodels[i],'Change of f (all):', (f1-f2_all)/f1)
        if f2_sprx!=0:
            print('Model:',testmodels[i],'Change of f (spx):', (f1-f2_sprx)/f1)
    ''' 
ax1.set_ylabel('$f$ [Hz]',fontsize=14)

ax2.set_title(r'Spike minima and maxima vs $I$',fontsize=16)
for i in range(Nmodels):
    color_these_all = color_pv_all[i]
    color_these_spx = color_pv_sprx[i]
    ### AP AMPL:
    I_AP_ampl_everywhere_this   = I_AP_ampl_everywhere_all[i]
    I_AP_ampl_sprx_this         = I_AP_ampl_sprx_all[i]
    avg_AP_ampl_everywhere_this = avg_AP_ampl_everywhere_all[i]
    avg_AP_ampl_sprx_this       = avg_AP_ampl_sprx_all[i]
    rms_AP_ampl_everywhere_this = rms_AP_ampl_everywhere_all[i]
    rms_AP_ampl_sprx_this       = rms_AP_ampl_sprx_all[i]
    # Assuming Cm=1.0 and Cm=1.25 is all we need: # We want stipled lines!
    I_AP_ampl_everywhere_this1 = I_AP_ampl_everywhere_this[0]
    I_AP_ampl_everywhere_this2 = I_AP_ampl_everywhere_this[1]
    I_AP_ampl_sprx_this1       = I_AP_ampl_sprx_this[0]
    I_AP_ampl_sprx_this2       = I_AP_ampl_sprx_this[1]
    avg_AP_ampl_everywhere_this1   = avg_AP_ampl_everywhere_this[0]
    avg_AP_ampl_everywhere_this2   = avg_AP_ampl_everywhere_this[1]
    avg_AP_ampl_sprx_this1         = avg_AP_ampl_sprx_this[0]
    avg_AP_ampl_sprx_this2         = avg_AP_ampl_sprx_this[1]
    rms_AP_ampl_everywhere_this1   = rms_AP_ampl_everywhere_this[0]
    rms_AP_ampl_everywhere_this2   = rms_AP_ampl_everywhere_this[1]
    rms_AP_ampl_sprx_this1         = rms_AP_ampl_sprx_this[0]
    rms_AP_ampl_sprx_this2         = rms_AP_ampl_sprx_this[1]
    #### AP MIN
    I_AP_mins_everywhere_this   = I_AP_mins_everywhere_all[i]
    I_AP_mins_sprx_this         = I_AP_mins_sprx_all[i]
    avg_AP_mins_everywhere_this = avg_AP_mins_everywhere_all[i]
    avg_AP_mins_sprx_this       = avg_AP_mins_sprx_all[i]
    rms_AP_mins_everywhere_this = rms_AP_mins_everywhere_all[i]
    rms_AP_mins_sprx_this       = rms_AP_mins_sprx_all[i]
    for j in range(Ng):
        ### Everywhere:
        ax2.errorbar(I_AP_ampl_everywhere_this[j], avg_AP_ampl_everywhere_this[j], yerr=rms_AP_ampl_everywhere_this[j], color=color_these_all[j], capsize=2, linewidth=mylinewidth)
        ax2.errorbar(I_AP_mins_everywhere_this[j], avg_AP_mins_everywhere_this[j], yerr=rms_AP_mins_everywhere_this[j], color=color_these_all[j],capsize=2, linewidth=mylinewidth)
        ### Somaprox:
        ax2.errorbar(I_AP_ampl_sprx_this[j], avg_AP_ampl_sprx_this[j], yerr=rms_AP_ampl_sprx_this[j], ls='--', color=color_these_spx[j],capsize=2, linewidth=mylinewidth)
        ax2.errorbar(I_AP_mins_sprx_this[j], avg_AP_mins_sprx_this[j], yerr=rms_AP_mins_sprx_this[j], ls='--', color=color_these_spx[j],capsize=2, linewidth=mylinewidth)


#I_AP_halfwidth_everywhere_all = I_AP_halfwidth_everywhere_all
#avg_AP_halfwidth_everywhere_all = avg_AP_halfwidth_everywhere_all
#rms_AP_halfwidth_everywhere_all = rms_AP_halfwidth_everywhere_all
#avg_AP_halfwidth_sprx_all = avg_AP_halfwidth_sprx_all
#rms_AP_halfwidth_sprx_all = rms_AP_halfwidth_sprx_all

ax3.set_title(r'Spike duration at %s mV vs $I$' % str(spikedurat),fontsize=16)
for i in range(Nmodels):
    color_these_all = color_pv_all[i]
    color_these_spx = color_pv_sprx[i]
    ### AP halfwidth:
    I_AP_halfwidth_everywhere_this   = I_AP_halfwidth_everywhere_all[i]
    I_AP_halfwidth_sprx_this         = I_AP_halfwidth_sprx_all[i]
    avg_AP_halfwidth_everywhere_this = avg_AP_halfwidth_everywhere_all[i]
    avg_AP_halfwidth_sprx_this       = avg_AP_halfwidth_sprx_all[i]
    rms_AP_halfwidth_everywhere_this = rms_AP_halfwidth_everywhere_all[i]
    rms_AP_halfwidth_sprx_this       = rms_AP_halfwidth_sprx_all[i]
    for j in range(Ng):
        ### Everywhere:
        ax3.errorbar(I_AP_halfwidth_everywhere_this[j], avg_AP_halfwidth_everywhere_this[j], yerr=rms_AP_halfwidth_everywhere_this[j], color=color_these_all[j],capsize=2, linewidth=mylinewidth,label=r'%i all: %s=%.2e' % (testmodels[i],plotlabel,varyg[j]))
        ### Somaprox:
        ax3.errorbar(I_AP_halfwidth_sprx_this[j], avg_AP_halfwidth_sprx_this[j], yerr=rms_AP_halfwidth_sprx_this[j],color=color_these_spx[j],capsize=2, ls='--', linewidth=mylinewidth,label=r'%i somaprox: %s=%.2e' % (testmodels[i],plotlabel,varyg[j]))
ax3.set_xlabel(r'$I$ [nA]',fontsize=14)
ax3.set_ylabel('Spike duration at %s mV [ms]' % str(spikedurat),fontsize=14)
ax3.legend(bbox_to_anchor =(0.5,-0.27), loc='lower center',ncol=3)

#I_ISI_sprx_all = I_ISI_sprx_all
#I_ISI_everywhere_all = I_ISI_everywhere_all
#avg_ISI_everywhere_all = avg_ISI_everywhere_all
#rms_ISI_everywhere_all = rms_ISI_everywhere_all
#avg_ISI_sprx_all = avg_ISI_sprx_all
#rms_ISI_sprx_all = rms_ISI_sprx_all

ax4.set_title(r'Interspike interval vs $I$',fontsize=16)
for i in range(Nmodels):
    color_these_all = color_pv_all[i]
    color_these_spx = color_pv_sprx[i]
    ### ISI:
    I_ISI_everywhere_this   = I_ISI_everywhere_all[i]
    I_ISI_sprx_this         = I_ISI_sprx_all[i]
    avg_ISI_everywhere_this = avg_ISI_everywhere_all[i]
    avg_ISI_sprx_this       = avg_ISI_sprx_all[i]
    rms_ISI_everywhere_this = rms_ISI_everywhere_all[i]
    rms_ISI_sprx_this       = rms_ISI_sprx_all[i]
    # Assuming Cm=1.0 and Cm=1.25 is all we need: # We want stipled lines!
    for j in range(Ng):
        ### Everywhere:
        ax4.errorbar(I_ISI_everywhere_this[j], avg_ISI_everywhere_this[j], yerr=rms_ISI_everywhere_this[j],color=color_these_all[j],capsize=2,label=r'%i all: %s=%.2e' % (testmodels[i],plotlabel,varyg[j]), linewidth=mylinewidth)
        ### Somaprox:
        ax4.errorbar(I_ISI_sprx_this[j], avg_ISI_sprx_this[j], yerr=rms_ISI_sprx_this[j], ls='--', color=color_these_spx[j],capsize=2,label=r'%i somaprox: %s=%.2e' % (testmodels[i],plotlabel,varyg[j]), linewidth=mylinewidth)
ax4.set_xlabel(r'$I$ [nA]',fontsize=14)
ax4.set_ylabel('Interspike interval [ms]',fontsize=14)
#ax4.legend(loc='upper right')
#fig.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.savefig(plotname)

'''
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all[-1]
I_Nspikes_sprx_this       = I_Nspikes_sprx_all[-1]
Nspikes_everywhere_this   = Nspikes_everywhere_all[-1]
Nspikes_sprx_this         = Nspikes_sprx_all[-1]
# Assuming Cm=1.0 and Cm=1.25 is all we need: # We want stipled lines!
I_Nspikes_everywhere_this1 = I_Nspikes_everywhere_this[0]
I_Nspikes_everywhere_this2 = I_Nspikes_everywhere_this[1]
I_Nspikes_sprx_this1       = I_Nspikes_sprx_this[0]
I_Nspikes_sprx_this2       = I_Nspikes_sprx_this[1]
Nspikes_everywhere_this1   = Nspikes_everywhere_this[0]
Nspikes_everywhere_this2   = Nspikes_everywhere_this[1]
Nspikes_sprx_this1         = Nspikes_sprx_this[0]
Nspikes_sprx_this2         = Nspikes_sprx_this[1]
f1 = 0
f2_all = 0
f2_sprx = 0
for j in range(len(I_Nspikes_everywhere_this1)):
    Ithis = I_Nspikes_everywhere_this1
    if Ithis[j]==0.4:
        f1 = Nspikes_everywhere_this1[j]
for j in range(len(I_Nspikes_everywhere_this2)):
    Ithis = I_Nspikes_everywhere_this2
    if Ithis[j]==0.4:
        f2_all = Nspikes_everywhere_this2[j]
for j in range(len(I_Nspikes_sprx_this2)):
    Ithis = I_Nspikes_sprx_this2
    if Ithis[j]==0.4:
        f2_sprx = Nspikes_sprx_this2[j]
if f1!=0:
    if f2_all!=0:
        print('Model:',testmodels[i],'Change of f (all):', (f1-f2_all)/f1)
    if f2_sprx!=0:
        print('Model:',testmodels[i],'Change of f (spx):', (f1-f2_sprx)/f1)
'''

# Should I have checked thresholds too? --> A lot of work, a lot of data.
#					--> Can always do afterwards and add
#					--> WILL probably need this afterwards
plt.show()


    