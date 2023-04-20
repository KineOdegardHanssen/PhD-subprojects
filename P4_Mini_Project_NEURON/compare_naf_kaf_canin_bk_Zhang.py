import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

# change the default font family
plt.rcParams.update({'font.family':'Arial'})

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

spikedurat = -40
idur       = 1000 #100 # ms
idelay     = 10
v_init     = -86.8 # mV
Ra         = 100
somasize   = 10 # 15 # 
dtexp      = -7

fig = plt.figure(figsize=(18,8),dpi=300)

gs = gridspec.GridSpec(2, 6)

ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
ax3 = plt.subplot(gs[0, 4:6])
ax4 = plt.subplot(gs[1, 0:2])
ax5 = plt.subplot(gs[1, 2:4])
ax6 = plt.subplot(gs[1, 4:6])

#fig.suptitle(r'Properties',fontsize=20)

ax1.set_title(r'A',loc='left',fontsize=18)
ax2.set_title(r'B',loc='left',fontsize=18)
ax3.set_title(r'C',loc='left',fontsize=18)
ax4.set_title(r'D',loc='left',fontsize=18)
ax5.set_title(r'E',loc='left',fontsize=18)
ax6.set_title(r'F',loc='left',fontsize=18)
ax1.set_title(r'$f$',fontsize=18)
ax2.set_title(r'$V_{\mathregular{max}}$',fontsize=18)
ax3.set_title(r'$V_{\mathregular{min}}$',fontsize=18)
ax4.set_title(r'Duration at %i mV' %spikedurat,fontsize=18)
ax5.set_title(r'ISI',fontsize=18)
ax6.set_title(r'$V_{\mathregular{rest}}$',fontsize=18)

    
model_folders = ['','canin_Konstantoudaki/','canin_Konstantoudaki/BK_I_Zhang/']
labels = ['base','canin','canin, BK']
Nmodels = len(model_folders)

outfolder = 'Compare/Soma%i/' % somasize
plotname = outfolder+'compare_properties_naf_kaf_canin_BK_Zhang.png'

'''# I don't quite know how to write this to file
outfilename_Nspikes = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I.txt'
outfilename_APampl  = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Vmax_vs_I.txt'
outfilename_APmins  = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Vmin_vs_I.txt'
outfilename_APdhw   = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_sdurat%s_vs_I.txt' % str(spikedurat)
outfilename_ISI     = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_ISI_vs_I.txt'

plotname_Nspikes = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I.png'
plotname_APampl  = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Vmax_vs_I.png'
plotname_APmins  = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Vmin_vs_I.png'
plotname_APdhw   = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_sdurat%s_vs_I.png' % str(spikedurat)
plotname_ISI     = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_ISI_vs_I.png'
'''# Should I have Vrest too? Would be nice


Vrest = np.zeros(Nmodels)
Vrest_rms = np.zeros(Nmodels)

cm = 1.0

for i in range(Nmodels):
    Nspikes      = []
    iamp_Nspikes = []
    avg_ISI      = []
    rms_ISI      = []
    iamp_ISI     = []
    iamp_AP_ampl = []
    avg_AP_ampl  = []
    rms_AP_ampl  = []
    iamp_AP_mins = []
    avg_AP_mins  = []
    rms_AP_mins  = []
    iamp_AP_halfwidth = []
    avg_AP_halfwidth  = []
    rms_AP_halfwidth  = []
    
    # Get names
    infolder = model_folders[i]+'Results/Soma%i/' % somasize
    infilename_Nspikes = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I.txt'
    infilename_APampl  = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Vmax_vs_I.txt'
    infilename_APmins  = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Vmin_vs_I.txt'
    infilename_APdhw   = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_sdurat%s_vs_I.txt' % str(spikedurat)
    infilename_ISI     = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_ISI_vs_I.txt'
    infilename_Vrest   = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Vrest.txt'
    
    # Open files
    infile_Nspikes = open(infilename_Nspikes,'r')
    infile_APampl  = open(infilename_APampl,'r')
    infile_APmins  = open(infilename_APmins,'r')
    infile_APdhw   = open(infilename_APdhw,'r')
    infile_ISI     = open(infilename_ISI,'r')
    infile_Vrest   = open(infilename_Vrest,'r')
    print('infilename_Vrest:',infilename_Vrest)
    
    lines_Nspikes = infile_Nspikes.readlines()
    lines_APampl  = infile_APampl.readlines()
    lines_APmins  = infile_APmins.readlines()
    lines_APdhw   = infile_APdhw.readlines()
    lines_ISI     = infile_ISI.readlines()
    line_Vrest    = infile_Vrest.readline()
    
    for line in lines_Nspikes:
        words = line.split()
        if len(words)>0:
            iamp_Nspikes.append(float(words[0]))
            Nspikes.append(float(words[1]))
    
    for line in lines_APampl:
        words = line.split()
        if len(words)>0:
            iamp_AP_ampl.append(float(words[0]))
            avg_AP_ampl.append(float(words[1]))
            rms_AP_ampl.append(float(words[2]))
    
    for line in lines_ISI:
        words = line.split()
        if len(words)>0:
            iamp_ISI.append(float(words[0]))
            avg_ISI.append(float(words[1]))
            rms_ISI.append(float(words[2]))
    
    for line in lines_APmins:
        words = line.split()
        if len(words)>0:
            iamp_AP_mins.append(float(words[0]))
            avg_AP_mins.append(float(words[1]))
            rms_AP_mins.append(float(words[2]))
    
    for line in lines_APdhw:
        words = line.split()
        if len(words)>0:
            iamp_AP_halfwidth.append(float(words[0]))
            avg_AP_halfwidth.append(float(words[1]))
            rms_AP_halfwidth.append(float(words[2]))
    
    ax1.plot(iamp_Nspikes,Nspikes,label=labels[i])
    ax2.errorbar(iamp_AP_ampl,avg_AP_ampl,yerr=rms_AP_ampl,label=labels[i])
    ax3.errorbar(iamp_AP_mins,avg_AP_mins,yerr=rms_AP_mins,label=labels[i])
    ax4.errorbar(iamp_AP_halfwidth,avg_AP_halfwidth,yerr=rms_AP_halfwidth,label=labels[i])
    ax5.errorbar(iamp_ISI,avg_ISI,yerr=rms_ISI,label=labels[i])
    
    words_Vrest  = line_Vrest.split()
    Vrest[i]     = float(words_Vrest[0])
    Vrest_rms[i] = float(words_Vrest[1])
    
    infile_Nspikes.close()
    infile_APampl.close()
    infile_APmins.close()
    infile_APdhw.close()
    infile_ISI.close()
    infile_Vrest.close()

ax1.set_xlabel('$I$ (nA)',fontsize=14)
ax1.set_ylabel('$f$ (Hz)',fontsize=14)
ax1.legend(loc='upper left')#,ncol=1)

ax2.set_xlabel('$I$ (nA)',fontsize=14)
ax2.set_ylabel('$V$ (mV)',fontsize=14)
ax2.legend(loc='upper left')#,ncol=1)

ax3.set_xlabel('$I$ (nA)',fontsize=14)
ax3.set_ylabel('$V$ (mV)',fontsize=14)
ax3.legend(loc='upper left')#,ncol=1)

ax4.set_xlabel('$I$ (nA)',fontsize=14)
ax4.set_ylabel('Duration (ms)',fontsize=14)
ax4.legend(loc='upper left')#,ncol=1)

ax5.set_xlabel('$I$ (nA)',fontsize=14)
ax5.set_ylabel('Interval (ms)',fontsize=14)
ax5.legend(loc='upper right')#,ncol=1)

br1 = np.arange(Nmodels)
ax6.set_xlabel('Model',fontsize=14)
ax6.set_ylabel('$V_{\mathregular{rest}}$',fontsize=14)
ax6.errorbar(br1,Vrest,yerr=Vrest_rms,label=labels[i])
plt.xticks(br1, [labels[0],labels[1],labels[2]])


print('Vrest:',Vrest)
print('Vrest_rms:',Vrest_rms)

fig.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.savefig(plotname)
plt.show()
