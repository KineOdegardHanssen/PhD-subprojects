import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

# change the default font family
plt.rcParams.update({'font.family':'Arial'})

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

cm         = 1.0
skiptime   = 500
spikedurat = -40
idur       = 1000 #100 # ms
idelay     = 10
v_init     = -86.8 # mV
Ra         = 100
somasize   = 10 # 15 # 
dtexp      = -7

fig = plt.figure(figsize=(18,8),dpi=300)

gs = gridspec.GridSpec(2, 8)

ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
ax3 = plt.subplot(gs[0, 4:6])
ax4 = plt.subplot(gs[0, 6:8])
ax5 = plt.subplot(gs[1, 0:2])
ax6 = plt.subplot(gs[1, 2:4])
ax7 = plt.subplot(gs[1, 4:6])
ax8 = plt.subplot(gs[1, 6:8])

#fig.suptitle(r'Properties',fontsize=20)

ax1.set_title(r'A',loc='left',pad=50,fontsize=18)
ax2.set_title(r'B',loc='left',pad=50,fontsize=18)
ax3.set_title(r'C',loc='left',pad=50,fontsize=18)
ax4.set_title(r'D',loc='left',pad=50,fontsize=18)
ax5.set_title(r'E',loc='left',pad=50,fontsize=18)
ax6.set_title(r'F',loc='left',pad=50,fontsize=18)
ax7.set_title(r'G',loc='left',pad=50,fontsize=18)
ax8.set_title(r'H',loc='left',pad=50,fontsize=18)

outfolder = 'Compare/Soma%i/' % somasize
plotname = outfolder+'compare_allmodels_fI.png'

#### Subplot 1 ####    
model_folders = ['','CaHVA_Allen/','CaHVA_Allen/bk/']
labels = ['base','CaHVA','CaHVA, bk']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
Nmodels = len(model_folders)

ax1.set_title(r'CaHVA (Allen); bk (Hjorth)',fontsize=12)
for i in range(Nmodels):
    Nspikes      = []
    iamp_Nspikes = []
    
    # Get names
    infolder = model_folders[i]+'Results/Soma%i/' % somasize
    infilename_Nspikes = infolder+infilenames[i]
    
    # Open files
    infile_Nspikes = open(infilename_Nspikes,'r')
    
    lines_Nspikes = infile_Nspikes.readlines()
    
    for line in lines_Nspikes:
        words = line.split()
        if len(words)>0:
            iamp_Nspikes.append(float(words[0]))
            Nspikes.append(float(words[1]))
    
    ax1.plot(iamp_Nspikes,Nspikes,label=labels[i])
    
    infile_Nspikes.close()

#### Subplot 2 ####
model_folders = ['','CaHVA_Allen/','CaHVA_Allen/BK_I_Zhang/']
labels = ['base','CaHVA','CaHVA, BK']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
Nmodels = len(model_folders)

cm = 1.0

ax2.set_title(r'CaHVA (Allen); BK (Zhang)',fontsize=12)
for i in range(Nmodels):
    Nspikes      = []
    iamp_Nspikes = []
    
    # Get names
    infolder = model_folders[i]+'Results/Soma%i/' % somasize
    infilename_Nspikes = infolder+infilenames[i]
    
    # Open files
    infile_Nspikes = open(infilename_Nspikes,'r')
    
    lines_Nspikes = infile_Nspikes.readlines()
    
    for line in lines_Nspikes:
        words = line.split()
        if len(words)>0:
            iamp_Nspikes.append(float(words[0]))
            Nspikes.append(float(words[1]))
    
    ax2.plot(iamp_Nspikes,Nspikes,label=labels[i])
    
    infile_Nspikes.close()

#### Subplot 3 ####
model_folders = ['','CaHVA_Allen/','CaHVA_Allen/BK_AitOuares/']
labels = ['base','CaHVA','CaHVA, BK']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
Nmodels = len(model_folders)

cm = 1.0

ax3.set_title(r'CaHVA (Allen); BK (Ait Ouares)',fontsize=12)
for i in range(Nmodels):
    Nspikes      = []
    iamp_Nspikes = []
    
    # Get names
    infolder = model_folders[i]+'Results/Soma%i/' % somasize
    infilename_Nspikes = infolder+infilenames[i]
    
    # Open files
    infile_Nspikes = open(infilename_Nspikes,'r')
    
    lines_Nspikes = infile_Nspikes.readlines()
    
    for line in lines_Nspikes:
        words = line.split()
        if len(words)>0:
            iamp_Nspikes.append(float(words[0]))
            Nspikes.append(float(words[1]))
    
    ax3.plot(iamp_Nspikes,Nspikes,label=labels[i])
    
    infile_Nspikes.close()


#### Subplot 4 ####
model_folders = ['','CaHVA_Allen/','CaHVA_Allen/SK_AitOuares/']
labels = ['base','CaHVA','CaHVA, SK']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
Nmodels = len(model_folders)

cm = 1.0

ax4.set_title(r'CaHVA (Allen); SK (Ait Ouares)',fontsize=12)
for i in range(Nmodels):
    Nspikes      = []
    iamp_Nspikes = []
    
    # Get names
    infolder = model_folders[i]+'Results/Soma%i/' % somasize
    infilename_Nspikes = infolder+infilenames[i]
    
    # Open files
    infile_Nspikes = open(infilename_Nspikes,'r')
    
    lines_Nspikes = infile_Nspikes.readlines()
    
    for line in lines_Nspikes:
        words = line.split()
        if len(words)>0:
            iamp_Nspikes.append(float(words[0]))
            Nspikes.append(float(words[1]))
    
    ax4.plot(iamp_Nspikes,Nspikes,label=labels[i])
    
    infile_Nspikes.close()

#### Subplot 5 ####
model_folders = ['','CaHVA_Allen/','CaHVA_Allen/SK_Allen/']
labels = ['base','CaHVA','CaHVA, SK']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gSK1.0p_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt']
Nmodels = len(model_folders)

cm = 1.0

ax5.set_title(r'CaHVA (Allen); SK (Allen)',fontsize=12)
for i in range(Nmodels):
    Nspikes      = []
    iamp_Nspikes = []
    
    # Get names
    infolder = model_folders[i]+'Results/Soma%i/' % somasize
    infilename_Nspikes = infolder+infilenames[i]
    
    # Open files
    infile_Nspikes = open(infilename_Nspikes,'r')
    
    lines_Nspikes = infile_Nspikes.readlines()
    
    for line in lines_Nspikes:
        words = line.split()
        if len(words)>0:
            iamp_Nspikes.append(float(words[0]))
            Nspikes.append(float(words[1]))
    
    ax5.plot(iamp_Nspikes,Nspikes,label=labels[i])
    
    infile_Nspikes.close()


#### Subplot 6 ####
model_folders = ['','canin_Konstantoudaki/','canin_Konstantoudaki/bk/']
labels = ['base','canin','canin, bk']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gcanin1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
Nmodels = len(model_folders)

cm = 1.0

ax6.set_title(r'canin (Konstantoudaki); bk (Hjorth)',fontsize=12)
for i in range(Nmodels):
    Nspikes      = []
    iamp_Nspikes = []
    
    # Get names
    infolder = model_folders[i]+'Results/Soma%i/' % somasize
    infilename_Nspikes = infolder+infilenames[i]
    
    # Open files
    infile_Nspikes = open(infilename_Nspikes,'r')
    
    lines_Nspikes = infile_Nspikes.readlines()
    
    for line in lines_Nspikes:
        words = line.split()
        if len(words)>0:
            iamp_Nspikes.append(float(words[0]))
            Nspikes.append(float(words[1]))
    
    ax6.plot(iamp_Nspikes,Nspikes,label=labels[i])
    
    infile_Nspikes.close()


#### Subplot 7 ####
model_folders = ['','canin_Konstantoudaki/','canin_Konstantoudaki/BK_I_Zhang/']
labels = ['base','canin','canin, BK']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gcanin1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
Nmodels = len(model_folders)

cm = 1.0

ax7.set_title(r'canin (Konstantoudaki); BK (Zhang)',fontsize=12)
for i in range(Nmodels):
    Nspikes      = []
    iamp_Nspikes = []
    
    # Get names
    infolder = model_folders[i]+'Results/Soma%i/' % somasize
    infilename_Nspikes = infolder+infilenames[i]
    
    # Open files
    infile_Nspikes = open(infilename_Nspikes,'r')
    
    lines_Nspikes = infile_Nspikes.readlines()
    
    for line in lines_Nspikes:
        words = line.split()
        if len(words)>0:
            iamp_Nspikes.append(float(words[0]))
            Nspikes.append(float(words[1]))
    
    ax7.plot(iamp_Nspikes,Nspikes,label=labels[i])
    
    infile_Nspikes.close()


#### Subplot 8 ####
model_folders = ['','canin_Konstantoudaki/','canin_Konstantoudaki/BK_AitOuares/']
labels = ['base','canin','canin, BK']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gcanin1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
Nmodels = len(model_folders)

cm = 1.0

ax8.set_title(r'canin (Konstantoudaki); BK (Ait Ouares)',loc='right',fontsize=12)
for i in range(Nmodels):
    Nspikes      = []
    iamp_Nspikes = []
    
    # Get names
    infolder = model_folders[i]+'Results/Soma%i/' % somasize
    infilename_Nspikes = infolder+infilenames[i]
    
    # Open files
    infile_Nspikes = open(infilename_Nspikes,'r')
    
    lines_Nspikes = infile_Nspikes.readlines()
    
    for line in lines_Nspikes:
        words = line.split()
        if len(words)>0:
            iamp_Nspikes.append(float(words[0]))
            Nspikes.append(float(words[1]))
    
    ax8.plot(iamp_Nspikes,Nspikes,label=labels[i])
    
    infile_Nspikes.close()



ax1.set_xlabel('$I$ (nA)',fontsize=14)
ax1.set_ylabel('$f$ (Hz)',fontsize=14)
ax1.legend(loc='upper left')#,ncol=1)

ax2.set_xlabel('$I$ (nA)',fontsize=14)
ax2.set_ylabel('$f$ (Hz)',fontsize=14)
ax2.legend(loc='upper left')#,ncol=1)

ax3.set_xlabel('$I$ (nA)',fontsize=14)
ax3.set_ylabel('$f$ (Hz)',fontsize=14)
ax3.legend(loc='upper left')#,ncol=1)

ax4.set_xlabel('$I$ (nA)',fontsize=14)
ax4.set_ylabel('$f$ (Hz)',fontsize=14)
ax4.legend(loc='upper left')#,ncol=1)

ax5.set_xlabel('$I$ (nA)',fontsize=14)
ax5.set_ylabel('$f$ (Hz)',fontsize=14)
ax5.legend(loc='upper left')#,ncol=1)

ax6.set_xlabel('$I$ (nA)',fontsize=14)
ax6.set_ylabel('$f$ (Hz)',fontsize=14)
ax6.legend(loc='upper left')#,ncol=1)

ax7.set_xlabel('$I$ (nA)',fontsize=14)
ax7.set_ylabel('$f$ (Hz)',fontsize=14)
ax7.legend(loc='upper left')#,ncol=1)

ax8.set_xlabel('$I$ (nA)',fontsize=14)
ax8.set_ylabel('$f$ (Hz)',fontsize=14)
ax8.legend(loc='upper left')#,ncol=1)


fig.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.savefig(plotname)
plt.show()
