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

fig = plt.figure(figsize=(22,12),dpi=300)
        
gs = gridspec.GridSpec(3, 11)

ax1  = plt.subplot(gs[0, 0:2])
ax2  = plt.subplot(gs[0, 2:4])
ax3  = plt.subplot(gs[0, 4:6])
ax4  = plt.subplot(gs[0, 6:8])
ax5  = plt.subplot(gs[0, 8:10])
axr1 = plt.subplot(gs[0, 10:11])
ax6  = plt.subplot(gs[1, 0:2])
ax7  = plt.subplot(gs[1, 2:4])
ax8  = plt.subplot(gs[1, 4:6])
ax9  = plt.subplot(gs[1, 6:8])
ax10 = plt.subplot(gs[1, 8:10])
axr2 = plt.subplot(gs[1, 10:11])
ax11 = plt.subplot(gs[2, 0:2])
ax12 = plt.subplot(gs[2, 2:4])
ax13 = plt.subplot(gs[2, 4:6])
ax14 = plt.subplot(gs[2, 6:8])
ax15 = plt.subplot(gs[2, 8:10])
axr3 = plt.subplot(gs[2, 10:11])

ax1.set_title(r'A',y=0,pad=-14,fontsize=18)
ax2.set_title(r'B',y=0,pad=-14,fontsize=18)
ax3.set_title(r'C',y=0,pad=-14,fontsize=18)
ax4.set_title(r'D',y=0,pad=-14,fontsize=18)
ax5.set_title(r'E',y=0,pad=-14,fontsize=18)
ax6.set_title(r'F',y=0,pad=-14,fontsize=18)
ax7.set_title(r'G',y=0.9,pad=-1,fontsize=18)
ax8.set_title(r'H',y=0.9,pad=-1,fontsize=18)
ax9.set_title(r'I',y=0,pad=-14,fontsize=18)
ax10.set_title(r'J',y=0,pad=-14,fontsize=18)
ax11.set_title(r'K',y=0,pad=-14,fontsize=18)
ax12.set_title(r'L',y=0,pad=-14,fontsize=18)
ax13.set_title(r'M',y=0,pad=-14,fontsize=18)
ax14.set_title(r'N',y=0,pad=-14,fontsize=18)
ax15.set_title(r'O',y=0,pad=-14,fontsize=18)

outfolder = 'Compare/Soma%i/' % somasize
plotname = outfolder+'compare_allmodels_fI_all_mainmodelgs_varyghva_getgSK.png'

### CaHVA_Allen, bk:
gcahvas = [1.0,1.25,1.5,2.0,2.5]
gbk     = 2 # Try this out 
Ng = len(gcahvas)
colors = ['tab:blue']
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))
colors.append('tab:orange')
#### Subplot 1 ####    
model_folders = ['']
labels = ['base']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']

for gcahva in gcahvas:
    labels.append(r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahva))
    model_folders.append('CaHVA_Allen/')
    infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+ '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
model_folders.append('CaHVA_Allen/bk/gSK0.0028/')
labels.append(r'1.0$\bar{g}_\mathregular{CaHVA}$, %s$\bar{g}_\mathregular{bk}$' % str(gbk))
Nmodels = len(model_folders)

ax1.set_title(r'CaHVA (Allen); bk (Hjorth)',loc='right',fontsize=12)
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
    
    if i==0 or i==Nmodels-1:
        ax1.plot(iamp_Nspikes,Nspikes,label=labels[i],color=colors[i])
    else:
        ax1.plot(iamp_Nspikes,Nspikes,color=colors[i])

    infile_Nspikes.close()

### CaHVA_Allen, BK_I_Zhang:
gcahvas = [1.0,1.25,1.5,2.0,2.5]
gbk     = 1000 
Ng = len(gcahvas)
colors = ['tab:blue']
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))
colors.append('tab:orange')
model_folders = ['']
labels = ['base']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']

for gcahva in gcahvas:   
    #### Subplot 2 ####
    model_folders.append('CaHVA_Allen/')
    labels.append(r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahva))
    infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
labels.append(r'1.0$\bar{g}_\mathregular{CaHVA}$, %s$\bar{g}_\mathregular{BK}$' % str(gbk))
model_folders.append('CaHVA_Allen/BK_I_Zhang/gSK0.0028/')
infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
Nmodels = len(model_folders)
    
ax2.set_title(r'CaHVA (Allen); BK (Zhang)',loc='right',fontsize=12)
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
    if i==0 or i==Nmodels-1:
        ax2.plot(iamp_Nspikes,Nspikes,label=labels[i],color=colors[i])
    else:
        ax2.plot(iamp_Nspikes,Nspikes,color=colors[i])
    
    infile_Nspikes.close()

### CaHVA_Allen, BK_AitOuares:
gcahvas = [1.0,1.25,1.5,2.0,2.5] # Change these every turn # !!!!!!!!!!!!!!!!!!
gbk     = 2000 
Ng = len(gcahvas)
colors = ['tab:blue']
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))
colors.append('tab:orange')
labels = ['base']
model_folders = ['']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
for gcahva in gcahvas: 
    #### Subplot 3 ####
    model_folders.append('CaHVA_Allen/')
    labels.append(r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahva))
    infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
labels.append(r'1.0$\bar{g}_\mathregular{CaHVA}$, %s$\bar{g}_\mathregular{BK}$' % str(gbk))
model_folders.append('CaHVA_Allen/BK_AitOuares/gSK0.0028/')
infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
Nmodels = len(model_folders)
    
ax3.set_title(r'CaHVA (Allen); BK (Ait Ouares)',loc='right',fontsize=12)
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
    
    if i==0 or i==Nmodels-1:
        ax3.plot(iamp_Nspikes,Nspikes,label=labels[i],color=colors[i])
    else:
        ax3.plot(iamp_Nspikes,Nspikes,color=colors[i])
    
    infile_Nspikes.close()

### CaHVA_Allen, SK_AitOuares:
gcahvas = [1.0,1.25,1.5,2.0,2.5] # Change these every turn # !!!!!!!!!!!!!!!!!!
gbk     = 10 # Try this out  # !!!!!!!!!!!!!!!!!!
Ng = len(gcahvas)
colors = ['tab:blue']
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))
colors.append('tab:orange')
model_folders = ['']
labels = ['base']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
for gcahva in gcahvas:         
    #### Subplot 4 ####
    model_folders.append('CaHVA_Allen/')
    labels.append(r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahva))
    infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm) + '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
model_folders.append('CaHVA_Allen/SK_AitOuares/gSK0.0028/')
labels.append(r'1.0$\bar{g}_\mathregular{CaHVA}$, %s$\bar{g}_\mathregular{SK}$' % str(gbk))
infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gsk'+str(gbk)+'p' + '_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
Nmodels = len(model_folders)
        
ax4.set_title(r'CaHVA (Allen); SK (Ait Ouares)',loc='right',fontsize=12)
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
    
    if i==0 or i==Nmodels-1:
        ax4.plot(iamp_Nspikes,Nspikes,label=labels[i],color=colors[i])
    else:
        ax4.plot(iamp_Nspikes,Nspikes,color=colors[i])
    
    infile_Nspikes.close()

### CaHVA_Allen, SK_Allen:
gcahvas = [1.0,1.25,1.5,2.0,2.5] # Change these every turn # !!!!!!!!!!!!!!!!!!
gbk     = 10 # Try this out  # !!!!!!!!!!!!!!!!!!
Ng = len(gcahvas)
colors = ['tab:blue']
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))
colors.append('tab:orange')
labels = ['base']
model_folders = ['']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
for gcahva in gcahvas:
    #### Subplot 5 ####
    model_folders.append('CaHVA_Allen/')
    labels.append(r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahva))
    infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
model_folders.append('CaHVA_Allen/SK_Allen/')
labels.append(r'1.0$\bar{g}_\mathregular{CaHVA}$, %s$\bar{g}_\mathregular{SK}$' % str(gbk))
infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gSK'+str(gbk)+'p' + '_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
Nmodels = len(model_folders)
        
ax5.set_title(r'CaHVA (Allen); SK (Allen)',loc='right',fontsize=12)
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
    
    if i==0 or i==Nmodels-1:
        ax5.plot(iamp_Nspikes,Nspikes,label=labels[i],color=colors[i])
    else:
        ax5.plot(iamp_Nspikes,Nspikes,color=colors[i])
    
    infile_Nspikes.close()

axr1.set_frame_on(False)
axr1.get_xaxis().set_visible(False)
axr1.get_yaxis().set_visible(False)
for i in range(1,Nmodels-1):
    axr1.plot('-',color=colors[i],label=labels[i])
axr1.legend(loc='center')
### canin, bk:
gcahvas = [1.0,2.0,3.0,5.0,7.0,10] # Change these every turn # !!!!!!!!!!!!!!!!!!
gbk     = 2.0 # Try this out  # !!!!!!!!!!!!!!!!!!
Ng = len(gcahvas)
colors = ['tab:blue']
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))
colors.append('tab:orange')
labels = ['base']
model_folders = ['']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
for gcahva in gcahvas:       
    #### Subplot 6 ####
    model_folders.append('canin_Konstantoudaki/gSK0.0028/')
    labels.append(r'%s$\bar{g}_\mathregular{canin}$' % str(gcahva))
    infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
model_folders.append('canin_Konstantoudaki/bk/gSK0.0028/')
labels.append(r'1.0$\bar{g}_\mathregular{canin}$, %s$\bar{g}_\mathregular{bk}$' % str(gbk))
infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
Nmodels = len(model_folders)
        
ax6.set_title(r'canin (Konstantoudaki); bk (Hjorth)',loc='right',fontsize=12)
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
    
    if i==0 or i==Nmodels-1:
        ax6.plot(iamp_Nspikes,Nspikes,label=labels[i],color=colors[i])
    else:
        ax6.plot(iamp_Nspikes,Nspikes,color=colors[i])
    
    infile_Nspikes.close()
            
### canin, BK_Zhang: ##################### yooooooo!
gcahvas = [1.0,2.0,3.0,5.0,7.0,10]
gbk     = 100000
Ng = len(gcahvas)
colors = ['tab:blue']
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))
colors.append('tab:orange')
labels = ['base']
model_folders = ['']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
for gcahva in gcahvas:
    #### Subplot 7 ####
    model_folders.append('canin_Konstantoudaki/gSK0.0028/')
    labels.append(r'%s$\bar{g}_\mathregular{canin}$' % str(gcahva))
    infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
model_folders.append('canin_Konstantoudaki/BK_I_Zhang/gSK0.0028/')
labels.append(r'1.0$\bar{g}_\mathregular{canin}$, %s$\bar{g}_\mathregular{BK}$' % str(gbk))
infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
Nmodels = len(model_folders)

ax7.set_title(r'canin (Konstantoudaki); BK (Zhang)',loc='right',fontsize=12)
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
    
    if i==0 or i==Nmodels-1:
        ax7.plot(iamp_Nspikes,Nspikes,label=labels[i],color=colors[i])
    else:
        ax7.plot(iamp_Nspikes,Nspikes,color=colors[i])
    
    infile_Nspikes.close()
     
### canin, BK_AO: ##################### yooooooo!
gcahvas = [1.0,2.0,3.0,5.0,7.0,10] # Change these every turn # !!!!!!!!!!!!!!!!!!
gbk     = 10000 # Try this out  # !!!!!!!!!!!!!!!!!!
Ng = len(gcahvas)
colors = ['tab:blue']#[(1.0,1.0,1.0,1.0)]#
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))
colors.append('tab:orange')
model_folders = ['']
labels = ['base']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
for gcahva in gcahvas:        
    #### Subplot 8 ####
    model_folders.append('canin_Konstantoudaki/gSK0.0028/')
    labels.append(r'%s$\bar{g}_\mathregular{canin}$' % str(gcahva))
    infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
model_folders.append('canin_Konstantoudaki/BK_AitOuares/gSK0.0028/')
labels.append(r'1.0$\bar{g}_\mathregular{canin}$, %s$\bar{g}_\mathregular{BK}$' % str(gbk))
infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
Nmodels = len(model_folders)
    
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
    
    if i==0 or i==Nmodels-1:
        ax8.plot(iamp_Nspikes,Nspikes,label=labels[i],color=colors[i])
    else:
        ax8.plot(iamp_Nspikes,Nspikes,color=colors[i])
    
    infile_Nspikes.close()

### canin, SK_AO: ##################### yooooooo!
gcahvas = [1.0,2.0,3.0,5.0,7.0,10] # Change these every turn # !!!!!!!!!!!!!!!!!!
gbk     = 10 # Try this out  # !!!!!!!!!!!!!!!!!!
Ng = len(gcahvas)
colors = ['tab:blue']
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))
colors.append('tab:orange')
labels = ['base']
model_folders = ['']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
for gcahva in gcahvas:        
    #### Subplot 9 ####
    model_folders.append('canin_Konstantoudaki/gSK0.0028/')
    labels.append(r'%s$\bar{g}_\mathregular{canin}$' % str(gcahva))
    infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
model_folders.append('canin_Konstantoudaki/SK_AitOuares/gSK0.0028/')
labels.append(r'1.0$\bar{g}_\mathregular{canin}$, %s$\bar{g}_\mathregular{SK}$' % str(gbk))
infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gsk'+str(gbk)+'p' + '_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
Nmodels = len(model_folders)
    
ax9.set_title(r'canin (Konstantoudaki); SK (Ait Ouares)',loc='right',fontsize=12) #'_gsk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p_
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
    
    if i==0 or i==Nmodels-1:
        ax9.plot(iamp_Nspikes,Nspikes,label=labels[i],color=colors[i])
    else:
        ax9.plot(iamp_Nspikes,Nspikes,color=colors[i])
    
    infile_Nspikes.close()
    
### canin, SK_Allen: ##################### yooooooo!
gcahvas = [1.0,2.0,3.0,5.0,7.0,10] # Change these every turn # !!!!!!!!!!!!!!!!!!
gbk     = 10 # Try this out  # !!!!!!!!!!!!!!!!!!
Ng = len(gcahvas)
colors = ['tab:blue']
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))
colors.append('tab:orange')
labels = ['base']
model_folders = ['']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
for gcahva in gcahvas:           
    #### Subplot 10 ####
    model_folders.append('canin_Konstantoudaki/gSK0.0028/')
    labels.append(r'%s$\bar{g}_\mathregular{canin}$' % str(gcahva))
    infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
model_folders.append('canin_Konstantoudaki/SK_Allen/gSK0.0028/')
labels.append(r'1.0$\bar{g}_\mathregular{canin}$, %s$\bar{g}_\mathregular{SK}$' % str(gbk))
infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gsk'+str(gbk)+'p' + '_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
Nmodels = len(model_folders)
    
ax10.set_title(r'canin (Konstantoudaki); SK (Allen)',loc='right',fontsize=12)
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
    
    if i==0 or i==Nmodels-1:
        ax10.plot(iamp_Nspikes,Nspikes,label=labels[i],color=colors[i])
    else:
        ax10.plot(iamp_Nspikes,Nspikes,color=colors[i])
    
    infile_Nspikes.close() 

axr2.set_frame_on(False)
axr2.get_xaxis().set_visible(False)
axr2.get_yaxis().set_visible(False)
for i in range(1,Nmodels-1):
    axr2.plot('-',color=colors[i],label=labels[i])
axr2.legend(loc='center')

### caq, bk: ##################### yooooooo!
gcahvas = [1.0,2.0,3.0,5.0,7.0] # Change these every turn # !!!!!!!!!!!!!!!!!!
gbk     = 5.0 # Try this out  # !!!!!!!!!!!!!!!!!!
Ng = len(gcahvas)
colors = ['tab:blue']
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))
colors.append('tab:orange')
labels = ['base']
model_folders = ['']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
for gcahva in gcahvas:   
    #### Subplot 11 ####
    model_folders.append('caq/gSK0.0028/')
    labels.append(r'%s$\bar{g}_\mathregular{caq}$' % str(gcahva))
    infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
model_folders.append('caq/bk/gSK0.0028/')
labels.append(r'1.0$\bar{g}_\mathregular{caq}$, %s$\bar{g}_\mathregular{bk}$' % str(gbk))
infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
Nmodels = len(model_folders)
        
ax11.set_title(r'caq (Hjorth); bk (Hjorth)',loc='right',fontsize=12)
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
    
    if i==0 or i==Nmodels-1:
        ax11.plot(iamp_Nspikes,Nspikes,label=labels[i],color=colors[i])
    else:
        ax11.plot(iamp_Nspikes,Nspikes,color=colors[i])
    
    infile_Nspikes.close()

### caq, BK_Zhang: ##################### yooooooo!
gcahvas = [1.0,2.0,3.0,5.0,7.0] # Change these every turn # !!!!!!!!!!!!!!!!!!
gbk     = 10 # Try this out  # !!!!!!!!!!!!!!!!!!
Ng = len(gcahvas)
colors = ['tab:blue']
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))
colors.append('tab:orange')
labels = ['base']
model_folders = ['']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
for gcahva in gcahvas:            
    #### Subplot 12 ####
    model_folders.append('caq/gSK0.0028/')
    labels.append(r'%s$\bar{g}_\mathregular{caq}$' % str(gcahva))
    infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
model_folders.append('caq/BK_I_Zhang/gSK0.0028/')
labels.append(r'1.0$\bar{g}_\mathregular{caq}$, %s$\bar{g}_\mathregular{BK}$' % str(gbk))
infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
Nmodels = len(model_folders)
        
ax12.set_title(r'caq (Hjorth); BK (Zhang)',loc='right',fontsize=12)
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
    
    if i==0 or i==Nmodels-1:
        ax12.plot(iamp_Nspikes,Nspikes,label=labels[i],color=colors[i])
    else:
        ax12.plot(iamp_Nspikes,Nspikes,color=colors[i])
    
    infile_Nspikes.close()
        
### caq, BK_A0: ##################### yooooooo!
gcahvas = [1.0,2.0,3.0,5.0,7.0] # Change these every turn # !!!!!!!!!!!!!!!!!!
gbk     = 5000 # Try this out  # !!!!!!!!!!!!!!!!!!
Ng = len(gcahvas)
colors = ['tab:blue']
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))
colors.append('tab:orange')
labels = ['base']
model_folders = ['']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
for gcahva in gcahvas:        
    #### Subplot 13 ####
    model_folders.append('caq/gSK0.0028/')
    labels.append(r'%s$\bar{g}_\mathregular{caq}$' % str(gcahva))
    infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
model_folders.append('caq/BK_AitOuares/gSK0.0028/')
labels.append(r'1.0$\bar{g}_\mathregular{caq}$, %s$\bar{g}_\mathregular{BK}$' % str(gbk))
infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
Nmodels = len(model_folders)
    
ax13.set_title(r'caq (Hjorth); BK (Ait Ouares)',loc='right',fontsize=12)
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
    
    if i==0 or i==Nmodels-1:
        ax13.plot(iamp_Nspikes,Nspikes,label=labels[i],color=colors[i])
    else:
        ax13.plot(iamp_Nspikes,Nspikes,color=colors[i])
        
    infile_Nspikes.close()

### caq, SK_AO: ##################### yooooooo!
gcahvas = [1.0,2.0,3.0,5.0,7.0] # Change these every turn # !!!!!!!!!!!!!!!!!!
gbk     = 15 # Try this out  # !!!!!!!!!!!!!!!!!!
Ng = len(gcahvas)
colors = ['tab:blue']
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))
colors.append('tab:orange')
labels = ['base']
model_folders = ['']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
for gcahva in gcahvas:        
    #### Subplot 14 ####
    model_folders.append('caq/gSK0.0028/')
    labels.append(r'%s$\bar{g}_\mathregular{caq}$' % str(gcahva))
    infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
model_folders.append('caq/SK_AitOuares/gSK0.0028/')
labels.append(r'1.0$\bar{g}_\mathregular{caq}$, %s$\bar{g}_\mathregular{SK}$' % str(gbk))
infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
Nmodels = len(model_folders)
    
ax14.set_title(r'caq (Hjorth); SK (Ait Ouares)',loc='right',fontsize=12) #'_gsk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p_
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
    
    if i==0 or i==Nmodels-1:
        ax14.plot(iamp_Nspikes,Nspikes,label=labels[i],color=colors[i])
    else:
        ax14.plot(iamp_Nspikes,Nspikes,color=colors[i])
    
    infile_Nspikes.close()
        
### caq, SK_Allen: ##################### yooooooo!
gcahvas = [1.0,2.0,3.0,5.0,7.0] # Change these every turn # !!!!!!!!!!!!!!!!!!
gbk     = 5.0 # Try this out  # !!!!!!!!!!!!!!!!!!
Ng = len(gcahvas)
colors = ['tab:blue']
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))
colors.append('tab:orange')
labels = ['base']
model_folders = ['']
infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt']
for gcahva in gcahvas:            
    #### Subplot 15 ####
    model_folders.append('caq/gSK0.0028/')
    labels.append(r'%s$\bar{g}_\mathregular{caq}$' % str(gcahva))
    infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
model_folders.append('caq/SK_Allen/gSK0.0028/')
labels.append(r'1.0$\bar{g}_\mathregular{caq}$, %s$\bar{g}_\mathregular{SK}$' % str(gbk))
infilenames.append('somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt')
Nmodels = len(model_folders)
   
ax15.set_title(r'caq (Hjorth); SK (Allen)',loc='right',fontsize=12)
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
    
    if i==0 or i==Nmodels-1:
        ax15.plot(iamp_Nspikes,Nspikes,label=labels[i],color=colors[i])
    else:
        ax15.plot(iamp_Nspikes,Nspikes,color=colors[i])
    
    infile_Nspikes.close()

axr3.set_frame_on(False)
axr3.get_xaxis().set_visible(False)
axr3.get_yaxis().set_visible(False)
for i in range(1,Nmodels-1):
    axr3.plot('-',color=colors[i],label=labels[i])
axr3.legend(loc='center')

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
ax7.legend()#(loc='center right')#,ncol=1)

ax8.set_xlabel('$I$ (nA)',fontsize=14)
ax8.set_ylabel('$f$ (Hz)',fontsize=14)
ax8.legend()#(loc='lower right')#,ncol=1)

ax9.set_xlabel('$I$ (nA)',fontsize=14)
ax9.set_ylabel('$f$ (Hz)',fontsize=14)
ax9.legend(loc='upper left')#,ncol=1)

ax10.set_xlabel('$I$ (nA)',fontsize=14)
ax10.set_ylabel('$f$ (Hz)',fontsize=14)
ax10.legend(loc='upper left')#,ncol=1)

ax11.set_xlabel('$I$ (nA)',fontsize=14)
ax11.set_ylabel('$f$ (Hz)',fontsize=14)
ax11.legend(loc='upper left')#,ncol=1)

ax12.set_xlabel('$I$ (nA)',fontsize=14)
ax12.set_ylabel('$f$ (Hz)',fontsize=14)
ax12.legend(loc='upper left')#,ncol=1)

ax13.set_xlabel('$I$ (nA)',fontsize=14)
ax13.set_ylabel('$f$ (Hz)',fontsize=14)
ax13.legend(loc='upper left')#,ncol=1)

ax14.set_xlabel('$I$ (nA)',fontsize=14)
ax14.set_ylabel('$f$ (Hz)',fontsize=14)
ax14.legend(loc='upper left')#,ncol=1)

ax15.set_xlabel('$I$ (nA)',fontsize=14)
ax15.set_ylabel('$f$ (Hz)',fontsize=14)
ax15.legend(loc='upper left')#,ncol=1)
        
fig.tight_layout()#(rect=[0, 0.03, 1, 0.95])
        
plt.savefig(plotname)
#plt.show()
