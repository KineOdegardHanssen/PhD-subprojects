import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

# change the default font family
plt.rcParams.update({'font.family':'Arial'})

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

cm         = 1.0
spikedurat = -40
idur       = 1000 #100 # ms
idelay     = 10
v_init     = -86 # mV
Ra         = 100
somasize   = 10 # 15 # 
dtexp      = -7
skiptime   = 500

fig = plt.figure(figsize=(6,5),dpi=300)

    
model_folders = ['','CaHVA_Allen/Shift_gCaHVA_V/Vshift_0/','CaHVA_Allen/SK_Allen/Shift_gCaHVA_V/Vshift_0/']
labels = ['base','CaHVA',r'0.2$\cdot$CaHVA, SK']
Nmodels = len(model_folders)

outfolder = 'Compare/Soma%i/' % somasize
plotname = outfolder+'fI_naf_kaf_CaHVA_SK_Allen_skiptime'+str(skiptime)+'_commonprops.png'


# Namestrings
namestring1 = '_gCaHVA1.0p_'   
namestring2 = '_gSK1.0p_gCaHVA0.2p_'   

# Nspikes
filename_Nspikes_base     = 'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'
filename_Nspikes_CaHVA    = 'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring1+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'
filename_Nspikes_CaHVA_SK = 'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring2+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'
# Gathering
filenames_Nspikes = [filename_Nspikes_base,filename_Nspikes_CaHVA,filename_Nspikes_CaHVA_SK]


Vrest = np.zeros(Nmodels)
Vrest_rms = np.zeros(Nmodels)

cm = 1.0

for i in range(Nmodels):
    Nspikes      = []
    iamp_Nspikes = []
    
    # Get names
    infolder = model_folders[i]+'Results/Soma%i/' % somasize
    infilename_Nspikes = infolder+filenames_Nspikes[i]
    
    # Open files
    infile_Nspikes = open(infilename_Nspikes,'r')
    
    lines_Nspikes = infile_Nspikes.readlines()
    
    for line in lines_Nspikes:
        words = line.split()
        if len(words)>0:
            iamp_Nspikes.append(float(words[0]))
            Nspikes.append(float(words[1]))
    
    plt.plot(iamp_Nspikes,Nspikes,label=labels[i])
    
    infile_Nspikes.close()

plt.xlabel('$I$ (nA)',fontsize=14)
plt.ylabel('$f$ (Hz)',fontsize=14)
plt.legend(loc='upper left')#,ncol=1)

fig.tight_layout()

plt.savefig(plotname)
plt.show()
