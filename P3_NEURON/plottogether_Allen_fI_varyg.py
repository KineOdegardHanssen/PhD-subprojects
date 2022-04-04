import numpy as np
import matplotlib.pyplot as plt

# change the default font family
plt.rcParams.update({'font.family':'Arial'})

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend', fontsize=13)

mylinewidth = 2

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

namestringfirst = ''
varygHVA = [0.1,0.3,0.5,0.7,0.9,1.0,1.1,1.5,2.0,3.0]
namestringfirstHVA = namestringfirst + '_gCaHVA'
plotlabelgHVA = r'$\bar{g}_{\mathregular{CaHVA}}$'
varygNaV = [0.8,1.0,2.0,3.0]
namestringfirstNaV = namestringfirst + '_gNaV'
plotlabelgNaV = r'$\bar{g}_{\mathregular{NaV}}$'
varygKv31 = [0.1,0.3,0.5,1.0,2.0,2.5,3.0]    
namestringfirstKv31 = namestringfirst + '_gKv31'
plotlabelgKv31 = r'$\bar{g}_{\mathregular{Kv31}}$'
varygKT = [0.3,0.5,1.0,1.5,2.0]    
namestringfirstKT = namestringfirst + '_gKT'
plotlabelgKT = r'$\bar{g}_{\mathregular{KT}}$'


plotstringHVA  = '_vary_CaHVA'
plotstringNaV  = '_vary_NaV'
plotstringKv31 = '_vary_Kv31'
plotstringKT   = '_vary_KT'

i_master = [0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5]
cm = 1.0 # Keep this a scalar for now
    
NI  = len(i_master)
NgHVA = len(varygHVA)

i_master_everywhere_all_gHVA   = []
Nspikes_everywhere_all_gHVA    = []
I_Nspikes_everywhere_all_gHVA  = []

i_master_sprx_all_gHVA         = []
Nspikes_sprx_all_gHVA          = []
I_Nspikes_sprx_all_gHVA        = []

for testmodel in testmodels:
    i_master_everywhere_gHVA   = []
    Nspikes_everywhere_gHVA    = []
    I_Nspikes_everywhere_gHVA  = []
    
    #i_master_sprx_Na         = []
    #Nspikes_sprx_Na          = []
    #I_Nspikes_sprx_Na        = []

    infolder      = 'figures/%i/' % (testmodel)
    vrestfolder   = infolder 
    for g in varygHVA:
        Nspikes   = []
        I_Nspikes = []
        
        namestringHVA = namestringfirstHVA+str(g)+'p'
        
        # Set names
        infilename_Nspikes = infolder+'%s_%i_cmfall'%(namestringHVA,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
        print('infilename_Nspikes:',infilename_Nspikes)
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
        
        Nspikes_everywhere_gHVA.append(Nspikes)
        I_Nspikes_everywhere_gHVA.append(I_Nspikes)
    
    
        ############# SOMAPROX ##################################
        #Nspikes   = []
        #I_Nspikes = []
        
        # Set names
        #infilename_Nspikes = infolder+'%s_%i_cmfsprx'%(namestringHVA,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
        # Read files
        #infile_Nspikes = open(infilename_Nspikes,'r')
    
        #lines_Nspikes = infile_Nspikes.readlines()
        #Nlines_Nspikes = len(lines_Nspikes)
    
        #for i in range(Nlines_Nspikes):
        #    words_Nspikes = lines_Nspikes[i].split()
        #    if len(words_Nspikes)>0:
        #        I_Nspikes.append(float(words_Nspikes[0]))
        #        Nspikes.append(float(words_Nspikes[1]))

        #infile_Nspikes.close()
        
        #Nspikes_sprx_Na.append(Nspikes)
        #I_Nspikes_sprx_Na.append(I_Nspikes)
    print('-------------------------------------')
    Nspikes_everywhere_all_gHVA.append(Nspikes_everywhere_gHVA)
    I_Nspikes_everywhere_all_gHVA.append(I_Nspikes_everywhere_gHVA)
    
    #Nspikes_sprx_all_Na.append(Nspikes_sprx_Na)
    #I_Nspikes_sprx_all_Na.append(I_Nspikes_sprx_Na)

# Plotting
# May choose to specify the color again: # (can skip some of the complicated stuff here since I plot cell for cell)
color_pv_all   = [['#1f77b4','darkblue','#7f7f7f','b','tab:blue'],['#d62728','#ff7f0e','#d62728','m','tab:pink'],['#2ca02c','#8c564b','#9467bd','g','tab:brown']]
color_pv_sprx  = [['rebeccapurple','#e377c2','#8c564b','c','tab:purple'],['#9467bd','#bcbd22','#17becf','xkcd:sky blue','tab:orange'],['olive','darkgreen','mediumseagreen','tab:green','tab:olive']]

plotfolder = 'figures/Comparemodels/'
plotname = plotfolder+'fI_varyg_AllenPV.png'

## avg and rms:
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9), (ax10, ax11, ax12)) = plt.subplots(4, 3, figsize=(20,25))
fig.suptitle(r'Frequency $f$ vs $I$',fontsize=16)

ax1.set_title(r'A',loc='left',fontsize=18)
ax2.set_title(r'B',loc='left',fontsize=18)
ax3.set_title(r'C',loc='left',fontsize=18)
ax4.set_title(r'D',loc='left',fontsize=18)
ax5.set_title(r'E',loc='left',fontsize=18)
ax6.set_title(r'F',loc='left',fontsize=18)
ax7.set_title(r'G',loc='left',fontsize=18)
ax8.set_title(r'H',loc='left',fontsize=18)
ax9.set_title(r'I',loc='left',fontsize=18)
ax10.set_title(r'J',loc='left',fontsize=18)
ax11.set_title(r'K',loc='left',fontsize=18)
ax12.set_title(r'L',loc='left',fontsize=18)

# Add legends! #I could probably have looped this... ax1., ax2. etc. makes it hard.
ax1.set_title(r'Vary $\bar{g}_{\mathregular{CaHVA}}$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gHVA[0]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gHVA[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gHVA[0]
#Nspikes_sprx_this         = Nspikes_sprx_all_gHVA[0]
for j in range(NgHVA):
    ### Everywhere:
    ax1.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygHVA[j],plotlabelgHVA), linewidth=mylinewidth)
    ### Somaprox:
    #ax1.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %.1f*%s' % (varygHVA[j],plotlabelgHVA), linewidth=mylinewidth)
ax1.set_xlabel('$I$ (nA)',fontsize=14)
ax1.set_ylabel('$f$ (Hz)',fontsize=14)
ax1.legend(loc='upper left',ncol=1)
#ax1.legend(bbox_to_anchor =(0.5,-0.27), loc='lower center',ncol=3)

ax2.set_title(r'Vary $\bar{g}_{\mathregular{CaHVA}}$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gHVA[1]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gHVA[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gHVA[1]
#Nspikes_sprx_this         = Nspikes_sprx_all_gHVA[1]
for j in range(NgHVA):
    ### Everywhere:
    ax2.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygHVA[j],plotlabelgHVA), linewidth=mylinewidth)
    ### Somaprox:
    #ax2.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %.1f*%s' % (varygHVA[j],plotlabelgHVA), linewidth=mylinewidth)
ax2.set_xlabel('$I$ (nA)',fontsize=14)
ax2.set_ylabel('$f$ (Hz)',fontsize=14)

ax3.set_title(r'Vary $\bar{g}_{\mathregular{CaHVA}}$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gHVA[2]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gHVA[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gHVA[2]
#Nspikes_sprx_this         = Nspikes_sprx_all_gHVA[2]
for j in range(NgHVA):
    ### Everywhere:
    ax3.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygHVA[j],plotlabelgHVA), linewidth=mylinewidth)
    ### Somaprox:
    #ax3.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %.1f*%s' % (varygHVA[j],plotlabelgHVA), linewidth=mylinewidth)
ax3.set_xlabel('$I$ (nA)',fontsize=14)
ax3.set_ylabel('$f$ (Hz)',fontsize=14)

########################### NaV ####################################################
NgNaV = len(varygNaV)

i_master_everywhere_all_gNaV   = []
Nspikes_everywhere_all_gNaV    = []
I_Nspikes_everywhere_all_gNaV  = []

i_master_sprx_all_gNaV         = []
Nspikes_sprx_all_gNaV          = []
I_Nspikes_sprx_all_gNaV        = []

for testmodel in testmodels:
    i_master_everywhere_gNaV   = []
    Nspikes_everywhere_gNaV    = []
    I_Nspikes_everywhere_gNaV  = []
    
    #i_master_sprx_gNaV         = []
    #Nspikes_sprx_gNaV          = []
    #I_Nspikes_sprx_gNaV        = []

    infolder      = 'figures/%i/' % (testmodel)
    vrestfolder   = infolder 
    for g in varygNaV:
        Nspikes   = []
        I_Nspikes = []
        
        namestringNaV = namestringfirstNaV+str(g)+'p'
        
        # Set names
        infilename_Nspikes = infolder+'%s_%i_cmfall'%(namestringNaV,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
        print('infilename_Nspikes:',infilename_Nspikes)
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
        
        Nspikes_everywhere_gNaV.append(Nspikes)
        I_Nspikes_everywhere_gNaV.append(I_Nspikes)
    
    
        ############# SOMAPROX ##################################
        #Nspikes   = []
        #I_Nspikes = []
        
        # Set names
        #infilename_Nspikes = infolder+'%s_%i_cmfsprx'%(namestringNaV,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
        # Read files
        #infile_Nspikes = open(infilename_Nspikes,'r')
    
        #lines_Nspikes = infile_Nspikes.readlines()
        #Nlines_Nspikes = len(lines_Nspikes)
    
        #for i in range(Nlines_Nspikes):
        #    words_Nspikes = lines_Nspikes[i].split()
        #    if len(words_Nspikes)>0:
        #        I_Nspikes.append(float(words_Nspikes[0]))
        #        Nspikes.append(float(words_Nspikes[1]))

        #infile_Nspikes.close()
        
        #Nspikes_sprx_gNaV.append(Nspikes)
        #I_Nspikes_sprx_gNaV.append(I_Nspikes)
    print('-------------------------------------')
    Nspikes_everywhere_all_gNaV.append(Nspikes_everywhere_gNaV)
    I_Nspikes_everywhere_all_gNaV.append(I_Nspikes_everywhere_gNaV)
    
    #Nspikes_sprx_all_gNaV.append(Nspikes_sprx_gNaV)
    #I_Nspikes_sprx_all_gNaV.append(I_Nspikes_sprx_gNaV)

# Plotting

# Add legends! #I could probably have looped this... ax1., ax2. etc. makes it hard.
ax4.set_title(r'Vary $\bar{g}_{\mathregular{NaV}}$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gNaV[0]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gNaV[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gNaV[0]
#Nspikes_sprx_this         = Nspikes_sprx_all_gNaV[0]
for j in range(NgNaV):
    ### Everywhere:
    ax4.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygNaV[j],plotlabelgNaV), linewidth=mylinewidth)
    ### Somaprox:
    #ax4.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %.1f*%s' % (varygNaV[j],plotlabelgNaV), linewidth=mylinewidth)
ax4.set_xlabel('$I$ (nA)',fontsize=14)
ax4.set_ylabel('$f$ (Hz)',fontsize=14)
ax4.legend(loc='upper left',ncol=1)
#ax1.legend(bbox_to_anchor =(0.5,-0.27), loc='lower center',ncol=3)

ax5.set_title(r'Vary $\bar{g}_{\mathregular{NaV}}$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gNaV[1]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gNaV[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gNaV[1]
#Nspikes_sprx_this         = Nspikes_sprx_all_gNaV[1]
for j in range(NgNaV):
    ### Everywhere:
    ax5.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygNaV[j],plotlabelgNaV), linewidth=mylinewidth)
    ### Somaprox:
    #ax5.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %.1f*%s' % (varygNaV[j],plotlabelgNaV), linewidth=mylinewidth)
ax5.set_xlabel('$I$ (nA)',fontsize=14)
ax5.set_ylabel('$f$ (Hz)',fontsize=14)

ax6.set_title(r'Vary $\bar{g}_{\mathregular{NaV}}$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gNaV[2]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gNaV[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gNaV[2]
#Nspikes_sprx_this         = Nspikes_sprx_all_gNaV[2]
for j in range(NgNaV):
    ### Everywhere:
    ax6.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygNaV[j],plotlabelgNaV), linewidth=mylinewidth)
    ### Somaprox:
    #ax6.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %.1f*%s' % (varygNaV[j],plotlabelgNaV), linewidth=mylinewidth)
ax6.set_xlabel('$I$ (nA)',fontsize=14)
ax6.set_ylabel('$f$ (Hz)',fontsize=14)



########################### Kv31 #################################
NgKv31 = len(varygKv31)

i_master_everywhere_all_gKv31   = []
Nspikes_everywhere_all_gKv31    = []
I_Nspikes_everywhere_all_gKv31  = []

i_master_sprx_all_gKv31         = []
Nspikes_sprx_all_gKv31          = []
I_Nspikes_sprx_all_gKv31        = []

for testmodel in testmodels:
    i_master_everywhere_gKv31   = []
    Nspikes_everywhere_gKv31    = []
    I_Nspikes_everywhere_gKv31  = []
    
    #i_master_sprx_gKv31         = []
    #Nspikes_sprx_gKv31          = []
    #I_Nspikes_sprx_gKv31        = []

    infolder      = 'figures/%i/' % (testmodel)
    vrestfolder   = infolder 
    for g in varygKv31:
        Nspikes   = []
        I_Nspikes = []
        
        namestringKv31 = namestringfirstKv31+str(g)+'p'
        
        # Set names
        infilename_Nspikes = infolder+'%s_%i_cmfall'%(namestringKv31,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
        print('infilename_Nspikes:',infilename_Nspikes)
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
        
        Nspikes_everywhere_gKv31.append(Nspikes)
        I_Nspikes_everywhere_gKv31.append(I_Nspikes)
    
    
        ############# SOMAPROX ##################################
        #Nspikes   = []
        #I_Nspikes = []
        
        # Set names
        #infilename_Nspikes = infolder+'%s_%i_cmfsprx'%(namestringNaV,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
        # Read files
        #infile_Nspikes = open(infilename_Nspikes,'r')
    
        #lines_Nspikes = infile_Nspikes.readlines()
        #Nlines_Nspikes = len(lines_Nspikes)
    
        #for i in range(Nlines_Nspikes):
        #    words_Nspikes = lines_Nspikes[i].split()
        #    if len(words_Nspikes)>0:
        #        I_Nspikes.append(float(words_Nspikes[0]))
        #        Nspikes.append(float(words_Nspikes[1]))

        #infile_Nspikes.close()
        
        #Nspikes_sprx_gKv31.append(Nspikes)
        #I_Nspikes_sprx_gKv31.append(I_Nspikes)
    print('-------------------------------------')
    Nspikes_everywhere_all_gKv31.append(Nspikes_everywhere_gKv31)
    I_Nspikes_everywhere_all_gKv31.append(I_Nspikes_everywhere_gKv31)
    
    #Nspikes_sprx_all_gKv31.append(Nspikes_sprx_gKv31)
    #I_Nspikes_sprx_all_gKv31.append(I_Nspikes_sprx_gKv31)

# Plotting

# Add legends! #I could probably have looped this... ax1., ax2. etc. makes it hard.
ax7.set_title(r'Vary $\bar{g}_{\mathregular{Kv3.1}}$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gKv31[0]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gKv31[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gKv31[0]
#Nspikes_sprx_this         = Nspikes_sprx_all_gKv31[0]
for j in range(NgKv31):
    ### Everywhere:
    ax7.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygKv31[j],plotlabelgKv31), linewidth=mylinewidth)
    ### Somaprox:
    #ax7.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %.1f*%s' % (varygKv31[j],plotlabelgKv31), linewidth=mylinewidth)
ax7.set_xlabel('$I$ (nA)',fontsize=14)
ax7.set_ylabel('$f$ (Hz)',fontsize=14)
ax7.legend(loc='upper left',ncol=1)

ax8.set_title(r'Vary $\bar{g}_{\mathregular{Kv3.1}}$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gKv31[1]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gKv31[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gKv31[1]
#Nspikes_sprx_this         = Nspikes_sprx_all_gKv31[1]
for j in range(NgKv31):
    ### Everywhere:
    ax8.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygKv31[j],plotlabelgKv31), linewidth=mylinewidth)
    ### Somaprox:
    #ax8.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %.1f*%s' % (varygKv31[j],plotlabelgKv31), linewidth=mylinewidth)
ax8.set_xlabel('$I$ (nA)',fontsize=14)
ax8.set_ylabel('$f$ (Hz)',fontsize=14)

ax9.set_title(r'Vary $\bar{g}_{\mathregular{Kv3.1}}$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gKv31[2]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gKv31[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gKv31[2]
#Nspikes_sprx_this         = Nspikes_sprx_all_gKv31[2]
for j in range(NgKv31):
    ### Everywhere:
    ax9.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygKv31[j],plotlabelgKv31), linewidth=mylinewidth)
    ### Somaprox:
    #ax9.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %.1f*%s' % (varygKv31[j],plotlabelgKv31), linewidth=mylinewidth)
ax9.set_xlabel('$I$ (nA)',fontsize=14)
ax9.set_ylabel('$f$ (Hz)',fontsize=14)

############################ K_T ############################################
NgKT = len(varygKT)

i_master_everywhere_all_gKT   = []
Nspikes_everywhere_all_gKT    = []
I_Nspikes_everywhere_all_gKT  = []

i_master_sprx_all_gKT         = []
Nspikes_sprx_all_gKT          = []
I_Nspikes_sprx_all_gKT        = []

for testmodel in testmodels:
    i_master_everywhere_gKT   = []
    Nspikes_everywhere_gKT    = []
    I_Nspikes_everywhere_gKT  = []
    
    #i_master_sprx_gKT         = []
    #Nspikes_sprx_gKT          = []
    #I_Nspikes_sprx_gKT        = []

    infolder      = 'figures/%i/' % (testmodel)
    vrestfolder   = infolder 
    for g in varygKT:
        Nspikes   = []
        I_Nspikes = []
        
        namestringKT = namestringfirstKT+str(g)+'p'
        
        # Set names
        infilename_Nspikes = infolder+'%s_%i_cmfall'%(namestringKT,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
        print('infilename_Nspikes:',infilename_Nspikes)
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
        
        Nspikes_everywhere_gKT.append(Nspikes)
        I_Nspikes_everywhere_gKT.append(I_Nspikes)
    
    
        ############# SOMAPROX ##################################
        #Nspikes   = []
        #I_Nspikes = []
        
        # Set names
        #infilename_Nspikes = infolder+'%s_%i_cmfsprx'%(namestringKT,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
        # Read files
        #infile_Nspikes = open(infilename_Nspikes,'r')
    
        #lines_Nspikes = infile_Nspikes.readlines()
        #Nlines_Nspikes = len(lines_Nspikes)
    
        #for i in range(Nlines_Nspikes):
        #    words_Nspikes = lines_Nspikes[i].split()
        #    if len(words_Nspikes)>0:
        #        I_Nspikes.append(float(words_Nspikes[0]))
        #        Nspikes.append(float(words_Nspikes[1]))

        #infile_Nspikes.close()
        
        #Nspikes_sprx_gKT.append(Nspikes)
        #I_Nspikes_sprx_gKT.append(I_Nspikes)
    print('-------------------------------------')
    Nspikes_everywhere_all_gKT.append(Nspikes_everywhere_gKT)
    I_Nspikes_everywhere_all_gKT.append(I_Nspikes_everywhere_gKT)
    
    #Nspikes_sprx_all_gKT.append(Nspikes_sprx_gKT)
    #I_Nspikes_sprx_all_gKT.append(I_Nspikes_sprx_gKT)

# Plotting

# Add legends! #I could probably have looped this... ax1., ax2. etc. makes it hard.
ax10.set_title(r'Vary $\bar{g}_{\mathregular{KT}}$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gKT[0]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gKT[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gKT[0]
#Nspikes_sprx_this         = Nspikes_sprx_all_gKT[0]
for j in range(NgKT):
    ### Everywhere:
    ax10.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygKT[j],plotlabelgKT), linewidth=mylinewidth)
    ### Somaprox:
    #ax10.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %.1f*%s' % (varygKT[j],plotlabelgKT), linewidth=mylinewidth)
ax10.set_xlabel('$I$ (nA)',fontsize=14)
ax10.set_ylabel('$f$ (Hz)',fontsize=14)
ax10.legend(loc='upper left',ncol=1)

ax5.set_title(r'Vary $\bar{g}_{\mathregular{KT}}$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gKT[1]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gKT[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gKT[1]
#Nspikes_sprx_this         = Nspikes_sprx_all_gKT[1]
for j in range(NgKT):
    ### Everywhere:
    ax11.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygKT[j],plotlabelgKT), linewidth=mylinewidth)
    ### Somaprox:
    #ax11.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %.1f*%s' % (varygKT[j],plotlabelgKT), linewidth=mylinewidth)
ax11.set_xlabel('$I$ (nA)',fontsize=14)
ax11.set_ylabel('$f$ (Hz)',fontsize=14)

ax12.set_title(r'Vary $\bar{g}_{\mathregular{KT}}$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gKT[2]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gKT[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gKT[2]
#Nspikes_sprx_this         = Nspikes_sprx_all_gKT[2]
for j in range(NgKT):
    ### Everywhere:
    ax12.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygKT[j],plotlabelgKT), linewidth=mylinewidth)
    ### Somaprox:
    #ax12.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %.1f*%s' % (varygKT[j],plotlabelgKT), linewidth=mylinewidth)
ax12.set_xlabel('$I$ (nA)',fontsize=14)
ax12.set_ylabel('$f$ (Hz)',fontsize=14)

fig.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.savefig(plotname)

plt.show()
