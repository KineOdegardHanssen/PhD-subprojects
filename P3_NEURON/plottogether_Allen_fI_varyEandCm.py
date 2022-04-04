import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

# change the default font family
plt.rcParams.update({'font.family':'Arial'})

plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('legend', fontsize=15)
fig = plt.figure(figsize=(20,20),dpi=300)
fig.suptitle(r'Frequency $f$ vs $I$',fontsize=20)

gs = gridspec.GridSpec(4, 8)

ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
ax3 = plt.subplot(gs[0, 4:6])
ax5 = plt.subplot(gs[1, 0:2])
ax6 = plt.subplot(gs[1, 2:4])
ax7 = plt.subplot(gs[1, 4:6])
ax9 = plt.subplot(gs[2, 0:2])
ax10 = plt.subplot(gs[2, 2:4])
ax11 = plt.subplot(gs[2, 4:6])
ax13 = plt.subplot(gs[3, 0:2])
ax14 = plt.subplot(gs[3, 2:4])
ax15 = plt.subplot(gs[3, 4:6])

ax1.set_title(r'A',loc='left',fontsize=18)
ax2.set_title(r'B',loc='left',fontsize=18)
ax3.set_title(r'C',loc='left',fontsize=18)
ax5.set_title(r'E',loc='left',fontsize=18)
ax6.set_title(r'F',loc='left',fontsize=18)
ax7.set_title(r'G',loc='left',fontsize=18)
ax9.set_title(r'I',loc='left',fontsize=18)
ax10.set_title(r'J',loc='left',fontsize=18)
ax11.set_title(r'K',loc='left',fontsize=18)
ax13.set_title(r'M',loc='left',fontsize=18)
ax14.set_title(r'N',loc='left',fontsize=18)
ax15.set_title(r'O',loc='left',fontsize=18)

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
varyENa = [53,63,73]
namestringfirstNa = namestringfirst + 'ENa'
plotlabelNa = r'$E_{\mathregular{Na}}$'
varyEK = [-127,-117,-107]
namestringfirstK = namestringfirst + 'EK'
plotlabelK = r'$E_{\mathregular{K}}$'
varyEpas = [-20,-10,0]   
namestringfirstpas = namestringfirst + 'Epasplus'
varycao = [2.0,20.0,200.0]
namestringfirstCa = namestringfirst + 'ECa'
namestringfirstcao = namestringfirst + '_cao'
plotlabelCa = r'$E_{\mathregular{Ca}}$'

plotstringNa  = '_vary_Na'
plotstringK   = '_vary_K'
plotstringCa  = '_vary_Ca'
plotstringpas = '_vary_pas'

i_master = [0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5]
cms = [1.0,1.5] 

Ncms = len(cms)
NI  = len(i_master)
NENa = len(varyENa)*Ncms # Do I need ALL, though?

all_pairs_Na = []

i_master_everywhere_all_Na   = []
Nspikes_everywhere_all_Na    = []
I_Nspikes_everywhere_all_Na  = []

i_master_sprx_all_Na         = []
Nspikes_sprx_all_Na          = []
I_Nspikes_sprx_all_Na        = []

for testmodel in testmodels:
    i_master_everywhere_Na   = []
    Nspikes_everywhere_Na    = []
    I_Nspikes_everywhere_Na  = []
    
    i_master_sprx_Na         = []
    Nspikes_sprx_Na          = []
    I_Nspikes_sprx_Na        = []
    
    pairs = []

    infolder      = 'figures/%i/' % (testmodel)
    vrestfolder   = infolder 
    for cm in cms:
        for E in varyENa:
            pairs.append([cm,E])
            Nspikes   = []
            I_Nspikes = []
            
            namestringNa = namestringfirstNa+str(E)
            
            # Set names
            infilename_Nspikes = infolder+'%s_%i_cmfall'%(namestringNa,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
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
            
            Nspikes_everywhere_Na.append(Nspikes)
            I_Nspikes_everywhere_Na.append(I_Nspikes)
    
            '''
            ############# SOMAPROX ##################################
            Nspikes   = []
            I_Nspikes = []
            
            # Set names
            infilename_Nspikes = infolder+'%s_%i_cmfsprx'%(namestringNa,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
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
                
            Nspikes_sprx_Na.append(Nspikes)
            I_Nspikes_sprx_Na.append(I_Nspikes)
            '''
    Nspikes_everywhere_all_Na.append(Nspikes_everywhere_Na)
    I_Nspikes_everywhere_all_Na.append(I_Nspikes_everywhere_Na)
    
    #Nspikes_sprx_all_Na.append(Nspikes_sprx_Na)
    #I_Nspikes_sprx_all_Na.append(I_Nspikes_sprx_Na)
    all_pairs_Na.append(pairs)

# Plotting
# May choose to specify the color again: # (can skip some of the complicated stuff here since I plot cell for cell)
color_pv_all   = [['#1f77b4','darkblue','#7f7f7f','b','tab:blue'],['#d62728','#ff7f0e','#d62728','m','tab:pink'],['#2ca02c','#8c564b','#9467bd','g','tab:brown']]
color_pv_sprx  = [['rebeccapurple','#e377c2','#8c564b','c','tab:purple'],['#9467bd','#bcbd22','#17becf','xkcd:sky blue','tab:orange'],['olive','darkgreen','mediumseagreen','tab:green','tab:olive']]

plotfolder = 'figures/Comparemodels/'
plotname = plotfolder+'fI_varyEandCm_AllenPV.png'

# Add legends! #I could probably have looped this... ax1., ax2. etc. makes it hard.
ax1.set_title(r'Vary $E_{\mathregular{Na}}$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Na[0]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Na[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Na[0]
#Nspikes_sprx_this         = Nspikes_sprx_all_Na[0]
thesepairs                = all_pairs_Na[0]
for j in range(NENa):
    pairs   = thesepairs[j]
    cmthis  = pairs[0]
    ENathis = pairs[1]
    ### Everywhere:
    ax1.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'$C_m$*%.1f, %s=%i mV' % (cmthis,plotlabelNa,ENathis), linewidth=mylinewidth)
    ### Somaprox:
    #ax1.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %s=%i' % (plotlabelNa,varyENa[j]), linewidth=mylinewidth)
ax1.set_xlabel('$I$ (nA)',fontsize=16)
ax1.set_ylabel('$f$ (Hz)',fontsize=16)
ax1.legend(loc='lower right',ncol=1)
#ax1.legend(bbox_to_anchor =(0.5,-0.27), loc='lower center',ncol=3)

ax2.set_title(r'Vary $E_{\mathregular{Na}}$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Na[1]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Na[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Na[1]
#Nspikes_sprx_this         = Nspikes_sprx_all_Na[1]
thesepairs                = all_pairs_Na[1]
for j in range(NENa):
    pairs   = thesepairs[j]
    cmthis  = pairs[0]
    ENathis = pairs[1]
    ### Everywhere:
    ax2.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'$C_m$*%.1f, %s=%i mV' % (cmthis,plotlabelNa,ENathis), linewidth=mylinewidth)
    ### Somaprox:
    #ax2.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.2f' % (plotlabelNa,varyENa[j]), linewidth=mylinewidth)
ax2.set_xlabel('$I$ (nA)',fontsize=16)
ax2.set_ylabel('$f$ (Hz)',fontsize=16)

ax3.set_title(r'Vary $E_{\mathregular{Na}}$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Na[2]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Na[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Na[2]
#Nspikes_sprx_this         = Nspikes_sprx_all_Na[2]
thesepairs                = all_pairs_Na[2]
for j in range(NENa):
    pairs   = thesepairs[j]
    cmthis  = pairs[0]
    ENathis = pairs[1]
    ### Everywhere:
    ax3.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'$C_m$*%.1f, %s=%i mV' % (cmthis,plotlabelNa,ENathis), linewidth=mylinewidth)
    ### Somaprox:
    #ax3.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.2f' % (plotlabelNa,varyENa[j]), linewidth=mylinewidth)
ax3.set_xlabel('$I$ (nA)',fontsize=16)
ax3.set_ylabel('$f$ (Hz)',fontsize=16)

### Do the same for K: #################################################################

NEK = len(varyEK)*Ncms

i_master_everywhere_all_K   = []
Nspikes_everywhere_all_K    = []
I_Nspikes_everywhere_all_K  = []

i_master_sprx_all_K         = []
Nspikes_sprx_all_K          = []
I_Nspikes_sprx_all_K        = []

all_pairs_K = []
for testmodel in testmodels:
    i_master_everywhere_K   = []
    Nspikes_everywhere_K    = []
    I_Nspikes_everywhere_K  = []
    
    i_master_sprx_K         = []
    Nspikes_sprx_K          = []
    I_Nspikes_sprx_K        = []
    
    pairs = []

    infolder      = 'figures/%i/' % (testmodel)
    vrestfolder   = infolder 
    for cm in cms:
        for E in varyEK:
            pairs.append([cm,E])
            Nspikes   = []
            I_Nspikes = []
            
            namestringK = namestringfirstK+str(E)
            
            # Set names
            infilename_Nspikes = infolder+'%s_%i_cmfall'%(namestringK,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
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
            
            Nspikes_everywhere_K.append(Nspikes)
            I_Nspikes_everywhere_K.append(I_Nspikes)
            
            '''
            ############# SOMAPROX ##################################
            Nspikes   = []
            I_Nspikes = []
            
            # Set names
            infilename_Nspikes = infolder+'%s_%i_cmfsprx'%(namestringK,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
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
            
            Nspikes_sprx_K.append(Nspikes)
            I_Nspikes_sprx_K.append(I_Nspikes)
            '''
    Nspikes_everywhere_all_K.append(Nspikes_everywhere_K)
    I_Nspikes_everywhere_all_K.append(I_Nspikes_everywhere_K)
    
    #Nspikes_sprx_all_K.append(Nspikes_sprx_K)
    #I_Nspikes_sprx_all_K.append(I_Nspikes_sprx_K)
    all_pairs_K.append(pairs)


ax5.set_title(r'Vary $E_{\mathregular{K}}$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_K[0]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_K[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_K[0]
#Nspikes_sprx_this         = Nspikes_sprx_all_K[0]
thesepairs                = all_pairs_K[0]
for j in range(NEK):
    pairs  = thesepairs[j]
    cmthis = pairs[0]
    EKthis = pairs[1]
    ### Everywhere:
    ax5.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'$C_m$*%.1f, %s=%i mV' % (cmthis,plotlabelK,EKthis), linewidth=mylinewidth)
    ### Somaprox:
    #ax5.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %s=%i' % (plotlabelK,varyEK[j]), linewidth=mylinewidth)
ax5.set_xlabel('$I$ (nA)',fontsize=16)
ax5.set_ylabel('$f$ (Hz)',fontsize=16)
ax5.legend(loc='lower right',ncol=1)
#ax5.legend(bbox_to_anchor =(0.5,-0.27), loc='lower center',ncol=3)

ax6.set_title(r'Vary $E_{\mathregular{K}}$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_K[1]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_K[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_K[1]
#Nspikes_sprx_this         = Nspikes_sprx_all_K[1]
thesepairs                = all_pairs_K[1]
for j in range(NEK):
    pairs  = thesepairs[j]
    cmthis = pairs[0]
    EKthis = pairs[1]
    ### Everywhere:
    ax6.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'$C_m$*%.1f, %s=%i mV' % (cmthis,plotlabelK,EKthis), linewidth=mylinewidth)
    ### Somaprox:
    #ax6.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.2f' % (plotlabelK,varyEK[j]), linewidth=mylinewidth)
ax6.set_xlabel('$I$ (nA)',fontsize=16)
ax6.set_ylabel('$f$ (Hz)',fontsize=16)

ax7.set_title(r'Vary $E_{\mathregular{K}}$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_K[2]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_K[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_K[2]
#Nspikes_sprx_this         = Nspikes_sprx_all_K[2]
thesepairs                = all_pairs_K[2]
for j in range(NEK):
    pairs  = thesepairs[j]
    cmthis = pairs[0]
    EKthis = pairs[1]
    ### Everywhere:
    ax7.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'$C_m$*%.1f, %s=%i mV' % (cmthis,plotlabelK,EKthis), linewidth=mylinewidth)
    ### Somaprox:
    #ax7.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.2f' % (plotlabelK,varyEK[j]), linewidth=mylinewidth)
ax7.set_xlabel('$I$ (nA)',fontsize=16)
ax7.set_ylabel('$f$ (Hz)',fontsize=16)

################ Do the same for Epas: #####################################
NEpas = len(varyEpas)*len(cms)

plotlabelEpas = []

i_master_everywhere_all_pas   = []
Nspikes_everywhere_all_pas    = []
I_Nspikes_everywhere_all_pas  = []

i_master_sprx_all_pas         = []
Nspikes_sprx_all_pas          = []
I_Nspikes_sprx_all_pas        = []

all_pairs_pas = []
for testmodel in testmodels:
    i_master_everywhere_pas   = []
    Nspikes_everywhere_pas    = []
    I_Nspikes_everywhere_pas  = []
    
    i_master_sprx_pas         = []
    Nspikes_sprx_pas          = []
    I_Nspikes_sprx_pas        = []
    
    pairs = []

    infolder      = 'figures/%i/' % (testmodel)
    vrestfolder   = infolder
    for cm in cms:
        for E in varyEpas:
            pairs.append([cm,E])
            if E>=0:
                Estring = '$E_{\mathregular{pas}}$ + %i' % E
            else:
                Estring = '$E_{\mathregular{pas}}$ - %i' % abs(E)
            plotlabelEpas.append(Estring)
            Nspikes   = []
            I_Nspikes = []
            
            namestringpas = namestringfirstpas+str(E)
            
            # Set names
            infilename_Nspikes = infolder+'%s_%i_cmfall'%(namestringpas,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
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
            
            Nspikes_everywhere_pas.append(Nspikes)
            I_Nspikes_everywhere_pas.append(I_Nspikes)
        
        
            ############# SOMAPROX ##################################
            '''
            Nspikes   = []
            I_Nspikes = []
            
            # Set names
            infilename_Nspikes = infolder+'%s_%i_cmfsprx'%(namestringpas,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
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
            
            Nspikes_sprx_pas.append(Nspikes)
            I_Nspikes_sprx_pas.append(I_Nspikes)
            '''
    Nspikes_everywhere_all_pas.append(Nspikes_everywhere_pas)
    I_Nspikes_everywhere_all_pas.append(I_Nspikes_everywhere_pas)
    all_pairs_pas.append(pairs)
    
    #Nspikes_sprx_all_pas.append(Nspikes_sprx_pas)
    #I_Nspikes_sprx_all_pas.append(I_Nspikes_sprx_pas)

# Plotting

ax9.set_title(r'Vary $E_{\mathregular{pas}}$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_pas[0]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_pas[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_pas[0]
#Nspikes_sprx_this         = Nspikes_sprx_all_pas[0]
thesepairs                = all_pairs_pas[0]
for j in range(NEpas):
    pairs    = thesepairs[j]
    cmthis   = pairs[0]
    Epasthis = pairs[1]
    if Epasthis==0:
        Epasthis = '+0'
    else:
        Epasthis=str(Epasthis)
    ### Everywhere:
    ax9.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'$C_m$*%.1f, %s%s mV' % (cmthis,plotlabelEpas[j],Epasthis), linewidth=mylinewidth)
    ### Somaprox:
    #ax9.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %s' % plotlabelEpas[j], linewidth=mylinewidth)
ax9.set_xlabel('$I$ (nA)',fontsize=14)
ax9.set_ylabel('$f$ (Hz)',fontsize=14)
ax9.legend(loc='upper left',ncol=1)
#ax9.legend(bbox_to_anchor =(0.5,-0.27), loc='lower center',ncol=3)

ax10.set_title(r'Vary $E_{\mathregular{pas}}$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_pas[1]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_pas[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_pas[1]
#Nspikes_sprx_this         = Nspikes_sprx_all_pas[1]
thesepairs                = all_pairs_pas[1]
for j in range(NEpas):
    pairs    = thesepairs[j]
    cmthis   = pairs[0]
    Epasthis = pairs[1]
    if Epasthis==0:
        Epasthis = '+0'
    else:
        Epasthis=str(Epasthis)
    ### Everywhere:
    ax10.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'$C_m$*%.1f, %s%s mV' % (cmthis,plotlabelEpas[j],Epasthis), linewidth=mylinewidth)
    ### Somaprox:
    #ax10.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s' % plotlabelEpas[j], linewidth=mylinewidth)
ax10.set_xlabel('$I$ (nA)',fontsize=14)
ax10.set_ylabel('$f$ (Hz)',fontsize=14)

ax11.set_title(r'Vary $E_{\mathregular{pas}}$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_pas[2]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_pas[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_pas[2]
#Nspikes_sprx_this         = Nspikes_sprx_all_pas[2]
thesepairs                = all_pairs_pas[2]
for j in range(NEpas):
    pairs    = thesepairs[j]
    cmthis   = pairs[0]
    Epasthis = pairs[1]
    if Epasthis==0:
        Epasthis = '+0'
    else:
        Epasthis=str(Epasthis)
    ### Everywhere:
    ax11.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'$C_m$*%.1f, %s%s mV' % (cmthis,plotlabelEpas[j],Epasthis), linewidth=mylinewidth)
    ### Somaprox:
    #ax11.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s' % plotlabelEpas[j], linewidth=mylinewidth)
ax11.set_xlabel('$I$ (nA)',fontsize=14)
ax11.set_ylabel('$f$ (Hz)',fontsize=14)

## CA:
### Do the same for Ca: #################################################################
R = 8.314   # JK-1mol-1
F = 9.648e4 # Cmol-1
T = 307.15  # K
prefactor = 1000*R*T/(2.0*F)
cai0 = 1e-4 # mM

NECa = len(varycao)*Ncms
ECa0 = np.zeros(NECa)

i_master_everywhere_all_Ca   = []
Nspikes_everywhere_all_Ca    = []
I_Nspikes_everywhere_all_Ca  = []

i_master_sprx_all_Ca         = []
Nspikes_sprx_all_Ca          = []
I_Nspikes_sprx_all_Ca        = []

all_pairs_Ca = []
for testmodel in testmodels:
    i_master_everywhere_Ca   = []
    Nspikes_everywhere_Ca    = []
    I_Nspikes_everywhere_Ca  = []
    
    i_master_sprx_Ca         = []
    Nspikes_sprx_Ca          = []
    I_Nspikes_sprx_Ca        = []
    
    pairs = []

    infolder      = 'figures/%i/' % (testmodel)
    vrestfolder   = infolder 
    for cm in cms:
        for cao in varycao:
            E = prefactor*np.log(cao/cai0)
            pairs.append([cm,E])
            Nspikes   = []
            I_Nspikes = []
            
            namestringCa = namestringfirstcao+str(cao)
            
            # Set names
            infilename_Nspikes = infolder+'%s_%i_cmfall'%(namestringCa,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
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
            
            Nspikes_everywhere_Ca.append(Nspikes)
            I_Nspikes_everywhere_Ca.append(I_Nspikes)
            
            '''
            ############# SOMAPROX ##################################
            Nspikes   = []
            I_Nspikes = []
            
            # Set names
            infilename_Nspikes = infolder+'%s_%i_cmfsprx'%(namestringCa,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
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
            
            Nspikes_sprx_Ca.append(Nspikes)
            I_Nspikes_sprx_Ca.append(I_Nspikes)
            '''
    Nspikes_everywhere_all_Ca.append(Nspikes_everywhere_Ca)
    I_Nspikes_everywhere_all_Ca.append(I_Nspikes_everywhere_Ca)
    
    #Nspikes_sprx_all_K.append(Nspikes_sprx_Ca)
    #I_Nspikes_sprx_all_K.append(I_Nspikes_sprx_Ca)
    all_pairs_Ca.append(pairs)


ax13.set_title(r'Vary $E_{\mathregular{Ca}}$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Ca[0]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Ca[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Ca[0]
#Nspikes_sprx_this         = Nspikes_sprx_all_Ca[0]
thesepairs                = all_pairs_Ca[0]
for j in range(NECa):
    pairs  = thesepairs[j]
    cmthis = pairs[0]
    ECathis = pairs[1]
    ### Everywhere:
    ax13.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'$C_m$*%.1f, %s=%.2f mV' % (cmthis,plotlabelK,ECathis), linewidth=mylinewidth)
    ### Somaprox:
    #ax13.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %s=%i' % (plotlabelK,varyECa[j]), linewidth=mylinewidth)
ax13.set_xlabel('$I$ (nA)',fontsize=16)
ax13.set_ylabel('$f$ (Hz)',fontsize=16)
ax13.legend(loc='lower right',ncol=1)
#ax13.legend(bbox_to_anchor =(0.5,-0.27), loc='lower center',ncol=3)

ax14.set_title(r'Vary $E_{\mathregular{Ca}}$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Ca[1]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Ca[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Ca[1]
#Nspikes_sprx_this         = Nspikes_sprx_all_Ca[1]
thesepairs                = all_pairs_Ca[1]
for j in range(NECa):
    pairs  = thesepairs[j]
    cmthis = pairs[0]
    ECathis = pairs[1]
    ### Everywhere:
    ax14.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'$C_m$*%.1f, %s=%i mV' % (cmthis,plotlabelK,ECathis), linewidth=mylinewidth)
    ### Somaprox:
    #ax14.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.2f' % (plotlabelK,varyECa[j]), linewidth=mylinewidth)
ax14.set_xlabel('$I$ (nA)',fontsize=16)
ax14.set_ylabel('$f$ (Hz)',fontsize=16)

ax15.set_title(r'Vary $E_{\mathregular{Ca}}$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Ca[2]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Ca[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Ca[2]
#Nspikes_sprx_this         = Nspikes_sprx_all_Ca[2]
thesepairs                = all_pairs_Ca[2]
for j in range(NECa):
    pairs  = thesepairs[j]
    cmthis = pairs[0]
    ECathis = pairs[1]
    ### Everywhere:
    ax15.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'$C_m$*%.1f, %s=%i mV' % (cmthis,plotlabelK,ECathis), linewidth=mylinewidth)
    ### Somaprox:
    #ax15.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.2f' % (plotlabelK,varyECa[j]), linewidth=mylinewidth)
ax15.set_xlabel('$I$ (nA)',fontsize=16)
ax15.set_ylabel('$f$ (Hz)',fontsize=16)

############################ Differences ################################################

idefNa       = 0 # Index that holds the defaul value of ENa
lastdiffs_Na = np.zeros((Nmodels,NENa))
for j in range(Nmodels):
    I_Nspikes_everywhere_all = I_Nspikes_everywhere_all_Na[j]
    Nspikes_everywhere_all = Nspikes_everywhere_all_Na[j]
    for i in range(NENa):
        im_Na, rd1_Na, rd2_Na = reldiffs(Nspikes_everywhere_all[i],Nspikes_everywhere_all[idefNa],I_Nspikes_everywhere_all[i],I_Nspikes_everywhere_all[idefNa]) 
        if len(rd1_Na)>0:
            lastdiffs_Na[j,i] = rd1_Na[-1]
        else:
            lastdiffs_Na[j,i] = 0

idefK = 2 # Index that holds the defaul value of EK
lastdiffs_K = np.zeros((Nmodels,NEK))
for j in range(Nmodels):
    I_Nspikes_everywhere_all = I_Nspikes_everywhere_all_K[j]
    Nspikes_everywhere_all = Nspikes_everywhere_all_K[j]
    for i in range(NEK):
        im_K, rd1_K, rd2_K = reldiffs(Nspikes_everywhere_all[i],Nspikes_everywhere_all[idefK],I_Nspikes_everywhere_all[i],I_Nspikes_everywhere_all[idefK]) 
        if len(rd1_K)>0:
            lastdiffs_K[j,i] = rd1_K[-1]
        else:
            lastdiffs_K[j,i] = 0

idefpas = 2 # Index that holds the defaul value of Epas
lastdiffs_pas = np.zeros((Nmodels,NEpas))
for j in range(Nmodels):
    I_Nspikes_everywhere_all = I_Nspikes_everywhere_all_pas[j]
    Nspikes_everywhere_all = Nspikes_everywhere_all_pas[j]
    for i in range(NEpas):
        im_pas, rd1_pas, rd2_pas = reldiffs(Nspikes_everywhere_all[i],Nspikes_everywhere_all[idefpas],I_Nspikes_everywhere_all[i],I_Nspikes_everywhere_all[idefpas]) 
        if len(rd1_pas)>0:
            lastdiffs_pas[j,i] = rd1_pas[-1]
        else:
            lastdiffs_pas[j,i] = 0

idefCa = 2 # Index that holds the defaul value of ECa
lastdiffs_Ca = np.zeros((Nmodels,NECa))
for j in range(Nmodels):
    I_Nspikes_everywhere_all = I_Nspikes_everywhere_all_Ca[j]
    Nspikes_everywhere_all = Nspikes_everywhere_all_Ca[j]
    for i in range(NECa):
        im_Ca, rd1_Ca, rd2_Ca = reldiffs(Nspikes_everywhere_all[i],Nspikes_everywhere_all[idefCa],I_Nspikes_everywhere_all[i],I_Nspikes_everywhere_all[idefCa]) 
        if len(rd1_Ca)>0:
            lastdiffs_Ca[j,i] = rd1_Ca[-1]
        else:
            lastdiffs_Ca[j,i] = 0

ax4 = plt.subplot(gs[0, 6:8])
ax4.set_title(r'G',loc='left',fontsize=18)
barWidth = 0.18
change_ena1   = [lastdiffs_Na[0,1],lastdiffs_Na[1,1],lastdiffs_Na[2,1]]
change_ena2   = [lastdiffs_Na[0,2],lastdiffs_Na[1,2],lastdiffs_Na[2,2]]
change_ena3   = [lastdiffs_Na[0,3],lastdiffs_Na[1,3],lastdiffs_Na[2,3]]
change_ena4 = [lastdiffs_Na[0,4],lastdiffs_Na[1,4],lastdiffs_Na[2,4]]
change_ena5 = [lastdiffs_Na[0,5],lastdiffs_Na[1,5],lastdiffs_Na[2,5]]

br1 = np.arange(len(change_ena2))
br2 = [x+barWidth for x in br1]
br3 = [x+barWidth for x in br2]
br4 = [x+barWidth for x in br3]
br5 = [x+barWidth for x in br4]
brcenter = br3
pairs_Na = all_pairs_Na[0]
pNa1 = pairs_Na[1]
pNa2 = pairs_Na[2]
pNa3 = pairs_Na[3]
pNa4 = pairs_Na[4]
pNa5 = pairs_Na[5]

# How do I fix the difference between values of Allen and HH? Run more sim.s for HH?
ax4.bar(br1, change_ena1, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{Na}}$=%i mV' % (pNa1[0],pNa1[1]))
ax4.bar(br2, change_ena2, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{Na}}$=%i mV' % (pNa2[0],pNa2[1]))
ax4.bar(br3, change_ena3, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{Na}}$=%i mV' % (pNa3[0],pNa3[1]))
ax4.bar(br4, change_ena4, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{Na}}$=%i mV' % (pNa4[0],pNa4[1]))
ax4.bar(br5, change_ena5, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{Na}}$=%i mV' % (pNa5[0],pNa5[1]))
ax4.axhline(y=0,color='k',linestyle='-',linewidth=1.0)
plt.xticks(brcenter, ['437','965','407'])
ax4.set_xlabel('Model',fontsize=16)
ax4.set_ylabel(r'Relative difference at max. current',fontsize=16)
ax4.set_title(r'Difference from original $E_{\mathregular{Na}}$',fontsize=18)
ax4.legend(loc='lower left')

ax8 = plt.subplot(gs[1, 6:8])
ax8.set_title(r'H',loc='left',fontsize=18)
change_K1   = [lastdiffs_K[0,0],lastdiffs_K[1,0],lastdiffs_K[2,0]]
change_K2   = [lastdiffs_K[0,1],lastdiffs_K[1,1],lastdiffs_K[2,1]]
change_K3   = [lastdiffs_K[0,3],lastdiffs_K[1,3],lastdiffs_K[2,3]]
change_K4   = [lastdiffs_K[0,4],lastdiffs_K[1,4],lastdiffs_K[2,4]]
change_K5   = [lastdiffs_K[0,5],lastdiffs_K[1,5],lastdiffs_K[2,5]]

br1 = np.arange(len(change_K2))
br2 = [x+barWidth for x in br1]
br3 = [x+barWidth for x in br2]
br4 = [x+barWidth for x in br3]
br5 = [x+barWidth for x in br4]
brcenter = br3
pairs_K = all_pairs_K[0]
pNa1 = pairs_K[0]
pNa2 = pairs_K[1]
pNa3 = pairs_K[3]
pNa4 = pairs_K[4]
pNa5 = pairs_K[5]

# How do I fix the difference between values of Allen and HH? Run more sim.s for HH?
ax8.bar(br1, change_K1, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{K}}$=%i mV' % (pNa1[0],pNa1[1]))
ax8.bar(br2, change_K2, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{K}}$=%i mV' % (pNa2[0],pNa2[1]))
ax8.bar(br3, change_K3, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{K}}$=%i mV' % (pNa4[0],pNa4[1]))
ax8.bar(br4, change_K4, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{K}}$=%i mV' % (pNa5[0],pNa5[1]))
ax8.bar(br5, change_K5, width=barWidth, label=r'$C_m$=%.1f , $E_{\mathregular{K}}$=%i mV' % (pNa5[0],pNa5[1]))
ax8.axhline(y=0,color='k',linestyle='-',linewidth=1.0)
plt.xticks(brcenter, ['437','965','407','OC','BAS'])
ax8.set_xlabel('Model',fontsize=16)
ax8.set_ylabel(r'Relative difference at max. current',fontsize=16)
ax8.set_title(r'Difference from original $E_{\mathregular{K}}$',fontsize=18)
ax8.legend(loc='lower left')

## pas
ax12 = plt.subplot(gs[2, 6:8])
ax12.set_title(r'I',loc='left',fontsize=18)
change_pas1   = [lastdiffs_pas[0,0],lastdiffs_pas[1,0],lastdiffs_pas[2,0]]
change_pas2   = [lastdiffs_pas[0,1],lastdiffs_pas[1,1],lastdiffs_pas[2,1]]
change_pas3   = [lastdiffs_pas[0,3],lastdiffs_pas[1,3],lastdiffs_pas[2,3]]
change_pas4   = [lastdiffs_pas[0,4],lastdiffs_pas[1,4],lastdiffs_pas[2,4]]
change_pas5   = [lastdiffs_pas[0,5],lastdiffs_pas[1,5],lastdiffs_pas[2,5]]

br1 = np.arange(len(change_pas2))
br2 = [x+barWidth for x in br1]
br3 = [x+barWidth for x in br2]
br4 = [x+barWidth for x in br3]
br5 = [x+barWidth for x in br4]
brcenter = br3
pairs_pas = all_pairs_pas[0]
pNa1 = pairs_pas[0]
pNa2 = pairs_pas[1]
pNa3 = pairs_pas[3]
pNa4 = pairs_pas[4]
pNa5 = pairs_pas[5]

# How do I fix the difference between values of Allen and HH? Run more sim.s for HH?
ax12.bar(br1, change_pas1, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{pas}}$%i mV' % (pNa1[0],pNa1[1]))
ax12.bar(br2, change_pas2, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{pas}}$%i mV' % (pNa2[0],pNa2[1]))
ax12.bar(br3, change_pas3, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{pas}}$%i mV' % (pNa3[0],pNa3[1]))
ax12.bar(br4, change_pas4, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{pas}}$%i mV' % (pNa4[0],pNa4[1]))
ax12.bar(br5, change_pas5, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{pas}}$%i' % (pNa5[0],pNa5[1]))
ax12.axhline(y=0,color='k',linestyle='-',linewidth=1.0)
plt.xticks(brcenter, ['437','965','407','OC','BAS'])
ax12.set_xlabel('Model',fontsize=16)
ax12.set_ylabel(r'Relative difference at max. current',fontsize=16)
ax12.set_title(r'Difference from original $E_{\mathregular{pas}}$',fontsize=18)
ax12.legend(loc='lower left')

## pas
ax16 = plt.subplot(gs[3, 6:8])
ax16.set_title(r'J',loc='left',fontsize=18)
change_Ca1   = [lastdiffs_Ca[0,1],lastdiffs_Ca[1,1],lastdiffs_Ca[2,1]]
change_Ca2   = [lastdiffs_Ca[0,2],lastdiffs_Ca[1,2],lastdiffs_Ca[2,2]]
change_Ca3   = [lastdiffs_Ca[0,3],lastdiffs_Ca[1,3],lastdiffs_Ca[2,3]]
change_Ca4   = [lastdiffs_Ca[0,4],lastdiffs_Ca[1,4],lastdiffs_Ca[2,4]]
change_Ca5   = [lastdiffs_Ca[0,5],lastdiffs_Ca[1,5],lastdiffs_Ca[2,5]]

br1 = np.arange(len(change_Ca2))
br2 = [x+barWidth for x in br1]
br3 = [x+barWidth for x in br2]
br4 = [x+barWidth for x in br3]
br5 = [x+barWidth for x in br4]
brcenter = br3
pairs_Ca = all_pairs_Ca[0]
pNa1 = pairs_Ca[0]
pNa2 = pairs_Ca[1]
pNa3 = pairs_Ca[3]
pNa4 = pairs_Ca[4]
pNa5 = pairs_Ca[5]

# How do I fix the difference between values of Allen and HH? Run more sim.s for HH?
ax16.bar(br1, change_Ca1, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{Ca}}$=%.2f mV' % (pNa1[0],pNa1[1]))
ax16.bar(br2, change_Ca2, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{Ca}}$=%.2f mV' % (pNa2[0],pNa2[1]))
ax16.bar(br3, change_Ca3, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{Ca}}$=%.2f mV' % (pNa3[0],pNa3[1]))
ax16.bar(br4, change_Ca4, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{Ca}}$=%.2f mV' % (pNa4[0],pNa4[1]))
ax16.bar(br5, change_Ca5, width=barWidth, label=r'$C_m$*%.1f , $E_{\mathregular{Ca}}$=%.2f mV' % (pNa5[0],pNa5[1]))
ax16.axhline(y=0,color='k',linestyle='-',linewidth=1.0)
plt.xticks(brcenter, ['437','965','407','OC','BAS'])
ax16.set_xlabel('Model',fontsize=16)
ax16.set_ylabel(r'Relative difference at max. current',fontsize=16)
ax16.set_title(r'Difference from original $E_{\mathregular{Ca}}$',fontsize=18)
ax16.legend(loc='lower left')

fig.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.savefig(plotname)

plt.show()
