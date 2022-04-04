import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

# change the default font family
plt.rcParams.update({'font.family':'Arial'})

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend', fontsize=13)

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
varyENa = [40,50,53,63,73]
namestringfirstNa = namestringfirst + 'ENa'
plotlabelNa = r'$E_{\mathregular{Na}}$'
varyEK = [-127,-117,-107,-100,-90,-80]
namestringfirstK = namestringfirst + 'EK'
plotlabelK = r'$E_{\mathregular{K}}$'
varyEpas = [-20,-10,0,10,20] # Ugh, should have zero too     
namestringfirstpas = namestringfirst + 'Epasplus'
plotlabelEpas = []
for i in range(len(varyEpas)):
    Eval = varyEpas[i]
    if Eval>=0:
        Estring = '$E_{\mathregular{pas}}$ + %i' % Eval
    else:
        Estring = '$E_{\mathregular{pas}}$ - %i' % abs(Eval)
    plotlabelEpas.append(Estring)
varycao = [0.02,0.2,2.0,20.0,200.0]
namestringfirstCa = namestringfirst + 'ECa'
namestringfirstcao = namestringfirst + '_cao'
plotlabelCa = r'$E_{\mathregular{Ca}}$'


plotstringNa  = '_vary_Na'
plotstringK   = '_vary_K'
plotstringCa  = '_vary_Ca'
plotstringpas = '_vary_pas'

i_master = [0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5]
cm = 1.0 # Keep this a scalar for now
    
NI  = len(i_master)
NENa = len(varyENa)

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

    infolder      = 'Allen_test_changecapacitance/figures/%i/' % (testmodel)
    vrestfolder   = infolder 
    for E in varyENa:
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
    
    
        ############# SOMAPROX ##################################
        #Nspikes   = []
        #I_Nspikes = []
        #
        ## Set names
        #infilename_Nspikes = infolder+'%s_%i_cmfsprx'%(namestringNa,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
        ## Read files
        #infile_Nspikes = open(infilename_Nspikes,'r')
        #
        #lines_Nspikes = infile_Nspikes.readlines()
        #Nlines_Nspikes = len(lines_Nspikes)
        #
        #for i in range(Nlines_Nspikes):
        #    words_Nspikes = lines_Nspikes[i].split()
        #    if len(words_Nspikes)>0:
        #        I_Nspikes.append(float(words_Nspikes[0]))
        #        Nspikes.append(float(words_Nspikes[1]))
        #
        #infile_Nspikes.close()
        #
        #Nspikes_sprx_Na.append(Nspikes)
        #I_Nspikes_sprx_Na.append(I_Nspikes)
    Nspikes_everywhere_all_Na.append(Nspikes_everywhere_Na)
    I_Nspikes_everywhere_all_Na.append(I_Nspikes_everywhere_Na)
    
    #Nspikes_sprx_all_Na.append(Nspikes_sprx_Na)
    #I_Nspikes_sprx_all_Na.append(I_Nspikes_sprx_Na)

# Plotting
# May choose to specify the color again: # (can skip some of the complicated stuff here since I plot cell for cell)
color_pv_all   = [['#1f77b4','darkblue','#7f7f7f','b','tab:blue'],['#d62728','#ff7f0e','#d62728','m','tab:pink'],['#2ca02c','#8c564b','#9467bd','g','tab:brown']]
color_pv_sprx  = [['rebeccapurple','#e377c2','#8c564b','c','tab:purple'],['#9467bd','#bcbd22','#17becf','xkcd:sky blue','tab:orange'],['olive','darkgreen','mediumseagreen','tab:green','tab:olive']]

plotfolder = 'Comparemodels/All/'
plotname = plotfolder+'fI_varyE_all_shifts.png'

## avg and rms:
#fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9), (ax10, ax11, ax12)) = plt.subplots(4, 3, figsize=(20,25))

fig = plt.figure(figsize=(25,25),dpi=300)
fig.suptitle(r'Frequency $f$ vs $I$',fontsize=16)

gs = gridspec.GridSpec(5, 10)

ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
ax3 = plt.subplot(gs[0, 4:6])
ax4 = plt.subplot(gs[0, 6:8])
ax5 = plt.subplot(gs[0, 8:10])
ax6 = plt.subplot(gs[1, 0:2])
ax7 = plt.subplot(gs[1, 2:4])
ax8 = plt.subplot(gs[1, 4:6])
ax9 = plt.subplot(gs[1, 6:8])
ax10 = plt.subplot(gs[1, 8:10])
ax11 = plt.subplot(gs[2, 0:2])
ax12 = plt.subplot(gs[2, 2:4])
ax13 = plt.subplot(gs[2, 4:6])
ax14 = plt.subplot(gs[2, 6:8])
ax15 = plt.subplot(gs[2, 8:10])
ax16 = plt.subplot(gs[3, 2:4])
ax17 = plt.subplot(gs[3, 4:6])
ax18 = plt.subplot(gs[3, 6:8])

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
ax13.set_title(r'M',loc='left',fontsize=18)
ax14.set_title(r'N',loc='left',fontsize=18)
ax15.set_title(r'O',loc='left',fontsize=18)
ax16.set_title(r'P',loc='left',fontsize=18)
ax17.set_title(r'Q',loc='left',fontsize=18)
ax18.set_title(r'R',loc='left',fontsize=18)

# Add legends! #I could probably have looped this... ax1., ax2. etc. makes it hard.
ax1.set_title(r'Vary $E_{\mathregular{Na}}$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Na[0]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Na[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Na[0]
#Nspikes_sprx_this         = Nspikes_sprx_all_Na[0]
for j in range(NENa):
    ### Everywhere:
    ax1.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s=%i mV' % (plotlabelNa,varyENa[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax1.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %s=%i mV' % (plotlabelNa,varyENa[j]), linewidth=mylinewidth)
ax1.set_xlabel('$I$ (nA)',fontsize=14)
ax1.set_ylabel('$f$ (Hz)',fontsize=14)
ax1.legend(loc='upper left',ncol=1)
#ax1.legend(bbox_to_anchor =(0.5,-0.27), loc='lower center',ncol=3)

ax2.set_title(r'Vary $E_{\mathregular{Na}}$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Na[1]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Na[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Na[1]
#Nspikes_sprx_this         = Nspikes_sprx_all_Na[1]
for j in range(NENa):
    ### Everywhere:
    ax2.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s=%.2f mV' % (plotlabelNa,varyENa[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax2.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.2f mV' % (plotlabelNa,varyENa[j]), linewidth=mylinewidth)
ax2.set_xlabel('$I$ (nA)',fontsize=14)
ax2.set_ylabel('$f$ (Hz)',fontsize=14)

ax3.set_title(r'Vary $E_{\mathregular{Na}}$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Na[2]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Na[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Na[2]
#Nspikes_sprx_this         = Nspikes_sprx_all_Na[2]
for j in range(NENa):
    ### Everywhere:
    ax3.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s=%.2f mV' % (plotlabelNa,varyENa[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax3.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.2f mV' % (plotlabelNa,varyENa[j]), linewidth=mylinewidth)
ax3.set_xlabel('$I$ (nA)',fontsize=14)
ax3.set_ylabel('$f$ (Hz)',fontsize=14)

### Do the same for K: #################################################################

NEK = len(varyEK)

i_master_everywhere_all_K   = []
Nspikes_everywhere_all_K    = []
I_Nspikes_everywhere_all_K  = []

i_master_sprx_all_K         = []
Nspikes_sprx_all_K          = []
I_Nspikes_sprx_all_K        = []

for testmodel in testmodels:
    i_master_everywhere_K   = []
    Nspikes_everywhere_K    = []
    I_Nspikes_everywhere_K  = []
    
    i_master_sprx_K         = []
    Nspikes_sprx_K          = []
    I_Nspikes_sprx_K        = []

    infolder      = 'Allen_test_changecapacitance/figures/%i/' % (testmodel)
    vrestfolder   = infolder 
    for E in varyEK:
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
    
    
        ############# SOMAPROX ##################################
        '''
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


ax6.set_title(r'Vary $E_{\mathregular{K}}$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_K[0]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_K[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_K[0]
#Nspikes_sprx_this         = Nspikes_sprx_all_K[0]
for j in range(NEK):
    ### Everywhere:
    ax6.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s=%i mV' % (plotlabelK,varyEK[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax6.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %s=%i mV' % (plotlabelK,varyEK[j]), linewidth=mylinewidth)
ax6.set_xlabel('$I$ (nA)',fontsize=14)
ax6.set_ylabel('$f$ (Hz)',fontsize=14)
ax6.legend(loc='lower right',ncol=1)
#ax6.legend(bbox_to_anchor =(0.5,-0.27), loc='lower center',ncol=3)

ax7.set_title(r'Vary $E_{\mathregular{K}}$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_K[1]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_K[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_K[1]
#Nspikes_sprx_this         = Nspikes_sprx_all_K[1]
for j in range(NEK):
    ### Everywhere:
    ax7.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s=%.2f mV' % (plotlabelK,varyEK[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax7.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.2f mV' % (plotlabelK,varyEK[j]), linewidth=mylinewidth)
ax7.set_xlabel('$I$ (nA)',fontsize=14)
ax7.set_ylabel('$f$ (Hz)',fontsize=14)

ax8.set_title(r'Vary $E_{\mathregular{K}}$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_K[2]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_K[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_K[2]
#Nspikes_sprx_this         = Nspikes_sprx_all_K[2]
for j in range(NEK):
    ### Everywhere:
    ax8.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s=%.2f mV' % (plotlabelK,varyEK[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax8.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.2f mV' % (plotlabelK,varyEK[j]), linewidth=mylinewidth)
ax8.set_xlabel('$I$ (nA)',fontsize=14)
ax8.set_ylabel('$f$ (Hz)',fontsize=14)

################ Do the same for Epas: #####################################
NEpas = len(varyEpas)

i_master_everywhere_all_pas   = []
Nspikes_everywhere_all_pas    = []
I_Nspikes_everywhere_all_pas  = []

i_master_sprx_all_pas         = []
Nspikes_sprx_all_pas          = []
I_Nspikes_sprx_all_pas        = []

for testmodel in testmodels:
    i_master_everywhere_pas   = []
    Nspikes_everywhere_pas    = []
    I_Nspikes_everywhere_pas  = []
    
    i_master_sprx_pas         = []
    Nspikes_sprx_pas          = []
    I_Nspikes_sprx_pas        = []

    infolder      = 'Allen_test_changecapacitance/figures/%i/' % (testmodel)
    vrestfolder   = infolder 
    for E in varyEpas:
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
    
    #Nspikes_sprx_all_pas.append(Nspikes_sprx_pas)
    #I_Nspikes_sprx_all_pas.append(I_Nspikes_sprx_pas)

# Plotting

ax11.set_title(r'Vary $E_{\mathregular{pas}}$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_pas[0]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_pas[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_pas[0]
#Nspikes_sprx_this         = Nspikes_sprx_all_pas[0]
for j in range(NEpas):
    ### Everywhere:
    ax11.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s mV' % plotlabelEpas[j], linewidth=mylinewidth)
    ### Somaprox:
    #ax11.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %s mV' % plotlabelEpas[j], linewidth=mylinewidth)
ax11.set_xlabel('$I$ (nA)',fontsize=14)
ax11.set_ylabel('$f$ (Hz)',fontsize=14)
ax11.legend(loc='upper left',ncol=1)
#ax11.legend(bbox_to_anchor =(0.5,-0.27), loc='lower center',ncol=3)

ax12.set_title(r'Vary $E_{\mathregular{pas}}$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_pas[1]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_pas[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_pas[1]
#Nspikes_sprx_this         = Nspikes_sprx_all_pas[1]
for j in range(NEpas):
    ### Everywhere:
    ax12.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s mV' % plotlabelEpas[j], linewidth=mylinewidth)
    ### Somaprox:
    #ax12.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s mV' % plotlabelEpas[j], linewidth=mylinewidth)
ax12.set_xlabel('$I$ (nA)',fontsize=14)
ax12.set_ylabel('$f$ (Hz)',fontsize=14)

ax13.set_title(r'Vary $E_{\mathregular{pas}}$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_pas[2]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_pas[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_pas[2]
#Nspikes_sprx_this         = Nspikes_sprx_all_pas[2]
for j in range(NEpas):
    ### Everywhere:
    ax13.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s mV' % plotlabelEpas[j], linewidth=mylinewidth)
    ### Somaprox:
    #ax13.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s mV' % plotlabelEpas[j], linewidth=mylinewidth)
ax13.set_xlabel('$I$ (nA)',fontsize=14)
ax13.set_ylabel('$f$ (Hz)',fontsize=14)

############################ Ca ############################
R = 8.314   # JK-1mol-1
F = 9.648e4 # Cmol-1
T = 307.15  # K
prefactor = 1000*R*T/(2.0*F)
cai0 = 1e-4 # mM
NECa = len(varycao)
ECa0 = np.zeros(NECa)

#print('prefactor:',prefactor)
#print('-----------------------------------')
for i in range(NECa):
    #print('varycao[i]:',varycao[i])
    #print('varycao[i]/cai0:',varycao[i]/cai0)
    ECa0[i] = prefactor*np.log(varycao[i]/cai0)
    #print('np.log(varycao[i]/cai0):',np.log(varycao[i]/cai0))
    #print('ECa0[i]:',ECa0[i])
    #print('-----------------------------------')

i_master_everywhere_all_Ca   = []
Nspikes_everywhere_all_Ca    = []
I_Nspikes_everywhere_all_Ca  = []

i_master_sprx_all_Ca         = []
Nspikes_sprx_all_Ca          = []
I_Nspikes_sprx_all_Ca        = []

for testmodel in testmodels:
    i_master_everywhere_Ca   = []
    Nspikes_everywhere_Ca    = []
    I_Nspikes_everywhere_Ca  = []
    
    i_master_sprx_Ca         = []
    Nspikes_sprx_Ca          = []
    I_Nspikes_sprx_Ca        = []

    infolder      = 'Allen_test_changecapacitance/figures/%i/' % (testmodel)
    vrestfolder   = infolder 
    for cao in varycao:
        Nspikes   = []
        I_Nspikes = []
        
        namestringcao = namestringfirstcao+str(cao)
        
        # Set names
        infilename_Nspikes = infolder+'%s_%i_cmfall'%(namestringcao,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
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
    
    
        ############# SOMAPROX ##################################
        #Nspikes   = []
        #I_Nspikes = []
        
        # Set names
        #infilename_Nspikes = infolder+'%s_%i_cmfsprx'%(namestringcao,testmodel)+str(cm)+'_idur%i_varyiamp'% idur+'_manual_Nspikes_vs_I.txt'
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
        
        #Nspikes_sprx_pas.append(Nspikes)
        #I_Nspikes_sprx_pas.append(I_Nspikes)
    Nspikes_everywhere_all_Ca.append(Nspikes_everywhere_Ca)
    I_Nspikes_everywhere_all_Ca.append(I_Nspikes_everywhere_Ca)
    
    #Nspikes_sprx_all_pas.append(Nspikes_sprx_pas)
    #I_Nspikes_sprx_all_pas.append(I_Nspikes_sprx_pas)

# Plotting

ax16.set_title(r'Vary $E_{\mathregular{Ca}}(t=0)$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Ca[0]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Ca[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Ca[0]
#Nspikes_sprx_this         = Nspikes_sprx_all_Ca[0]
for j in range(NECa):
    ### Everywhere:
    ax16.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s=%.1f mV' % (plotlabelCa,ECa0[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax16.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %s=%.1f mV' % (plotlabelCa,ECa0[j]), linewidth=mylinewidth)
ax16.set_xlabel('$I$ (nA)',fontsize=14)
ax16.set_ylabel('$f$ (Hz)',fontsize=14)
ax16.legend(loc='upper left',ncol=1)

ax17.set_title(r'Vary $E_{\mathregular{Ca}}(t=0)$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Ca[1]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Ca[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Ca[1]
#Nspikes_sprx_this         = Nspikes_sprx_all_Ca[1]
for j in range(NECa):
    ### Everywhere:
    ax17.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s=%.1f mV' % (plotlabelCa,ECa0[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax17.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.1f mV % (plotlabelCa,ECa0[j]), linewidth=mylinewidth)
ax17.set_xlabel('$I$ (nA)',fontsize=14)
ax17.set_ylabel('$f$ (Hz)',fontsize=14)

ax18.set_title(r'Vary $E_{\mathregular{Ca}}(t=0)$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Ca[2]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Ca[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Ca[2]
#Nspikes_sprx_this         = Nspikes_sprx_all_Ca[2]
for j in range(NECa):
    ### Everywhere:
    ax18.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s=%.1f mV' % (plotlabelCa,ECa0[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax18.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.1f mV' % (plotlabelCa,ECa0[j]), linewidth=mylinewidth)
ax18.set_xlabel('$I$ (nA)',fontsize=14)
ax18.set_ylabel('$f$ (Hz)',fontsize=14)
################ SOMA AND BAS! ##############################################################

mylinewidth = 2

cm         = 1.0
idur       = 1000 #100 # ms
idelay     = 100
v_init     = -86.5 # mV
Ra         = 150
somasize   = 10 # 15 # 
dendlen    = 1000
denddiam   = 1
nsegments  = 200 
spikedurat = -20

varymech = 'Na'
varyE_bool = True
varyE = [45,50,60,70,80]
folderstring = 'VaryNa/'
plotstring = '_vary' 
plotstring = plotstring + '_Na'
labelstring = r'$E_{\mathregular{Na}}$'

NENa_BAS = len(varyE)
NENa_OC  = NENa_BAS

i_master = [0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5]
NI  = len(i_master)

i_master_everywhere  = []
Nspikes_everywhere   = []
I_Nspikes_everywhere = []
    
i_master_sprx  = []
Nspikes_sprx   = []
I_Nspikes_sprx = []
    
i_master_onecomp  = []
Nspikes_onecomp   = []
I_Nspikes_onecomp = []

for E in varyE:
    changestring =''
    changestring = changestring+'_E'+str(E)+'_gdf'

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
    infolder = 'Ball-and-stick models/BAS_somaHH_dendpassive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+folderstring
    infilename_Nspikes = infolder+'basHHdpas_cmall'+str(cm)+'_idur%i_varyiamp'% (idur) +'_E'+str(E)+'_manual_Nspikes_vs_I.txt'
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
    '''
    Nspikes         = []
    I_Nspikes       = []
        
    # Set names # FIX BAS NAMES!
    infolder = 'Ball-and-stick models/BAS_somaHH_dendpassive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+folderstring
    infilename_Nspikes = infolder+'basHHdpas_csprx'+str(cm)+'_idur%i_varyiamp'% (idur) +'_E'+str(E)+'_manual_Nspikes_vs_I.txt'
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

    print('i_master:',i_master) 
    print('Nspikes:',Nspikes)
    '''

    #################### Soma only, Hodgkin-Huxley ########################################
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

    # Default HH values:
    ena = E
    ek = -77
    el_hh = -54.3
    gnabar_hh = 0.12
    gkbar_hh = 0.036
    gl_hh = 0.0003

    infolder_shh = 'Somaonly/Results/IStim/Soma%i/' % somasize
    hhstring = '_ena'+str(ena)+'_ek'+str(ek)+'_el'+str(el_hh)+'_gnabar'+str(gnabar_hh)+'_gkbar'+str(gkbar_hh)+'_gl'+str(gl_hh)
    infilename_Nspikes = infolder_shh+'somaonlyHH_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+hhstring+'_Nspikes_vs_I.txt'
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

# Save in order to find differences later (neater this way):
Nspikes_everywhere_Na   = Nspikes_everywhere
I_Nspikes_everywhere_Na = I_Nspikes_everywhere

Nspikes_onecomp_Na   = Nspikes_onecomp
I_Nspikes_onecomp_Na = I_Nspikes_onecomp

# Plotting

color_bashhall  = ['#1f77b4','#d62728']
color_bashhsprx = ['#ff7f0e','#9467bd']
color_somahh    = ['#2ca02c','#8c564b']

## avg and rms:

plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('legend', fontsize=16)

print('varyE:', varyE)
print('len(varyE):',len(varyE))
fig.suptitle(r'Frequency $f$ vs $I$',fontsize=18)
ax4.set_title('Ball-and-stick model, changing %s' % labelstring,fontsize=16)
ax5.set_title('One compartment model, changing %s' % labelstring,fontsize=16)
for i in range(len(varyE)):
    print('i:',i)
    ax4.plot(I_Nspikes_everywhere[i], Nspikes_everywhere[i],label=r'$E_{\mathregular{Na}}$ = %s mV' % (str(varyE[i])), linewidth=mylinewidth)
    #ax4.plot(I_Nspikes_sprx[i], Nspikes_sprx[i],label=r'$E_{\mathregular{Na}}$ = %s, sprx' % (str(varyE[i])), linewidth=mylinewidth,ls='--')
    ax5.plot(I_Nspikes_onecomp[i], Nspikes_onecomp[i],label=r'$E_{\mathregular{Na}}$ = %s mV' % (str(varyE[i])), linewidth=mylinewidth)
    
ax1.set_xlabel('$I$ (nA)',fontsize=16)
ax1.set_ylabel('$f$ (Hz)',fontsize=16)
ax2.set_xlabel('$I$ (nA)',fontsize=16)
ax2.set_ylabel('$f$ (Hz)',fontsize=16)
ax3.set_xlabel('$I$ (nA)',fontsize=16)
ax3.set_ylabel('$f$ (Hz)',fontsize=16)
ax4.set_xlabel('$I$ (nA)',fontsize=16)
ax4.set_ylabel('$f$ (Hz)',fontsize=16)
ax5.set_xlabel('$I$ (nA)',fontsize=16)
ax5.set_ylabel('$f$ (Hz)',fontsize=16)
ax6.set_xlabel('$I$ (nA)',fontsize=16)
ax6.set_ylabel('$f$ (Hz)',fontsize=16)
ax4.legend(loc='lower right')
ax5.legend(loc='lower right')

#################### K ########################
varymech = 'K'
varyE_bool = True
varyE = [-70,-73,-77,-87,-97]
folderstring = 'VaryK/'
plotstring = '_vary' 
plotstring = plotstring + '_K'
labelstring = r'$E_{\mathregular{K}}$'
NEK_BAS = len(varyE)
NEK_OC  = NEK_BAS

i_master_everywhere  = []
Nspikes_everywhere   = []
I_Nspikes_everywhere = []
    
i_master_sprx  = []
Nspikes_sprx   = []
I_Nspikes_sprx = []
    
i_master_onecomp  = []
Nspikes_onecomp   = []
I_Nspikes_onecomp = []

for E in varyE:
    changestring =''
    changestring = changestring+'_E'+str(E)+'_gdf'

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
    infolder = 'Ball-and-stick models/BAS_somaHH_dendpassive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+folderstring
    infilename_Nspikes = infolder+'basHHdpas_cmall'+str(cm)+'_idur%i_varyiamp'% (idur) +'_E'+str(E)+'_manual_Nspikes_vs_I.txt'
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
    '''
    Nspikes         = []
    I_Nspikes       = []
        
    # Set names
    infolder = 'Ball-and-stick models/BAS_somaHH_dendpassive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+folderstring
    infilename_Nspikes = infolder+'basHHdpas_csprx'+str(cm)+'_idur%i_varyiamp'% (idur) +'_E'+str(E)+'_manual_Nspikes_vs_I.txt'
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

    print('i_master:',i_master) 
    print('Nspikes:',Nspikes)
    '''

    #################### Soma only, Hodgkin-Huxley ########################################
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

    # Default HH values:
    ena = 50
    ek = E
    el_hh = -54.3
    gnabar_hh = 0.12
    gkbar_hh = 0.036
    gl_hh = 0.0003

    infolder_shh = 'Somaonly/Results/IStim/Soma%i/' % somasize
    hhstring = '_ena'+str(ena)+'_ek'+str(ek)+'_el'+str(el_hh)+'_gnabar'+str(gnabar_hh)+'_gkbar'+str(gkbar_hh)+'_gl'+str(gl_hh)
    infilename_Nspikes = infolder_shh+'somaonlyHH_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+hhstring+'_Nspikes_vs_I.txt'
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

# Save in order to find differences later:
Nspikes_everywhere_K   = Nspikes_everywhere
I_Nspikes_everywhere_K = I_Nspikes_everywhere

Nspikes_onecomp_K   = Nspikes_onecomp
I_Nspikes_onecomp_K = I_Nspikes_onecomp

# Plotting

ax9.set_title('Ball-and-stick model, changing %s' % labelstring,fontsize=16)
ax10.set_title('One compartment model, changing %s' % labelstring,fontsize=16)
for i in range(len(varyE)):
    ax9.plot(I_Nspikes_everywhere[i], Nspikes_everywhere[i],label=r'$E_{\mathregular{K}}$ = %s mV' % (str(varyE[i])), linewidth=mylinewidth)
    #ax9.plot(I_Nspikes_sprx[i], Nspikes_sprx[i],label=r'$E_{\mathregular{K}}$ = %s mV, sprx' % (str(varyE[i])), linewidth=mylinewidth,ls='--')
    ax10.plot(I_Nspikes_onecomp[i], Nspikes_onecomp[i],label=r'$E_{\mathregular{K}}$ = %s mV' % (str(varyE[i])), linewidth=mylinewidth)
    
ax9.legend(loc='lower right')
ax10.legend(loc='lower right')

##################### Leak #############################
varymech = 'leak'
varyE_bool = True
varyE = [-34.3,-44.3,-54.3,-64.3,-74.3]
folderstring = 'VaryL/'
plotstring = '_vary' 
plotstring = plotstring + '_L'
labelstring = r'$E_L$'
NEl_BAS = len(varyE)
NEl_OC  = NEl_BAS

i_master_everywhere  = []
Nspikes_everywhere   = []
I_Nspikes_everywhere = []
    
i_master_sprx  = []
Nspikes_sprx   = []
I_Nspikes_sprx = []
    
i_master_onecomp  = []
Nspikes_onecomp   = []
I_Nspikes_onecomp = []

for E in varyE:
    changestring =''
    changestring = changestring+'_E'+str(E)+'_gdf'

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
    infolder = 'Ball-and-stick models/BAS_somaHH_dendpassive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+folderstring
    infilename_Nspikes = infolder+'basHHdpas_cmall'+str(cm)+'_idur%i_varyiamp'% (idur) +'_E'+str(E)+'_manual_Nspikes_vs_I.txt'
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
    '''
    Nspikes         = []
    I_Nspikes       = []
        
    # Set names
    infolder = 'Ball-and-stick models/BAS_somaHH_dendpassive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+folderstring
    infilename_Nspikes = infolder+'basHHdpas_csprx'+str(cm)+'_idur%i_varyiamp'% (idur) +'_E'+str(E)+'_manual_Nspikes_vs_I.txt'
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

    print('i_master:',i_master) 
    print('Nspikes:',Nspikes)
    '''
    
    #################### Soma only, Hodgkin-Huxley ########################################
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

    # Default HH values:
    ena = 50
    ek = -77
    el_hh = E
    gnabar_hh = 0.12
    gkbar_hh = 0.036
    gl_hh = 0.0003

    infolder_shh = 'Somaonly/Results/IStim/Soma%i/' % somasize
    hhstring = '_ena'+str(ena)+'_ek'+str(ek)+'_el'+str(el_hh)+'_gnabar'+str(gnabar_hh)+'_gkbar'+str(gkbar_hh)+'_gl'+str(gl_hh)
    infilename_Nspikes = infolder_shh+'somaonlyHH_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+hhstring+'_Nspikes_vs_I.txt'
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

Nspikes_everywhere_l   = Nspikes_everywhere
I_Nspikes_everywhere_l = I_Nspikes_everywhere

Nspikes_onecomp_l   = Nspikes_onecomp
I_Nspikes_onecomp_l = I_Nspikes_onecomp

# Plotting

ax14.set_title('Ball-and-stick model, changing %s' % labelstring,fontsize=16)
ax15.set_title('One compartment model, changing %s' % labelstring,fontsize=16)
for i in range(len(varyE)):
    ax14.plot(I_Nspikes_everywhere[i], Nspikes_everywhere[i],label=r'$E_L$ = %s mV' % (str(varyE[i])), linewidth=mylinewidth)
    #ax14.plot(I_Nspikes_sprx[i], Nspikes_sprx[i],label=r'$E_L$ = %s mV, sprx' % (str(varyE[i])), linewidth=mylinewidth,ls='--')
    ax15.plot(I_Nspikes_onecomp[i], Nspikes_onecomp[i],label=r'$E_L$ = %s mV' % (str(varyE[i])), linewidth=mylinewidth)

ax14.legend(loc='lower right')
ax15.legend(loc='lower right')

############################ Differences ################################################

idefNa       = 2 # Index that holds the defaul value of ENa
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

idefCa       = 2 # Index that holds the defaul value of ECa
lastdiffs_Ca = np.zeros((Nmodels,NECa))
for j in range(Nmodels):
    I_Nspikes_everywhere_all = I_Nspikes_everywhere_all_Ca[j]
    Nspikes_everywhere_all = Nspikes_everywhere_all_Ca[j]
    for i in range(NECa):
        im_Ca, rd1_Ca, rd2_Ca = reldiffs(Nspikes_everywhere_all[i],Nspikes_everywhere_all[idefNa],I_Nspikes_everywhere_all[i],I_Nspikes_everywhere_all[idefCa]) 
        if len(rd1_Ca)>0:
            lastdiffs_Ca[j,i] = rd1_Ca[-1]
        else:
            lastdiffs_Ca[j,i] = 0

idefpas       = 2 # Index that holds the defaul value of Epas
lastdiffs_pas = np.zeros((Nmodels,NEpas))
for j in range(Nmodels):
    I_Nspikes_everywhere_all = I_Nspikes_everywhere_all_pas[j]
    Nspikes_everywhere_all = Nspikes_everywhere_all_pas[j]
    for i in range(NEpas):
        im_pas, rd1_pas, rd2_pas = reldiffs(Nspikes_everywhere_all[i],Nspikes_everywhere_all[idefNa],I_Nspikes_everywhere_all[i],I_Nspikes_everywhere_all[idefpas]) 
        if len(rd1_pas)>0:
            lastdiffs_pas[j,i] = rd1_pas[-1]
        else:
            lastdiffs_pas[j,i] = 0

#### Onecomp & BAS ######
NEl_BAS  = len(Nspikes_everywhere_l)
NEl_OC   = len(Nspikes_onecomp_l)

# Na

idefNa_BAS       = 1 # Index that holds the defaul value of ENa
lastdiffs_Na_BAS = np.zeros(NENa_BAS)
for i in range(NENa_BAS):
    im_Na_BAS, rd1_Na_BAS, rd2_Na_BAS = reldiffs(Nspikes_everywhere_Na[i],Nspikes_everywhere_Na[idefNa_BAS],I_Nspikes_everywhere_Na[i],I_Nspikes_everywhere_Na[idefNa_BAS]) 
    if len(rd1_Na_BAS)>0:
        lastdiffs_Na_BAS[i] = rd1_Na_BAS[-1]
    else:
        lastdiffs_Na_BAS[i] = 0

idefNa_OC       = 1 # Index that holds the defaul value of ENa
lastdiffs_Na_OC = np.zeros(NENa_OC)
for i in range(NENa_OC):
    im_Na_OC, rd1_Na_OC, rd2_Na_OC = reldiffs(Nspikes_onecomp_Na[i],Nspikes_onecomp_Na[idefNa_OC],I_Nspikes_onecomp_Na[i],I_Nspikes_onecomp_Na[idefNa_OC]) 
    if len(rd1_Na_OC)>0:
        lastdiffs_Na_OC[i] = rd1_Na_OC[-1]
    else:
        lastdiffs_Na_OC[i] = 0

# K

idefK_BAS       = 2 # Index that holds the defaul value of EK
lastdiffs_K_BAS = np.zeros(NEK_BAS)
for i in range(NEK_BAS):
    im_K_BAS, rd1_K_BAS, rd2_K_BAS = reldiffs(Nspikes_everywhere_K[i],Nspikes_everywhere_K[idefK_BAS],I_Nspikes_everywhere_K[i],I_Nspikes_everywhere_K[idefK_BAS]) 
    if len(rd1_K_BAS)>0:
        lastdiffs_K_BAS[i] = rd1_K_BAS[-1]
    else:
        lastdiffs_K_BAS[i] = 0

idefK_OC       = 2 # Index that holds the defaul value of EK
lastdiffs_K_OC = np.zeros(NEK_OC)
for i in range(NEK_OC):
    im_K_OC, rd1_K_OC, rd2_K_OC = reldiffs(Nspikes_onecomp_K[i],Nspikes_onecomp_K[idefK_BAS],I_Nspikes_onecomp_K[i],I_Nspikes_onecomp_K[idefK_BAS]) 
    if len(rd1_K_OC)>0:
        lastdiffs_K_OC[i] = rd1_K_OC[-1]
    else:
        lastdiffs_K_OC[i] = 0

# Leak

idefl_BAS       = 2 # Index that holds the defaul value of El
lastdiffs_l_BAS = np.zeros(NEl_BAS)
for i in range(NEl_BAS):
    im_l_BAS, rd1_l_BAS, rd2_l_BAS = reldiffs(Nspikes_everywhere_l[i],Nspikes_everywhere_l[idefl_BAS],I_Nspikes_everywhere_l[i],I_Nspikes_everywhere_l[idefl_BAS]) 
    if len(rd1_l_BAS)>0:
        lastdiffs_l_BAS[i] = rd1_l_BAS[-1]
    else:
        lastdiffs_l_BAS[i] = 0

idefl_OC       = 2 # Index that holds the defaul value of El
lastdiffs_l_OC = np.zeros(NEl_OC)
for i in range(NEl_OC):
    im_l_OC, rd1_l_OC, rd2_l_OC = reldiffs(Nspikes_onecomp_l[i],Nspikes_onecomp_l[idefl_BAS],I_Nspikes_onecomp_l[i],I_Nspikes_onecomp_l[idefl_BAS]) 
    if len(rd1_l_OC)>0:
        lastdiffs_l_OC[i] = rd1_l_OC[-1]
    else:
        lastdiffs_l_OC[i] = 0


ax19 = plt.subplot(gs[4, 1:3])
ax19.set_title(r'S',loc='left',fontsize=18)

############ FIND MORE INDICES AND ADD MORE STUFF!!! ##############
barWidth = 0.25
change_ena2   = [lastdiffs_Na[0,3],lastdiffs_Na[1,3],lastdiffs_Na[2,3],lastdiffs_Na_BAS[2],lastdiffs_Na_OC[2]]
change_ena3 = [lastdiffs_Na[0,4],lastdiffs_Na[1,4],lastdiffs_Na[2,4],lastdiffs_Na_BAS[3],lastdiffs_Na_OC[3]]

br1 = np.arange(len(change_ena2))
br2 = [x+barWidth for x in br1]
brcenter = [x+0.5*barWidth for x in br1]

# How do I fix the difference between values of Allen and HH? Run more sim.s for HH?
ax19.bar(br1, change_ena2, width=barWidth, label=r'$E_{\mathregular{Na}}$+10 mV')#=%s'%str(varyENa[3]))
ax19.bar(br2, change_ena3, width=barWidth, label=r'$E_{\mathregular{Na}}$+20 mV')#=%s'%str(varyENa[4]))
ax19.axhline(y=0,color='k',linestyle='-',linewidth=1.0)
plt.xticks(brcenter, ['437','965','407','OC','BAS'])
ax19.set_xlabel('Model',fontsize=14)
ax19.set_ylabel(r'Relative difference at max. current',fontsize=14)
ax19.set_title(r'Difference from original $E_{\mathregular{Na}}$',fontsize=15)
ax19.legend(loc='lower left')

ax20 = plt.subplot(gs[4, 3:5])
ax20.set_title(r'T',loc='left',fontsize=18)
change_K2   = [lastdiffs_K[0,1],lastdiffs_K[1,1],lastdiffs_K[2,1],lastdiffs_K_BAS[3],lastdiffs_K_OC[3]]
change_K3   = [lastdiffs_K[0,0],lastdiffs_K[1,0],lastdiffs_K[2,0],lastdiffs_K_BAS[4],lastdiffs_K_OC[4]]

br1 = np.arange(len(change_K2))
br2 = [x+barWidth for x in br1]
brcenter = [x+0.5*barWidth for x in br1]

# How do I fix the difference between values of Allen and HH? Run more sim.s for HH?
ax20.bar(br1, change_K2, width=barWidth, label=r'$E_{\mathregular{K}}$-10 mV')#'%str(varyEK[3]))
ax20.bar(br2, change_K3, width=barWidth, label=r'$E_{\mathregular{K}}$-20 mV')#=%s'%str(varyEK[4]))
ax20.axhline(y=0,color='k',linestyle='-',linewidth=1.0)
plt.xticks(brcenter, ['437','965','407','OC','BAS'])
ax20.set_xlabel('Model',fontsize=14)
ax20.set_ylabel(r'Relative difference at max. current',fontsize=14)
ax20.set_title(r'Difference from original $E_{\mathregular{K}}$',fontsize=15)
ax20.legend(loc='upper right')

ax21 = plt.subplot(gs[4, 5:7])
ax21.set_title(r'U',loc='left',fontsize=18)
change_pas2   = [lastdiffs_pas[0,1],lastdiffs_pas[1,1],lastdiffs_pas[2,1],lastdiffs_l_BAS[3],lastdiffs_l_OC[3]]
change_pas3   = [lastdiffs_pas[0,0],lastdiffs_pas[1,0],lastdiffs_pas[2,0],lastdiffs_l_BAS[4],lastdiffs_l_OC[4]]

br1 = np.arange(len(change_pas2))
br2 = [x+barWidth for x in br1]
brcenter = [x+0.5*barWidth for x in br1]

# How do I fix the difference between values of Allen and HH? Run more sim.s for HH?
ax21.bar(br1, change_pas2, width=barWidth, label=r'$E_{\mathregular{L}}$-10 mV')#'%str(varyK[3]))
ax21.bar(br2, change_pas3, width=barWidth, label=r'$E_{\mathregular{L}}$-20 mV')#=%s'%str(varyK[4]))
ax21.axhline(y=0,color='k',linestyle='-',linewidth=1.0)
plt.xticks(brcenter, ['437','965','407','OC','BAS'])
ax21.set_xlabel('Model',fontsize=14)
ax21.set_ylabel(r'Relative difference at max. current',fontsize=14)
ax21.set_title(r'Difference from original $E_{\mathregular{L}}$',fontsize=15)
ax21.legend(loc='upper left')


ax22 = plt.subplot(gs[4, 7:9])
ax22.set_title(r'V',loc='left',fontsize=18)
change_ca2   = [lastdiffs_Ca[0,3],lastdiffs_Ca[1,3],lastdiffs_Ca[2,3]]
change_ca3   = [lastdiffs_Ca[0,4],lastdiffs_Ca[1,4],lastdiffs_Ca[2,4]]

br1 = np.arange(len(change_ca2))
br2 = [x+barWidth for x in br1]
brcenter = [x+0.5*barWidth for x in br1]

# How do I fix the difference between values of Allen and HH? Run more sim.s for HH?
ax22.bar(br1, change_ca2, width=barWidth, label=r'$E_{\mathregular{Ca}}$=%.1f mV'%ECa0[3])
ax22.bar(br2, change_ca3, width=barWidth, label=r'$E_{\mathregular{Ca}}$=%.1f mV'%ECa0[4])
ax22.axhline(y=0,color='k',linestyle='-',linewidth=1.0)
plt.xticks(brcenter, ['437','965','407'])
ax22.set_xlabel('Model',fontsize=14)
ax22.set_ylabel(r'Relative difference at max. current',fontsize=14)
ax22.set_title(r'Difference from original $E_{\mathregular{Ca}}$',fontsize=15)
ax22.legend(loc='upper right')

fig.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.savefig(plotname)

plt.show()
