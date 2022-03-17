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
varyENa = [40,50,53,60,70]
namestringfirstNa = namestringfirst + 'ENa'
plotlabelNa = r'$E_{\mathregular{Na}}$'
varyEK = [-120,-107,-100,-90,-80]
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

    infolder      = 'figures/%i/' % (testmodel)
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
    Nspikes_everywhere_all_Na.append(Nspikes_everywhere_Na)
    I_Nspikes_everywhere_all_Na.append(I_Nspikes_everywhere_Na)
    
    Nspikes_sprx_all_Na.append(Nspikes_sprx_Na)
    I_Nspikes_sprx_all_Na.append(I_Nspikes_sprx_Na)

# Plotting
# May choose to specify the color again: # (can skip some of the complicated stuff here since I plot cell for cell)
color_pv_all   = [['#1f77b4','darkblue','#7f7f7f','b','tab:blue'],['#d62728','#ff7f0e','#d62728','m','tab:pink'],['#2ca02c','#8c564b','#9467bd','g','tab:brown']]
color_pv_sprx  = [['rebeccapurple','#e377c2','#8c564b','c','tab:purple'],['#9467bd','#bcbd22','#17becf','xkcd:sky blue','tab:orange'],['olive','darkgreen','mediumseagreen','tab:green','tab:olive']]

plotfolder = 'figures/Comparemodels/'
plotname = plotfolder+'fI_varyE_AllenPV_withCa.png'

## avg and rms:
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9), (ax10, ax11, ax12)) = plt.subplots(4, 3, figsize=(22,20))
fig.suptitle(r'Frequency $f$ vs $I$',fontsize=16)

# Add legends! #I could probably have looped this... ax1., ax2. etc. makes it hard.
ax1.set_title(r'Vary $E_{\mathregular{Na}}$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Na[0]
I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Na[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Na[0]
Nspikes_sprx_this         = Nspikes_sprx_all_Na[0]
for j in range(NENa):
    ### Everywhere:
    ax1.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s=%i' % (plotlabelNa,varyENa[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax1.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %s=%i' % (plotlabelNa,varyENa[j]), linewidth=mylinewidth)
ax1.set_xlabel('$I$ (nA)',fontsize=14)
ax1.set_ylabel('$f$ (Hz)',fontsize=14)
ax1.legend(loc='upper left',ncol=1)
#ax1.legend(bbox_to_anchor =(0.5,-0.27), loc='lower center',ncol=3)

ax2.set_title(r'Vary $E_{\mathregular{Na}}$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Na[1]
I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Na[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Na[1]
Nspikes_sprx_this         = Nspikes_sprx_all_Na[1]
for j in range(NENa):
    ### Everywhere:
    ax2.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s=%.2f' % (plotlabelNa,varyENa[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax2.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.2f' % (plotlabelNa,varyENa[j]), linewidth=mylinewidth)
ax2.set_xlabel('$I$ (nA)',fontsize=14)
ax2.set_ylabel('$f$ (Hz)',fontsize=14)

ax3.set_title(r'Vary $E_{\mathregular{Na}}$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Na[2]
I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Na[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Na[2]
Nspikes_sprx_this         = Nspikes_sprx_all_Na[2]
for j in range(NENa):
    ### Everywhere:
    ax3.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s=%.2f' % (plotlabelNa,varyENa[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax3.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.2f' % (plotlabelNa,varyENa[j]), linewidth=mylinewidth)
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

    infolder      = 'figures/%i/' % (testmodel)
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
    Nspikes_everywhere_all_K.append(Nspikes_everywhere_K)
    I_Nspikes_everywhere_all_K.append(I_Nspikes_everywhere_K)
    
    Nspikes_sprx_all_K.append(Nspikes_sprx_K)
    I_Nspikes_sprx_all_K.append(I_Nspikes_sprx_K)


ax4.set_title(r'Vary $E_{\mathregular{K}}$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_K[0]
I_Nspikes_sprx_this       = I_Nspikes_sprx_all_K[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_K[0]
Nspikes_sprx_this         = Nspikes_sprx_all_K[0]
for j in range(NEK):
    ### Everywhere:
    ax4.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s=%i' % (plotlabelK,varyEK[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax4.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %s=%i' % (plotlabelK,varyEK[j]), linewidth=mylinewidth)
ax4.set_xlabel('$I$ (nA)',fontsize=14)
ax4.set_ylabel('$f$ (Hz)',fontsize=14)
ax4.legend(loc='lower right',ncol=1)
#ax4.legend(bbox_to_anchor =(0.5,-0.27), loc='lower center',ncol=3)

ax5.set_title(r'Vary $E_{\mathregular{K}}$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_K[1]
I_Nspikes_sprx_this       = I_Nspikes_sprx_all_K[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_K[1]
Nspikes_sprx_this         = Nspikes_sprx_all_K[1]
for j in range(NEK):
    ### Everywhere:
    ax5.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s=%.2f' % (plotlabelK,varyEK[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax5.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.2f' % (plotlabelK,varyEK[j]), linewidth=mylinewidth)
ax5.set_xlabel('$I$ (nA)',fontsize=14)
ax5.set_ylabel('$f$ (Hz)',fontsize=14)

ax6.set_title(r'Vary $E_{\mathregular{K}}$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_K[2]
I_Nspikes_sprx_this       = I_Nspikes_sprx_all_K[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_K[2]
Nspikes_sprx_this         = Nspikes_sprx_all_K[2]
for j in range(NEK):
    ### Everywhere:
    ax6.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s=%.2f' % (plotlabelK,varyEK[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax6.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.2f' % (plotlabelK,varyEK[j]), linewidth=mylinewidth)
ax6.set_xlabel('$I$ (nA)',fontsize=14)
ax6.set_ylabel('$f$ (Hz)',fontsize=14)

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

    infolder      = 'figures/%i/' % (testmodel)
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
    Nspikes_everywhere_all_pas.append(Nspikes_everywhere_pas)
    I_Nspikes_everywhere_all_pas.append(I_Nspikes_everywhere_pas)
    
    Nspikes_sprx_all_pas.append(Nspikes_sprx_pas)
    I_Nspikes_sprx_all_pas.append(I_Nspikes_sprx_pas)

# Plotting

ax7.set_title(r'Vary $E_{\mathregular{pas}}$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_pas[0]
I_Nspikes_sprx_this       = I_Nspikes_sprx_all_pas[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_pas[0]
Nspikes_sprx_this         = Nspikes_sprx_all_pas[0]
for j in range(NEpas):
    ### Everywhere:
    ax7.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s' % plotlabelEpas[j], linewidth=mylinewidth)
    ### Somaprox:
    #ax7.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %s' % plotlabelEpas[j], linewidth=mylinewidth)
ax7.set_xlabel('$I$ (nA)',fontsize=14)
ax7.set_ylabel('$f$ (Hz)',fontsize=14)
ax7.legend(loc='upper left',ncol=1)
#ax7.legend(bbox_to_anchor =(0.5,-0.27), loc='lower center',ncol=3)

ax8.set_title(r'Vary $E_{\mathregular{pas}}$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_pas[1]
I_Nspikes_sprx_this       = I_Nspikes_sprx_all_pas[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_pas[1]
Nspikes_sprx_this         = Nspikes_sprx_all_pas[1]
for j in range(NEpas):
    ### Everywhere:
    ax8.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s' % plotlabelEpas[j], linewidth=mylinewidth)
    ### Somaprox:
    #ax8.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s' % plotlabelEpas[j], linewidth=mylinewidth)
ax8.set_xlabel('$I$ (nA)',fontsize=14)
ax8.set_ylabel('$f$ (Hz)',fontsize=14)

ax9.set_title(r'Vary $E_{\mathregular{pas}}$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_pas[2]
I_Nspikes_sprx_this       = I_Nspikes_sprx_all_pas[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_pas[2]
Nspikes_sprx_this         = Nspikes_sprx_all_pas[2]
for j in range(NEpas):
    ### Everywhere:
    ax9.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s' % plotlabelEpas[j], linewidth=mylinewidth)
    ### Somaprox:
    #ax9.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s' % plotlabelEpas[j], linewidth=mylinewidth)
ax9.set_xlabel('$I$ (nA)',fontsize=14)
ax9.set_ylabel('$f$ (Hz)',fontsize=14)

############################ Ca ############################
R = 8.314   # JK-1mol-1
F = 9.648e4 # Cmol-1
T = 307.15  # K
prefactor = R*T/(2.0*F)
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

    infolder      = 'figures/%i/' % (testmodel)
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

ax10.set_title(r'Vary $E_{\mathregular{Ca}}(t=0)$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Ca[0]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Ca[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Ca[0]
#Nspikes_sprx_this         = Nspikes_sprx_all_Ca[0]
for j in range(NECa):
    ### Everywhere:
    ax10.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s==%.2f' % (plotlabelCa,ECa0[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax10.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %s=%.2f' % (plotlabelCa,ECa0[j]), linewidth=mylinewidth)
ax10.set_xlabel('$I$ (nA)',fontsize=14)
ax10.set_ylabel('$f$ (Hz)',fontsize=14)
ax10.legend(loc='upper left',ncol=1)

ax11.set_title(r'Vary $E_{\mathregular{Ca}}(t=0)$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Ca[1]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Ca[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Ca[1]
#Nspikes_sprx_this         = Nspikes_sprx_all_Ca[1]
for j in range(NECa):
    ### Everywhere:
    ax11.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s=%.2f' % (plotlabelCa,ECa0[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax11.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.2f' % (plotlabelCa,ECa0[j]), linewidth=mylinewidth)
ax11.set_xlabel('$I$ (nA)',fontsize=14)
ax11.set_ylabel('$f$ (Hz)',fontsize=14)

ax12.set_title(r'Vary $E_{\mathregular{Ca}}(t=0)$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_Ca[2]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_Ca[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_Ca[2]
#Nspikes_sprx_this         = Nspikes_sprx_all_Ca[2]
for j in range(NECa):
    ### Everywhere:
    ax12.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s=%.2f' % (plotlabelCa,ECa0[j]), linewidth=mylinewidth)
    ### Somaprox:
    #ax12.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s=%.2f' % (plotlabelCa,ECa0[j]), linewidth=mylinewidth)
ax12.set_xlabel('$I$ (nA)',fontsize=14)
ax12.set_ylabel('$f$ (Hz)',fontsize=14)

fig.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.savefig(plotname)

plt.show()
