import matplotlib.pyplot as plt
import numpy as np
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


fig = plt.figure(figsize=(30,15),dpi=300)
fig.suptitle(r'Frequency $f$ vs $I$',fontsize=16)

gs = gridspec.GridSpec(2, 10)

ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
ax3 = plt.subplot(gs[0, 4:6])
ax4 = plt.subplot(gs[0, 6:8])
ax5 = plt.subplot(gs[0, 8:10])
ax6 = plt.subplot(gs[1, 0:2])
ax7 = plt.subplot(gs[1, 2:4])
ax8 = plt.subplot(gs[1, 4:6])


ax1.set_title(r'A',loc='left',fontsize=18)
ax2.set_title(r'B',loc='left',fontsize=18)
ax3.set_title(r'C',loc='left',fontsize=18)
ax4.set_title(r'D',loc='left',fontsize=18)
ax5.set_title(r'E',loc='left',fontsize=18)
ax6.set_title(r'F',loc='left',fontsize=18)
ax7.set_title(r'G',loc='left',fontsize=18)
ax8.set_title(r'H',loc='left',fontsize=18)

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

plotstringHVA  = '_vary_CaHVA'
plotstringNaV  = '_vary_NaV'

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

    infolder      = 'Allen_test_changecapacitance/figures/%i/' % (testmodel)
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

plotfolder = 'Comparemodels/All/'

plotname = plotfolder+'fI_varyg_allmodels.png'

ax6.set_title(r'Vary $\bar{g}_{\mathregular{CaHVA}}$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gHVA[0]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gHVA[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gHVA[0]
#Nspikes_sprx_this         = Nspikes_sprx_all_gHVA[0]
for j in range(NgHVA):
    ### Everywhere:
    ax6.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygHVA[j],plotlabelgHVA), linewidth=mylinewidth)
    ### Somaprox:
    #ax6.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %.1f*%s' % (varygHVA[j],plotlabelgHVA), linewidth=mylinewidth)
ax6.set_xlabel('$I$ (nA)',fontsize=14)
ax6.set_ylabel('$f$ (Hz)',fontsize=14)
ax6.legend(loc='upper left',ncol=1)

ax7.set_title(r'Vary $\bar{g}_{\mathregular{CaHVA}}$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gHVA[1]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gHVA[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gHVA[1]
#Nspikes_sprx_this         = Nspikes_sprx_all_gHVA[1]
for j in range(NgHVA):
    ### Everywhere:
    ax7.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygHVA[j],plotlabelgHVA), linewidth=mylinewidth)
    ### Somaprox:
    #ax7.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %.1f*%s' % (varygHVA[j],plotlabelgHVA), linewidth=mylinewidth)
ax7.set_xlabel('$I$ (nA)',fontsize=14)
ax7.set_ylabel('$f$ (Hz)',fontsize=14)

ax8.set_title(r'Vary $\bar{g}_{\mathregular{CaHVA}}$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gHVA[2]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gHVA[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gHVA[2]
#Nspikes_sprx_this         = Nspikes_sprx_all_gHVA[2]
for j in range(NgHVA):
    ### Everywhere:
    ax8.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygHVA[j],plotlabelgHVA), linewidth=mylinewidth)
    ### Somaprox:
    #ax8.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %.1f*%s' % (varygHVA[j],plotlabelgHVA), linewidth=mylinewidth)
ax8.set_xlabel('$I$ (nA)',fontsize=14)
ax8.set_ylabel('$f$ (Hz)',fontsize=14)

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

    infolder      = 'Allen_test_changecapacitance/figures/%i/' % (testmodel)
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
ax1.set_title(r'Vary $\bar{g}_{\mathregular{NaV}}$, Allen model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gNaV[0]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gNaV[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gNaV[0]
#Nspikes_sprx_this         = Nspikes_sprx_all_gNaV[0]
for j in range(NgNaV):
    ### Everywhere:
    ax1.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygNaV[j],plotlabelgNaV), linewidth=mylinewidth)
    ### Somaprox:
    #ax1.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %.1f*%s' % (varygNaV[j],plotlabelgNaV), linewidth=mylinewidth)
ax1.set_xlabel('$I$ (nA)',fontsize=14)
ax1.set_ylabel('$f$ (Hz)',fontsize=14)
ax1.legend(loc='upper left',ncol=1)

ax2.set_title(r'Vary $\bar{g}_{\mathregular{NaV}}$, Allen model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gNaV[1]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gNaV[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gNaV[1]
#Nspikes_sprx_this         = Nspikes_sprx_all_gNaV[1]
for j in range(NgNaV):
    ### Everywhere:
    ax2.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygNaV[j],plotlabelgNaV), linewidth=mylinewidth)
    ### Somaprox:
    #ax2.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %.1f*%s' % (varygNaV[j],plotlabelgNaV), linewidth=mylinewidth)
ax2.set_xlabel('$I$ (nA)',fontsize=14)
ax2.set_ylabel('$f$ (Hz)',fontsize=14)

ax3.set_title(r'Vary $\bar{g}_{\mathregular{NaV}}$, Allen model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gNaV[2]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gNaV[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gNaV[2]
#Nspikes_sprx_this         = Nspikes_sprx_all_gNaV[2]
for j in range(NgNaV):
    ### Everywhere:
    ax3.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*%s' % (varygNaV[j],plotlabelgNaV), linewidth=mylinewidth)
    ### Somaprox:
    #ax3.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %.1f*%s' % (varygNaV[j],plotlabelgNaV), linewidth=mylinewidth)
ax3.set_xlabel('$I$ (nA)',fontsize=14)
ax3.set_ylabel('$f$ (Hz)',fontsize=14)

################################# SOMAONLY ##################################################
cm         = 1.0
iamps      = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5]
idur       = 1000
idelay     = 10
tstop      = idur+idelay+10. #2.
dtexp      = -8 # Soma only, so this should be fast regardless
somasize   = 10
Ra         = 100.
v_init     = -70 
Ni         = len(iamps)

# Default HH values:
ena = 50
ek = -77
el_hh = -54.3
gnabar_hh = 0.12
gkbar_hh = 0.036
gl_hh = 0.0003

gnabar_hh_factors = [0.7,0.75,0.8,0.9,1.0,1.1,1.2,1.5,2.0,2.5,3.0]#[0.7,0.8,0.9,1.0,1.1,1.2,1.5,2,3,4,5,10]
Ngs = len(gnabar_hh_factors)
gnabars_hh = np.zeros(Ngs)
for i in range(Ngs):
    gnabars_hh[i] = gnabar_hh*gnabar_hh_factors[i]

NgNa_oc = Ngs
plotcolors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan','b','k','g','r','m','indigo','maroon','darkolivegreen','dodgerblue','palevioletred','darkseagreen','midnightblue','yellow','darkgoldenrod']

Nspikes_somaonly = []
Iamps_somaonly   = []
for i in range(Ngs): # This is going to be slow
    gnabar_hh = gnabars_hh[i]
    hhstring = '_ena'+str(ena)+'_ek'+str(ek)+'_el'+str(el_hh)+'_gnabar'+str(gnabar_hh)+'_gkbar'+str(gkbar_hh)+'_gl'+str(gl_hh)
    Npeaksbefore = 0
    theseiamps   = []
    theseNspikes = []
    for j in range(Ni):
        iamp = iamps[j]
        folder = 'Somaonly/Results/IStim/Soma%i/current_idur%i_iamp'%(somasize,idur)+str(iamp)+'/'
        
        filename = folder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+hhstring+'_Ra%i_vinit' %Ra+str(v_init)+'_Nspikes_vs_I.txt'
        try:
            infile = open(filename,'r')
        except:
            break
        line = infile.readline()
        words = line.split()
        Npeaks = int(words[0]) 
        infile.close()
        
        theseiamps.append(iamp)
        theseNspikes.append(Npeaks)
        
    ax4.plot(theseiamps,theseNspikes,color=plotcolors[i],label=r'%s*$\bar{g}_{\mathregular{Na}}$' % str(gnabar_hh_factors[i]))
    Nspikes_somaonly.append(theseNspikes)
    Iamps_somaonly.append(theseiamps)
ax4.set_title(r'Vary $\bar{g}_{\mathregular{Na}}$, One comp.',fontsize=16)
ax4.set_xlabel('$I$ (nA)',fontsize=14)
ax4.set_ylabel('$f$ (Hz)',fontsize=14)
ax4.legend(loc='lower right')


############################### Ball-and-stick #####################################
cm         = 1.0
iamps      = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5]
cm         = 1.0
spikedurat = -40
idur       = 1000 #100 # ms
idelay     = 100
v_init     = -65 # mV
Ra         = 100
somasize   = 10 # 15 # 
dendlen    = 1000
denddiam   = 1
nsegments  = 200 
Ni         = len(iamps)

# Default HH values:
ena = 50
ek = -77
el_hh = -54.3
gnabar_hh = 0.12
gkbar_hh = 0.036
gl_hh = 0.0003

folderstring = 'VaryNa/'
namestring = 'na'
gnabar_hh_factors = [1.0,1.2,1.5,2.0,2.5,3.0]
Ngs = len(gnabar_hh_factors)
varygs = np.zeros(Ngs)
for i in range(Ngs):
    varygs[i] = gnabar_hh*gnabar_hh_factors[i]

NgNa_bas = Ngs
plotcolors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan','b','k','g','r','m','indigo','maroon','darkolivegreen','dodgerblue','palevioletred','darkseagreen','midnightblue','yellow','darkgoldenrod']


Nspikes_bas = []
Iamps_bas   = []
for i in range(Ngs):
    infolder = 'Ball-and-stick models/BAS_somaHH_dendpassive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+folderstring
    infilename = infolder+'basHHdpas_manual_cm'+str(cm)+'_g'+namestring+str(varygs[i])+'_Nspikes_vs_Iamp.txt'
    Npeaksbefore = 0
    theseiamps   = []
    theseNspikes = []
    
    infile = open(infilename,'r')
    lines = infile.readlines()
    for line in lines:
        words = line.split()
        if len(words)>0:
            theseiamps.append(float(words[0]))
            theseNspikes.append(float(words[1]))
    infile.close()
        
    ax5.plot(theseiamps,theseNspikes,color=plotcolors[i],label=r'%.1f*$\bar{g}_{\mathregular{Na}}$' % gnabar_hh_factors[i])
    Nspikes_bas.append(theseNspikes)
    Iamps_bas.append(theseiamps)
ax5.set_title(r'Vary $\bar{g}_{\mathregular{Na}}$, Ball-and-stick',fontsize=16)
ax5.set_xlabel('$I$ (nA)',fontsize=14)
ax5.set_ylabel('$f$ (Hz)',fontsize=14)
ax5.legend(loc='lower right')


###################### Relative differences #################################
########### CaHVA ############
currat_allmodels_cahva    = []
maxdiff_allmodels_cahva   = []
maxdiff1_allmodels_cahva  = []
maxdiff2_allmodels_cahva  = []
lastdiff1_allmodels_cahva = []
cm_cahva_allmodels_store  = []
im_cahva_allmodels_store  = []
rd1_cahva_allmodels_store = []
lastdiffs_sorted_cahva    = np.zeros((Nmodels,NgHVA))
for j in range(Nmodels):
    currat_cahva    = []
    maxdiff_cahva   = []
    maxdiff1_cahva  = []
    maxdiff2_cahva  = []
    lastdiff1_cahva = []
    cm_cahva_store  = []
    im_cahva_store  = []
    rd1_cahva_store = []
    I_Nspikes_everywhere_all = I_Nspikes_everywhere_all_gHVA[j]
    Nspikes_everywhere_all = Nspikes_everywhere_all_gHVA[j]
    for i in range(NgHVA):
        im_cahva, rd1_cahva, rd2_cahva = reldiffs(Nspikes_everywhere_all[i],Nspikes_everywhere_all[5],I_Nspikes_everywhere_all[i],I_Nspikes_everywhere_all[5]) 
        if len(rd1_cahva)>0:
            maxrd1 = max(rd1_cahva)
            minrd1 = min(rd1_cahva)
            if maxrd1<abs(minrd1):
                maxrd1 = minrd1
            maxdiff1_cahva.append(maxrd1)
        if len(rd2_cahva)>0:
            maxrd2 = max(rd2_cahva)
            minrd2 = min(rd2_cahva)
            if maxrd2<abs(minrd2):
                maxrd2 = minrd2
            maxdiff2_cahva.append(maxrd2)
        if maxrd1<maxrd2:
            maxrd = maxrd2
        else:
            maxrd = maxrd1
        maxdiff_cahva.append(maxrd)
        if len(rd1_cahva)>0 or len(rd2_cahva)>0:
            currat_cahva.append(im_cahva[np.argmax(abs(rd1_cahva))])
        if len(rd1_cahva)>0:
            im_cahva_store.append(im_cahva)
            lastdiff1_cahva.append(rd1_cahva[-1])
            lastdiffs_sorted_cahva[j,i] = rd1_cahva[-1]
        else:
            lastdiffs_sorted_cahva[j,i] = 0
                
    currat_allmodels_cahva.append(currat_cahva)
    maxdiff_allmodels_cahva.append(maxdiff_cahva)
    maxdiff1_allmodels_cahva.append(maxdiff1_cahva)
    maxdiff2_allmodels_cahva.append(maxdiff2_cahva)
    lastdiff1_allmodels_cahva.append(lastdiff1_cahva)
    cm_cahva_allmodels_store.append(cm_cahva_store)
    im_cahva_allmodels_store.append(im_cahva_store)
    rd1_cahva_allmodels_store.append(rd1_cahva_store)

ax9 = plt.subplot(gs[1, 6:8])
ax9.set_title(r'I',loc='left',fontsize=18)

barWidth = 0.25
change_gcahva2   = [lastdiffs_sorted_cahva[0,8],lastdiffs_sorted_cahva[1,8],lastdiffs_sorted_cahva[2,8]]
change_gcahva3 = [lastdiffs_sorted_cahva[0,9],lastdiffs_sorted_cahva[1,9],lastdiffs_sorted_cahva[2,9]]

br1 = np.arange(len(change_gcahva2))
br2 = [x+barWidth for x in br1]

ax9.bar(br1, change_gcahva2, width=barWidth, label=r'2*$\bar{g}_{\mathregular{CaHVA}}$')
ax9.bar(br2, change_gcahva3, width=barWidth, label=r'3*$\bar{g}_{\mathregular{CaHVA}}$')
ax9.axhline(y=0,color='k',linestyle='-',linewidth=1.0)
plt.xticks(br1, ['437','965','407'])
ax9.set_xlabel('Model',fontsize=14)
ax9.set_ylabel(r'Relative difference at max. current',fontsize=14)
ax9.set_title(r'Difference from original $C_m$ and $\bar{g}_{\mathregular{CaHVA}}$',fontsize=15)
ax9.legend(loc='upper right')

############ NaV ###############
currat_allmodels_nav    = []
maxdiff_allmodels_nav   = []
maxdiff1_allmodels_nav  = []
maxdiff2_allmodels_nav  = []
lastdiff1_allmodels_nav = []
cm_nav_allmodels_store  = []
im_nav_allmodels_store  = []
rd1_nav_allmodels_store = []
lastdiffs_sorted_nav    = np.zeros((Nmodels,NgNaV))
for j in range(Nmodels):
    currat_nav    = []
    maxdiff_nav   = []
    maxdiff1_nav  = []
    maxdiff2_nav  = []
    lastdiff1_nav = []
    cm_nav_store  = []
    im_nav_store  = []
    rd1_nav_store = []
    I_Nspikes_everywhere_all = I_Nspikes_everywhere_all_gNaV[j]
    Nspikes_everywhere_all = Nspikes_everywhere_all_gNaV[j]
    for i in range(NgNaV):
        im_nav, rd1_nav, rd2_nav = reldiffs(Nspikes_everywhere_all[i],Nspikes_everywhere_all[1],I_Nspikes_everywhere_all[i],I_Nspikes_everywhere_all[1]) 
        if len(rd1_nav)>0:
            maxrd1 = max(rd1_nav)
            minrd1 = min(rd1_nav)
            if maxrd1<abs(minrd1):
                maxrd1 = minrd1
            maxdiff1_nav.append(maxrd1)
        if len(rd2_nav)>0:
            maxrd2 = max(rd2_nav)
            minrd2 = min(rd2_nav)
            if maxrd2<abs(minrd2):
                maxrd2 = minrd2
            maxdiff2_nav.append(maxrd2)
        if maxrd1<maxrd2:
            maxrd = maxrd2
        else:
            maxrd = maxrd1
        maxdiff_nav.append(maxrd)
        if len(rd1_nav)>0 or len(rd2_nav)>0:
            currat_nav.append(im_nav[np.argmax(abs(rd1_nav))])
        if len(rd1_nav)>0:
            im_nav_store.append(im_nav)
            lastdiff1_nav.append(rd1_nav[-1])
            lastdiffs_sorted_nav[j,i] = rd1_nav[-1]
        else:
            lastdiffs_sorted_nav[j,i] = 0
                
    currat_allmodels_nav.append(currat_nav)
    maxdiff_allmodels_nav.append(maxdiff_nav)
    maxdiff1_allmodels_nav.append(maxdiff1_nav)
    maxdiff2_allmodels_nav.append(maxdiff2_nav)
    lastdiff1_allmodels_nav.append(lastdiff1_nav)
    cm_nav_allmodels_store.append(cm_nav_store)
    im_nav_allmodels_store.append(im_nav_store)
    rd1_nav_allmodels_store.append(rd1_nav_store)

## Somaonly
lastdiffs_sorted_na_oc = np.zeros(NgNa_oc)

currat_na_oc    = []
maxdiff_na_oc   = []
maxdiff1_na_oc  = []
maxdiff2_na_oc  = []
lastdiff1_na_oc = []
cm_na_oc_store  = []
im_na_oc_store  = []
rd1_na_oc_store = []
for i in range(NgNa_oc):
    im_na_oc, rd1_na_oc, rd2_na_oc = reldiffs(Nspikes_somaonly[i],Nspikes_somaonly[4],Iamps_somaonly[i],Iamps_somaonly[4]) 
    if len(rd1_na_oc)>0:
        maxrd1 = max(rd1_na_oc)
        minrd1 = min(rd1_na_oc)
        if maxrd1<abs(minrd1):
            maxrd1 = minrd1
        maxdiff1_na_oc.append(maxrd1)
    if len(rd2_na_oc)>0:
        maxrd2 = max(rd2_na_oc)
        minrd2 = min(rd2_na_oc)
        if maxrd2<abs(minrd2):
            maxrd2 = minrd2
        maxdiff2_na_oc.append(maxrd2)
    if maxrd1<maxrd2:
        maxrd = maxrd2
    else:
        maxrd = maxrd1
    maxdiff_na_oc.append(maxrd)
    if len(rd1_na_oc)>0 or len(rd2_na_oc)>0:
        currat_na_oc.append(im_na_oc[np.argmax(abs(rd1_na_oc))])
    if len(rd1_na_oc)>0:
        im_nav_store.append(im_na_oc)
        lastdiff1_na_oc.append(rd1_na_oc[-1])
        lastdiffs_sorted_na_oc[i] = rd1_na_oc[-1]
    else:
        lastdiffs_sorted_na_oc[i] = 0
## BAS
lastdiffs_sorted_na_bas = np.zeros(NgNa_bas)

currat_na_bas    = []
maxdiff_na_bas   = []
maxdiff1_na_bas  = []
maxdiff2_na_bas  = []
lastdiff1_na_bas = []
cm_na_bas_store  = []
im_na_bas_store  = []
rd1_na_bas_store = []
for i in range(NgNa_bas):
    im_na_bas, rd1_na_bas, rd2_na_bas = reldiffs(Nspikes_bas[i],Nspikes_bas[0],Iamps_bas[i],Iamps_bas[0]) 
    if len(rd1_na_bas)>0:
        maxrd1 = max(rd1_na_bas)
        minrd1 = min(rd1_na_bas)
        if maxrd1<abs(minrd1):
            maxrd1 = minrd1
        maxdiff1_na_bas.append(maxrd1)
    if len(rd2_na_bas)>0:
        maxrd2 = max(rd2_na_bas)
        minrd2 = min(rd2_na_bas)
        if maxrd2<abs(minrd2):
            maxrd2 = minrd2
        maxdiff2_na_bas.append(maxrd2)
    if maxrd1<maxrd2:
        maxrd = maxrd2
    else:
        maxrd = maxrd1
    maxdiff_na_bas.append(maxrd)
    if len(rd1_na_bas)>0 or len(rd2_na_bas)>0:
        currat_na_bas.append(im_na_bas[np.argmax(abs(rd1_na_bas))])
    if len(rd1_na_bas)>0:
        im_nav_store.append(im_na_bas)
        lastdiff1_na_bas.append(rd1_na_bas[-1])
        lastdiffs_sorted_na_bas[i] = rd1_na_bas[-1]
    else:
        lastdiffs_sorted_na_bas[i] = 0

ax10 = plt.subplot(gs[1, 8:10])
ax10.set_title(r'J',loc='left',fontsize=18)

barWidth = 0.25
change_gna2   = [lastdiffs_sorted_nav[0,2],lastdiffs_sorted_nav[1,2],lastdiffs_sorted_nav[2,2],lastdiffs_sorted_na_oc[8],lastdiffs_sorted_na_bas[3]]
change_gna3 = [lastdiffs_sorted_nav[0,3],lastdiffs_sorted_nav[1,3],lastdiffs_sorted_nav[2,3],lastdiffs_sorted_na_oc[10],lastdiffs_sorted_na_bas[5]]

br1 = np.arange(len(change_gna2))
br2 = [x+barWidth for x in br1]

ax10.bar(br1, change_gna2, width=barWidth, label=r'2*$\bar{g}_{\mathregular{Na}}$')
ax10.bar(br2, change_gna3, width=barWidth, label=r'3*$\bar{g}_{\mathregular{Na}}$')
ax10.axhline(y=0,color='k',linestyle='-',linewidth=1.0)
plt.xticks(br1, ['437','965','407','One comp.','BAS'])
ax10.set_xlabel('Model',fontsize=14)
ax10.set_ylabel(r'Relative difference at max. current',fontsize=14)
ax10.set_title(r'Difference from original $C_m$ and $\bar{g}_{\mathregular{Na}}$',fontsize=15)
ax10.legend(loc='upper right')

fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(plotname)
plt.show()

