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
cms        = [1.0,1.5]

namestringfirst = ''
varygHVA = [1.0,1.5]
namestringfirstHVA = namestringfirst + '_gCaHVA'
plotlabelgHVA = r'$\bar{g}_{\mathregular{CaHVA}}$'
varygNaV = [0.5,0.8,1.0,1.5,2.0,3.0,4.0,5.0]
namestringfirstNaV = namestringfirst + '_gNaV'
plotlabelgNaV = r'$\bar{g}_{\mathregular{NaV}}$'

plotstringHVA  = '_vary_CaHVA'
plotstringNaV  = '_vary_NaV'

i_master = [0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5]
    
NI  = len(i_master)
NgHVA = len(varygHVA)*len(cms)-1

gHVA_legends = []
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
    for cm in cms:
        for g in varygHVA:
            if cm==1.0 and g==1.5:
                continue
            else:
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
                thisstring = r'%.1f*$C_m$, %.1f*$\bar{g}_{\mathregular{CaHVA}}$' % (cm,g)
                gHVA_legends.append(thisstring)
    Nspikes_everywhere_all_gHVA.append(Nspikes_everywhere_gHVA)
    I_Nspikes_everywhere_all_gHVA.append(I_Nspikes_everywhere_gHVA)
        
    #Nspikes_sprx_all_gHVA.append(Nspikes_sprx_gHVA)
    #I_Nspikes_sprx_all_gHVA.append(I_Nspikes_sprx_gHVA)

# Plotting
# May choose to specify the color again: # (can skip some of the complicated stuff here since I plot cell for cell)
color_pv_all   = [['#1f77b4','darkblue','#7f7f7f','b','tab:blue'],['#d62728','#ff7f0e','#d62728','m','tab:pink'],['#2ca02c','#8c564b','#9467bd','g','tab:brown']]
color_pv_sprx  = [['rebeccapurple','#e377c2','#8c564b','c','tab:purple'],['#9467bd','#bcbd22','#17becf','xkcd:sky blue','tab:orange'],['olive','darkgreen','mediumseagreen','tab:green','tab:olive']]

print('gHVA_legends:',gHVA_legends)

plotfolder = 'figures/Comparemodels/'
plotname = plotfolder+'fI_varyCmandg_AllenPV_withlowg.png'

## avg and rms:
#fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9), (ax10, ax11, ax12)) = plt.subplots(4, 3, figsize=(20,25))

fig = plt.figure(figsize=(15,15),dpi=300)

gs = gridspec.GridSpec(3, 6)

ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
ax3 = plt.subplot(gs[0, 4:6])
ax4 = plt.subplot(gs[1, 0:2])
ax5 = plt.subplot(gs[1, 2:4])
ax6 = plt.subplot(gs[1, 4:6])

fig.suptitle(r'Frequency $f$ vs $I$',fontsize=16)

ax1.set_title(r'A',loc='left',fontsize=18)
ax2.set_title(r'B',loc='left',fontsize=18)
ax3.set_title(r'C',loc='left',fontsize=18)
ax4.set_title(r'D',loc='left',fontsize=18)
ax5.set_title(r'E',loc='left',fontsize=18)
ax6.set_title(r'F',loc='left',fontsize=18)

# Add legends! #I could probably have looped this... ax1., ax2. etc. makes it hard.
ax1.set_title(r'Vary $\bar{g}_{\mathregular{CaHVA}}$, model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gHVA[0]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gHVA[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gHVA[0]
#Nspikes_sprx_this         = Nspikes_sprx_all_gHVA[0]
for j in range(NgHVA):
    ### Everywhere:
    ax1.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s' % gHVA_legends[j], linewidth=mylinewidth)
    ### Somaprox:
    #ax1.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'sprx: %s' % gHVA_legends[j], linewidth=mylinewidth)
ax1.set_xlabel('$I$ (nA)',fontsize=14)
ax1.set_ylabel('$f$ (Hz)',fontsize=14)
ax1.legend(loc='upper left',ncol=1)
#ax1.legend(bbox_to_anchor =(0.5,-0.27), loc='lower center',ncol=3)

ax2.set_title(r'Vary $\bar{g}_{\mathregular{CaHVA}}$, model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gHVA[1]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gHVA[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gHVA[1]
#Nspikes_sprx_this         = Nspikes_sprx_all_gHVA[1]
for j in range(NgHVA):
    ### Everywhere:
    ax2.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s' % gHVA_legends[j], linewidth=mylinewidth)
    ### Somaprox:
    #ax2.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s' % gHVA_legends[j], linewidth=mylinewidth)
ax2.set_xlabel('$I$ (nA)',fontsize=14)
ax2.set_ylabel('$f$ (Hz)',fontsize=14)

ax3.set_title(r'Vary $\bar{g}_{\mathregular{CaHVA}}$, model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gHVA[2]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gHVA[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gHVA[2]
#Nspikes_sprx_this         = Nspikes_sprx_all_gHVA[2]
for j in range(NgHVA):
    ### Everywhere:
    ax3.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%s' % gHVA_legends[j], linewidth=mylinewidth)
    ### Somaprox:
    #ax3.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %.1f*%s' % (varygHVA[j],plotlabelgHVA), linewidth=mylinewidth)
ax3.set_xlabel('$I$ (nA)',fontsize=14)
ax3.set_ylabel('$f$ (Hz)',fontsize=14)

########################### NaV ####################################################
NgNaV = len(varygNaV)*len(cms)-len(varygNaV)+1

i_master_everywhere_all_gNaV   = []
Nspikes_everywhere_all_gNaV    = []
I_Nspikes_everywhere_all_gNaV  = []

i_master_sprx_all_gNaV         = []
Nspikes_sprx_all_gNaV          = []
I_Nspikes_sprx_all_gNaV        = []

gNaV_legends = []
for testmodel in testmodels:
    i_master_everywhere_gNaV   = []
    Nspikes_everywhere_gNaV    = []
    I_Nspikes_everywhere_gNaV  = []
    
    #i_master_sprx_Na         = []
    #Nspikes_sprx_Na          = []
    #I_Nspikes_sprx_Na        = []

    infolder      = 'figures/%i/' % (testmodel)
    vrestfolder   = infolder 
    for cm in cms:
        for g in varygNaV:
            if (cm==1.0 and g!=1.0):
                continue
            else:
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
            
                #Nspikes_sprx_NaV.append(Nspikes)
                #I_Nspikes_sprx_NaV.append(I_Nspikes)
                gNaV_legends.append([cm,g])
    print('-------------------------------------')
    Nspikes_everywhere_all_gNaV.append(Nspikes_everywhere_gNaV)
    I_Nspikes_everywhere_all_gNaV.append(I_Nspikes_everywhere_gNaV)
    
    #Nspikes_sprx_all_Na.append(Nspikes_sprx_Na)
    #I_Nspikes_sprx_all_Na.append(I_Nspikes_sprx_Na)

# Plotting

# Add legends! #I could probably have looped this... ax1., ax2. etc. makes it hard.
ax4.set_title(r'Vary $\bar{g}_{\mathregular{NaV}}$, model %i' % testmodels[0],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gNaV[0]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gNaV[0]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gNaV[0]
#Nspikes_sprx_this         = Nspikes_sprx_all_gNaV[0]
for j in range(NgNaV):
    theselegends_nav = gNaV_legends[j]
    ### Everywhere:
    ax4.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*$C_m$, %.1f*$\bar{g}_{NaV}$' % (theselegends_nav[0],theselegends_nav[1]), linewidth=mylinewidth)
    ### Somaprox:
    #ax4.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'%.1f$C_m$, %.1f*$\bar{g}_{NaV}$' % (theselegends_nav[0],theselegends_nav[1]), linewidth=mylinewidth)
ax4.set_xlabel('$I$ (nA)',fontsize=14)
ax4.set_ylabel('$f$ (Hz)',fontsize=14)
#ax4.legend(loc='upper left',ncol=1)
#ax1.legend(bbox_to_anchor =(0.5,-0.27), loc='lower center',ncol=3)

ax5.set_title(r'Vary $\bar{g}_{\mathregular{NaV}}$, model %i' % testmodels[1],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gNaV[1]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gNaV[1]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gNaV[1]
#Nspikes_sprx_this         = Nspikes_sprx_all_gNaV[1]
for j in range(NgNaV):
    ### Everywhere:
    theselegends_nav = gNaV_legends[j]
    ### Everywhere:
    ax5.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f$C_m$, %.1f*$\bar{g}_{NaV}$' % (theselegends_nav[0],theselegends_nav[1]), linewidth=mylinewidth)
    ### Somaprox:
    #ax5.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s' % gNaV_legends, linewidth=mylinewidth)
ax5.set_xlabel('$I$ (nA)',fontsize=14)
ax5.set_ylabel('$f$ (Hz)',fontsize=14)

ax6.set_title(r'Vary $\bar{g}_{\mathregular{NaV}}$, model %i' % testmodels[2],fontsize=16)
I_Nspikes_everywhere_this = I_Nspikes_everywhere_all_gNaV[2]
#I_Nspikes_sprx_this       = I_Nspikes_sprx_all_gNaV[2]
Nspikes_everywhere_this   = Nspikes_everywhere_all_gNaV[2]
#Nspikes_sprx_this         = Nspikes_sprx_all_gNaV[2]
for j in range(NgNaV):
    theselegends_nav = gNaV_legends[j]
    ### Everywhere:
    ax6.plot(I_Nspikes_everywhere_this[j], Nspikes_everywhere_this[j],label=r'%.1f*$C_m$, %.1f*$\bar{g}_{NaV}$' % (theselegends_nav[0],theselegends_nav[1]), linewidth=mylinewidth)
    ### Somaprox:
    #ax6.plot(I_Nspikes_sprx_this[j], Nspikes_sprx_this[j],'--',label=r'somaprox: %s' % gNaV_legends, linewidth=mylinewidth)
ax6.set_xlabel('$I$ (nA)',fontsize=14)
ax6.set_ylabel('$f$ (Hz)',fontsize=14)
ax6.legend(loc='upper left',ncol=1)

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
lastdiffs_sorted_cahva    = np.zeros((Nmodels,NgHVA-1))
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
    for i in range(1,NgHVA):
        im_cahva, rd1_cahva, rd2_cahva = reldiffs(Nspikes_everywhere_all[i],Nspikes_everywhere_all[0],I_Nspikes_everywhere_all[i],I_Nspikes_everywhere_all[0]) 
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
            lastdiffs_sorted_cahva[j,i-1] = rd1_cahva[-1]
        else:
            lastdiffs_sorted_cahva[j,i-1] = 0
                
    currat_allmodels_cahva.append(currat_cahva)
    maxdiff_allmodels_cahva.append(maxdiff_cahva)
    maxdiff1_allmodels_cahva.append(maxdiff1_cahva)
    maxdiff2_allmodels_cahva.append(maxdiff2_cahva)
    lastdiff1_allmodels_cahva.append(lastdiff1_cahva)
    cm_cahva_allmodels_store.append(cm_cahva_store)
    im_cahva_allmodels_store.append(im_cahva_store)
    rd1_cahva_allmodels_store.append(rd1_cahva_store)

ax7 = plt.subplot(gs[2, 1:3])
ax7.set_title(r'G',loc='left',fontsize=18)

barWidth = 0.25
changeCm1p5_cahva   = [lastdiffs_sorted_cahva[0,0],lastdiffs_sorted_cahva[1,0],lastdiffs_sorted_cahva[2,0]]
changeboth1p5_cahva = [lastdiffs_sorted_cahva[0,1],lastdiffs_sorted_cahva[1,1],lastdiffs_sorted_cahva[2,1]]

br1 = np.arange(len(changeCm1p5_cahva))
br2 = [x+barWidth for x in br1]

ax7.bar(br1, changeCm1p5_cahva, width=barWidth, label=r'1.5*$C_m$, 1.0*$\bar{g}_{\mathregular{CaHVA}}$')
ax7.bar(br2, changeboth1p5_cahva, width=barWidth, label=r'1.5*$C_m$, 1.5*$\bar{g}_{\mathregular{CaHVA}}$')
ax7.axhline(y=0,color='k',linestyle='-',linewidth=1.0)
plt.xticks(br1, ['437','965','407']) # Hmmmmm
ax7.set_xlabel('Model',fontsize=14)
ax7.set_ylabel(r'Relative difference at max. current',fontsize=14)
ax7.set_title(r'Difference from original $C_m$ and $\bar{g}_{\mathregular{CaHVA}}$',fontsize=15,loc='right')
ax7.legend(loc='upper right')

############ NaV ###############
currat_allmodels_nav    = []
maxdiff_allmodels_nav   = []
maxdiff1_allmodels_nav  = []
maxdiff2_allmodels_nav  = []
lastdiff1_allmodels_nav = []
cm_nav_allmodels_store  = []
im_nav_allmodels_store  = []
rd1_nav_allmodels_store = []
lastdiffs_sorted_nav    = np.zeros((Nmodels,NgNaV-1))
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
    for i in range(1,NgNaV):
        im_nav, rd1_nav, rd2_nav = reldiffs(Nspikes_everywhere_all[i],Nspikes_everywhere_all[0],I_Nspikes_everywhere_all[i],I_Nspikes_everywhere_all[0]) 
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
            lastdiffs_sorted_nav[j,i-1] = rd1_nav[-1]
        else:
            lastdiffs_sorted_nav[j,i-1] = 0
                
    currat_allmodels_nav.append(currat_nav)
    maxdiff_allmodels_nav.append(maxdiff_nav)
    maxdiff1_allmodels_nav.append(maxdiff1_nav)
    maxdiff2_allmodels_nav.append(maxdiff2_nav)
    lastdiff1_allmodels_nav.append(lastdiff1_nav)
    cm_nav_allmodels_store.append(cm_nav_store)
    im_nav_allmodels_store.append(im_nav_store)
    rd1_nav_allmodels_store.append(rd1_nav_store)

ax8 = plt.subplot(gs[2, 3:5])
ax8.set_title(r'H',loc='left',fontsize=18)

barWidth = 0.125
changeCm1p5_g05_nav   = [lastdiffs_sorted_nav[0,0],lastdiffs_sorted_nav[1,0],lastdiffs_sorted_nav[2,0]]
changeCm1p5_g08_nav   = [lastdiffs_sorted_nav[0,1],lastdiffs_sorted_nav[1,1],lastdiffs_sorted_nav[2,1]]
changeCm1p5_nav   = [lastdiffs_sorted_nav[0,2],lastdiffs_sorted_nav[1,2],lastdiffs_sorted_nav[2,2]]
changeboth1p5_nav = [lastdiffs_sorted_nav[0,3],lastdiffs_sorted_nav[1,3],lastdiffs_sorted_nav[2,3]]
changeCm1p5_g2_nav = [lastdiffs_sorted_nav[0,4],lastdiffs_sorted_nav[1,4],lastdiffs_sorted_nav[2,4]]
changeCm1p5_g3_nav = [lastdiffs_sorted_nav[0,5],lastdiffs_sorted_nav[1,5],lastdiffs_sorted_nav[2,5]]
changeCm1p5_g4_nav = [lastdiffs_sorted_nav[0,6],lastdiffs_sorted_nav[1,6],lastdiffs_sorted_nav[2,6]]
changeCm1p5_g5_nav = [lastdiffs_sorted_nav[0,7],lastdiffs_sorted_nav[1,7],lastdiffs_sorted_nav[2,7]]

br1 = np.arange(len(changeCm1p5_nav))
br2 = [x+barWidth for x in br1]
br3 = [x+2*barWidth for x in br1]
br4 = [x+3*barWidth for x in br1]
br5 = [x+4*barWidth for x in br1]
br6 = [x+5*barWidth for x in br1]
br7 = [x+6*barWidth for x in br1]
br8 = [x+7*barWidth for x in br1]
brcenter = [x+0.5*barWidth for x in br4]

#1.5*$C_m$, 
ax8.bar(br1, changeCm1p5_g05_nav, width=0.5*barWidth, label=r'0.5*$\bar{g}_{\mathregular{NaV}}$')
ax8.bar(br2, changeCm1p5_g08_nav, width=0.5*barWidth, label=r'0.8*$\bar{g}_{\mathregular{NaV}}$')
ax8.bar(br3, changeCm1p5_nav, width=0.5*barWidth, label=r'1.0*$\bar{g}_{\mathregular{NaV}}$')
ax8.bar(br4, changeboth1p5_nav, width=0.5*barWidth, label=r'1.5*$\bar{g}_{\mathregular{NaV}}$')
ax8.bar(br5, changeCm1p5_g2_nav, width=0.5*barWidth, label=r'2.0*$\bar{g}_{\mathregular{NaV}}$')
ax8.bar(br6, changeCm1p5_g3_nav, width=0.5*barWidth, label=r'3.0*$\bar{g}_{\mathregular{NaV}}$')
ax8.bar(br7, changeCm1p5_g4_nav, width=0.5*barWidth, label=r'4.0*$\bar{g}_{\mathregular{NaV}}$')
ax8.bar(br8, changeCm1p5_g5_nav, width=0.5*barWidth, label=r'5.0*$\bar{g}_{\mathregular{NaV}}$')
ax8.axhline(y=0,color='k',linestyle='-',linewidth=1.0)
plt.xticks(br1, ['437','965','407'])
ax8.set_xlabel('Model',fontsize=14)
ax8.set_ylabel(r'Relative difference at max. current',fontsize=14)
ax8.set_title(r'Difference from original $C_m$ and $\bar{g}_{\mathregular{NaV}}$',fontsize=15)
ax8.legend(loc='upper center',fontsize=12)#,ncol=2)

fig.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.savefig(plotname)

plt.show()