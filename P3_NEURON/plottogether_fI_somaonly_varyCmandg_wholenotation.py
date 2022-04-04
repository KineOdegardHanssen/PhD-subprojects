import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

# change the default font family
plt.rcParams.update({'font.family':'Arial'})

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend', fontsize=13)

mylinewidth = 2

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


cms        = [1.0,1.5]
spikedurat = -40
idur       = 1000 #100 # ms
idelay     = 100
v_init     = -70 # mV
Ra         = 100
somasize   = 10 # 15 # 
dendlen    = 1000
denddiam   = 1
nsegments  = 200 
gna        = 0.12
gfactors   = [0.7,0.8,0.9,1.0]
Ngoc = len(gfactors)*len(cms)
    
varymech = 'Na' # 'K' # 'leak'
varyE_bool = False # True # 
varyg_bool = True # False # 
varyE = 50 #[30,40,50,60,70] #[30,40,70]# Every Cmf has 50 and 60. Need to run again for the other values
    
varylist = [] # Should be redundant
plotstring = '_vary'
if varyE_bool==True:
    varylist = varyE
    plotstring = plotstring + 'E'
else:
    varylist = gfactors
    plotstring = plotstring + 'g'
      
if varymech=='Na':
    plotstring = plotstring + '_Na'
    namestring = 'na'
elif varymech=='leak':
    plotstring = plotstring + '_leak'
    namestring = 'l'
elif varymech=='K':
    plotstring = plotstring + '_K'
    namestring = 'k'

plotfolder = 'Results/IStim/Soma%i/' % somasize
plotname = plotfolder + 'somaonly_fI_varyCmandg_gfrange'+str(gfactors[0])+'to'+str(gfactors[-1])+'.png'

# Default HH values:
ena = 50
ek = -77
el_hh = -54.3
gnabar_hh = 0.12
gkbar_hh = 0.036
gl_hh = 0.0003

Nspikes_OC_all = []
Iamps_OC_all   = []

gNa_oc_legends = []
for cm in cms:
    for gfactor in gfactors:  
        gnabar_hh_new = gnabar_hh*gfactor
        iamps_these   = []
        Nspikes_these = []
        gNa_oc_legends.append([cm,gfactor])
        # Set names
        infolder   = 'Results/IStim/Soma%i/' % somasize
        hhstring = '_ena'+str(ena)+'_ek'+str(ek)+'_el'+str(el_hh)+'_gnabar'+str(gnabar_hh_new)+'_gkbar'+str(gkbar_hh)+'_gl'+str(gl_hh)
        infilename = infolder+'somaonlyHH_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+hhstring+'_Nspikes_vs_I.txt'
        infile = open(infilename,'r')
        lines = infile.readlines()
        
        for line in lines:
            words = line.split()
            if len(words)>0:
                iamps_these.append(float(words[0]))
                Nspikes_these.append(int(words[1]))
        infile.close()
        
        Iamps_OC_all.append(iamps_these)
        Nspikes_OC_all.append(Nspikes_these)

fig, (ax1, ax2) = plt.subplots(1,2,figsize=(20,10))
ax1.set_title('A',loc='left',fontsize=18)
ax1.set_title(r'Vary $\bar{g}_{\mathregular{Na}}$, One comp. HH',fontsize=16)
for j in range(Ngoc):
    theselegends_oc_na = gNa_oc_legends[j]
    ### Everywhere:
    ax1.plot(Iamps_OC_all[j], Nspikes_OC_all[j],label=r'%.1f*$C_m$, %.1f*$\bar{g}_{Na}$' % (theselegends_oc_na[0],theselegends_oc_na[1]), linewidth=mylinewidth)
ax1.set_xlabel('$I$ (nA)',fontsize=14)
ax1.set_ylabel('$f$ (Hz)',fontsize=14)
ax1.legend(loc='lower right',ncol=1)
#ax1.legend(bbox_to_anchor =(0.5,-0.27), loc='lower center',ncol=3)

################## Relative differences ##############################
lastdiffs_sorted_na_oc = np.zeros(Ngoc)
currat_na_oc    = []
maxdiff_na_oc   = []
maxdiff1_na_oc  = []
maxdiff2_na_oc  = []
lastdiff1_na_oc = []
cm_na_oc_store  = []
im_na_oc_store  = []
rd1_na_oc_store = []
for i in range(Ngoc):
    im_na_oc, rd1_na_oc, rd2_na_oc = reldiffs(Nspikes_OC_all[i],Nspikes_OC_all[3],Iamps_OC_all[i],Iamps_OC_all[3]) 
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
        im_na_oc_store.append(im_na_oc)
        lastdiff1_na_oc.append(rd1_na_oc[-1])
        lastdiffs_sorted_na_oc[i] = rd1_na_oc[-1]
    else:
        lastdiffs_sorted_na_oc[i] = 0

ax2.set_title(r'B',loc='left',fontsize=18)

barlabels = []
for cm in cms:
    for gfactor in gfactors:
        barlabel = r'(%s,%s)' %  (str(cm),str(gfactor))
        barlabels.append(barlabel)

barWidth = 0.25
'''
changeCm1p5_nav   = [lastdiffs_sorted_nav[0,0],lastdiffs_sorted_nav[1,0],lastdiffs_sorted_nav[2,0],lastdiffs_sorted_na_oc[1],lastdiffs_sorted_na_bas[1]]
changeboth1p5_nav = [lastdiffs_sorted_nav[0,1],lastdiffs_sorted_nav[1,1],lastdiffs_sorted_nav[2,1],lastdiffs_sorted_na_oc[2],lastdiffs_sorted_na_bas[2]]
changeCm1p5_g2_nav = [lastdiffs_sorted_nav[0,2],lastdiffs_sorted_nav[1,2],lastdiffs_sorted_nav[2,2],lastdiffs_sorted_na_oc[3],lastdiffs_sorted_na_bas[3]]
changeCm1p5_g3_nav = [lastdiffs_sorted_nav[0,3],lastdiffs_sorted_nav[1,3],lastdiffs_sorted_nav[2,3],lastdiffs_sorted_na_oc[4],lastdiffs_sorted_na_bas[4]]
changeCm1p5_g4_nav = [lastdiffs_sorted_nav[0,4],lastdiffs_sorted_nav[1,4],lastdiffs_sorted_nav[2,4],lastdiffs_sorted_na_oc[5],lastdiffs_sorted_na_bas[5]]
changeCm1p5_g5_nav = [lastdiffs_sorted_nav[0,5],lastdiffs_sorted_nav[1,5],lastdiffs_sorted_nav[2,5],lastdiffs_sorted_na_oc[6],lastdiffs_sorted_na_bas[6]]

br1 = np.arange(len(changeCm1p5_nav))
br2 = [x+barWidth for x in br1]
br3 = [x+2*barWidth for x in br1]
br4 = [x+3*barWidth for x in br1]
br5 = [x+4*barWidth for x in br1]
br6 = [x+5*barWidth for x in br1]
brcenter = [x+0.5*barWidth for x in br3]
'''

br1 = np.arange(len(lastdiffs_sorted_na_oc))

#1.5*$C_m$, 
ax2.bar(br1, lastdiffs_sorted_na_oc, width=barWidth)
#ax2.bar(br1, changeCm1p5_nav, width=0.5*barWidth, label=r'1.0*$\bar{g}_{\mathregular{Na}}$')
#ax2.bar(br2, changeboth1p5_nav, width=0.5*barWidth, label=r'1.5*$\bar{g}_{\mathregular{Na}}$')
#ax2.bar(br3, changeCm1p5_g2_nav, width=0.5*barWidth, label=r'2.0*$\bar{g}_{\mathregular{Na}}$')
#ax2.bar(br4, changeCm1p5_g3_nav, width=0.5*barWidth, label=r'3.0*$\bar{g}_{\mathregular{Na}}$')
#ax2.bar(br5, changeCm1p5_g4_nav, width=0.5*barWidth, label=r'4.0*$\bar{g}_{\mathregular{Na}}$')
#ax2.bar(br6, changeCm1p5_g5_nav, width=0.5*barWidth, label=r'5.0*$\bar{g}_{\mathregular{Na}}$')
ax2.axhline(y=0,color='k',linestyle='-',linewidth=1.0)
plt.xticks(br1, barlabels)#['437','965','407','OC','BAS'])
ax2.set_xlabel('Model',fontsize=14)
ax2.set_ylabel(r'Relative difference at max. current',fontsize=14)
ax2.set_title(r'Difference from original $C_m$ and $\bar{g}_{\mathregular{Na}}$',fontsize=15)
#ax2.legend(loc='upper right',fontsize=12)#,ncol=2)

fig.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.savefig(plotname)
