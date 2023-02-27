import numpy
import matplotlib.pyplot as plt
from matplotlib import gridspec

def avg_and_rms(x):
    N = len(x)
    avgx = numpy.mean(x)
    rmsx = 0
    for i in range(N):
        rmsx += (avgx-x[i])**2
    rmsx = numpy.sqrt(rmsx/(N-1))
    return avgx,rmsx

skiptime   = 500
spikedurat = -40
idur       = 1000 #100 # ms
idelay     = 10
v_init     = -86 # mV
Ra         = 100
somasize   = 10 # 15 # 
dtexp      = -7
t_before_rec = -600.


# change the default font family
#plt.rcParams.update({'font.family':'Arial'})

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

fig = plt.figure(figsize=(18,18),dpi=300)

gs = gridspec.GridSpec(2, 4)

ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
ax3 = plt.subplot(gs[1, 0:2])
ax4 = plt.subplot(gs[1, 2:4])

#fig.suptitle(r'Properties',fontsize=20)

ax1.set_title(r'A',loc='left',fontsize=18)
ax2.set_title(r'B',loc='left',fontsize=18)
ax3.set_title(r'C',loc='left',fontsize=18)
ax4.set_title(r'D',loc='left',fontsize=18)
ax1.set_title(r'$f$',fontsize=18)
ax2.set_title(r'$V_{\mathregular{max}}$',fontsize=18)
ax3.set_title(r'$V_{\mathregular{min}}$',fontsize=18)
ax4.set_title(r'Duration at %i mV' %spikedurat,fontsize=18)

ax1.set_xlabel(r'$I$ [nA]',fontsize=14)
ax1.set_ylabel(r'$N_{spikes}$',fontsize=14)
ax2.set_xlabel(r'$I$ [nA]',fontsize=14)
ax2.set_ylabel(r'Peak voltage [mV]',fontsize=14)
ax3.set_xlabel(r'$I$ [nA]',fontsize=14)
ax3.set_ylabel(r'Peak minima [mV]',fontsize=14)
ax4.set_xlabel(r'$I$ [nA]',fontsize=14)
ax4.set_ylabel(r'Spike duration at %s mV [ms]' % str(spikedurat),fontsize=14)


model_folder = '' # We are in the model folder
    
varyE = 0
varymech = 'None' # 'Epas' # 'EK' # 'ENa' # 
outnamestring = 'vary'
if varymech=='ENa':
    varyE = 63 #[40,50,60,70]
    outnamestring = outnamestring + '_ENa'+str(varyE)
elif varymech=='EK':
    varyE = -97#-107
    outnamestring = outnamestring + '_EK'+str(varyE)
elif varymech=='Epas':
    varyE = -20 # Vary by shifts
    outnamestring = outnamestring + '_Epasshift'+str(varyE)

vary_naf     = False
vary_kaf     = False
vary_SK      = True # False # 
vary_Ca_HVA  = True # False
vary_gpas    = False # 

gnafs   = 1.0
gkafs   = 1.0
gsks    = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,2.0]#,2.5,5.0,10.0]
gcahva  = 0.2
gpases  = 1.0
Ng      = len(gsks)

colors = []
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))

if vary_naf==True:
    outnamestring = outnamestring + '_gnaf'
if vary_kaf==True:
    outnamestring = outnamestring + '_gkaf'
if vary_SK==True:
    outnamestring = outnamestring + '_gSK'
if vary_Ca_HVA==True:
    outnamestring = outnamestring + '_gCaHVA0p2'
if vary_gpas==True: 
    outnamestring = outnamestring + '_gpas'
outnamestring = outnamestring +'_'
  
figname =   'Results/Soma%i/' % somasize+outnamestring+'_s'+str(skiptime)+'.png'

i = 0
cm = 1.0
for gsk in gsks: # Ugh, this was not that general after all...
    namestring = ''
    thelabel = ''
    if vary_naf==True:
        namestring = namestring + '_gnaf'+str(gnaf)+'p'
        thelabel   = thelabel + r'%.2f$\bar{g}_\mathregular{naf}$' % gnaf
    if vary_kaf==True:
        namestring = namestring + '_gkaf'+str(gkaf)+'p'
        thelabel   = thelabel + r'%.2f$\bar{g}_\mathregular{kaf}$' % gkaf
    if vary_SK==True:
        namestring = namestring + '_gSK'+str(gsk)+'p'
        thelabel   = thelabel + r'%.2f$\bar{g}_\mathregular{SK}$' % gsk
    if vary_Ca_HVA==True:
        namestring = namestring + '_gCaHVA'+str(gcahva)+'p'
        #thelabel   = thelabel + r'$\bar{g}_\mathregular{CaHVA}$, ' % gcahva
    if vary_gpas==True: 
        namestring = namestring + '_gpas'+str(gpas)+'p'
        thelabel   = thelabel + r'%.2f$\bar{g}_\mathregular{L}$' % gpas
    namestring = namestring +'_'
    
    # Set names
    infolder = 'Results/Soma%i/' % somasize
    infilename_Nspikes = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'
    infilename_APampl  = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Vmax_vs_I_s'+str(skiptime)+'.txt'
    infilename_APmins  = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Vmin_vs_I_s'+str(skiptime)+'.txt'
    infilename_APdhw   = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_sdurat%s_vs_I_s'% str(spikedurat)+str(skiptime)+'.txt' 
    # make files
    infile_Nspikes = open(infilename_Nspikes,'r')
    infile_APampl  = open(infilename_APampl,'r')
    infile_APmins  = open(infilename_APmins,'r')
    infile_APdhw   = open(infilename_APdhw,'r')
    lines_Nspikes = infile_Nspikes.readlines()
    lines_APampl  = infile_APampl.readlines()
    lines_APmins  = infile_APmins.readlines()
    lines_APdhw   = infile_APdhw.readlines()
    
    I_Nspikes  = []
    Nspikes    = []
    I_maxima   = []
    maxima     = []
    maxima_rms = []
    I_minima   = []
    minima     = []
    minima_rms = []
    I_dur      = []
    dur        = []
    dur_rms    = []
    
    for line in lines_Nspikes:
        words = line.split()
        if len(words)>0:
            I_Nspikes.append(float(words[0]))
            Nspikes.append(float(words[1]))
    
    for line in lines_APampl:
        words = line.split()
        if len(words)>0:
            I_maxima.append(float(words[0]))
            maxima.append(float(words[1]))
            maxima_rms.append(float(words[2]))
    
    for line in lines_APmins:
        words = line.split()
        if len(words)>0:
            I_minima.append(float(words[0]))
            minima.append(float(words[1]))
            minima_rms.append(float(words[2]))
    
    for line in lines_APdhw:
        words = line.split()
        if len(words)>0:
            I_dur.append(float(words[0]))
            dur.append(float(words[1]))
            dur_rms.append(float(words[2]))
     
    infile_Nspikes.close()
    infile_APampl.close()
    infile_APmins.close()
    infile_APdhw.close()
    
    if gcahva==1.0:
        thislinewidth = 3
    else:
        thislinewidth = 2
    
    ax1.plot(I_Nspikes,Nspikes, color=colors[i],linewidth=thislinewidth, label=thelabel)
    ax2.errorbar(I_maxima,maxima, yerr=maxima_rms, capsize=2, color=colors[i],linewidth=thislinewidth, label=thelabel)
    ax3.errorbar(I_minima,minima, yerr=minima_rms, capsize=2, color=colors[i],linewidth=thislinewidth, label=thelabel)
    ax4.errorbar(I_dur,dur, yerr=dur_rms, capsize=2, color=colors[i],linewidth=thislinewidth, label=thelabel)
    i+=1

plt.suptitle(r'%.2f$\bar{g}_\mathregular{CaHVA}$' % gcahva,fontsize=18)
ax1.legend(loc='upper left')#,ncol=1)
ax2.legend(loc='upper left')#,ncol=1)
ax3.legend(loc='upper left')#,ncol=1)
ax4.legend(loc='upper left')#,ncol=1)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(figname)
