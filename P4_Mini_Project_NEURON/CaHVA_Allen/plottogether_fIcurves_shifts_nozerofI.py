import numpy
import matplotlib.pyplot as plt
from matplotlib import gridspec

live = False # True # 

def avg_and_rms(x):
    N = len(x)
    avgx = numpy.mean(x)
    rmsx = 0
    for i in range(N):
        rmsx += (avgx-x[i])**2
    rmsx = numpy.sqrt(rmsx/(N-1))
    return avgx,rmsx

spikedurat = -40
idur       = 1000 #100 # ms
idelay     = 10
v_init     = -86 # mV
Ra         = 100
somasize   = 10 # 15 # 
dtexp      = -7
t_before_rec = -600.
skiptime   = 500


# change the default font family
#plt.rcParams.update({'font.family':'Arial'})

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

# I'm keeping this in case I need it:
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
vary_Ca_HVA  = True # False
vary_gpas    = False # 

gnafs         = 1.0
gkafs         = 1.0
gcahvas       = [0.02,0.1,0.2,0.4]#[0.1,0.5,1.0,2.0] # Possibly change these. A lot of simulations need to be run first
gpases        = 1.0
vshifts       = [50,45,40,35,30,25,20,15,10,5,0,-5,-10,-15,-20,-25,-30,-35,-40,-45,-50]
vshifts_label = [-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50]
Ng            = len(gcahvas)
NV            = len(vshifts)

colors = []
for i in range(NV):
    colors.append((1.0-i/float(NV),0,i/float(NV),1.0-abs(0.5-i/float(NV))))

if live==False:
    fig = plt.figure(figsize=(18,10),dpi=300)
else:
    fig = plt.figure(figsize=(18,10))

gs = gridspec.GridSpec(3, 8, height_ratios=[4, 4, 1])

ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
ax3 = plt.subplot(gs[0, 4:6])
ax4 = plt.subplot(gs[0, 6:8])
ax5 = plt.subplot(gs[1, 0:2])
ax6 = plt.subplot(gs[1, 2:4])
ax7 = plt.subplot(gs[1, 4:6])
ax8 = plt.subplot(gs[1, 6:8])
ax9 = plt.subplot(gs[2, 0:8])


#fig.suptitle(r'Properties',fontsize=20)

ax1.set_title(r'A',loc='left',fontsize=18)
ax2.set_title(r'B',loc='left',fontsize=18)
ax3.set_title(r'C',loc='left',fontsize=18)
ax4.set_title(r'D',loc='left',fontsize=18)
ax5.set_title(r'E',loc='left',fontsize=18)
ax6.set_title(r'F',loc='left',fontsize=18)
ax7.set_title(r'G',loc='left',fontsize=18)
ax8.set_title(r'H',loc='left',fontsize=18)
ax1.set_title(r'$f-I$-plot, %.1f$g_{\mathregular{CaHVA}}$' % (5.0*gcahvas[0]),fontsize=18)
ax2.set_title(r'$f-I$-plot, %.1f$g_{\mathregular{CaHVA}}$' % (5.0*gcahvas[1]),fontsize=18)
ax3.set_title(r'$f-I$-plot, %.1f$g_{\mathregular{CaHVA}}$' % (5.0*gcahvas[2]),fontsize=18)
ax4.set_title(r'$f-I$-plot, %.1f$g_{\mathregular{CaHVA}}$' % (5.0*gcahvas[3]),fontsize=18)

ax1.set_xlabel(r'$I$ [nA]',fontsize=14)
ax1.set_ylabel(r'$f$ [Hz]',fontsize=14)
ax2.set_xlabel(r'$I$ [nA]',fontsize=14)
ax2.set_ylabel(r'$f$ [Hz]',fontsize=14)
ax3.set_xlabel(r'$I$ [nA]',fontsize=14)
ax3.set_ylabel(r'$f$ [Hz]',fontsize=14)
ax4.set_xlabel(r'$I$ [nA]',fontsize=14)
ax4.set_ylabel(r'$f$ [Hz]',fontsize=14)

outnamestring
if vary_naf==True:
    outnamestring = outnamestring + '_gnaf'
if vary_kaf==True:
    outnamestring = outnamestring + '_gkaf'
if vary_Ca_HVA==True:
    outnamestring = outnamestring + '_gCaHVA'
if vary_gpas==True: 
    outnamestring = outnamestring + '_gpas'
outnamestring = outnamestring +'_'

outfolder = 'Results/Soma10/Compare/' 
figname   = outfolder+'shift_bothsystems_s'+str(skiptime)+'.png'
infolder_Ca_base = 'Shift_gCaHVA_V/'

cm = 1.0
I_Nspikes_all  = []
Nspikes_all    = []
for gcahva in gcahvas:
    I_Nspikes_thisg  = []
    Nspikes_thisg    = []
    namestring = ''
    thelabel = ''
    if vary_naf==True:
        namestring = namestring + '_gnaf'+str(gnaf)+'p'
    if vary_kaf==True:
        namestring = namestring + '_gkaf'+str(gkaf)+'p'
    if vary_Ca_HVA==True:
        namestring = namestring + '_gCaHVA'+str(gcahva)+'p'
    if vary_gpas==True: 
        namestring = namestring + '_gpas'+str(gpas)+'p'
    namestring = namestring +'_'
    for vshift in vshifts:
        folder  = infolder_Ca_base+'Vshift_'+str(vshift)+'/Results/Soma%i/' % somasize
        
        # Set names
        infilename_Nspikes = folder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'
        # make files
        infile_Nspikes = open(infilename_Nspikes,'r')
        lines_Nspikes = infile_Nspikes.readlines()
        
        I_Nspikes  = []
        Nspikes    = []
        
        for line in lines_Nspikes:
            words = line.split()
            if len(words)>0:
                I_Nspikes.append(float(words[0]))
                Nspikes.append(float(words[1]))
         
        infile_Nspikes.close()
        
        I_Nspikes_thisg.append(I_Nspikes)
        Nspikes_thisg.append(Nspikes)

    I_Nspikes_all.append(I_Nspikes_thisg)
    Nspikes_all.append(Nspikes_thisg)

I_Nspikes_this = I_Nspikes_all[0]
Nspikes_this   = Nspikes_all[0]

for i in range(NV):
    I_Nspikes  = I_Nspikes_this[i]
    Nspikes    = Nspikes_this[i]
    if numpy.sum(Nspikes)>0:
        if vshifts[i]==0:
            ax1.plot(I_Nspikes, Nspikes,linewidth=3,color=colors[i])    
        else:
            ax1.plot(I_Nspikes, Nspikes,color=colors[i])

I_Nspikes_this = I_Nspikes_all[1]
Nspikes_this   = Nspikes_all[1]
for i in range(NV):
    I_Nspikes  = I_Nspikes_this[i]
    Nspikes    = Nspikes_this[i]
    if numpy.sum(Nspikes)>0:
        if vshifts[i]==0:
            ax2.plot(I_Nspikes, Nspikes,linewidth=3,color=colors[i])
        else:
            ax2.plot(I_Nspikes, Nspikes,linewidth=2,color=colors[i])

I_Nspikes_this = I_Nspikes_all[2]
Nspikes_this   = Nspikes_all[2]
for i in range(NV):
    I_Nspikes  = I_Nspikes_this[i]
    Nspikes    = Nspikes_this[i]
    if numpy.sum(Nspikes)>0:
        if vshifts[i]==0:
            ax3.plot(I_Nspikes, Nspikes,linewidth=3,color=colors[i])
        else:
            ax3.plot(I_Nspikes, Nspikes,linewidth=2,color=colors[i])

I_Nspikes_this = I_Nspikes_all[3]
Nspikes_this   = Nspikes_all[3]
for i in range(NV):
    I_Nspikes  = I_Nspikes_this[i]
    Nspikes    = Nspikes_this[i]
    if numpy.sum(Nspikes)>0:
        if vshifts[i]==0:
            ax4.plot(I_Nspikes, Nspikes,linewidth=3,color=colors[i])
        else:
            ax4.plot(I_Nspikes, Nspikes,linewidth=2,color=colors[i])


varyE = 0
varymech = 'None'
outnamestring = 'vary'
if varymech=='ENa':
    varyE = 63 
    outnamestring = outnamestring + '_ENa'+str(varyE)
elif varymech=='EK':
    varyE = -97
    outnamestring = outnamestring + '_EK'+str(varyE)
elif varymech=='Epas':
    varyE = -20
    outnamestring = outnamestring + '_Epasshift'+str(varyE)

vary_naf     = False
vary_kaf     = False
vary_SK      = True # False # 
vary_Ca_HVA  = True # False
vary_gpas    = False # 

gnafs         = 1.0
gkafs         = 1.0
gsks          = [0.1,0.5,1.0,2.0]
gcahva        = 0.2
gpases        = 1.0
vshifts       = [50,45,40,35,30,25,20,15,10,5,0,-5,-10,-15,-20,-25,-30,-35,-40,-45,-50]
vshifts_label = [-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50]
Ng            = len(gsks)
NV            = len(vshifts)

colors = []
for i in range(NV):
    colors.append((1.0-i/float(NV),0,i/float(NV),1.0-abs(0.5-i/float(NV))))


ax5.set_title(r'$f-I$-plot, %.1f$g_{\mathregular{SK}}$' % gsks[0],fontsize=18)
ax6.set_title(r'$f-I$-plot, %.1f$g_{\mathregular{SK}}$' % gsks[1],fontsize=18)
ax7.set_title(r'$f-I$-plot, %.1f$g_{\mathregular{SK}}$' % gsks[2],fontsize=18)
ax8.set_title(r'$f-I$-plot, %.1f$g_{\mathregular{SK}}$' % gsks[3],fontsize=18)

ax5.set_xlabel(r'$I$ [nA]',fontsize=14)
ax5.set_ylabel(r'$f$ [Hz]',fontsize=14)
ax6.set_xlabel(r'$I$ [nA]',fontsize=14)
ax6.set_ylabel(r'$f$ [Hz]',fontsize=14)
ax7.set_xlabel(r'$I$ [nA]',fontsize=14)
ax7.set_ylabel(r'$f$ [Hz]',fontsize=14)
ax8.set_xlabel(r'$I$ [nA]',fontsize=14)
ax8.set_ylabel(r'$f$ [Hz]',fontsize=14)

outnamestring
if vary_naf==True:
    outnamestring = outnamestring + '_gnaf'
if vary_kaf==True:
    outnamestring = outnamestring + '_gkaf'
if vary_SK==True:
    outnamestring = outnamestring + '_gSK'
if vary_Ca_HVA==True:
    outnamestring = outnamestring + '_gCaHVA'
if vary_gpas==True: 
    outnamestring = outnamestring + '_gpas'
outnamestring = outnamestring +'_'


infolder_SK_base = 'SK_Allen/Shift_gCaHVA_V/'

cm = 1.0
I_Nspikes_all  = []
Nspikes_all    = []
for gsk in gsks:
    I_Nspikes_thisg  = []
    Nspikes_thisg    = []
    namestring = ''
    thelabel = ''
    if vary_naf==True:
        namestring = namestring + '_gnaf'+str(gnaf)+'p'
    if vary_kaf==True:
        namestring = namestring + '_gkaf'+str(gkaf)+'p'
    if vary_SK==True:
        namestring = namestring + '_gSK'+str(gsk)+'p'
    if vary_Ca_HVA==True:
        namestring = namestring + '_gCaHVA'+str(gcahva)+'p'
    if vary_gpas==True: 
        namestring = namestring + '_gpas'+str(gpas)+'p'
    namestring = namestring +'_'
    for vshift in vshifts:
        folder  = infolder_SK_base+'Vshift_'+str(vshift)+'/Results/Soma%i/' % somasize
        
        # Set names
        infilename_Nspikes = folder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'
        # make files
        infile_Nspikes = open(infilename_Nspikes,'r')
        lines_Nspikes = infile_Nspikes.readlines()
        
        I_Nspikes  = []
        Nspikes    = []
        
        for line in lines_Nspikes:
            words = line.split()
            if len(words)>0:
                I_Nspikes.append(float(words[0]))
                Nspikes.append(float(words[1]))
         
        infile_Nspikes.close()
        
        I_Nspikes_thisg.append(I_Nspikes)
        Nspikes_thisg.append(Nspikes)

    I_Nspikes_all.append(I_Nspikes_thisg)
    Nspikes_all.append(Nspikes_thisg)

# First plot

I_Nspikes_this = I_Nspikes_all[0]
Nspikes_this   = Nspikes_all[0]

for i in range(NV):
    I_Nspikes  = I_Nspikes_this[i]
    Nspikes    = Nspikes_this[i]
    if numpy.sum(Nspikes)!=0:
        if vshifts[i]==0:
            ax5.plot(I_Nspikes, Nspikes,linewidth=3,color=colors[i])
        else:
            ax5.plot(I_Nspikes, Nspikes,linewidth=2,color=colors[i])

I_Nspikes_this = I_Nspikes_all[1]
Nspikes_this   = Nspikes_all[1]
for i in range(NV):
    I_Nspikes  = I_Nspikes_this[i]
    Nspikes    = Nspikes_this[i]
    if numpy.sum(Nspikes)!=0:
        if vshifts[i]==0:
            ax6.plot(I_Nspikes, Nspikes,linewidth=3,color=colors[i])
        else:
            ax6.plot(I_Nspikes, Nspikes,linewidth=2,color=colors[i])

I_Nspikes_this = I_Nspikes_all[2]
Nspikes_this   = Nspikes_all[2]
for i in range(NV):
    I_Nspikes  = I_Nspikes_this[i]
    Nspikes    = Nspikes_this[i]
    if numpy.sum(Nspikes)!=0:
        if vshifts[i]==0:
            ax7.plot(I_Nspikes, Nspikes,linewidth=3,color=colors[i])
        else:
            ax7.plot(I_Nspikes, Nspikes,linewidth=2,color=colors[i])

I_Nspikes_this = I_Nspikes_all[3]
Nspikes_this   = Nspikes_all[3]
for i in range(NV):
    I_Nspikes  = I_Nspikes_this[i]
    Nspikes    = Nspikes_this[i]
    if numpy.sum(Nspikes)!=0:
        if vshifts[i]==0:
            ax8.plot(I_Nspikes, Nspikes,linewidth=3,color=colors[i])
        else:
            ax8.plot(I_Nspikes, Nspikes,linewidth=2,color=colors[i])

ax9.set_frame_on(False)
ax9.get_xaxis().set_visible(False)
ax9.get_yaxis().set_visible(False)
for i in range(NV):
    ax9.plot('-',color=colors[i],label=r'$V_{\mathregular{shift}}$ = %s mV' % str(vshifts_label[i]))
ax9.legend(loc='upper center',ncol=7, fontsize=12.5)

fig.tight_layout()
if live==False:
    plt.savefig(figname)
else:
    plt.show()
