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

gnafs   = 1.0
gkafs   = 1.0
gcahvas = [0.02,0.1,0.2,0.4]#[0.1,0.5,1.0,2.0] # Possibly change these. A lot of simulations need to be run first
gpases  = 1.0
vshifts = [-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50]
Ng      = len(gcahvas)
NV      = len(vshifts)

colors = []
for i in range(NV):
    colors.append((1.0-i/float(NV),0,i/float(NV),1.0-abs(0.5-i/float(NV))))

if live==False:
    fig = plt.figure(figsize=(18,18),dpi=300)
else:
    fig = plt.figure(figsize=(18,18))

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

outfolder = 'Compare/' 
figname   = outfolder+outnamestring+'_s'+str(skiptime)+'.png'

cm = 1.0
I_Nspikes_all  = []
Nspikes_all    = []
#I_maxima_all   = []
#maxima_all     = []
#maxima_rms_all = []
#I_minima_all   = []
#minima_all     = []
#minima_rms_all = []
#I_dur_all      = []
#dur_all        = []
#dur_rms_all    = []
for gcahva in gcahvas:
    I_Nspikes_thisg  = []
    Nspikes_thisg    = []
    #I_maxima_thisg   = []
    #maxima_thisg     = []
    #maxima_rms_thisg = []
    #I_minima_thisg   = []
    #minima_thisg     = []
    #minima_rms_thisg = []
    #I_dur_thisg      = []
    #dur_thisg        = []
    #dur_rms_thisg    = []
    namestring = ''
    thelabel = ''
    if vary_naf==True:
        namestring = namestring + '_gnaf'+str(gnaf)+'p'
        #thelabel   = thelabel + r'$\bar{g}_\mathregular{naf}$=%.1f, ' % gnaf
    if vary_kaf==True:
        namestring = namestring + '_gkaf'+str(gkaf)+'p'
        #thelabel   = thelabel + r'$\bar{g}_\mathregular{kaf}$=%.1f, ' % gkaf
    if vary_Ca_HVA==True:
        namestring = namestring + '_gCaHVA'+str(gcahva)+'p'
        #thelabel   = thelabel + r'$\bar{g}_\mathregular{CaHVA}$=%.1f, ' % gcahva
    if vary_gpas==True: 
        namestring = namestring + '_gpas'+str(gpas)+'p'
        #thelabel   = thelabel + r'$\bar{g}_\mathregular{L}$=%.1f, ' % gpas
    namestring = namestring +'_'
    for vshift in vshifts:
        folder  = 'Vshift_'+str(vshift)+'/Results/Soma%i/' % somasize
        
        # Set names
        infilename_Nspikes = folder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'
        #infilename_APampl  = folder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Vmax_vs_I.txt'
        #infilename_APmins  = folder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Vmin_vs_I.txt'
        #infilename_APdhw   = folder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_sdurat%s_vs_I.txt' % str(spikedurat)
        # make files
        infile_Nspikes = open(infilename_Nspikes,'r')
        #infile_APampl  = open(infilename_APampl,'r')
        #infile_APmins  = open(infilename_APmins,'r')
        #infile_APdhw   = open(infilename_APdhw,'r')
        lines_Nspikes = infile_Nspikes.readlines()
        #lines_APampl  = infile_APampl.readlines()
        #lines_APmins  = infile_APmins.readlines()
        #lines_APdhw   = infile_APdhw.readlines()
        
        I_Nspikes  = []
        Nspikes    = []
        #I_maxima   = []
        #maxima     = []
        #maxima_rms = []
        #I_minima   = []
        #minima     = []
        #minima_rms = []
        #I_dur      = []
        #dur        = []
        #dur_rms    = []
        
        for line in lines_Nspikes:
            words = line.split()
            if len(words)>0:
                I_Nspikes.append(float(words[0]))
                Nspikes.append(float(words[1]))
        
        #for line in lines_APampl:
        #    words = line.split()
        #    if len(words)>0:
        #        I_maxima.append(float(words[0]))
        #        maxima.append(float(words[1]))
        #        maxima_rms.append(float(words[2]))
        
        #for line in lines_APmins:
        #    words = line.split()
        #    if len(words)>0:
        #        I_minima.append(float(words[0]))
        #        minima.append(float(words[1]))
        #        minima_rms.append(float(words[2]))
        
        #for line in lines_APdhw:
        #    words = line.split()
        #    if len(words)>0:
        #        I_dur.append(float(words[0]))
        #        dur.append(float(words[1]))
        #        dur_rms.append(float(words[2]))
         
        infile_Nspikes.close()
        #infile_APampl.close()
        #infile_APmins.close()
        #infile_APdhw.close()
        
        I_Nspikes_thisg.append(I_Nspikes)
        Nspikes_thisg.append(Nspikes)
        #I_maxima_thisg.append(I_maxima)
        #maxima_thisg.append(maxima)
        #maxima_rms_thisg.append(maxima_rms)
        #I_minima_thisg.append(I_minima)
        #minima_thisg.append(minima)
        #minima_rms_thisg.append(minima_rms)
        #I_dur_thisg.append(I_dur)
        #dur_thisg.append(dur)
        #dur_rms_thisg.append(dur_rms)

    I_Nspikes_all.append(I_Nspikes_thisg)
    Nspikes_all.append(Nspikes_thisg)
    #I_maxima_all.append(I_maxima_thisg)
    #maxima_all.append(maxima_thisg)
    #maxima_rms_all.append(maxima_rms_thisg)
    #I_minima_all.append(I_minima_thisg)
    #minima_all.append(minima_thisg)
    #minima_rms_all.append(minima_rms_thisg)
    #I_dur_all.append(I_dur_thisg)
    #dur_all.append(dur_thisg)
    #dur_rms_all.append(dur_rms_thisg)

# First plot

I_Nspikes_this = I_Nspikes_all[0]
Nspikes_this   = Nspikes_all[0]

for i in range(NV):
    I_Nspikes  = I_Nspikes_this[i]
    Nspikes    = Nspikes_this[i]
    if numpy.sum(Nspikes)>0:
        if vshifts[i]==0:
            ax1.plot(I_Nspikes, Nspikes, label=r'$V_{\mathregular{shift}}$ = %s mV' % str(vshifts[i]),linewidth=3,color=colors[i])    
        else:
            ax1.plot(I_Nspikes, Nspikes, label=r'$V_{\mathregular{shift}}$ = %s mV' % str(vshifts[i]),linewidth=2,color=colors[i])

I_Nspikes_this = I_Nspikes_all[1]
Nspikes_this   = Nspikes_all[1]
for i in range(NV):
    I_Nspikes  = I_Nspikes_this[i]
    Nspikes    = Nspikes_this[i]
    if numpy.sum(Nspikes)>0:
        if vshifts[i]==0:
            ax2.plot(I_Nspikes, Nspikes, label=r'$V_{\mathregular{shift}}$ = %s mV' % str(vshifts[i]),linewidth=3,color=colors[i])
        else:
            ax2.plot(I_Nspikes, Nspikes, label=r'$V_{\mathregular{shift}}$ = %s mV' % str(vshifts[i]),linewidth=2,color=colors[i])

I_Nspikes_this = I_Nspikes_all[2]
Nspikes_this   = Nspikes_all[2]
for i in range(NV):
    I_Nspikes  = I_Nspikes_this[i]
    Nspikes    = Nspikes_this[i]
    if numpy.sum(Nspikes)>0:
        if vshifts[i]==0:
            ax3.plot(I_Nspikes, Nspikes, label=r'$V_{\mathregular{shift}}$ = %s mV' % str(vshifts[i]),linewidth=3,color=colors[i])
        else:
            ax3.plot(I_Nspikes, Nspikes, label=r'$V_{\mathregular{shift}}$ = %s mV' % str(vshifts[i]),linewidth=2,color=colors[i])

I_Nspikes_this = I_Nspikes_all[3]
Nspikes_this   = Nspikes_all[3]
for i in range(NV):
    I_Nspikes  = I_Nspikes_this[i]
    Nspikes    = Nspikes_this[i]
    if numpy.sum(Nspikes)>0:
        if vshifts[i]==0:
            ax4.plot(I_Nspikes, Nspikes, label=r'$V_{\mathregular{shift}}$ = %s mV' % str(vshifts[i]),linewidth=3,color=colors[i])
        else:
            ax4.plot(I_Nspikes, Nspikes, label=r'$V_{\mathregular{shift}}$ = %s mV' % str(vshifts[i]),linewidth=2,color=colors[i])


ax1.legend(loc='lower right')#,ncol=1)
ax2.legend(loc='lower right')#,ncol=1)
ax3.legend(loc='lower right')#,ncol=1)
ax4.legend(loc='lower right')#,ncol=1)
plt.tight_layout()
if live==False:
    plt.savefig(figname)
else:
    plt.show()
