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

vary_naf     = False
vary_kaf     = False
vary_SK      = True # False # 
vary_Ca_HVA  = True # False
vary_gpas    = False # 

gnafs   = 1.0
gkafs   = 1.0
gsks    = [1.0,3.337]#[0.1,1.0,2.5,10.0]#
gcahva  = 0.2
gpases  = 1.0
vshifts = [0,-14.5]
N       = len(gsks)

plt.figure(figsize=(6,4),dpi=300)
plt.xlabel(r'$I$ (nA)',fontsize=14)
plt.ylabel(r'$f$ (Hz)',fontsize=14)

outfolder = 'Compare/' 
figname   = outfolder+'realistic_params_s'+str(skiptime)+'.png'

cm = 1.0
I_Nspikes_all  = []
Nspikes_all    = []
for i in range(N):
    gsk = gsks[i]
    vshift = vshifts[i]
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
    folder  = 'Vshift_'+str(vshift)+'/Results/Soma%i/' % somasize
        
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
    
    plt.plot(I_Nspikes,Nspikes, label=r'$V_{\mathregular{shift}}$ = %s mV, %s$\bar{g}_\mathregular{SK}$' % (str(vshift),str(gsk)))

plt.legend(loc='lower right')#,ncol=1)
plt.tight_layout()
plt.savefig(figname)
plt.show()
