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

vshift      = 0
spikedurat = -40
idur       = 1000 #100 # ms
idelay     = 10
v_init     = -86 # mV
Ra         = 100
somasize   = 10 # 15 # 
dtexp      = -7
skiptime   = 500
t_before_rec = -600.

cm = 1.0


# change the default font family
#plt.rcParams.update({'font.family':'Arial'})

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
    
fig = plt.figure(figsize=(18,6),dpi=300)
    
gs = gridspec.GridSpec(1, 6)
    
ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
ax3 = plt.subplot(gs[0, 4:6])
    
#fig.suptitle(r'Properties',fontsize=20)

ax1.set_title(r'A',loc='left',fontsize=18)
ax2.set_title(r'B',loc='left',fontsize=18)
ax3.set_title(r'C',loc='left',fontsize=18)
ax1.set_title(r'Models',fontsize=18)
ax2.set_title(r'Altering $\bar{g}_\mathregular{CaHVA}$, no SK',fontsize=18)
ax3.set_title(r'Altering $\bar{g}_\mathregular{SK}$',fontsize=18)
    
ax1.set_xlabel(r'$I$ [nA]',fontsize=14)
ax1.set_ylabel(r'$f$ [Hz]',fontsize=14)
ax2.set_xlabel(r'$I$ [nA]',fontsize=14)
ax2.set_ylabel(r'$f$ [Hz]',fontsize=14)
ax3.set_xlabel(r'$I$ [nA]',fontsize=14)
ax3.set_ylabel(r'$f$ [Hz]',fontsize=14)

model_folders = ['','CaHVA_Allen/Shift_gCaHVA_V/Vshift_0/','CaHVA_Allen/SK_Allen/Shift_gCaHVA_V/Vshift_0/']
labels = ['base','CaHVA',r'0.2$\cdot$CaHVA, SK']
Nmodels = len(model_folders)

outfolder = 'Compare/Soma%i/' % somasize
plotname = outfolder+'fI_naf_kaf_CaHVA_SK_Allen_skiptime'+str(skiptime)+'_commonprops.png'


# Namestrings
namestring1 = '_gCaHVA1.0p_'   
namestring2 = '_gSK1.0p_gCaHVA0.2p_'   

# Nspikes
filename_Nspikes_base     = 'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'
filename_Nspikes_CaHVA    = 'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring1+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'
filename_Nspikes_CaHVA_SK = 'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring2+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'
# Gathering
filenames_Nspikes = [filename_Nspikes_base,filename_Nspikes_CaHVA,filename_Nspikes_CaHVA_SK]


Vrest = numpy.zeros(Nmodels)
Vrest_rms = numpy.zeros(Nmodels)

for i in range(Nmodels):
    Nspikes      = []
    iamp_Nspikes = []
    
    # Get names
    infolder = model_folders[i]+'Results/Soma%i/' % somasize
    infilename_Nspikes = infolder+filenames_Nspikes[i]
    
    # Open files
    infile_Nspikes = open(infilename_Nspikes,'r')
    
    lines_Nspikes = infile_Nspikes.readlines()
    
    for line in lines_Nspikes:
        words = line.split()
        if len(words)>0:
            iamp_Nspikes.append(float(words[0]))
            Nspikes.append(float(words[1]))
    
    ax1.plot(iamp_Nspikes,Nspikes,label=labels[i])
    
    infile_Nspikes.close()

ax1.legend(loc='upper left')#,ncol=1)

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
gcahvas = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,2.0,2.5]
gpases  = 1.0
Ng      = len(gcahvas)

colors = []
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))

if vary_naf==True:
    outnamestring = outnamestring + '_gnaf'
if vary_kaf==True:
    outnamestring = outnamestring + '_gkaf'
if vary_Ca_HVA==True:
    outnamestring = outnamestring + '_gCaHVA'
if vary_gpas==True: 
    outnamestring = outnamestring + '_gpas'
outnamestring = outnamestring +'_'

folder    = 'CaHVA_Allen/Shift_gCaHVA_V/Vshift_'+str(vshift)+'/Results/Soma%i/' % somasize
outfolder = 'Compare/Soma10/'
figname   = outfolder+'mainmodel_many_fI.png'

i = 0
for gcahva in gcahvas:
    namestring = ''
    thelabel = ''
    if vary_naf==True:
        namestring = namestring + '_gnaf'+str(gnaf)+'p'
        #thelabel   = thelabel + r'%.1f$\bar{g}_\mathregular{naf}$, ' % gnaf
    if vary_kaf==True:
        namestring = namestring + '_gkaf'+str(gkaf)+'p'
        #thelabel   = thelabel + r'%.1f$\bar{g}_\mathregular{kaf}$, ' % gkaf 
    if vary_Ca_HVA==True:
        namestring = namestring + '_gCaHVA'+str(gcahva)+'p'    
        thelabel   = thelabel + r'%.1f$\bar{g}_\mathregular{CaHVA}$, ' % gcahva
    if vary_gpas==True: 
        namestring = namestring + '_gpas'+str(gpas)+'p'
        #thelabel   = thelabel + r'%.1f$\bar{g}_\mathregular{L}$, ' % gpas
    namestring = namestring +'_'
    
    # Set names
    infilename_Nspikes = folder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'
    # make files    
    infile_Nspikes = open(infilename_Nspikes,'r')
    lines_Nspikes = infile_Nspikes.readlines()
    
    I_Nspikes = []
    Nspikes   = []
    
    for line in lines_Nspikes:
        words = line.split()
        if len(words)>0:
            I_Nspikes.append(float(words[0]))
            Nspikes.append(float(words[1]))
    
    infile_Nspikes.close()
    
    if gcahva==1.0:
        thislinewidth = 3
    else:
        thislinewidth = 2
    
    if numpy.sum(Nspikes)!=0:
        ax2.plot(I_Nspikes,Nspikes, label=thelabel,linewidth=thislinewidth,color=colors[i])
    i+=1
    
ax2.legend(loc='lower right',ncol=2)


vary_SK      = True # False # 
vary_Ca_HVA  = True # False
    
gnafs   = 1.0
gkafs   = 1.0
gsks    = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,2.0,2.5]
gcahva  = 0.2
gpases  = 1.0
Ng      = len(gsks)
    
colors = []
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))

folder    = 'CaHVA_Allen/SK_Allen/Shift_gCaHVA_V/Vshift_'+str(vshift)+'/Results/Soma%i/' % somasize    
i = 0
for gsk in gsks:
    namestring = ''
    thelabel = ''
    if vary_naf==True:
        namestring = namestring + '_gnaf'+str(gnaf)+'p'
        #thelabel   = thelabel + r'%.1f$\bar{g}_\mathregular{naf}$, ' % gnaf
    if vary_kaf==True:
        namestring = namestring + '_gkaf'+str(gkaf)+'p'
        #thelabel   = thelabel + r'%.1f$\bar{g}_\mathregular{kaf}$, ' % gkaf    
    if vary_SK==True:
        namestring = namestring + '_gSK'+str(gsk)+'p'
        thelabel   = thelabel + r'%.1f$\bar{g}_\mathregular{SK}$' % gsk
    if vary_Ca_HVA==True:
        namestring = namestring + '_gCaHVA'+str(gcahva)+'p'    
        #thelabel   = thelabel + r'%.1f$\bar{g}_\mathregular{CaHVA}$, ' % gcahva
    if vary_gpas==True: 
        namestring = namestring + '_gpas'+str(gpas)+'p'
        #thelabel   = thelabel + r'%.1f$\bar{g}_\mathregular{L}$, ' % gpas
    namestring = namestring +'_'
    
    # Set names
    infilename_Nspikes = folder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'
    # make files    
    infile_Nspikes = open(infilename_Nspikes,'r')
    lines_Nspikes = infile_Nspikes.readlines()
    
    I_Nspikes = []
    Nspikes   = []
    
    for line in lines_Nspikes:
        words = line.split()
        if len(words)>0:
            I_Nspikes.append(float(words[0]))
            Nspikes.append(float(words[1]))
    
    infile_Nspikes.close()
    
    if gsk==1.0:
        thislinewidth = 3
    else:
        thislinewidth = 2
    
    if numpy.sum(Nspikes)!=0:
        ax3.plot(I_Nspikes,Nspikes, label=thelabel,linewidth=thislinewidth,color=colors[i])
    i+=1

ax3.legend(loc='lower right',ncol=2)
plt.tight_layout()
plt.savefig(figname)
#plt.show()
