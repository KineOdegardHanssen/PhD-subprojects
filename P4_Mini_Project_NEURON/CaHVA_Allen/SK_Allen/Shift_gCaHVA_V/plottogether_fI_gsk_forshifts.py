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

vshifts     = [-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50]
spikedurat = -40
idur       = 1000 #100 # ms
idelay     = 10
v_init     = -86 # mV
Ra         = 100
somasize   = 10 # 15 # 
dtexp      = -7
skiptime   = 500
t_before_rec = -600.


# change the default font family
#plt.rcParams.update({'font.family':'Arial'})


for vshift in vshifts:
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=14)
    
    plt.figure(figsize=(6,5),dpi=300)
    
    plt.xlabel(r'$I$ [nA]',fontsize=14)
    plt.ylabel(r'$f$ [Hz]',fontsize=14)
    
    
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
    gsks    = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,2.0,2.5]
    gcahvas = [0.2]
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
        outnamestring = outnamestring + '_gCaHVA'
    if vary_gpas==True: 
        outnamestring = outnamestring + '_gpas'
    outnamestring = outnamestring +'_'
    
    folder    = 'Vshift_'+str(vshift)+'/Results/Soma%i/' % somasize
    outfolder = 'Compare/'
    figname   = outfolder+outnamestring+'_Vshift'+str(vshift)+'Ng'+str(Ng)+'_fI.png'
    
    cm = 1.0
    for gcahva in gcahvas:
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
                plt.plot(I_Nspikes,Nspikes, label=thelabel,linewidth=thislinewidth,color=colors[i])
            i+=1

    plt.legend(loc='lower right',ncol=2)
    plt.tight_layout()
    plt.savefig(figname)
    #plt.show()
