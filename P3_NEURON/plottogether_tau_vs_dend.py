import numpy as np
import matplotlib.pyplot as plt
import math

varydiam = True # False # 

zoomed = False
somasize = 10
dendlen  = 100 # 1 # 2 # 5 # 10 # 20 # 50 # 
denddiam = 1 # 2 # 5 # 10 # 20 # 
if varydiam==True:
    denddiams = [0,1,2,5,10,20]
    N = len(denddiams)
else:
    dendlens = [0,1,2,5,10,20,50,100]
    N = len(dendlens)
    
Ra     = 100
v_init = -70

idur = 1000   # ms
iamp = -0.1   # nA 
currentfolder = 'current_idur%i_iamp'%idur+str(iamp)+'/'

taus_at_1 = []
taus_at_2 = []
taus_at_5 = []
filenames = []
somaonly_folder = 'Somaonly/pas/Results/IStim/Soma10/'+currentfolder
filename_nodend = somaonly_folder+'somaonly_cms_idur%i_iamp'%idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas_tau_highCms.txt'
filenames.append(filename_nodend)
plotfolder = 'Comparemodels/BAS_vs_somaonly_passive/'

if varydiam==True:
    endsnippet = '_varydenddiam'
    otherparamsn = '_len%i' % dendlen
else:
    endsnippet = '_varydendlen'
    otherparamsn = '_diam%i' % denddiam

plotname   = plotfolder + 'BAS_vs_onecomp_cms_idur%i_iamp'%idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas'+otherparamsn+'_tau'+endsnippet
plotname_2 = plotfolder + 'BAS_vs_onecomp_cm1_idur%i_iamp'%idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas'+otherparamsn+'_tau'+endsnippet+'_percentagediff.png'
plotname_3 = plotfolder + 'BAS_cm1and2and5_idur%i_iamp'%idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas'+otherparamsn+'_tau'+endsnippet+'_percentagediff.png'
if zoomed==True:
    plotname=plotname+'_zoomed.png'
else:
    plotname=plotname+'.png'

for i in range(1,N):
    if varydiam==True:
        denddiam = denddiams[i]
    else:
        dendlen = dendlens[i]
    folder = 'Ball-and-stick models/BAS_passive/Results/IStim/Soma%i/dendlen%i/denddiam%i/' % (somasize,dendlen,denddiam)+currentfolder
    filename = folder + 'baspass_cms_idur%i_iamp' % idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas_tau_highCms.txt'
    filenames.append(filename)

plt.figure(figsize=(6,5))
for i in range(N):
    file = open(filenames[i],'r')
    lines = file.readlines()
    
    cms = []
    taus = []
    for line in lines:
        words = line.split()
        if len(words)>1:
            cm = float(words[0])
            tau = float(words[1])
            cms.append(cm)
            taus.append(tau)
            if cm==1:
                taus_at_1.append(tau)
            elif cm==2:
                taus_at_2.append(tau)
            elif cm==5:
                taus_at_5.append(tau)
    if zoomed==False:
        if varydiam==True:
            plt.plot(cms,taus,label='Dendrite diam: %i' % denddiams[i])
        else:
            plt.plot(cms,taus,label='Dendrite len: %i' % dendlens[i])
    else:
        if varydiam==True:
            plt.plot(cms[2:7],taus[2:7],label='Dendrite diam: %i' % denddiams[i])
        else:
            plt.plot(cms[2:7],taus[2:7],label='Dendrite len: %i' % dendlens[i])
    file.close()
plt.xlabel(r'$C_m$ [$\mu$F/cm$^2$] (Cell parameter)')
plt.ylabel(r'$\tau_m$ [ms] (Measured value)')
plt.title(r'$\tau_m$ vs $C_m$, passive cell, ball-and-stick')
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig(plotname)

fracs = np.zeros(N-1)
fracs12 = np.zeros(N)
fracs15 = np.zeros(N)
plt.figure(figsize=(6,5))
plt.ylabel('% difference')
if varydiam==True:
    print('dendlen=', dendlen, ', varydiam:')
    print('Difference, one compartment and thickest dendrite (tau0/tau_thick):', taus_at_1[0]/taus_at_1[-1])
    print('Relative difference:', (taus_at_1[0]-taus_at_1[-1])/taus_at_1[0])
    print('thickness 1:', denddiams[0], 'thickness 2:', denddiams[-1])
    for i in range(1,N):
        fracs[i-1] = 100*(taus_at_1[0]-taus_at_1[i])/taus_at_1[0] # Will be 0 at i=0
    plt.plot(denddiams[1:],fracs)
    plt.xlabel('Dendrite diameter (nm)')
    plt.title(r'Difference, $\tau$, one comp. and BAS, dend.len. %i' % dendlen)
    plt.tight_layout()
    plt.savefig(plotname_2)
    # New plot
    for i in range(N):
        fracs12[i] = 100*(taus_at_1[i]-taus_at_2[i])/taus_at_1[i]
        fracs15[i] = 100*(taus_at_1[i]-taus_at_5[i])/taus_at_1[i]
    plt.figure(figsize=(6,5))
    plt.plot(denddiams,fracs12,label=r'$C_m=1$ vs $C_m=2$')
    plt.plot(denddiams,fracs15,label=r'$C_m=1$ vs $C_m=5$')
    plt.xlabel('Dendrite diameter (nm)')
    plt.ylabel('% difference')
    plt.title(r'Difference, $\tau$, one comp. and BAS, dend.len. %i' % dendlen)
    plt.tight_layout()
    plt.legend()
    plt.savefig(plotname_3)
else:
    print('denddiam=', denddiam, ', varylen:')
    print('Difference, one compartment and longest dendrite (tau0/tau_long):', taus_at_1[0]/taus_at_1[-1])
    print('Relative difference:', (taus_at_1[0]-taus_at_1[-1])/taus_at_1[0])
    print('length 1:', dendlens[0], 'length 2:', dendlens[-1])
    for i in range(1,N):
        fracs[i-1] = 100*(taus_at_1[0]-taus_at_1[i])/taus_at_1[0] # Will be 0 at i=0
    plt.plot(dendlens[1:],fracs)
    plt.xlabel('Dendrite length (nm)')
    plt.title(r'Difference, $\tau$, one comp. and BAS, dend.diam. %i' % denddiam)
    plt.tight_layout()
    plt.savefig(plotname_2)
    # New plot
    for i in range(N):
        fracs12[i] = 100*(taus_at_1[i]-taus_at_2[i])/taus_at_1[i]
        fracs15[i] = 100*(taus_at_1[i]-taus_at_5[i])/taus_at_1[i]
    plt.figure(figsize=(6,5))
    plt.plot(dendlens,fracs12,label=r'$C_m=1$ vs $C_m=2$')
    plt.plot(dendlens,fracs15,label=r'$C_m=1$ vs $C_m=5$')
    plt.xlabel('Dendrite length (nm)')
    plt.ylabel('% difference')
    plt.title(r'Difference, $\tau$, one comp. and BAS, dend.diam. %i' % denddiam)
    plt.tight_layout()
    plt.legend()
    plt.savefig(plotname_3)
plt.show()
