import numpy as np
import matplotlib.pyplot as plt
import math

varydiam = False

zoomed = False
somasize = 10
dendlen  = 10
denddiam = 10
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

plotname = plotfolder + 'BAS_vs_onecomp_cms_idur%i_iamp'%idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas'+otherparamsn+'_tau'+endsnippet
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
plt.show()

if varydiam==True:
    print('Difference, one compartment and thickest dendrite (tau0/tau_thick):', taus_at_1[0]/taus_at_1[-1])
    print('thickness 1:', denddiams[0], 'thickness 2:', denddiams[-1])
else:
    print('Difference, one compartment and longest dendrite (tau0/tau_long):', taus_at_1[0]/taus_at_1[-1])
    print('length 1:', dendlens[0], 'length 2:', dendlens[-1])
