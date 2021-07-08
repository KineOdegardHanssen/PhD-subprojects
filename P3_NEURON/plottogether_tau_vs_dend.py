import numpy as np
import matplotlib.pyplot as plt
import math

# Later on: Maybe include option to vary the dendrite length too. Possibly do that in another script.

zoomed = False
somasize = 10
dendlen  = 10
denddiams = [0,1,2,5,10,20]
Ra     = 100
v_init = -70
N = len(denddiams)

idur = 1000   # ms
iamp = -0.1   # nA 
currentfolder = 'current_idur%i_iamp'%idur+str(iamp)+'/'

filenames = []
somaonly_folder = 'Somaonly/pas/Results/IStim/Soma10/'+currentfolder
filename_nodend = somaonly_folder+'somaonly_cms_idur%i_iamp'%idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas_tau_highCms.txt'
filenames.append(filename_nodend)
plotfolder = 'Comparemodels/BAS_vs_somaonly_passive/'
plotname = plotfolder + 'BAS_vs_onecomp_cms_idur%i_iamp'%idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas_tau'
if zoomed==True:
    plotname=plotname+'_zoomed.png'
else:
    plotname=plotname+'.png'

for i in range(1,N):
    denddiam = denddiams[i]
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
            cms.append(float(words[0]))
            taus.append(float(words[1]))
    if zoomed==False:
        plt.plot(cms,taus,label='Dendrite diam: %i' % denddiams[i])
    else:
        plt.plot(cms[2:7],taus[2:7],label='Dendrite diam: %i' % denddiams[i])
    file.close()
plt.xlabel(r'$C_m$ [$\mu$F/cm$^2$] (Cell parameter)')
plt.ylabel(r'$\tau_m$ [ms] (Measured value)')
plt.title(r'$\tau_m$ vs $C_m$, passive cell, ball-and-stick')
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig(plotname)
plt.show()
    