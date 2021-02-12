from pylab import *
import matplotlib.pyplot as plt  # To plot
import numpy as np
import math

hitprobs = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9] # Probability of hitting obstacle (and be sent back)
Nprobs   = len(hitprobs)
Ds       = np.zeros(Nprobs)
ds       = np.zeros(Nprobs)
ds2      = np.zeros(Nprobs)
Nsteps   = 100   # Increase #Maybe for other values of hitprob
Nreal    = 10000

if hitprobs[0]==0:
    outfilename = 'D_vs_hitprobs_Nsteps%i_Nreal%i_2D_with0.txt' %(Nsteps,Nreal)
    plotname = 'D_vs_hitprobs_Nsteps%i_Nreal%i_2D_with0.png' %(Nsteps,Nreal)
else:
    outfilename = 'D_vs_hitprobs_Nsteps%i_Nreal%i_2D.txt' %(Nsteps,Nreal)
    plotname = 'D_vs_hitprobs_Nsteps%i_Nreal%i_2D.png' %(Nsteps,Nreal)
plotname2 = 'D_vs_1dsqrthitprobs_Nsteps%i_Nreal%i_2D.png' %(Nsteps,Nreal)
plotname3 = 'D_vs_1dhitprobs_Nsteps%i_Nreal%i_2D.png' %(Nsteps,Nreal)
outfile = open(outfilename,'w')

for i in range(Nprobs):
    infilename = '2D_hitprob'+str(hitprobs[i])+'_Nsteps%i_Nreal%i_D.txt' %(Nsteps,Nreal)
    infile = open(infilename,'r')
    line = infile.readline()
    Ds[i] = float(line.split()[0])
    if hitprobs[0]!=0:
        ds[i] = 1/np.sqrt(hitprobs[i])
        ds2[i] = 1/hitprobs[i]
    outfile.write('%i %.16f\n' % (hitprobs[i],Ds[i]))
outfile.close()

plt.figure(figsize=(6,5))
plt.plot(hitprobs,Ds)
plt.xlabel(r'Probability of hitting obstacle')
plt.ylabel(r'Diffusion constant $D$')
plt.title(r'Probability of hitting obstacle vs $D$')
plt.tight_layout()
plt.savefig(plotname)

if hitprobs[0]!=0:
    plt.figure(figsize=(6,5))
    plt.plot(ds,Ds)
    plt.xlabel(r'$1/\sqrt{hitprob}$')
    plt.ylabel(r'Diffusion constant $D$')
    plt.title(r'$1/\sqrt{hitprob}$ vs $D$')
    plt.tight_layout()
    plt.savefig(plotname2)
    
    plt.figure(figsize=(6,5))
    plt.plot(ds2,Ds)
    plt.xlabel(r'$1/hitprob$')
    plt.ylabel(r'Diffusion constant $D$')
    plt.title(r'$1/hitprob$ vs $D$')
    plt.tight_layout()
    plt.savefig(plotname3)
plt.show()