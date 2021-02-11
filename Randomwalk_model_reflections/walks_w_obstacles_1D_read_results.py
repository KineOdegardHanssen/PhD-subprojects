from pylab import *
import matplotlib.pyplot as plt  # To plot
import numpy as np
import math

hitprobs = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9] # Probability of hitting obstacle (and be sent back)
Nprobs   = len(hitprobs)
Ds       = np.zeros(Nprobs)
Nsteps   = 100   # Increase #Maybe for other values of hitprob
Nreal    = 100000

outfilename = 'D_vs_hitprobs_Nsteps%i_Nreal%i.txt' %(Nsteps,Nreal)
plotname = 'D_vs_hitprobs_Nsteps%i_Nreal%i.png' %(Nsteps,Nreal)
outfile = open(outfilename,'w')

for i in range(Nprobs):
    infilename = 'hitprob'+str(hitprobs[i])+'_Nsteps%i_Nreal%i_D.txt' %(Nsteps,Nreal)
    infile = open(infilename,'r')
    line = infile.readline()
    Ds[i] = float(line.split()[0])
    outfile.write('%i %.16f\n' % (hitprobs[i],Ds[i]))
outfile.close()

plt.figure(figsize=(6,5))
plt.plot(hitprobs,Ds)
plt.xlabel(r'Probability of hitting obstacle')
plt.ylabel(r'Diffusion constant $D$')
plt.title(r'Probability of hitting obstacle vs $D$')
plt.tight_layout()
plt.show()
plt.savefig(plotname)