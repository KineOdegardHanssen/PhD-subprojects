from pylab import *
import matplotlib.pyplot as plt  # To plot
import numpy as np
import math

r_hardsphere = 0.8
reflfac      = 1.2
Nsteps       = 1000
Nreal        = 10000

ds = np.array([2,3,4,5,6,7,8,10,15,25,50,75,100])
Nd = len(ds)
Ds = np.zeros(Nd)
Ds_rms = np.zeros(Nd)

outfilename = 'D_vs_d_Nsteps%i_Nreal%i_rsphere'+str(r_hardsphere)+'_norefl.txt'
plotname    = 'D_vs_d_Nsteps%i_Nreal%i_rsphere'+str(r_hardsphere)+'_norefl.png'
outfile     = open(outfilename,'w')

for i in range(Nd):
    d = ds[i]
    infilename_D = '2D_d'+str(d)+'_rsphere'+str(r_hardsphere)+'_Nsteps%i_Nreal%i_norefl_D.txt' %(Nsteps,Nreal)
    infile   = open(infilename_D,'r')
    line     = infile.readline()
    words    = line.split()
    Ds[i]    = float(words[0])
    Ds_rms[i]= float(words[1])
    outfile.write('%i %.16f %.16f\n' % (d,Ds[i],Ds_rms[i]))
    infile.close()
outfile.close()

plt.figure(figsize=(6,5))
plt.plot(ds,Ds,'-o')
plt.xlabel(r'Spacing $d$')
plt.ylabel(r'Diffusion constant $D$')
plt.title(r'$D$ vs $d$')
plt.tight_layout()
plt.savefig(plotname)
plt.show()