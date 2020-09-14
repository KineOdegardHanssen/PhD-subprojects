import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob

def rmsd(x,y):
    Nx = len(x)
    Ny = len(y)
    if Nx!=Ny:
        print('WARNING! Nx!=Ny. Could not calculate rmsd value')
        return 'WARNING! Nx!=Ny. Could not calculate rmsd value'
    delta = 0
    for i in range(Nx):
        delta += (x[i]-y[i])*(x[i]-y[i])
    delta = np.sqrt(delta/(Nx-1))
    return delta

# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
psigma   = 1 # For instance 
pmass    = 1 # I don't use this anymore, but it got stuck in the file names. Have constant mass density now.
damp     = 10
# Input booleans for file selection:
bulkdiffusion = False
substrate     = False
bulk_cut      = False # hmmm...
confignrs     = np.arange(1,1001)

basepath_base      = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/'
endlocation        = basepath_base+'D_vs_d/Brush/Sigma_bead_' +str(psigma) +'/Nocut/'

## File to read
brushfilename = endlocation + 'D_vs_d.txt'
#
brushfile = open(brushfilename, 'r')
## Files to write to
outfilename  = endlocation+'Dpar_div_Dz_vs_d_dynamic.txt'
plotname     = endlocation+'Dpar_div_Dz_vs_d_dynamic.png'
plotname_smalld = endlocation+'Dpar_div_Dz_vs_d_dynamic_smalld.png'

#
outfile   = open(outfilename, 'w')

lines = brushfile.readlines()
N = len(lines)-1

# ds
spacings = []
# Ds
Dparallel_div_Dzs = np.zeros(N)
# Ds, stdv
Dparallel_div_Dzs_stdv = np.zeros(N)


for i in range(1,N+1):
    words = lines[i].split()
    j = i-1
    
    spacing = float(words[0]) 
    spacings.append(spacing)
    Dz = float(words[5])
    Dparallel = float(words[9])
    
    Dparallel_div_Dzs[j] = Dparallel/Dz
    Dparallel_div_Dzs_stdv[j] = abs(Dparallel_div_Dzs[j])*np.sqrt((float(words[10])/Dparallel)**2+(float(words[6])/Dz)**2)
    
    outfile.write('%.5e %.5e %.5e\n' % (spacing, Dparallel_div_Dzs[j], Dparallel_div_Dzs_stdv[j]))


outfile.close()
brushfile.close()

plt.figure(figsize=(6,5))
plt.errorbar(spacings, Dparallel_div_Dzs, yerr=Dparallel_div_Dzs_stdv, capsize=2)
plt.xlabel(r'$d$')
plt.ylabel(r'$D_{par}/D_z$')
plt.title(r'$D_{par}/D_z$ vs $d$')
plt.tight_layout()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.savefig(plotname)

N10 = spacings.index(10)

plt.figure(figsize=(6,5))
plt.errorbar(spacings[0:N10], Dparallel_div_Dzs[0:N10], yerr=Dparallel_div_Dzs_stdv[0:N10], capsize=2)
plt.xlabel(r'$d$')
plt.ylabel(r'$D_{par}/D_z$')
plt.title(r'$D_{par}/D_z$ vs $d$')
plt.tight_layout()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.savefig(plotname_smalld)
