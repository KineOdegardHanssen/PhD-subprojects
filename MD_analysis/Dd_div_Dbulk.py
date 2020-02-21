from skimage import morphology, measure, io, util   # For performing morphology operations (should maybe import more from skimage...)
from mpl_toolkits.mplot3d import Axes3D             # Plotting in 3D
import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
from scipy.ndimage import measurements, convolve    # Should have used the command below, but changing now will lead to confusion
from scipy import ndimage                           # For Euclidean distance measurement
import maintools_percolation as perctools
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
confignrs     = np.arange(1,1001)
filestext     = '_seed'+str(confignrs[0])+'to'+str(confignrs[-1])

endlocation = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/'

## Files to read
bulkfilename  = endlocation + 'diffusion_bulk'+filestext+'.txt'
brushfilename = endlocation + 'D_vs_d.txt'
#
bulkfile  = open(bulkfilename, 'r')
brushfile = open(brushfilename, 'r')
## Files to write to
outfilename  = endlocation+'Dd_div_Dbulk_vs_d.txt'
plotname     = endlocation+'Dd_div_Dbulk_vs_d.png'
#
outfile = open(outfilename, 'w')

# Write header
outfile.write('d   D_R2/Dbulk   sigmaD_R2/Dbulk; D_z2/Dbulk  sigmaD_z2/Dbulk  sigmaD_z2/Dbulk; D_par2/Dbulk sigmaD_par2/Dbulk\n')


# D_R2  sigmaD_R2 b_R2 sigmab_R2; D_z2  sigmaD_z2  b_z2  sigmaD_z2; D_par2 sigmaD_par2  b_par2  sigmab_par2
bulklines = bulkfile.readlines()
bulkline  = bulklines[1]
words     = bulkline.split()

# Ds
DRs_bulk = float(words[0])
Dzs_bulk = float(words[4])
Dparallel_bulk = float(words[8])

# Ds, stdv
DRs_stdv_bulk = float(words[1])
Dzs_stdv_bulk = float(words[5])
Dparallel_stdv_bulk = float(words[9])

bulkfile.close()

lines = brushfile.readlines()
N = len(lines)

# ds
spacings = np.zeros(N)
# Ds
DRs = np.zeros(N)
Dxs = np.zeros(N)
Dys = np.zeros(N)
Dzs = np.zeros(N)
Dparallel = np.zeros(N)
# Ds, stdv
DRs_stdv = np.zeros(N)
Dxs_stdv = np.zeros(N)
Dys_stdv = np.zeros(N)
Dzs_stdv = np.zeros(N)
Dparallel_stdv = np.zeros(N)


for i in range(1,N):
    words = lines[i].split()
    
    spacings[i] = float(words[0])
    DRs[i] = float(words[1])/DRs_bulk
    Dzs[i] = float(words[5])/Dzs_bulk
    Dparallel[i] = float(words[9])/Dparallel_bulk
    # Ds, stdv
    DRs_stdv[i] = abs(DRs_bulk)*np.sqrt((float(words[2])/DRs[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dzs_stdv[i] = abs(Dzs_bulk)*np.sqrt((float(words[6])/Dzs[i])**2+(Dzs_stdv_bulk/Dzs_bulk)**2)
    Dparallel_stdv[i] = abs(Dparallel_bulk)*np.sqrt((float(words[10])/Dparallel[i])**2+(Dparallel_stdv_bulk/Dparallel_bulk)**2)
    
    outfile.write('%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n' % (spacings[i], DRs[i], DRs_stdv[i], Dzs[i], Dzs_stdv[i], Dparallel[i], Dparallel_stdv[i]))


brushfile.close()


plt.figure(figsize=(6,5))
plt.errorbar(spacings, DRs, yerr=DRs_stdv, capsize=2, label=r'$D_R/D_{bulk,R}$')
plt.errorbar(spacings, Dzs, yerr=Dzs_stdv, capsize=2, label=r'$D_\perp/D_{bulk,\perp}$')
plt.errorbar(spacings, Dparallel, yerr=Dparallel_stdv, capsize=2, label=r'$D_\parallel/D_{bulk,\parallel}$')
plt.xlabel(r'$d$')
plt.ylabel(r'Diffusion constant $D/D_{bulk}$')
plt.title('Diffusion constant $D/D_{bulk}$ vs $d$')
plt.tight_layout()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='upper left')
plt.savefig(plotname)

