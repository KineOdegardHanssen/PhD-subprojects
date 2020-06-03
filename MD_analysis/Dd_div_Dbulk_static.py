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
bulk_cut      = True
confignrs     = np.arange(1,1001)
filestext     = '_seed'+str(confignrs[0])+'to'+str(confignrs[-1])

endlocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/'

basepath_base      = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/'
endlocation_static = basepath_base+'D_vs_d/'

## Files to read
bulkfilename  = endlocation + 'diffusion_bulk'+filestext
if bulk_cut==True:
    bulkfilename = bulkfilename +'_cut.txt'
else:
    bulkfilename = bulkfilename +'_uncut.txt'

brushfilename = endlocation_static + 'D_vs_d.txt'
#
bulkfile  = open(bulkfilename, 'r')
brushfile = open(brushfilename, 'r')
## Files to write to
outfilename  = endlocation_static+'Dd_div_Dbulk_vs_d_static'
plotname     = endlocation_static+'Dd_div_Dbulk_vs_d_static'
outfilename_2  = endlocation_static+'Dd_div_Dbulk_vs_d_static_2'
plotname_2     = endlocation_static+'Dd_div_Dbulk_vs_d_static_2'
outfilename_3  = endlocation_static+'D_vs_d_static.txt'
plotname_3     = endlocation_static+'D_vs_d_static.png'
if bulk_cut==True:
    outfilename   = outfilename   + '_cut.txt'
    outfilename_2 = outfilename_2 + '_cut.txt'
    plotname      = plotname      + '_cut.png'
    plotname_2    = plotname_2    + '_cut.png'
else:
    outfilename   = outfilename   + '_uncut.txt'
    outfilename_2 = outfilename_2 + '_uncut.txt'
    plotname      = plotname      + '_uncut.png'
    plotname_2    = plotname_2    + '_uncut.png'
#
outfile   = open(outfilename, 'w')
outfile_2 = open(outfilename_2, 'w')
outfile_3 = open(outfilename_3, 'w')

# Write header
outfile.write('d   D_R2/DR_bulk   sigmaD_R2/DR_bulk; D_z2/Dz_Dbulk  sigmaD_z2/Dz_Dbulk  sigmaD_z2/Dz_Dbulk; D_par2/Dpar_bulk sigmaD_par2/Dpar_bulk\n')
outfile_2.write('d   D_R2/Dbulk   sigmaD_R2/Dbulk; D_z2/Dbulk  sigmaD_z2/Dbulk  sigmaD_z2/Dbulk; D_par2/Dbulk sigmaD_par2/Dbulk\n')
outfile_3.write('d   D_R2   sigmaD_R2; D_z2  sigmaD_z2  sigmaD_z2; D_par2 sigmaD_par2\n')


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
N = len(lines)-1

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

# Ds, 2
DRs_v2 = np.zeros(N)
Dxs_v2 = np.zeros(N)
Dys_v2 = np.zeros(N)
Dzs_v2 = np.zeros(N)
Dparallel_v2 = np.zeros(N)
# Ds, stdv
DRs_stdv_v2 = np.zeros(N)
Dxs_stdv_v2 = np.zeros(N)
Dys_stdv_v2 = np.zeros(N)
Dzs_stdv_v2 = np.zeros(N)
Dparallel_stdv_v2 = np.zeros(N)


# Ds, undivided
DRs_ud = np.zeros(N)
Dxs_ud = np.zeros(N)
Dys_ud = np.zeros(N)
Dzs_ud = np.zeros(N)
Dparallel_ud = np.zeros(N)
# Ds, stdv
DRs_stdv_ud = np.zeros(N)
Dxs_stdv_ud = np.zeros(N)
Dys_stdv_ud = np.zeros(N)
Dzs_stdv_ud = np.zeros(N)
Dparallel_stdv_ud = np.zeros(N)

for i in range(1,N+1):
    words = lines[i].split()
    j = i-1
    
    spacings[i-1] = float(words[0])
    DRs[i-1] = float(words[1])/DRs_bulk
    Dzs[i-1] = float(words[5])/Dzs_bulk
    Dparallel[i-1] = float(words[9])/Dparallel_bulk
    # Ds, stdv
    DRs_stdv[i-1] = abs(DRs_bulk)*np.sqrt((float(words[2])/DRs[i-1])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dzs_stdv[i-1] = abs(Dzs_bulk)*np.sqrt((float(words[6])/Dzs[i-1])**2+(Dzs_stdv_bulk/Dzs_bulk)**2)
    Dparallel_stdv[i-1] = abs(Dparallel_bulk)*np.sqrt((float(words[10])/Dparallel[i-1])**2+(Dparallel_stdv_bulk/Dparallel_bulk)**2)
    
    # v2
    DRs_v2[i-1] = float(words[1])/DRs_bulk
    Dzs_v2[i-1] = float(words[5])/DRs_bulk
    Dparallel_v2[i-1] = float(words[9])/DRs_bulk
    # Ds, stdv
    DRs_stdv_v2[i-1] = abs(DRs_bulk)*np.sqrt((float(words[2])/DRs_v2[i-1])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dzs_stdv_v2[i-1] = abs(DRs_bulk)*np.sqrt((float(words[6])/Dzs_v2[i-1])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dparallel_stdv_v2[i-1] = abs(DRs_bulk)*np.sqrt((float(words[10])/Dparallel_v2[i-1])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    
    # Undivided
    DRs_ud[i-1] = float(words[1])
    Dzs_ud[i-1] = float(words[5])
    Dparallel_ud[i-1] = float(words[9])
    # Ds, stdv
    DRs_stdv_ud[i-1] = float(words[2])
    Dzs_stdv_ud[i-1] = float(words[6])
    Dparallel_stdv_ud[i-1] = float(words[10])
    
    outfile.write('%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n' % (spacings[j], DRs[j], DRs_stdv[j], Dzs[j], Dzs_stdv[j], Dparallel[j], Dparallel_stdv[j]))
    outfile_2.write('%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n' % (spacings[j], DRs_v2[j], DRs_stdv_v2[j], Dzs_v2[j], Dzs_stdv_v2[j], Dparallel_v2[j], Dparallel_stdv_v2[j]))
    outfile_3.write('%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n' % (spacings[j], DRs_ud[j], DRs_stdv_ud[j], Dzs_ud[j], Dzs_stdv_ud[j], Dparallel_ud[j], Dparallel_stdv_ud[j]))


outfile.close()
outfile_2.close()
outfile_3.close()
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


plt.figure(figsize=(6,5))
plt.errorbar(spacings, DRs_v2, yerr=DRs_stdv_v2, capsize=2, label=r'$D_R/D_{bulk,R}$')
plt.errorbar(spacings, Dzs_v2, yerr=Dzs_stdv_v2, capsize=2, label=r'$D_\perp/D_{bulk,R}$')
plt.errorbar(spacings, Dparallel_v2, yerr=Dparallel_stdv_v2, capsize=2, label=r'$D_\parallel/D_{bulk,R}$')
plt.xlabel(r'$d$')
plt.ylabel(r'Diffusion constant $D/D_{bulk}$')
plt.title('Diffusion constant $D/D_{bulk}$ vs $d$')
plt.tight_layout()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='upper left')
plt.savefig(plotname_2)


plt.figure(figsize=(6,5))
plt.errorbar(spacings, DRs_ud, yerr=DRs_stdv_ud, capsize=2, label=r'$D_R$')
plt.errorbar(spacings, Dzs_ud, yerr=Dzs_stdv_ud, capsize=2, label=r'$D_\perp$')
plt.errorbar(spacings, Dparallel_ud, yerr=Dparallel_stdv_ud, capsize=2, label=r'$D_\parallel$')
plt.xlabel(r'$d$')
plt.ylabel(r'Diffusion constant $D$')
plt.title('Diffusion constant $D$ vs $d$')
plt.tight_layout()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='upper left')
plt.savefig(plotname_3)
