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
moresigmas    = True
confignrs     = np.arange(1,1001)
filestext     = '_seed'+str(confignrs[0])+'to'+str(confignrs[-1])

endlocation = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/'
if moresigmas==True:
    endlocation = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Moresigmas/'

basepath_base      = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/'
endlocation_static = basepath_base+'D_vs_d/'

## Files to read
brushfilename_dyn  = endlocation + 'D_vs_d.txt'
brushfilename_stat = endlocation_static + 'D_vs_d_static.txt'
## Files to write to
plotname     = endlocation_static+'Dd_dyn_vs_stat.png'
plotname_cut = endlocation_static+'Dd_dyn_vs_stat_cut.png'


# Dynamic sims:
brushfile_dyn = open(brushfilename_dyn, 'r')

lines = brushfile_dyn.readlines()
N_dyn = len(lines)-1

# ds
spacings_dyn = np.zeros(N_dyn)
# Ds
DRs_dyn = np.zeros(N_dyn)
Dxs_dyn = np.zeros(N_dyn)
Dys_dyn = np.zeros(N_dyn)
Dzs_dyn = np.zeros(N_dyn)
Dparallel_dyn = np.zeros(N_dyn)
# Ds, stdv
DRs_stdv_dyn = np.zeros(N_dyn)
Dxs_stdv_dyn = np.zeros(N_dyn)
Dys_stdv_dyn = np.zeros(N_dyn)
Dzs_stdv_dyn = np.zeros(N_dyn)
Dparallel_stdv_dyn = np.zeros(N_dyn)

for i in range(1,N_dyn+1):
    words = lines[i].split()
    j = i-1
    
    spacings_dyn[j] = float(words[0])
    DRs_dyn[j] = float(words[1])
    Dzs_dyn[j] = float(words[5])
    Dparallel_dyn[j] = float(words[9])
    # Ds, stdv
    DRs_stdv_dyn[j] = float(words[2])
    Dzs_stdv_dyn[j] = float(words[6])
    Dparallel_stdv_dyn[j] = float(words[10])
    
brushfile_dyn.close()

##########

# Static sims:
brushfile_stat = open(brushfilename_stat, 'r')

lines = brushfile_stat.readlines()
N_stat = len(lines)-1

# ds
spacings_stat = np.zeros(N_stat)
# Ds
DRs_stat = np.zeros(N_stat)
Dxs_stat = np.zeros(N_stat)
Dys_stat = np.zeros(N_stat)
Dzs_stat = np.zeros(N_stat)
Dparallel_stat = np.zeros(N_stat)
# Ds, stdv
DRs_stdv_stat = np.zeros(N_stat)
Dxs_stdv_stat = np.zeros(N_stat)
Dys_stdv_stat = np.zeros(N_stat)
Dzs_stdv_stat = np.zeros(N_stat)
Dparallel_stdv_stat = np.zeros(N_stat)

for i in range(1,N_stat+1):
    words = lines[i].split()
    j = i-1
    
    spacings_stat[j] = float(words[0])
    DRs_stat[j] = float(words[1])
    Dzs_stat[j] = float(words[3])
    Dparallel_stat[j] = float(words[5])
    # Ds, stdv
    DRs_stdv_stat[j] = float(words[2])
    Dzs_stdv_stat[j] = float(words[4])
    Dparallel_stdv_stat[j] = float(words[6])
    
brushfile_stat.close()


plt.figure(figsize=(6,5))
plt.errorbar(spacings_dyn, DRs_dyn, yerr=DRs_stdv_dyn, capsize=2, label=r'$D_R$, dyn.')
plt.errorbar(spacings_dyn, Dzs_dyn, yerr=Dzs_stdv_dyn, capsize=2, label=r'$D_\perp$, dyn.')
plt.errorbar(spacings_dyn, Dparallel_dyn, yerr=Dparallel_stdv_dyn, capsize=2, label=r'$D_\parallel$, dyn.')
plt.errorbar(spacings_stat, DRs_stat, yerr=DRs_stdv_stat, capsize=2, label=r'$D_R$, stat.')
plt.errorbar(spacings_stat, Dzs_stat, yerr=Dzs_stdv_stat, capsize=2, label=r'$D_\perp$, stat.')
plt.errorbar(spacings_stat, Dparallel_stat, yerr=Dparallel_stdv_stat, capsize=2, label=r'$D_\parallel$, stat.')
if moresigmas==True:
    plt.xlabel(r'$d/\sigma_b$')
else:
    plt.xlabel(r'$d$')
plt.ylabel(r'Diffusion constant $D$')

if moresigmas==True:
    plt.title('Diffusion constant $D$ vs $d/\sigma_b$ for dynamic and static brushes')
else:
    plt.title('Diffusion constant $D$ vs $d$ for dynamic and static brushes')
plt.tight_layout()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='upper left')
plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.errorbar(spacings_dyn, DRs_dyn, yerr=DRs_stdv_dyn, capsize=2, label=r'$D_R$, dyn.')
plt.errorbar(spacings_dyn, Dzs_dyn, yerr=Dzs_stdv_dyn, capsize=2, label=r'$D_\perp$, dyn.')
plt.errorbar(spacings_dyn, Dparallel_dyn, yerr=Dparallel_stdv_dyn, capsize=2, label=r'$D_\parallel$, dyn.')
plt.errorbar(spacings_stat, DRs_stat, yerr=DRs_stdv_stat, capsize=2, label=r'$D_R$, stat.')
plt.errorbar(spacings_stat, Dzs_stat, yerr=Dzs_stdv_stat, capsize=2, label=r'$D_\perp$, stat.')
#plt.errorbar(spacings_stat, Dparallel_stat, yerr=Dparallel_stdv_stat, capsize=2, label=r'$D_\parallel$, stat.')
if moresigmas==True:
    plt.xlabel(r'$d/\sigma_b$')
else:
    plt.xlabel(r'$d$')
plt.ylabel(r'Diffusion constant $D$')

if moresigmas==True:
    plt.title('Diffusion constant $D$ vs $d/\sigma_b$ for dynamic and static brushes')
else:
    plt.title('Diffusion constant $D$ vs $d$ for dynamic and static brushes')
plt.axis([0,10,0,6e-7])
plt.tight_layout()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='upper left')
plt.savefig(plotname_cut)

