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
bulk_cut      = False
moresigmas    = False
beadplot      = True
confignrs     = np.arange(1,1001)
filestext     = '_seed'+str(confignrs[0])+'to'+str(confignrs[-1])

endlocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_forest/D_vs_d/Nocut/'
if moresigmas==True:
    endlocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_forest/D_vs_d/Brush/Moresigmas/Nocut/'

## Files to read
bulkfilename  = endlocation + 'diffusion_bulk'+filestext
if bulk_cut==True:
    bulkfilename = bulkfilename +'_cut.txt'
else:
    bulkfilename = bulkfilename +'_uncut.txt'

brushfilename = endlocation + 'D_vs_d_forest.txt'
#
bulkfile  = open(bulkfilename, 'r')
brushfile = open(brushfilename, 'r')
## Files to write to
outfilename_2  = endlocation+'Dd_div_Dbulk_vs_d_2_forest'
plotname_2     = endlocation+'Dd_div_Dbulk_vs_d_2_forest'
outfilename_3  = endlocation+'D_vs_d_dynamic_forest.txt'
plotname_3     = endlocation+'D_vs_d_dynamic_forest.png'
plotname_small = endlocation+'Dd_div_Dbulk_vs_d_2_smalld_forest.png'
if bulk_cut==True:
    outfilename_2 = outfilename_2 + '_cut.txt'
    plotname_2    = plotname_2    + '_cut.png'
else:
    outfilename_2 = outfilename_2 + '_uncut.txt'
    plotname_2    = plotname_2    + '_uncut.png'
#
outfile_2 = open(outfilename_2, 'w')
outfile_3 = open(outfilename_3, 'w')

# Write header
if moresigmas==False:
    outfile_2.write('d   D_R2/Dbulk   sigmaD_R2/Dbulk; D_z2/Dbulk  sigmaD_z2/Dbulk  sigmaD_z2/Dbulk; D_par2/Dbulk sigmaD_par2/Dbulk\n')
    outfile_3.write('d   D_R2   sigmaD_R2; D_z2  sigmaD_z2  sigmaD_z2; D_par2 sigmaD_par2\n')
else:
    outfile_2.write('d/sigma_b   D_R2/Dbulk   sigmaD_R2/Dbulk; D_z2/Dbulk  sigmaD_z2/Dbulk  sigmaD_z2/Dbulk; D_par2/Dbulk sigmaD_par2/Dbulk\n')
    outfile_3.write('d/sigma_b   D_R2   sigmaD_R2; D_z2  sigmaD_z2  sigmaD_z2; D_par2 sigmaD_par2\n')

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
    
    spacings[j] = float(words[0])
    
    # v2
    DRs_v2[j] = float(words[1])/DRs_bulk
    Dzs_v2[j] = float(words[5])/DRs_bulk
    Dparallel_v2[j] = float(words[9])/DRs_bulk
    # Ds, stdv
    DRs_stdv_v2[j] = abs(DRs_bulk)*np.sqrt((float(words[2])/DRs_v2[j])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dzs_stdv_v2[j] = abs(DRs_bulk)*np.sqrt((float(words[6])/Dzs_v2[j])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dparallel_stdv_v2[j] = abs(DRs_bulk)*np.sqrt((float(words[10])/Dparallel_v2[j])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    
    # Undivided
    DRs_ud[j] = float(words[1])
    Dzs_ud[j] = float(words[5])
    Dparallel_ud[j] = float(words[9])
    # Ds, stdv
    DRs_stdv_ud[j] = float(words[2])
    Dzs_stdv_ud[j] = float(words[6])
    Dparallel_stdv_ud[j] = float(words[10])
    
    outfile_2.write('%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n' % (spacings[j], DRs_v2[j], DRs_stdv_v2[j], Dzs_v2[j], Dzs_stdv_v2[j], Dparallel_v2[j], Dparallel_stdv_v2[j]))
    outfile_3.write('%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n' % (spacings[j], DRs_ud[j], DRs_stdv_ud[j], Dzs_ud[j], Dzs_stdv_ud[j], Dparallel_ud[j], Dparallel_stdv_ud[j]))

outfile_2.close()
outfile_3.close()
brushfile.close()



plt.figure(figsize=(6,5))
plt.errorbar(spacings, DRs_v2, yerr=DRs_stdv_v2, capsize=2, label=r'$D_R/D_{bulk,R}$')
plt.errorbar(spacings, Dzs_v2, yerr=Dzs_stdv_v2, capsize=2, label=r'$D_\perp/D_{bulk,R}$')
plt.errorbar(spacings, Dparallel_v2, yerr=Dparallel_stdv_v2, capsize=2, label=r'$D_\parallel/D_{bulk,R}$')
if moresigmas==False:
    plt.xlabel(r'$d$')
else:
    plt.xlabel(r'$d/\sigma_b$')
plt.ylabel(r'Diffusion constant $D/D_{bulk}$')
if moresigmas==False:
    plt.title('Diffusion constant $D/D_{bulk}$ vs $d$')
else:
    plt.title('Diffusion constant $D/D_{bulk}$ vs $d/\sigma_b$')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='lower right')
###box = ax.get_position()
###ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
###ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.tight_layout()
plt.savefig(plotname_2)


plt.figure(figsize=(6,5))
if beadplot==True:
    plt.errorbar(spacings, DRs_v2, yerr=DRs_stdv_v2, capsize=2, fmt='none', label=r'$D_R/D_{bulk,R}$')
    plt.errorbar(spacings, Dzs_v2, yerr=Dzs_stdv_v2, capsize=2, fmt='none', label=r'$D_\perp/D_{bulk,R}$')
    plt.errorbar(spacings, Dparallel_v2, yerr=Dparallel_stdv_v2, capsize=2, fmt='none', label=r'$D_\parallel/D_{bulk,R}$')
else:
    plt.errorbar(spacings, DRs_v2, yerr=DRs_stdv_v2, capsize=2, label=r'$D_R/D_{bulk,R}$')
    plt.errorbar(spacings, Dzs_v2, yerr=Dzs_stdv_v2, capsize=2, label=r'$D_\perp/D_{bulk,R}$')
    plt.errorbar(spacings, Dparallel_v2, yerr=Dparallel_stdv_v2, capsize=2, label=r'$D_\parallel/D_{bulk,R}$')
if moresigmas==False:
    plt.xlabel(r'$d$')
else:
    plt.xlabel(r'$d/\sigma_b$')
plt.ylabel(r'Diffusion constant $D/D_{bulk}$')
if moresigmas==False:
    plt.title('Diffusion constant $D/D_{bulk}$ vs $d$')
else:
    plt.title('Diffusion constant $D/D_{bulk}$ vs $d/\sigma_b$')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='lower right')
#plt.axis([0,10,0,6e-1])
plt.axis([1.1,1.6,0,1e-1])
###box = ax.get_position()
###ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
###ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.tight_layout()
plt.savefig(plotname_small)

 

plt.figure(figsize=(6,5))
plt.errorbar(spacings, DRs_ud, yerr=DRs_stdv_ud, capsize=2, label=r'$D_R$')
plt.errorbar(spacings, Dzs_ud, yerr=Dzs_stdv_ud, capsize=2, label=r'$D_\perp$')
plt.errorbar(spacings, Dparallel_ud, yerr=Dparallel_stdv_ud, capsize=2, label=r'$D_\parallel$')
if moresigmas==False:
    plt.xlabel(r'$d$')
else:
    plt.xlabel(r'$d/\sigma_b$')
plt.ylabel(r'Diffusion constant $D$')
if moresigmas==False:
    plt.title('Diffusion constant $D$ vs $d$')
else:
    plt.title('Diffusion constant $D$ vs $d/\sigma_b$')
plt.tight_layout()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='lower right')
###box = ax.get_position()
###ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
###ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.tight_layout()
plt.savefig(plotname_3)
plt.show()
print('Done')