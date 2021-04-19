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

Nintervals = 10

# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
psigma   = 1 # For instance 
spacings = [1,1.25,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,10,15,25,50,75,100] 
pmass    = 1 # I don't use this anymore, but it got stuck in the file names. Have constant mass density now.
damp     = 10
N        = len(spacings)

# Input booleans for file selection:
bulkdiffusion = False
substrate     = False
confignrs     = np.arange(1,1001)

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

# bs
bRs = np.zeros(N)
bxs = np.zeros(N)
bys = np.zeros(N)
bzs = np.zeros(N)
bparallel = np.zeros(N)
# bs, stdv
bRs_stdv = np.zeros(N)
bxs_stdv = np.zeros(N)
bys_stdv = np.zeros(N)
bzs_stdv = np.zeros(N)
bparallel_stdv = np.zeros(N)


if bulkdiffusion==True:
    parentfolder = 'Pure_bulk/'
    filestext    = '_seed'+str(confignrs[0])+'to'+str(confignrs[-1])
    systemtype   = 'bulk'
    if substrate==True:
        parentfolder = 'Bulk_substrate/'
        systemtype   = 'substrate'
else:
    parentfolder = 'Brush/'
    systemtype   = 'brush'
    filestext    = '_config'+str(confignrs[0])+'to'+str(confignrs[-1])

endlocation_out = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/'+parentfolder+'Sigma_bead_' +str(psigma) + '/Nocut/'
outfilename  = endlocation_out+'D_vs_d_better_rms_Nintervals%i.txt' % Nintervals
plotname     = endlocation_out+'D_vs_d_better_rms_Nintervals%i.png' % Nintervals
plotname_fit = endlocation_out+'D_vs_d_fit_better_rms_Nintervals%i.png' % Nintervals
indfilename  = endlocation_out+'D_vs_d_fitindices_better_rms_Nintervals%i.txt' % Nintervals

outfile = open(outfilename, 'w')
outfile.write('d   D_R2   sigmaD_R2  b_R2 sigmab_R2; D_z2  sigmaD_z2 b_z2  sigmaD_z2; D_par2 sigmaD_par2  b_par2  sigmab_par2\n')

indexfile = open(indfilename, 'w')
indexfile.write('Start_index_R     end_index_R     Start_index_ort     end_index_ort     Start_index_par     end_index_par\n')

for i in range(N):
    spacing = spacings[i]
    endlocation_in   = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/' % damp +parentfolder+ 'Sigma_bead_' +str(psigma) + '/Nocut/'
    infilename = endlocation_in+'diffusion'+filestext+'_nocut_better_rms_Nintervals%i.txt' % Nintervals
    metaname   = endlocation_in+'diffusion_metadata'+filestext+'_nocut_better_rms_Nintervals%i.txt' % Nintervals

    #print('infilename_all:',infilename_all)
    
    # Read in:
    #### Automatic part
    ## Find the extent of the polymers: Max z-coord of beads in the chains
    infile = open(infilename, "r")
    lines = infile.readlines() # This takes some time
    # Getting the number of lines, etc.
    line = lines[1]
    words = line.split()
    
    # Ds
    DRs[i] = float(words[0])
    Dzs[i] = float(words[2])
    Dparallel[i] = float(words[4])
    # Ds, stdv
    DRs_stdv[i] = float(words[1])
    Dzs_stdv[i] = float(words[3])
    Dparallel_stdv[i] = float(words[5])
    
    infile.close()
    
    outfile.write('%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n' % (spacing, DRs[i], DRs_stdv[i],  Dzs[i], Dzs_stdv[i], Dparallel[i], Dparallel_stdv[i]))
    
    metafile   = open(metaname, 'r')
    mlines     = metafile.readlines()
    startindex = int(mlines[0].split()[1])
    endindex   = int(mlines[1].split()[1])
    metafile.close()
    indexfile.write('%.2f %i %i\n' % (spacing, startindex, endindex))

outfile.close()

spacings = np.array(spacings)
    
plt.figure(figsize=(6,5))
#plt.errorbar(psigmas, DRs, yerr=DRs_stdv, capsize=2, fmt='none', label=r'$D_R$')
plt.errorbar(spacings, DRs, yerr=DRs_stdv, capsize=2, label=r'$D_R$')
plt.errorbar(spacings, Dzs, yerr=Dzs_stdv, capsize=2, label=r'$D_\perp$')
plt.errorbar(spacings, Dparallel, yerr=Dparallel_stdv, capsize=2, label=r'$D_\parallel$')
plt.xlabel(r'$d$')
plt.ylabel(r'Diffusion constant $D$')
plt.title('Diffusion constant $D$, dynamic %s' % systemtype)
plt.tight_layout()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='upper left')
plt.savefig(plotname)
plt.show() # For testing purposes

