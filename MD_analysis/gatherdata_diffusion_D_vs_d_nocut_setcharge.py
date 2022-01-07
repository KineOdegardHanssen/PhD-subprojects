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
charge   = -0.25
spacings = [2,2.5,3,3.5,4,4.5,5,6,7,8,10] #np.array([5, 10]) # Will add d=100 to this list soon. Need even more points too.
pmass    = 1 # I don't use this anymore, but it got stuck in the file names. Have constant mass density now.
damp     = 10
N        = len(spacings)
Nrms     = 10 
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
    filestext    = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])

endlocation_out = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/'+parentfolder+'Sigma_bead_' +str(psigma) + '/Nocut/'+ 'Charge' + str(charge)+'/'
outfilename  = endlocation_out+'D_vs_d_charge'+str(charge)+'.txt'
plotname     = endlocation_out+'D_vs_d_charge'+str(charge)+'.png'
plotname_fit = endlocation_out+'D_vs_d_charge'+str(charge)+'_fit.png'
indfilename  = endlocation_out+'D_vs_d_charge'+str(charge)+'_fitindices.txt'

outfile = open(outfilename, 'w')
outfile.write('d   D_R2   sigmaD_R2  b_R2 sigmab_R2; D_z2  sigmaD_z2 b_z2  sigmaD_z2; D_par2 sigmaD_par2  b_par2  sigmab_par2\n')

indexfile = open(indfilename, 'w')
indexfile.write('Start_index_R     end_index_R     Start_index_ort     end_index_ort     Start_index_par     end_index_par\n')

for i in range(N):
    spacing = spacings[i]
    endlocation_in   = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/' % damp +parentfolder+ 'Sigma_bead_' +str(psigma) +'/Charge'+str(charge)+'/Nocut/'
    infilename = endlocation_in+'diffusion'+filestext+'_better_rms_Nestimates%i_charge' % Nrms +str(charge)+'.txt'
    metaname   = endlocation_in+'indices_for_fit_Nestimates%i_charge' % Nrms +str(charge)+'.txt'
    
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
    
    outfile.write('%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n' % (spacing, DRs[i], DRs_stdv[i], Dzs[i], Dzs_stdv[i], Dparallel[i], Dparallel_stdv[i]))
    
    metafile = open(metaname, 'r')
    mlines   = metafile.readlines()
    startindex_R   = int(mlines[0].split()[1])
    endindex_R     = int(mlines[1].split()[1])
    startindex_ort = int(mlines[2].split()[1])
    endindex_ort   = int(mlines[3].split()[1])
    startindex_par = int(mlines[4].split()[1])
    endindex_par   = int(mlines[5].split()[1])
    metafile.close()
    indexfile.write('%i %i %i %i %i %i\n' % (startindex_R, endindex_R, startindex_ort, endindex_ort, startindex_par, endindex_par))

outfile.close()

spacings = np.array(spacings)
    
plt.figure(figsize=(6,5))
#plt.errorbar(psigmas, DRs, yerr=DRs_stdv, capsize=2, fmt='none', label=r'$D_R$')
plt.errorbar(spacings, DRs, yerr=DRs_stdv, capsize=2, label=r'$D_R$')
plt.errorbar(spacings, Dzs, yerr=Dzs_stdv, capsize=2, label=r'$D_\perp$')
plt.errorbar(spacings, Dparallel, yerr=Dparallel_stdv, capsize=2, label=r'$D_\parallel$')
plt.xlabel(r'$d$')
plt.ylabel(r'Diffusion constant $D$')
plt.title('Diffusion constant $D$, d = %i nm, dynamic %s, charge %.2f' % (spacing, systemtype,charge))
plt.tight_layout()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='upper left')
plt.savefig(plotname)

