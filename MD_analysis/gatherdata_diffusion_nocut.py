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

### No changes made here yet.

# Function for curve_fit
def exponential_complete(s,A,P):
    return A*np.exp(-s/P)


# Don't think I use this anymore
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
damp    = 10
spacing = 10
if spacing==10:
    psigmas = [0.25, 0.5, 1, 1.5, 2, 3, 5, 10, 25]
elif spacing==5:
    psigmas = [0.25, 0.5, 1, 1.5, 2]
pmass   = 1.5
N       = len(psigmas)
# Input booleans for file selection:
bulkdiffusion = True
substrate     = False

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
    systemtype = 'bulk'
    startpart = '_'
    parentfolder = 'Pure_bulk/'
    if substrate==True:
        systemtype = 'substrate'
        parentfolder = 'Bulk_substrate/'
    namebase = 'bulkdiffusion'+startpart+'ljunits_spacing%i_Langevin_scaled_T3_pmass' % spacing
else:
    systemtype = 'brush'
    parentfolder = 'Brush/'

endlocation   = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing%i/damp%i_diffseedLgv/' % (spacing,damp) +parentfolder

outfilename  = endlocation+'D_vs_psigma.txt'
plotname     = endlocation+'D_vs_psigma.png'
plotname_log = endlocation+'D_vs_psigma_log.png'
plotname_fit = endlocation+'D_vs_psigma_fit.png'
indfilename  = endlocation+'D_vs_psigma_fitindices.txt'

outfile = open(outfilename, 'w')
outfile.write('psigma   D_R2   sigmaD_R2  b_R2 sigmab_R2; D_z2  sigmaD_z2 b_z2  sigmaD_z2; D_par2 sigmaD_par2  b_par2  sigmab_par2\n')

indexfile = open(indfilename, 'w')
indexfile.write('Start_index_R     end_index_R     Start_index_ort     end_index_ort     Start_index_par     end_index_par\n')

for i in range(N):
    psigma = psigmas[i]
    endlocation_in = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing%i/damp%i_diffseedLgv/' % (spacing,damp) +parentfolder+ 'Sigma_bead_' +str(psigma) + '/Nocut/'
    infilename = endlocation_in+'diffusion_seed1to1000.txt' #_timestep'+str(startindex)+'to'+str(endindex)+'.txt'
    metaname   = endlocation_in+'diffusion_metadata_seed1to1000.txt' # In original file set as:  endlocation+'lammpsdiffusion_qdrgr'+namebase_short+'_metadata.txt'

    #print('infilename_all:',infilename_all)
    
    # 0		1	2	3	4	5	6	7	8	9		10	11
    #D_R2   sigmaD_R2  b_R2 sigmab_R2; D_z2  sigmaD_z2 b_z2  sigmaD_z2; D_par2 sigmaD_par2  b_par2  sigmab_par2
    #-1.99660e-08 1.52954e-09 3.20458e-15 1.37937e-18 -1.29622e-07 1.49293e-09 3.25575e-15 1.37937e-18 1.09656e-07 2.81101e-10 -5.11720e-17 1.37937e-18
    
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
    Dzs[i] = float(words[4])
    Dparallel[i] = float(words[8])
    # Ds, stdv
    DRs_stdv[i] = float(words[1])
    Dzs_stdv[i] = float(words[5])
    Dparallel_stdv[i] = float(words[9])
    
    # bs
    bRs[i] = float(words[2])
    bzs[i] = float(words[6])
    bparallel[i] = float(words[10])
    
    # bs, stdv
    bRs_stdv[i] = float(words[3])
    bzs_stdv[i] = float(words[7])
    bparallel_stdv[i] = float(words[11])
    
    infile.close()
    
    outfile.write('%.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e \n' % (psigma, DRs[i], DRs_stdv[i], bRs[i], bRs_stdv[i], Dzs[i], Dzs_stdv[i], bzs[i], bzs_stdv[i], Dparallel[i], Dparallel_stdv[i], bparallel[i], bparallel[i]))
    
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

psigmas = np.array(psigmas)

poptR, pcovR  = curve_fit(exponential_complete, psigmas, DRs)#, maxfev=thismaxfev) #?
AR       = poptR[0]
expR     = poptR[1]

fitR = exponential_complete(psigmas, AR, expR)

poptz, pcovz  = curve_fit(exponential_complete, psigmas, Dzs)#, maxfev=thismaxfev) #?
Az       = poptz[0]
expz     = poptz[1]

fitz = exponential_complete(psigmas, Az, expz)

poptpar, pcovpar  = curve_fit(exponential_complete, psigmas, Dparallel)#, maxfev=thismaxfev) #?
Apar       = poptpar[0]
exppar     = poptpar[1]

fitpar = exponential_complete(psigmas, Apar, exppar)
    
plt.figure(figsize=(6,5))
#plt.errorbar(psigmas, DRs, yerr=DRs_stdv, capsize=2, fmt='none', label=r'$D_R$')
plt.errorbar(psigmas, DRs, yerr=DRs_stdv, capsize=2, label=r'$D_R$')
plt.errorbar(psigmas, Dzs, yerr=Dzs_stdv, capsize=2, label=r'$D_\perp$')
plt.errorbar(psigmas, Dparallel, yerr=Dparallel_stdv, capsize=2, label=r'$D_\parallel$')
plt.xlabel(r'$\sigma_p$')
plt.ylabel(r'Diffusion constant $D$')
plt.title('Diffusion constant $D$, d = %i nm, system %s' % (spacing, systemtype))
plt.tight_layout()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='upper right')
plt.savefig(plotname)
#plt.show()


plt.figure(figsize=(6,5))
plt.errorbar(psigmas, DRs, yerr=DRs_stdv, capsize=2, label=r'$D_R$')
plt.plot(psigmas, fitR, '--', label=r'$D_R$, fit')
plt.errorbar(psigmas, Dzs, yerr=Dzs_stdv, capsize=2, label=r'$D_\perp$')
plt.plot(psigmas, fitz, '--', label=r'$D_\perp$, fit')
plt.errorbar(psigmas, Dparallel, yerr=Dparallel_stdv, capsize=2, label=r'$D_\parallel$')
plt.plot(psigmas, fitpar, '--', label=r'$D_\parallel$, fit')
plt.xlabel(r'$\sigma_p$')
plt.ylabel(r'Diffusion constant $D$')
plt.title('Diffusion constant $D$, d = %i nm, system %s' % (spacing, systemtype))
plt.tight_layout()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='upper right')
plt.savefig(plotname_fit)


plt.figure(figsize=(6,5))
#plt.errorbar(psigmas, DRs, yerr=DRs_stdv, capsize=2, fmt='none', label=r'$D_R$')
plt.loglog(psigmas, DRs, '-o', label=r'$D_R$')
plt.loglog(psigmas, Dzs, '-o', label=r'$D_\perp$')
plt.loglog(psigmas, Dparallel, '-o', label=r'$D_\parallel$')
plt.xlabel(r'$\sigma_p$')
plt.ylabel(r'Diffusion constant $D$')
plt.title('Diffusion constant $D$, d = %i nm, system %s' % (spacing, systemtype))
plt.tight_layout()
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='upper right')
plt.savefig(plotname_log)


