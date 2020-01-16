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
spacings = np.array([5, 10]) # Will add d=100 to this list soon. Need even more points too.
pmass    = 1.5 # I don't use this anymore, but it got stuck in the file names. Have constant mass density now.
N        = len(spacings)
# Input booleans for file selection:
bulkdiffusion = False
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
    systemtype   = 'bulk'
    startpart    = '_'
    parentfolder = 'Pure_bulk/'
    if substrate==True:
        systemtype   = 'substrate'
        startpart    = '_withsubstrate_'
        parentfolder = 'Bulk_substrate/'
    namebase_out   = 'bulkdiffusion'+startpart+'ljunits_Langevin_scaled_T3_pmass'
    name_end       = '_sect_placeexact_ljepsilon0.730372054992096_ljcut1p122'
    folderbase     = 'Bulkdiffusion'+startpart+'ljunits_Langevin_scaled_T3_pmass1.5_sect_placeexact_ljepsilon0.730372054992096_ljcut1p122'
    namebase_short = namebase # In case I ever need to shorten it
    
    name_start = 'lammpsdiffusion_'
    
else:
    systemtype = 'brush'
    parentfolder = 'Brush/'
    #namebase = '_quadr_M9N101_ljunits_spacing%i_Langevin_scaled_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_pmass'% spacing + str(pmass)
    namebase_out = '_quadr_M9N101_ljunits_Langevin_scaled_Kangle14_Kbond140_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_pmass' + str(pmass)
    folderbase   = 'Quadr_M9N101_ljunits_Langevin_scaled_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_ljcut1p122'    
    name_start   = 'lammpsdiffusion_qdrgr'
    name_end     = '_sect_placeexact_ljcut1p122'

endlocation_gathered = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/'+parentfolder+folderbase+'/'
outfilename  = endlocation_gathered+'D_vs_d_'+name_start+namebase_out+name_end+'.txt'
plotname     = endlocation_gathered+'D_vs_d_'+name_start+namebase_out+name_end+'.png'
plotname_fit = endlocation_gathered+'D_vs_d_'+name_start+namebase_out+name_end+'_fit.png'
indfilename  = endlocation_gathered+'D_vs_d_'+name_start+namebase_out+name_end+'_fitindices.txt'

outfile = open(outfilename, 'w')
outfile.write('psigma   D_R2   sigmaD_R2  b_R2 sigmab_R2; D_z2  sigmaD_z2 b_z2  sigmaD_z2; D_par2 sigmaD_par2  b_par2  sigmab_par2\n')

indexfile = open(indfilename, 'w')
indexfile.write('Start_index_R     end_index_R     Start_index_ort     end_index_ort     Start_index_par     end_index_par\n')

for i in range(N):
    spacing = spacings[i]
    if bulkdiffusion==True:
        namebase = namebase = 'bulkdiffusion'+startpart+'ljunits_spacing%i_Langevin_scaled_T3_pmass' % spacing
    else:
        namebase = '_quadr_M9N101_ljunits_spacing%i_Langevin_scaled_Kangle14_Kbond140_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_pmass'% spacing+ str(pmass)
    endlocation   = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/'+parentfolder+folderbase+'/Spacing'+str(spacing)+'/Sigma_bead_'+str(psigma)+'/'
    namebase_start = '_quadr_M9N101_ljunits_'
    folderbase_mid = 'Langevin_scaled_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_pmass1.5'
    namebase_short = namebase+'_psigma'+str(psigma)+name_end
    infilename = endlocation+name_start+namebase_short+'_diffusion.txt' #_timestep'+str(startindex)+'to'+str(endindex)+'.txt'
    metaname   = endlocation+name_start+namebase_short+'_metadata.txt' # In original file set as:  endlocation+'lammpsdiffusion_qdrgr'+namebase_short+'_metadata.txt'

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

    
plt.figure(figsize=(6,5))
#plt.errorbar(psigmas, DRs, yerr=DRs_stdv, capsize=2, fmt='none', label=r'$D_R$')
plt.errorbar(spacings, DRs, yerr=DRs_stdv, capsize=2, label=r'$D_R$')
plt.errorbar(spacings, Dzs, yerr=Dzs_stdv, capsize=2, label=r'$D_\perp$')
plt.errorbar(spacings, Dparallel, yerr=Dparallel_stdv, capsize=2, label=r'$D_\parallel$')
plt.xlabel(r'$\sigma_p$')
plt.ylabel(r'Diffusion constant $D$')
plt.title('Diffusion constant $D$, d = %i nm, system %s' % (spacing, systemtype))
plt.tight_layout()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='upper left')
plt.savefig(plotname)

