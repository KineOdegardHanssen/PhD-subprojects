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

# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
bulkdiffusion = False
substrate     = False

spacing = 5
psigma  = 1.5

# Choosing which part to fit
startindex = 5000
endindex   = 10000

# I need to set the file name in an easier way, but for now I just use this:  ## Might want to add a loop too, if I have more files...

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
plotseed = 0
plotdirs = False
test_sectioned = False
seeds  = [23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]
Nseeds = len(seeds)
Nsteps = 80001
unitlength   = 1e-9
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime
print('timestepsize:', timestepsize)
Npartitions = 5 # For extracting more walks from one file (but is it really such a random walk here...?)

# Make a loop?

# Make file names
## Weird cutoff (bead):
#namebase    = '_quadr_M9N101_ljunits_spacing%i_Langevin_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_theta0is180_pmass1.5_sect_placeexact_ljcut1p122' %spacing
#folderbase  = 'Part_in_chgr_subst_all_quadr_M9N101_ljunits_Langevin_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_theta0is180_pmass1.5_sect_placeexact_ljcut1p122'
# Usual cutoff (bead):
if bulkdiffusion==True:
    startpart = '_'
    parentfolder = 'Pure_bulk/'
    if substrate==True:
        startpart = '_withsubstrate_'
        parentfolder = 'Bulk_substrate/'
    namebase = 'bulkdiffusion'+startpart+'ljunits_spacing%i_Langevin_scaled_T3_pmass1.5_sect_placeexact_ljepsilon0.730372054992096_ljcut1p122' % spacing
    folderbase = 'Bulkdiffusion'+startpart+'ljunits_Langevin_scaled_T3_pmass1.5_sect_placeexact_ljepsilon0.730372054992096_ljcut1p122'
    namebase_short = namebase # In case I ever need to shorten it
    
    endlocation   = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/'+parentfolder+folderbase+'/Spacing'+str(spacing)+'/Sigma_bead_'+str(psigma)+'/'
    infilename    = endlocation+'lammpsdiffusion_'+namebase_short+'_av_ds.txt'
    outfilename   = endlocation+'lammpsdiffusion_'+namebase_short+'_diffusion.txt'
    metaname      = endlocation+'lammpsdiffusion_'+namebase_short+'_diffusion_metadata.txt'
    plotname_R    = endlocation+'lammpsdiffusion_'+namebase_short+'_diffusion_R.png'
    plotname_dz   = endlocation+'lammpsdiffusion_'+namebase_short+'_diffusion_dz.png'
    plotname_dpar = endlocation+'lammpsdiffusion_'+namebase_short+'_diffusion_dpar.png'
else:
    namebase = '_quadr_M9N101_ljunits_spacing%i_Langevin_scaled_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_pmass1.5_psigma' % spacing +str(psigma)+'_sect_placeexact_ljcut1p122'
    namebase_short = '_quadr_M9N101_ljunits_spacing%i_Langevin_scaled_Kangle14_Kbond140_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_pmass1.5_psigma' % spacing +str(psigma)+'_sect_placeexact_ljcut1p122'
    folderbase  = 'Quadr_M9N101_ljunits_Langevin_scaled_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_ljcut1p122'
    
    endlocation   = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Brush/'+folderbase+'/Spacing'+str(spacing)+'/Sigma_bead_'+str(psigma)+'/'
    infilename    = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'_av_ds.txt'
    outfilename   = endlocation+'lammpsdiffusion_qdrgr'+namebase_short+'_diffusion.txt'
    metaname      = endlocation+'lammpsdiffusion_qdrgr'+namebase_short+'_metadata.txt'
    plotname_R    = endlocation+'lammpsdiffusion_qdrgr'+namebase_short+'_diffusion_R.png'
    plotname_dz   = endlocation+'lammpsdiffusion_qdrgr'+namebase_short+'_diffusion_dz.png'
    plotname_dpar = endlocation+'lammpsdiffusion_qdrgr'+namebase_short+'_diffusion_dpar.png'

# Set up arrays
Nsteps = endindex-startindex 
times  = np.zeros(Nsteps) 
dR2s   = np.zeros(Nsteps) 
dz2s   = np.zeros(Nsteps) 
dpar2s = np.zeros(Nsteps)

# Read file
infile = open(infilename,'r')
lines = infile.readlines()

#    	0			1		2		3		4			5		6
# (times_single[i], times_single_real[i], averageRs_SI[i], averagedxs_SI[i], averagedys_SI[i], averagedzs_SI[i], averagedparallel_SI[i]))
for i in range(startindex,endindex):
    j     = i-startindex
    line  = lines[i]
    words = line.split()
    times[j]  = float(words[1])
    dR2s[j]   = float(words[2])
    dz2s[j]   = float(words[5])
    dpar2s[j] = float(words[6])

infile.close()


# line fit, all
coeffs, covs = polyfit(times, dR2s, 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a_R2 = coeffs[0]
b_R2 = coeffs[1]
D_R2 = a_R2/6.
rms_D_R2 = np.sqrt(covs[0,0])/6.
rms_b_R2 = np.sqrt(covs[1,1])
 
fit_R2 = a_R2*times+b_R2

# line fit, dz
coeffs, covs = polyfit(times, dz2s, 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a_z2 = coeffs[0]
b_z2 = coeffs[1]
D_z2 = a_z2/6.
rms_D_z2 = np.sqrt(covs[0,0])/6.
rms_b_z2 = np.sqrt(covs[1,1])

fit_z2 = a_z2*times+b_z2

# line fit, parallel
coeffs, covs = polyfit(times, dpar2s, 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a_par2     = coeffs[0]
b_par2     = coeffs[1]
D_par2     = a_par2/6.
rms_D_par2 = np.sqrt(covs[0,0])/6.
rms_b_par2 = np.sqrt(covs[1,1])

fit_par2 = a_par2*times+b_par2

outfile = open(outfilename,'w')
outfile.write('D_R2  sigmaD_R2 b_R2 sigmab_R2; D_z2  sigmaD_z2  b_z2  sigmaD_z2; D_par2 sigmaD_par2  b_par2  sigmab_par2\n')
outfile.write('%.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n' % (D_R2, rms_D_R2, b_R2, rms_b_par2, D_z2, rms_D_z2, b_z2, rms_b_par2, D_par2, rms_D_par2, b_par2, rms_b_par2))
outfile.close()

metafile = open(metaname, 'w')
metafile.write('startindex: %i\n' % startindex)
metafile.write('endindex: %i\n' % endindex)
metafile.close()

plt.figure(figsize=(6,5))
plt.plot(times, dR2s, ',', label='Average, brush')
plt.plot(times, fit_R2, '--', label='Fit, average, brush')
plt.xlabel(r'Time (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title('RMSD in bulk, d = %i nm, SI' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_R)

plt.figure(figsize=(6,5))
plt.plot(times, dz2s, ',', label='Average, brush')
plt.plot(times, fit_z2, '--', label='Fit, average, brush')
plt.xlabel(r'Time (s)')
plt.ylabel(r'$dz^2$ [in unit length]')
plt.title('$dz^2$ in bulk, d = %i nm, SI' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_dz)

plt.figure(figsize=(6,5))
plt.plot(times, dpar2s, ',', label='Average, brush')
plt.plot(times, fit_par2, '--', label='Fit, average, brush')
plt.xlabel(r'Time (s)')
plt.ylabel(r'$dx^2+dy^2$ [in unit length]')
plt.title('$dx^2+dy^2$ in bulk, d = %i nm, SI' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_dpar)


