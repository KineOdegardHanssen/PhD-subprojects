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
spacing = 5

# I need to set the file name in an easier way, but for now I just use this:  ## Might want to add a loop too, if I have more files...

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them
plotseed = 0
plotdirs = True
seeds  = ['23', '29', '31', '37', '41', '43', '47', '53', '59', '61', '67', '71', '73', '79', '83', '89', '97', '101', '103', '107', '109', '113']
Nseeds = len(seeds)
Nsteps = 80001
namebase        = '_quadr_M9N101_ljunits_spacing5_Langevin_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_theta0is180_pmass1.5_sect_placeexact_ljcut1p122'
folderbase      = 'Part_in_chgr_subst_all_quadr_M9N101_ljunits_Langevin_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_theta0is180_pmass1.5_sect_placeexact_ljcut1p122'
endlocation     = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/'+folderbase+'/Spacing'+str(spacing)+'/'
outfilename     = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'.txt'
plotname        = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'.png'
plotname_dirs   = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'_dxdydzR2_seed'+str(seeds[plotseed])+'.png'
# Setting arrays
# These are all squared:
# All together:
allRs     = []
alldxs    = []
alldys    = []
alldzs    = []
# Separated by seed:
Rs_byseed    = []
dxs_byseed   = []
dys_byseed   = []
dzs_byseed   = []
times_byseed = []
# This is not squared, obviously:
alltimes  = []
   
for seed in seeds:
    print('On seed', seed)
    infilename_all  = 'part_in_chgr_subst_all'+namebase+'_seed'+seed+'.lammpstrj'
    infilename_free = 'part_in_chgr_subst_freeatom'+namebase+'_seed'+seed+'.lammpstrj'
    
    # Read in:
    #### Automatic part
    ## Find the extent of the polymers: Max z-coord of beads in the chains
    infile_all = open(infilename_all, "r")
    lines = infile_all.readlines() # This takes some time
    # Getting the number of lines, etc.
    totlines = len(lines)         # Total number of lines
    lineend = totlines-1          # Index of last element
    
    # Extracting the number of atoms:
    words = lines[3].split()
    Nall = int(words[0])
    N    = Nall
    
    skiplines   = 9             # If we hit 'ITEM:', skip this many steps...
    skipelem    = 0
    sampleevery = 0
    i           = int(math.ceil(skipelem*(Nall+9)))
    skiplines  += (Nall+skiplines)*sampleevery # Check!
    
    # For now: Only find the largest z-position among the beads in the chain. # Otherwise I must read off the atom number (will do that with the free bead anyways.)
    maxz = -1000 # No z-position is this small.
    
    time_start = time.process_time()
    counter = 0
    while i<totlines:
        words = lines[i].split()
        if (words[0]=='ITEM:' and words[1]=='TIMESTEP'): # Some double testing going on...
            if words[1]=='TIMESTEP':
                i+=skiplines
            elif words[1]=='NUMBER': # These will never kick in. 
                i+=7
            elif words[1]=='BOX':
                i+=5
            elif words[1]=='ATOMS':
                i+=1
        elif len(words)<8:
            i+=1
        else:
            # Find properties
            # Order:  id  type mol ux  uy  uz  vx  vy   vz
            #         [0] [1]  [2] [3] [4] [5] [6] [7]  [8]
            ind      = int(words[0])-1 # Atom ids go from zero to N-1.
            atomtype = int(words[1]) 
            #molID    = int(words[2])
            z        = float(words[5])
            if atomtype==2: # Moving polymer bead. Test if this is larger than maxz:
                if z>maxz:
                    maxz = z
            i+=1
    extent_polymers = maxz
    infile_all.close()
    
    ## Find the position of the free bead: # I reuse quite a bit of code here...
    infile_free = open(infilename_free, "r")
    lines = infile_free.readlines() # This takes some time
    # Getting the number of lines, etc.
    totlines = len(lines)         # Total number of lines
    lineend = totlines-1          # Index of last element
    
    # Extracting the number of atoms:
    words = lines[3].split()
    Nall = int(words[0])
    N    = Nall
    
    skiplines   = 9             # If we hit 'ITEM:', skip this many steps...
    skipelem    = 0
    sampleevery = 0
    i           = int(math.ceil(skipelem*(Nall+9)))
    skiplines  += (Nall+skiplines)*sampleevery # Check!
    
    # Setting arrays for treatment:
    freeatom_positions = [] # Fill with xu,yu,zu-coordinates. One position for each time step. # We have the same time step for the whole of the simulation. Can bother with time later
    times              = [] # Not useful anymore?
    
    # For now: Only find the largest z-position among the beads in the chain. # Otherwise I must read off the atom number (will do that with the free bead anyways.)
    maxz = -1000 # No z-position is this small.
    
    time_start = time.process_time()
    counter = 0
    while i<totlines:
        words = lines[i].split()
        if (words[0]=='ITEM:' and words[1]=='TIMESTEP'): # Some double testing going on...
            if words[1]=='TIMESTEP':
                words2 = lines[i+1].split() # The time step is on the next line
                t = float(words2[0])
                times.append(t)
                i+=skiplines
            elif words[1]=='NUMBER': # These will never kick in. 
                i+=7
            elif words[1]=='BOX':
                i+=5
            elif words[1]=='ATOMS':
                i+=1
        elif len(words)<8:
            i+=1
        else:
            # Find properties
            # Order:  id  type mol ux  uy  uz  vx  vy   vz
            #         [0] [1]  [2] [3] [4] [5] [6] [7]  [8]
            ind      = int(words[0])-1 # Atom ids go from zero to N-1.
            #atomtype = int(words[1]) 
            #molID    = int(words[2])
            x        = float(words[3])
            y        = float(words[4])
            z        = float(words[5])
            freeatom_positions.append(np.array([x,y,z]))
            counter+=1
            i+=1
    infile_free.close()
    dt = times[1]-times[0] # This might be handy
    
    time_end = time.process_time()
    
    pos_inpolymer = []
    
    #before_in
    
    # I do not take into account that the free bead can enter the polymer grid and then exit again. If that happens, there might be discontinuities or weird kinks in the data (if you look at the graphs...) # I do not have this in my test-dataset. Maybe make the bead lighter and see what happens?
    
    Nin   = 0
    maxzpol = 0
    for i in range(counter):
        thesepos = freeatom_positions[i]
        z        = thesepos[2]
        if z>extent_polymers: # If the polymer is in bulk # We don't want it to go back and forth between brush and bulk # That will cause discontinuities in our data
            break
        else:
            pos_inpolymer.append(thesepos)
            if z>maxzpol:
                maxzpol = z
            Nin+=1
    
    startpos_in   = pos_inpolymer[0]
        
    #######
    # Will divide into several walks with different starting points later on
    # Should store for RMS too? ... Then I need to have more starting points or more time frames.
    
    # Finding R2s and corresponding times
    # Maybe... Not plot the first ones since it takes some time to equilibrate.
    # But how much do I need to cut? Varies from plot to plot how much makes sense to cut (if it is even possible to say...)
    # Go for some percentage of all points?
    R_temp  = []
    dx_temp = []
    dy_temp = []
    dz_temp = []
    step_temp = []
    
    R_temp.append(0)
    dx_temp.append(0)
    dy_temp.append(0)
    dz_temp.append(0)
    step_temp.append(0)
    allRs.append(0)        # We will set this here since we know the value
    alltimes.append(0)
    for i in range(1,Nin):       
        this_in = pos_inpolymer[i]
        dist = this_in-startpos_in
        R2   = np.dot(dist,dist)
        dx2  = dist[0]*dist[0]
        dy2  = dist[1]*dist[1]
        dz2  = dist[2]*dist[2]
        # All together:
        allRs.append(R2)
        alldxs.append(dx2)
        alldys.append(dy2)
        alldzs.append(dz2)
        alltimes.append(i)
        # Separated by seed:
        R_temp.append(R2)
        dx_temp.append(dx2)
        dy_temp.append(dy2)
        dz_temp.append(dz2)
        step_temp.append(i)
    Rs_byseed.append(R_temp)
    dxs_byseed.append(dx_temp)
    dys_byseed.append(dy_temp)
    dzs_byseed.append(dz_temp)
    times_byseed.append(step_temp)

allRs     = np.array(allRs)
alltimes  = np.array(alltimes)

## Finding the diffusion coefficient (in polymer)
# Finding the last 25% of the data:               # Do this later! In loop! # But do now just for testing purposes.
# Performing the line fit:
coeffs_poly = polyfit(alltimes, allRs, 1)
a_poly = coeffs_poly[0]
b_poly = coeffs_poly[1]
D_poly = a_poly/6.
    
fit_poly = a_poly*alltimes+b_poly

print('D_poly:', D_poly)

outfile = open(outfilename, 'w')
outfile.write('D_poly: %.16f\n' % D_poly)
    
plt.figure(figsize=(6,5))
plt.plot(alltimes, allRs, ',', label='Data, brush')
plt.plot(alltimes, fit_poly, '--', label='Fit, brush')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title('Random walk in the vicinity of the polymer chains')
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.plot(times_byseed[plotseed], Rs_byseed[plotseed], label=r'$<R^2>$')
plt.plot(times_byseed[plotseed], dxs_byseed[plotseed], label=r'$<dx^2>$')
plt.plot(times_byseed[plotseed], dys_byseed[plotseed], label=r'$<dy^2>$')
plt.plot(times_byseed[plotseed], dzs_byseed[plotseed], label=r'$<dz^2>$')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title('Random walk in the vicinity of the polymer chains, seed %s' % seeds[plotseed])
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_dirs)

