from skimage import morphology, measure, io, util   # For performing morphology operations (should maybe import more from skimage...)
from mpl_toolkits.mplot3d import Axes3D             # Plotting in 3D
import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
from scipy.ndimage import measurements, convolve    # Should have used the command below, but changing now will lead to confusion
from scipy import ndimage                           # For Euclidean distance measurement
import numpy as np
import random
import math
import time
import os
import glob
import copy

# Treats all the _all-files regardless of whether the diffusing bead displays an unphysical trajectory.

#
damp = 10
M = 9
Ninch = 101
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
popup_plots = False
long    = False
spacing = 50
psigma  = 1
density = 0.238732414637843 # Yields mass 1 for bead of radius 1 nm
print('spacing:', spacing)
print('psigma:', psigma)
# I need to set the file name in an easier way, but for now I just use this:  ## Might want to add a loop too, if I have more files...

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T        = 3
plotseed = 0
plotdirs = False
zhigh        = 290    # For omitting unphysical trajectories
zlow         = -50    
confignrs    = np.arange(1,1001)
Nseeds       = len(confignrs)      # So that I don't have to change that much
maxz_av      = 0
filescounter = 0
if long==True:
    Nsteps   = 10001
else:
    Nsteps   = 2001 # 20001 before
#writeevery   = 10                  # I'm writing to file every this many time steps
unitlength   = 1e-9
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime#*writeevery
Npartitions  = 5 # For extracting more walks from one file (but is it really such a random walk here...?)
minlength    = int(floor(Nsteps/Npartitions)) # For sectioning the data
print('timestepsize:', timestepsize)

endlocation          = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp + str(psigma) + '/'
filestext            = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])
# Text files
outfilename_reez_frac_all = endlocation+filestext+'_reez_frac_all'
outfilename_reez_frac_avg = endlocation+filestext+'_reez_frac_avg_rms'
outfilename_reez_base_avg = endlocation+filestext+'_reez_base_avg_rms'
outfilename_reepar_avg = endlocation+filestext+'_reepar_avg_rms'

if long==True:
    outfilename_reez_frac_all = outfilename_reez_frac_all+'_long.txt'
    outfilename_reez_frac_avg = outfilename_reez_frac_avg +'_long.txt'
    outfilename_reez_base_avg = outfilename_reez_base_avg +'_long.txt'
    outfilename_reepar_avg = outfilename_reepar_avg +'_long.txt'
else:
    outfilename_reez_frac_all = outfilename_reez_frac_all+'.txt'
    outfilename_reez_frac_avg = outfilename_reez_frac_avg +'.txt'
    outfilename_reez_base_avg = outfilename_reez_base_avg +'.txt'
    outfilename_reepar_avg = outfilename_reepar_avg +'.txt'

outfile_reez_frac_all = open(outfilename_reez_frac_all,'w')

reezs_frac_all = []
reezs_frac_avg = 0
reezs_frac_rms = 0
reezs_base_all = []
reezs_base_avg = 0
reezs_base_rms = 0
reepars_all = []
reepar_avg = 0
reepar_rms = 0

Nread        = 0
skippedfiles = 0

for confignr in confignrs:
    print('On config number:', confignr)
    if long==True:
        infilename_all  = endlocation+'long/'+'all_confignr'+str(confignr)+'_long.lammpstrj'
    else:
        infilename_all  = endlocation+'all_confignr'+str(confignr)+'.lammpstrj'

    
    #print('infilename_all:',infilename_all)
    
    # Read in:
    #### Automatic part
    ## Find the extent of the polymers: Max z-coord of beads in the chains
    try:
        infile_all = open(infilename_all, "r")
    except:
        try:
            infilename_all = endlocation+'long/'+'all_confignr'+str(confignr)+'.lammpstrj'
            infile_all = open(infilename_all, "r")
        except:
            print('Oh, lammpstrj-file! Where art thou?')
            skippedfiles += 1
            continue # Skipping this file if it does not exist
    # Moving on, if the file
    filescounter += 1
    lines = infile_all.readlines() # This takes some time
    # Getting the number of lines, etc.
    totlines = len(lines)         # Total number of lines
    lineend = totlines-1          # Index of last element
    
    # Extracting the number of atoms:
    words = lines[3].split()
    Nall = int(words[0])
    N    = Nall

    chain_pos = np.zeros((M,Ninch,3)) # All atoms in chains. One time step because the first is the same for all configs
    
    skiplines   = 9             # If we hit 'ITEM:', skip this many steps...
    skipelem    = 0
    sampleevery = 0
    #i           = int(math.ceil(skipelem*(Nall+9)))
    skiplines  += (Nall+skiplines)*sampleevery # Check!
    i           = 0
    t_index     = -1
    
    time_start = time.process_time()
    counter = 0
    while i<totlines:
        words = lines[i].split()
        if words[0]=='ITEM:':
            if words[1]=='TIMESTEP':
                t_index+=1
                if t_index==2:
                    break
                i_array     = np.zeros(9)
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
            atomtype = int(words[1]) 
            if atomtype==1 or atomtype==2: # Chain bead. Treat
                molID    = int(words[2])-1
                x        = float(words[3])
                y        = float(words[4])
                z        = float(words[5])
                # Overwriting the first time step. Not that efficient, but won't use this script much
                chain_pos[molID,int(i_array[molID]),0] = x 
                chain_pos[molID,int(i_array[molID]),1] = y
                chain_pos[molID,int(i_array[molID]),2] = z
                i_array[molID] += 1
            i+=1
    infile_all.close()
    
    # Find radius of gyration:
    for ic in range(M): # Looping over chains
        # Finding average r:
        reex = chain_pos[ic,-1,0]-chain_pos[ic,0,0]
        reey = chain_pos[ic,-1,1]-chain_pos[ic,0,1]
        reez = chain_pos[ic,-1,2]-chain_pos[ic,0,2]
        reepar = np.sqrt(reex*reex+reey*reey)
        reefrac = reez/reepar
        reezs_frac_all.append(reefrac)
        reepars_all.append(reefrac)
        reezs_frac_avg += reez
        reepar_avg += reepar
        reezs_base_all.append(reez)
        reezs_base_avg += reefrac
        print('reezs_frac_avg (accumulating):',reezs_frac_avg)
        outfile_reez_frac_all.write('%.5f ' % reefrac)
    outfile_reez_frac_all.write('\n') # New line for each time step

'''
reezs_base_all = []
reezs_base_avg = 0
reezs_base_rms = 0
'''
            
outfile_reez_frac_all.close()
print('reez_frac_all:',reezs_frac_all)

## Base
NRz = len(reezs_base_all) # Should be filescounter*9, but this works too
reezs_base_avg/=NRz
reezs_base_rms=0
for i in range(NRz):
    reezs_base_rms += (reezs_base_all[i]-reezs_base_avg)**2
reezs_base_rms = np.sqrt(reezs_base_rms/(NRz-1))

outfile_reez_base_avg = open(outfilename_reez_base_avg,'w')
outfile_reez_base_avg.write('%.16e %.16e' % (reezs_base_avg,reezs_base_rms))
outfile_reez_base_avg.close()

## par
NRz = len(reepars_all) # Should be filescounter*9, but this works too
reepar_avg/=NRz
reepar_rms=0
for i in range(NRz):
    reepar_rms += (reepars_all[i]-reezs_base_avg)**2
reepar_rms = np.sqrt(reepar_rms/(NRz-1))

outfile_reepar_avg = open(outfilename_reepar_avg,'w')
outfile_reepar_avg.write('%.16e %.16e' % (reepar_avg,reepar_rms))
outfile_reepar_avg.close()

## Frac
NRz = len(reezs_frac_all) # Should be filescounter*9, but this works too
reezs_frac_avg/=NRz
reezs_frac_rms=0
for i in range(NRz):
    reezs_frac_rms += (reezs_frac_all[i]-reezs_frac_avg)**2
reezs_frac_rms = np.sqrt(reezs_frac_rms/(NRz-1))

outfile_reez_frac_avg = open(outfilename_reez_frac_avg,'w')
outfile_reez_frac_avg.write('%.16e %.16e' % (reezs_frac_avg,reezs_frac_rms))
outfile_reez_frac_avg.close()


print('d:',spacing)
print('Reez_frac:', reezs_frac_avg, ' +/-', reezs_frac_rms)
print('mean(reezs_frac_all):',np.mean(reezs_frac_all))
print('min(reezs_frac_all):',min(reezs_frac_all))
print('max(reezs_frac_all):',max(reezs_frac_all))
print('Reez_base:', reezs_base_avg, ' +/-', reezs_base_rms)
print('mean(reezs_base_all):',np.mean(reezs_base_all))
print('min(reezs_base_all):',min(reezs_base_all))
print('max(reezs_base_all):',max(reezs_base_all))


print('NRz:', NRz, '; filescounter*9:',filescounter*9)