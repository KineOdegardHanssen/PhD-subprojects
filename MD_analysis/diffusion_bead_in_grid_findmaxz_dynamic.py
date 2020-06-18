import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob
import copy

# For a given spacing, find the maximal z-value of the beads or the walker.

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
#
damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
popup_plots = False
spacings = [1,1.25,1.5,2,3,4,5,6,7,10,15,25,50,75,100]
psigma  = 1#.5
density = 0.238732414637843 # Yields mass 1 for bead of radius 1 nm
#pmass   = 1.5
print('spacing:', spacing)
print('psigma:', psigma)
# I need to set the file name in an easier way, but for now I just use this:  ## Might want to add a loop too, if I have more files...

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T        = 3
plotseed = 0
plotdirs = False
test_sectioned = False
zhigh          = 250
zlow           = -50
confignrs    = np.arange(1,1001)
Nseeds       = len(confignrs)      # So that I don't have to change that much
maxz_av      = 0
filescounter = 0
Nsteps       = 2001 # 20001 before
#writeevery   = 10                  # I'm writing to file every this many time steps
unitlength   = 1e-9
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime#*writeevery
Npartitions  = 5 # For extracting more walks from one file (but is it really such a random walk here...?)
minlength    = int(floor(Nsteps/Npartitions)) # For sectioning the data
print('timestepsize:', timestepsize)

endlocation_out  = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' + str(psigma) + '/Nocut/' # This is the new main folder, so I'm going to store the results there.
outfilename      = endlocation_out+'maxzs_vs_d.txt'
outfile          = open(outfilename,'w')
maxzs            =  []

for spacing in spacings:
    endlocation          = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp + str(psigma) + '/'
    filestext            = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])
    # Text files
    maxz = -1000 # No z-position is this small.
    for confignr in confignrs:
        print('d =', spacing , ' On config number:', confignr)
        infilename_all  = endlocation+'all_confignr'+str(confignr)+'.lammpstrj'      #endlocation+'all_density'+str(density)+'_confignr'+str(confignr)+'.lammpstrj'
        infilename_free = endlocation+'freeatom_confignr'+str(confignr)+'.lammpstrj' #endlocation+'freeatom_density'+str(density)+'_confignr'+str(confignr)+'.lammpstrj'
    
        # Read in:
        #### Automatic part
        ## Find the extent of the polymers: Max z-coord of beads in the chains
        try:
            infile_all = open(infilename_all, "r")
        except:
            print('Oh, lammpstrj-file! Where art thou?')
            continue # Skipping this file if it does not exist
        # Moving on, if the file
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
        maxz_temp = maxz
        
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
                    if z>maxz_temp:
                        maxz_temp = z
                i+=1
        ##maxz_av += maxz # Maybe I should use this instead? # No.
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
        
        # Setting arrays for treatment:
        positions = [] # Fill with xu,yu,zu-coordinates. One position for each time step. # We have the same time step for the whole of the simulation. Can bother with time later
        times     = [] # Not useful anymore?
        
        # For now: Only find the largest z-position among the beads in the chain. # Otherwise I must read off the atom number (will do that with the free bead anyways.)
        
        time_start = time.process_time()
        counter = 0
        zs_fortesting = []
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
                z        = float(words[5])
            zs_fortesting.append(z)
            counter+=1
            i+=1
        infile_free.close()
        dt = (times[1]-times[0])*timestepsize # This might be handy
        times_single = np.arange(Nsteps)#*dt
        times_single_real = np.arange(Nsteps)*dt
        
        if max(zs_fortesting)>zhigh:
            continue
        if min(zs_fortesting)<zlow:
            continue


        maxwalker = max(zs_fortesting)
        # Some testing for max value of z here
        if maxwalker>maxz_temp:
            maxz = maxwalker
        else:
            maxz = maxz_temp
        time_end = time.process_time()
    outfile.write('%.2f %.16e\n')
    maxzs.append(maxz)

