from skimage import morphology, measure, io, util   # For performing morphology operations (should maybe import more from skimage...)
from mpl_toolkits.mplot3d import Axes3D             # Plotting in 3D
import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
from scipy.ndimage import measurements, convolve    # Should have used the command below, but changing now will lead to confusion
from scipy import ndimage                           # For Euclidean distance measurement
import data_treatment as datr
import numpy as np
import random
import math
import time
import os
import glob
import copy

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
spacings  = [1,1.25,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,10,15,25,50,75,100]
#              1    1.25  1.5    2    2.5   3    3.5   4    4.5   5     6     7     8     10
longbools =[False,False,False,False,True,False,True,False,True,False,False,False,False,False, False,False,False,False,False]
Nd = len(spacings)

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
test_sectioned = False
zhigh          = 290    # For omitting unphysical trajectories
zlow           = -50    
testh          = 20     # For transport measurements (default 50, might change for tests)
confignrs    = np.arange(1,1001)#300)#20)#101)#1001)#22) 
Nseeds       = len(confignrs)      # So that I don't have to change that much
maxz_av      = 0
filescounter = 0
#if long==True:
#    Nsteps   = 10001
Nsteps   = 2001                     # We want to study the exits in a shorter time
#writeevery   = 10                  # I'm writing to file every this many time steps
unitlength   = 1e-9
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime#*writeevery
Npartitions  = 5 # For extracting more walks from one file (but is it really such a random walk here...?)
minlength    = int(floor(Nsteps/Npartitions)) # For sectioning the data
print('timestepsize:', timestepsize)

endlocation  = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Exitfrac_vs_d/'
filestext    = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])
# Text files
outfilename  = endlocation+'Nesc_div_Nread_'+filestext+'.txt'

# Plots
plotname        = endlocation+'Nesc_div_Nread_'+filestext+'.png'
plotname_smalld = endlocation+'Nesc_div_Nread_'+filestext+'_smalld.png'

## Setting arrays
# Prepare for sectioning distance data:
time_walks_SI, steps, partition_walks, numberofsamples, len_all, lengths, startpoints = datr.partition_holders_averaged(Nsteps,minlength)


escfraction = np.zeros(Nd)

outfile = open(outfilename, 'w')

for l in range(Nd):
    spacing = spacings[l]
    long =longbools[l]
    escapedcounter = 0
    Nread          = 0
    endlocation          = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp + str(psigma) + '/'
    for confignr in confignrs:
        print('d =', spacing,'; On config number:', confignr)
        if long==True:
            infilename_all  = endlocation+'long/'+'all_confignr'+str(confignr)+'_long.lammpstrj'
            infilename_free = endlocation+'long/'+'freeatom_confignr'+str(confignr)+'_long.lammpstrj'
        else:
            infilename_all  = endlocation+'all_confignr'+str(confignr)+'.lammpstrj'
            infilename_free = endlocation+'freeatom_confignr'+str(confignr)+'.lammpstrj'
        
        # Read in:
        #### Automatic part
        ## Find the extent of the polymers: Max z-coord of beads in the chains
        try:
            infile_all = open(infilename_all, "r")
        except:
            try:
                infilename_all = endlocation+'long/'+'all_confignr'+str(confignr)+'.lammpstrj'
                infilename_free = endlocation+'long/'+'freeatom_confignr'+str(confignr)+'.lammpstrj'
                infile_all = open(infilename_all, "r")
            except:
                print('Oh, lammpstrj-file! Where art thou?')
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
            if words[0]=='ITEM:':
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
        maxz_av += maxz
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
        vx = [] # x-component of the velocity. 
        vy = []
        vz = []
        positions = [] # Fill with xu,yu,zu-coordinates. One position for each time step. # We have the same time step for the whole of the simulation. Can bother with time later
        times     = [] # Not useful anymore?
        
        # For now: Only find the largest z-position among the beads in the chain. # Otherwise I must read off the atom number (will do that with the free bead anyways.)
        maxz = -1000 # No z-position is this small.
        
        time_start = time.process_time()
        counter = 0
        zs_fortesting = []
        while i<totlines:
            words = lines[i].split()
            if words[0]=='ITEM:': # Some double testing going on...
                if words[1]=='TIMESTEP':
                    words2 = lines[i+1].split() # The time step is on the next line
                    t = float(words2[0])
                    if t>200000:
                        break # Want to compare short and long simulations on equal footing
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
                positions.append(np.array([x,y,z]))
                vx.append(float(words[6]))
                vy.append(float(words[7]))
                vz.append(float(words[8]))
                zs_fortesting.append(z)
                counter+=1
                i+=1
        infile_free.close()
        dt = (times[1]-times[0])*timestepsize # This might be handy
        times_single = np.arange(Nsteps)#*dt
        times_single_real = np.arange(Nsteps)*dt
        
        time_end = time.process_time()
        
        if max(zs_fortesting)>zhigh:
            continue
        if min(zs_fortesting)<zlow:
            continue
        
        Nread += 1 # Counting the number of files we analyze
        pos_inpolymer = []
        
        Nin   = 0
        maxzpol = 0
        brokenthis = 0
        for i in range(counter):
            thesepos = positions[i]
            z        = thesepos[2]
            if z>extent_polymers: # If the polymer is in bulk # We don't want it to go back and forth between brush and bulk # That will cause discontinuities in our data
                brokenthis = 1
                escapedcounter +=1
                break
            else:
                pos_inpolymer.append(thesepos)
                if z>maxzpol:
                    maxzpol = z
                Nin+=1
    escfrac_this = escapedcounter/float(Nread)
    outfile.write('%.2f %.5e\n' % (spacing,escfrac_this))
    escfraction[l] = escfrac_this
outfile.close()

plt.figure(figsize=(6,5))
plt.plot(spacings,escfraction)
plt.xlabel('$d$ [nm]')
plt.ylabel(r'$N_{esc}/N_{sims}$')
plt.title(r'Fraction of escaped beads vs spacing $d$')
plt.tight_layout()
plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.plot(spacings[:13],escfraction[:13])
plt.xlabel('$d$ [nm]')
plt.ylabel(r'$N_{esc}/N_{sims}$')
plt.title(r'Fraction of escaped beads vs spacing $d$')
plt.tight_layout()
plt.savefig(plotname_smalld)
