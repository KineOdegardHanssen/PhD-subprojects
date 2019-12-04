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



# I need to set the file name in an easier way, but for now I just use this:  ## Might want to add a loop too, if I have more files...

nameend  = '_ppf_sect_minim'
namebase = 'M9N101_ljunits_gridspacing5_Langevin_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debyecutoff3_chargeelementary-1_effectivedielectric0.00881819074717447_T3_theta0is180_pmass1.5' + nameend
infilename = 'particle_in_chaingrid_all_quadratic_'+namebase+'.lammpstrj'
outfilename = 'lammpsdiffusion_qdrgr_'+namebase+'.txt'

# Read in:
#### Automatic part
infile = open(infilename, "r")
lines = infile.readlines() # This takes some time
# Getting the number of lines, etc.
totlines = len(lines)         # Total number of lines
lineend = totlines-1          # Index of last element

# Extracting the number of atoms:
words = lines[3].split()
print("words:", words)
print("words[0]:", words[0])
Nall = int(words[0])
N    = Nall

words = lines[5].split()
xmin = float(words[0])
xmax = float(words[1])

words = lines[6].split()
ymin = float(words[0])
ymax = float(words[1])

words = lines[7].split()
zmin = float(words[0])
zmax = float(words[1])

Lx = xmax-xmin
Ly = ymax-ymin
Lz = zmax-zmin

skiplines   = 9             # If we hit 'ITEM:', skip this many steps...
skipelem    = 0
sampleevery = 0
i           = int(math.ceil(skipelem*(Nall+9)))
skiplines  += (Nall+skiplines)*sampleevery # Check!

# Setting arrays for treatment:
freeatom_positions = [] # Fill with xu,yu,zu-coordinates. One position for each time step. # We have the same time step for the whole of the simulation. Can bother with time later
#polybeads_temp     = [] # Fill with the position of the polymer beads. Overwrite every time step (... could increase the risk of a superbug)
extent_polymers    = []
times              = [] # Not useful anymore?

# For now: Only find the largest z-position among the beads in the chain. # Otherwise I must read off the atom number (will do that with the free bead anyways.)
maxz = -1000 # No z-position is this small.

time_start = time.process_time()
counter = 0
while i<totlines:
    #print('i:',i)
    words = lines[i].split()
    if i==int(math.ceil(skipelem*(N+9))):
        print("In loop, words[0]=", words[0])
        print("First line I read:", words)
    if (words[0]=='ITEM:' and words[1]=='TIMESTEP'): # Some double testing going on...
        if words[1]=='TIMESTEP':
            words2 = lines[i+1].split() # The time step is on the next line
            t = float(words2[0])
            times.append(t)
            #print("New time step")
            #print("Loop time step:", j)
            if t!=0:
                extent_polymers.append(maxz)
                maxz = -1000
                #print('Is this what takes time?')
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
        elif atomtype==3:
            x        = float(words[3])
            y        = float(words[4])
            freeatom_positions.append([x,y,z])
            print('i:',i,'; counter:', counter)
            counter+=1
            if z<2:
                print('Breaking now.')
                break # Stop the treatment if the atom is close enough to the membrane
        i+=1
time_end = time.process_time()
print('Time to read through atoms:', time_end-time_start,'s')

pos_bulk = []
pos_inpolymer = []

print('len(extent_polymers):',len(extent_polymers))
print('len(freeatom_positions):',len(freeatom_positions))

print('max(extent_polymers):',max(extent_polymers))

#before_in

for i in range(counter):
    thesepos = freeatom_positions[i]
    z        = thesepos[2]
    if z>extent_polymers[i]: # If the polymer is in bulk
        pos_bulk.append(thesepos)      # Need some test to see if it has gone from bulk to polymer or vice versa
    else:
        pos_inpolymer.append(thesepos)
    if z<2:
        break # Just in case the test did not work last time

#print('pos_bulk:',pos_bulk)
#print('pos_inpolymer:',pos_inpolymer)

print('len(pos_bulk):',len(pos_bulk))
print('len(pos_inpolymer):',len(pos_inpolymer))
