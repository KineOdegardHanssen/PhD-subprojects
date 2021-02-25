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

# Resolution
lx = 0.5  # nm. Want this smaller than sigma.
ly = lx
thr = 0.002 # Revisit this after first run
autothr = True # Determines the implementation

#
damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
spacings = [1,2,3,4,5,6,7,8,10,15,25,50,75,100]
psigma  = 1
print('spacing:', spacing)
print('psigma:', psigma)
# I need to set the file name in an easier way, but for now I just use this:  ## Might want to add a loop too, if I have more files...

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T        = 3
plotseed = 0
plotdirs = False
test_sectioned = False
confignrs    = np.arange(1,1001)
Nseeds       = len(confignrs)      # So that I don't have to change that much
maxz_av      = 0
filescounter = 0
Nsteps       = 2001 # 20001 before
zhigh        = 250
zlow         = -50
#writeevery   = 10                  # I'm writing to file every this many time steps
unitlength   = 1e-9
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime#*writeevery
Npartitions  = 5 # For extracting more walks from one file (but is it really such a random walk here...?)
minlength    = int(floor(Nsteps/Npartitions)) # For sectioning the data
print('timestepsize:', timestepsize)
extent_polymers = 100
outcounter      = 0


filestext            = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])
outlocation          = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_forest/Projected_trajectories/'
if autothr==True:
    plotfilename         = outlocation +'relativefreearea_lx'+str(lx)+'_autothr.png'
    plotfilename_smalld  = outlocation +'relativefreearea_lx'+str(lx)+'_autothr.png'
else:
    plotfilename         = outlocation +'relativefreearea_lx'+str(lx)+'_thr'+str(thr)+'_nocut.png'
    plotfilename_smalld  = outlocation +'relativefreearea_lx'+str(lx)+'_thr'+str(thr)+'_nocut_smalld.png'

relfreeareas = []
for spacing in spacings:
    endlocation_in = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_forest/Spacing'+str(spacing)+'/'
    endlocation = endlocation_in + 'Results/' # Nocut is default
    if autothr==True:
        infilename  = endlocation+'relativefreearea_'+filestext+'_d'+str(spacing)+'_lx'+str(lx)+'_autothr.txt'
    else:
        infilename  = endlocation+'relativefreearea_'+filestext+'_d'+str(spacing)+'_lx'+str(lx)+'_thr'+str(thr)+'_nocut.txt'
    infile = open(infilename,'r')
    lines  = infile.readlines()
    for line in lines:
        words = line.split()
        if len(words)!=0:
            relfreeareas.append(float(words[0]))
    infile.close()

plt.figure(figsize=(6,5))
plt.plot(spacings,relfreeareas,'-o')
plt.xlabel(r'Chain spacing $d$')
plt.ylabel(r'Relative free area for bead')
plt.title('Relative free area for bead vs chain spacing')
plt.tight_layout()
plt.savefig(plotfilename)

plt.figure(figsize=(6,5))
plt.plot(spacings[:9],relfreeareas[:9],'-o')
plt.xlabel(r'Chain spacing $d$')
plt.ylabel(r'Relative free area for bead')
plt.title('Relative free area for bead vs chain spacing')
plt.tight_layout()
plt.savefig(plotfilename_smalld)
plt.show()
