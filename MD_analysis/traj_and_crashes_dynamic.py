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
long       = True
spacing    = 4.5
psigma     = 1
sigma_atom = psigma
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
#seeds  = [23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]
confignrs    = [1,200,500,750,999]
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
plotlocation          = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Singletrajectories/Spacing'+str(spacing)+'/'
filestext            = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])


########## Collision files:
basefolder_coll = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp+str(sigma_atom) + '/'
folder_coll = basefolder_coll + 'log-files/'
if long==True:
    folder_coll = folder_coll + 'long/'
outfolder_coll = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Singletrajectories/Spacing' + str(spacing) + '/'
ncollfilename = outfolder_coll + 'Ncoll_selectedconfigs.txt'
ncollfile     = open(ncollfilename,'w')
##########

## Setting arrays
# Prepare for sectioning distance data:
time_walks_SI, steps, partition_walks, numberofsamples, len_all, lengths, startpoints = datr.partition_holders_averaged(Nsteps,minlength)


Nread        = 0
skippedfiles = 0


for confignr in confignrs:
    print('On config number:', confignr)
    if long==True:
        infilename_all  = endlocation+'long/'+'all_confignr'+str(confignr)+'_long.lammpstrj'
        infilename_free = endlocation+'long/'+'freeatom_confignr'+str(confignr)+'_long.lammpstrj'
    else:
        infilename_all  = endlocation+'all_confignr'+str(confignr)+'.lammpstrj'
        infilename_free = endlocation+'freeatom_confignr'+str(confignr)+'.lammpstrj'
    plotname     = plotlocation+'traj_d'+str(spacing)+'_config%i.png' % confignr
    plotname_ind = plotlocation+'traj_d'+str(spacing)+'_config%i_index.png' % confignr
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
    
    # Mark point where the bead exits the brush (if it does)
    Nin   = 0
    maxzpol = 0
    brokenthis = 0
    exitt = False
    x_inpolymer  = []
    y_inpolymer  = []
    z_inpolymer  = []
    xy_inpolymer = []
    time_inpolymer = []
    time_real_inpolymer = []
    for i in range(counter):
        thesepos = positions[i]
        z        = thesepos[2]
        if z>extent_polymers: # If the polymer is in bulk # We don't want it to go back and forth between brush and bulk # That will cause discontinuities in our data
            # Saving info on exit times:
            exitt = i*dt
            break
        else:
            x = thesepos[0]
            y = thesepos[1]
            xy = np.sqrt(x*x+y*y)
            x_inpolymer.append(x)
            y_inpolymer.append(y)
            z_inpolymer.append(z)
            xy_inpolymer.append(xy)
            time_inpolymer.append(times_single[i])
            time_real_inpolymer.append(times_single_real[i])
            if z>maxzpol:
                maxzpol = z
            Nin+=1
    # Finding the collision times:
    filename_infolder = 'log.confignr%i_printevery10' % confignr
    if long==True:
        filename_infolder = filename_infolder + '_long'
    infilename = folder_coll + filename_infolder
    #print('infilename:', infilename)
    
    try:
        infile_coll = open(infilename, "r")
    except:
        try:
            filename_infolder = 'log.confignr%i' % confignr
            infilename = folder_coll + filename_infolder
            infile_coll = open(infilename, "r")
        except:
            print('Oh, log-file! Where art thou?')
            continue # Skipping this file if it does not exist
    # Moving on, if the file
    lines = infile_coll.readlines()
    
    startit = 0
    times_coll   = []
    crashtimes   = []
    pot_energies = []
    for line in lines:
        words = line.split()
        if len(words)>1:
            if words[0]=='Loop':
                break
            if startit==1:
                times_coll.append(int(words[0]))
                pot_energies.append(float(words[1]))
            if words[0]=='Step':
                startit = 1
    infile_coll.close()
    
    N = len(times_coll)
    crashbefore = 0
    thistime = 0
    i = 0
    while thistime<=200001:
        pot_en   = pot_energies[i]
        thistime = times_coll[i]
        if (pot_en!=0 and crashbefore==0):
            # Append collition times
            crashtimes.append(thistime/100.)
            crashbefore = 1
        elif pot_en==0:
            crashbefore = 0
        i+=1
    
    # Plotting the crashes
    zmax_traj  = max(z_inpolymer)
    xymax_traj = max(xy_inpolymer)
    maxtraj    = 0
    if zmax_traj>=xymax_traj:
        maxtraj  = zmax_traj
    else:
        maxtraj  = xymax_traj
    yvals    = np.zeros(2)
    yvals[1] = maxtraj
    Ncoll    = len(crashtimes)
    ncollfile.write('%i %i\n' % (confignr,Ncoll))
    
    #print('crashtimes:', crashtimes)
    # Plot
    plt.figure(figsize=(10,5))
    ax = plt.subplot(111)
    ax.plot(time_real_inpolymer,z_inpolymer,label=r'$z$')
    ax.plot(time_real_inpolymer,xy_inpolymer,label=r'$\sqrt{x^2+y^2}$')
    for l in range(Ncoll):
        tcoll = np.zeros(2) + crashtimes[l]*dt
        ax.plot(tcoll,yvals, '--', color='silver')
    plt.xlabel('Time [s]')
    plt.ylabel('Distance travelled [nm]')
    plt.title(r'$z(t)$ and $\sqrt{x^2(t)+y^2(t)}$, d=%.2f, config %i' % (spacing,confignr))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.savefig(plotname)
    
    plt.figure(figsize=(10,5))
    ax = plt.subplot(111)
    ax.plot(time_inpolymer,z_inpolymer,label=r'$z$')
    ax.plot(time_inpolymer,xy_inpolymer,label=r'$\sqrt{x^2+y^2}$')
    for l in range(Ncoll):
        tcoll = np.zeros(2) + crashtimes[l]
        ax.plot(tcoll,yvals, '--', color='silver')
    plt.xlabel('Time [s]')
    plt.ylabel('Distance travelled [nm]')
    plt.title(r'$z(t)$ and $\sqrt{x^2(t)+y^2(t)}$, d=%.2f, config %i' % (spacing,confignr))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.savefig(plotname_ind)
ncollfile.close()