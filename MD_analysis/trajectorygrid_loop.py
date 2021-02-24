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
#thr = 0.002 # Revisit this after first run

#
damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
spacings = [1,2,3,4,5,6,7,8,10,15,25,50,75,100]
psigma   = 1

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T        = 3
plotseed = 0
plotdirs = False
test_sectioned = False
confignrs    = np.arange(1,1001)
Nseeds       = len(confignrs)      # So that I don't have to change that much
maxz_av      = 0
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

outloc   = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_forest/Projected_trajectories/'
plotloc  = outloc + 'Spacings/'
metafilename = outloc + 'metafile.txt' 
filestext = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])

metafile = open(metafilename,'w')
metafile.write('Spacing; minval; maxval; Nx; Ny; filescounter\n')
for spacing in spacings:
    endlocation_in         = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_forest/Spacing'+str(spacing)+'/'
    endlocation           = endlocation_in + 'Results/' # Nocut is default
    outfilename_matrix    = endlocation+'matrix_'+filestext+'_d'+str(spacing)+'_lx'+str(lx)+'_nocut_newnormalization'
    outfilename_midline   = endlocation+'midline_'+filestext+'_d'+str(spacing)+'_lx'+str(lx)+'_nocut.txt'
    outfilename_relfrarea = endlocation+'relativefreearea_'+filestext+'_d'+str(spacing)+'_lx'+str(lx)+'_autothr_nocut.txt'
    plotname_matrix       = plotloc+'matrix_'+filestext+'_d'+str(spacing)+'_lx'+str(lx)+'_nocut_newnormalization.png'
    plotname_midline      = plotloc+'midline_'+filestext+'_d'+str(spacing)+'_lx'+str(lx)+'_nocut.png'
    
    ### Make matrix ###
    infilename = endlocation_in+'freeatom_confignr1.lammpstrj' # For extracting box size
    infile = open(infilename, 'r')
    lines = infile.readlines()
    xlines = lines[5]
    xmin = float(xlines.split()[0])
    xmax = float(xlines.split()[1])
    ylines = lines[6]               # Should be the same as x
    ymin = float(ylines.split()[0])
    ymax = float(ylines.split()[1])
    Lx   = xmax-xmin
    Ly   = ymax-ymin
    Nx   = int(Lx//lx) #int(Lx/lx)
    Ny   = int(Ly//ly) #int(Ly/ly)
    
    xvals = np.linspace(xmin,xmax,Nx)
    yvals = np.linspace(ymin,ymax,Ny) 
    xgrid, ygrid = np.meshgrid(xvals,yvals)
    grid = np.zeros((Nx,Ny))
    
    filescounter = 0
    lentime = 0
    for confignr in confignrs:
        print('d =', spacing, ', on config number:', confignr)
        infilename_free = endlocation_in+'freeatom_confignr'+str(confignr)+'.lammpstrj'
        
        # Read in:
        #### Automatic part
        ## Find the extent of the polymers: This is 100 for forest
        try:
            infile_free = open(infilename_free, "r")
        except:
            print('Oh, lammpstrj-file! Where art thou?')
            continue # Skipping this file if it does not exist
        # Moving on, if the file
        
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
        xes   = [] # NB! These are unwrapped!
        ys    = [] 
        times = []
        
        # For now: Only find the largest z-position among the beads in the chain. # Otherwise I must read off the atom number (will do that with the free bead anyways.)
        maxz = -1000 # No z-position is this small.
        
        time_start = time.process_time()
        counter = 0
        zs_fortesting = []
        while i<totlines:
            words = lines[i].split()
            if (words[0]=='ITEM:'):
                if words[1]=='TIMESTEP':
                    words2 = lines[i+1].split() # The time step is on the next line
                    t = float(words2[0])
                    times.append(t)
                    i+=skiplines
                elif words[1]=='NUMBER': # These will never kick in. 
                    i+=1
                elif words[1]=='BOX':
                    i+=1
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
                if z>extent_polymers:
                     break # Don't record trajectory if the bead is outside it # Bad to mix cut and nocut like this? Since I use the results from nocut?
                xes.append(x)
                ys.append(y)
                zs_fortesting.append(z)
                counter+=1
                i+=1
        infile_free.close()
        dt = (times[1]-times[0])*timestepsize # This might be handy
        time_end = time.process_time()
        lentime += len(times)
        
        if max(zs_fortesting)>zhigh:
            continue
        if min(zs_fortesting)<zlow:
           continue   
        filescounter += 1
    
        for i in range(len(xes)):
            x = xes[i]
            y = ys[i]
            # Rewrapping:
            while x>xmax: # May be several box lengths away.
                x-=Lx
            while x<xmin:
                x+=Lx
            while y>ymax:
                y-=Ly
            while y<ymin:
                y+=Ly
            xind = int((x-xmin)//lx) # Integer division to get on index format. Shift to start at 0.
            yind = int((y-ymin)//ly)
            grid[xind-1,yind-1]+=1 # Switch indexation? # Hope I don't get an overflow. # Normalize after each file?
    
    #grid*=Nx*Ny/lentime # Is this right? # Does not work, somehow.
    maxgrid = grid.max()
    mingrid = grid.min()
    grid/=maxgrid
    
    plt.figure(figsize=(6,5))
    #plt.contourf(xgrid,ygrid,grid)
    plt.imshow(grid, extent=[xmin/lx,xmax/lx,ymin/ly,ymax/ly])
    plt.colorbar()
    plt.title(r'Intensity')
    plt.savefig(plotname_matrix)
    
    np.save(outfilename_matrix,grid)

    # Write to meta file:
    #metafile.write('Spacing; minval; maxval; Nx; Ny; filescounter')
    metafile.write('%i %.16f %.16f %i %i %i\n' % (spacing,mingrid,maxgrid,Nx,Ny,filescounter))
    
    # For relative free area
    maxgrid = grid.max()
    mingrid = grid.min()
    # Skip this and do it another way?
    outfile_relfrarea = open(outfilename_relfrarea,'w')
    ##normgrid = grid/grid.max()
    thr = (maxgrid+mingrid)/2.
    gridmat_binary = grid>thr
    rel_area = np.sum(gridmat_binary)/(Nx*Ny)
    print('gridmat_binary:',gridmat_binary)
    print('rel_area:',rel_area)
    outfile_relfrarea.write('%.16f' % rel_area)
    print('Spacing:',spacing)
    #print('Min, normgrid:',normgrid.min())

    #Normalize first? # Yes.
    outfile_midline = open(outfilename_midline,'w')
    midline = np.zeros(Nx)
    for i in range(Nx):
        midline[i] = grid[i,int(Ny//2)] # Do I actually need an array?
        outfile_midline.write('%.16f %.16f\n' % (xvals[i]/spacing,midline[i]))
    outfile_midline.close()

    plt.figure(figsize=(6,5))
    #plt.contourf(xgrid,ygrid,grid)
    plt.imshow(grid, extent=[xmin/lx,xmax/lx,ymin/ly,ymax/ly]) 
    plt.colorbar()
    plt.title(r'Intensity')

metafile.close()
 