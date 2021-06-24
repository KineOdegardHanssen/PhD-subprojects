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
long    = False
spacing = 3
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

endlocation_in           = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_nostiffness/Spacing'+str(spacing)+'/'
endlocation              = endlocation_in+'Cut/'
filestext                = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])
# Text files
outfilename_ds           = endlocation+'av_ds_'+filestext
outfilename_gamma        = endlocation+'zimportance_'+filestext
outfilename_sections     = endlocation+'sections_'+filestext
outfilename_maxz         = endlocation+'maxz_az_'+filestext
outfilename_dt           = endlocation+'dt'
outfilename_cuts         = endlocation+'cuts'
outfilename_noexit       = endlocation+'noexitzs'
outfilename_skippedfiles = endlocation+'skippedfiles'
outfilename_exittimes    = endlocation+'exittimes'
outfilename_th           = endlocation+'ths_h%i' % testh +filestext

# Plots
plotname             = endlocation+filestext
plotname_all         = endlocation+'all_'+filestext
plotname_gamma       = endlocation+'zimportance_'+filestext
plotname_SI          = endlocation+'SI_'+filestext
plotname_parallel_orthogonal = endlocation+'par_ort_'+filestext
plotname_dx_dy_dz    = endlocation+'dx_dy_dz_'+filestext
plotname_parallel    = endlocation+'par_'+filestext
plotname_orthogonal  = endlocation+'ort_'+filestext
plotname_short_all   = endlocation+'short_all_'+filestext
plotname_velocity    = endlocation+'velocity_'+filestext
plotname_velocity_SI = endlocation+'velocity_SI_'+filestext
plotname_velocity_sq = endlocation+'velocity_sq_'+filestext
plotname_sectioned_average = endlocation+'sections_'+filestext
plotname_sectioned_average_vs_steps = endlocation+'sections_steps_'+filestext
plotname_traj_xy    = endlocation+'traj_xy_config'+str(confignrs[-1])
plotname_traj_xz    = endlocation+'traj_xz_config'+str(confignrs[-1])
plotname_traj_yz    = endlocation+'traj_yz_config'+str(confignrs[-1])
plotname_traj_xt    = endlocation+'traj_xt_config'+str(confignrs[-1])
plotname_traj_yt    = endlocation+'traj_yt_config'+str(confignrs[-1])
plotname_traj_zt    = endlocation+'traj_zt_config'+str(confignrs[-1])
plotname_th_hist    = endlocation+'th_hist_h%i' % testh +filestext
plotname_exittimes  = endlocation+'exittimes_'+filestext
plotname_exittimes_binned = endlocation+'exittimes_hist_'+filestext

if long==True:
    outfilename_ds           = outfilename_ds+'_long.txt'
    outfilename_gamma        = outfilename_gamma +'_long.txt'
    outfilename_sections     = outfilename_sections+'_long.txt'
    outfilename_maxz         = outfilename_maxz+'_long.txt'
    outfilename_dt           = outfilename_dt+'_long.txt'
    outfilename_cuts         = outfilename_cuts+'_long.txt'
    outfilename_noexit       = outfilename_noexit+'_long.txt'
    outfilename_skippedfiles = outfilename_skippedfiles+'_long.txt'
    outfilename_exittimes    = outfilename_exittimes+'_long.txt'
    outfilename_th           = outfilename_th+'_long.txt'
    
    # Plots
    plotname             = plotname+'_long.png'
    plotname_all         = plotname_all+'_long.png'
    plotname_gamma       = plotname_gamma+'_long.png'
    plotname_SI          = plotname_SI+'_long.png'
    plotname_parallel_orthogonal = plotname_parallel_orthogonal+'_long.png'
    plotname_dx_dy_dz    = plotname_dx_dy_dz+'_long.png'
    plotname_parallel    = plotname_parallel+'_long.png'
    plotname_orthogonal  = plotname_orthogonal+'_long.png'
    plotname_short_all   = plotname_short_all+'_long.png'
    plotname_velocity    = plotname_velocity+'_long.png'
    plotname_velocity_SI = plotname_velocity_SI+'_long.png'
    plotname_velocity_sq = plotname_velocity_sq+'_long.png'
    plotname_sectioned_average = plotname_sectioned_average+'_long.png'
    plotname_sectioned_average_vs_steps = plotname_sectioned_average_vs_steps+'_long.png'
    plotname_traj_xy    = plotname_traj_xy+'_long.png'
    plotname_traj_xz    = plotname_traj_xz+'_long.png'
    plotname_traj_yz    = plotname_traj_yz+'_long.png'
    plotname_traj_xt    = plotname_traj_xt+'_long.png'
    plotname_traj_yt    = plotname_traj_yt+'_long.png'
    plotname_traj_zt    = plotname_traj_zt+'_long.png'
    plotname_th_hist    = plotname_th_hist+'_long.png'
    plotname_exittimes  = plotname_exittimes+'_long.png'
    plotname_exittimes_binned = plotname_exittimes_binned+'_long.png'
else:
    outfilename_ds           = outfilename_ds+'.txt'
    outfilename_gamma        = outfilename_gamma +'.txt'
    outfilename_sections     = outfilename_sections+'.txt'
    outfilename_maxz         = outfilename_maxz+'.txt'
    outfilename_dt           = outfilename_dt+'.txt'
    outfilename_cuts         = outfilename_cuts+'.txt'
    outfilename_noexit       = outfilename_noexit+'.txt'
    outfilename_skippedfiles = outfilename_skippedfiles+'.txt'
    outfilename_exittimes    = outfilename_exittimes+'.txt'
    outfilename_th           = outfilename_th+'.txt'
    
    # Plots
    plotname             = plotname+'.png'
    plotname_all         = plotname_all+'.png'
    plotname_gamma       = plotname_gamma+'.png'
    plotname_SI          = plotname_SI+'.png'
    plotname_parallel_orthogonal = plotname_parallel_orthogonal+'.png'
    plotname_dx_dy_dz    = plotname_dx_dy_dz+'.png'
    plotname_parallel    = plotname_parallel+'.png'
    plotname_orthogonal  = plotname_orthogonal+'.png'
    plotname_short_all   = plotname_short_all+'.png'
    plotname_velocity    = plotname_velocity+'.png'
    plotname_velocity_SI = plotname_velocity_SI+'.png'
    plotname_velocity_sq = plotname_velocity_sq+'.png'
    plotname_sectioned_average = plotname_sectioned_average+'.png'
    plotname_sectioned_average_vs_steps = plotname_sectioned_average_vs_steps+'.png'
    plotname_traj_xy    = plotname_traj_xy+'.png'
    plotname_traj_xz    = plotname_traj_xz+'.png'
    plotname_traj_yz    = plotname_traj_yz+'.png'
    plotname_traj_xt    = plotname_traj_xt+'.png'
    plotname_traj_yt    = plotname_traj_yt+'.png'
    plotname_traj_zt    = plotname_traj_zt+'.png'
    plotname_th_hist    = plotname_th_hist+'.png'
    plotname_exittimes  = plotname_exittimes+'.png'
    plotname_exittimes_binned = plotname_exittimes_binned+'.png'

## Setting arrays
# Prepare for sectioning distance data:
time_walks_SI, steps, partition_walks, numberofsamples, len_all, lengths, startpoints = datr.partition_holders_averaged(Nsteps,minlength)

print('Sections, walks:')
print('minlength:', minlength)
print('lenghts:', lengths)
# These are all squared:
# All together:
allRs     = []
alldxs    = []
alldys    = []
alldzs    = []
# Averages:
averageR2s  = np.zeros(Nsteps)   # Distances, squared
averagedx2s = np.zeros(Nsteps)
averagedy2s = np.zeros(Nsteps)
averagedz2s = np.zeros(Nsteps)
averagedxs  = np.zeros(Nsteps)   # Distances, non-squared
averagedys  = np.zeros(Nsteps)
averagedzs  = np.zeros(Nsteps)
averagev2s  = np.zeros(Nsteps)  # Velocities, squared
averagevx2s = np.zeros(Nsteps)
averagevy2s = np.zeros(Nsteps)
averagevz2s = np.zeros(Nsteps)
averagevs   = np.zeros(Nsteps)  # Velocities
averagevxs  = np.zeros(Nsteps)
averagevys  = np.zeros(Nsteps)
averagevzs  = np.zeros(Nsteps)
averagevparallel  = np.zeros(Nsteps)
averagevparallel2 = np.zeros(Nsteps) # Velocities, squared
averagedparallel2 = np.zeros(Nsteps) # Distances, squared
average_counter   = np.zeros(Nsteps)
average_walks     = copy.deepcopy(partition_walks)
average_walks_SI  = copy.deepcopy(partition_walks)
average_walks_z2     = copy.deepcopy(partition_walks)
average_walks_z2_SI  = copy.deepcopy(partition_walks)
average_walks_p2     = copy.deepcopy(partition_walks)
average_walks_p2_SI  = copy.deepcopy(partition_walks)
average_counters     = copy.deepcopy(partition_walks) # Okay, this name is confusing.


# Separated by seed:
Rs_byseed    = []
dxs_byseed   = []
dys_byseed   = []
dzs_byseed   = []
gamma_byseed = []
times_byseed = []
gamma_avgs   = []
single_slopes = []
rmsds         = []
Nins          = []
#
sections_walk  = []
sections_steps = []
# This is not squared, obviously:
alltimes    = []
exittimes   = []
exitzs      = [] # A safeguard
testh_times = []
exiths      = [] # A safeguard


Nread        = 0
skippedfiles = 0

outfile_cuts = open(outfilename_cuts, 'w')
outfile_cuts.write('confignr time_cut zs\n')

outfile_noexit = open(outfilename_noexit,'w')
outfile_noexit.write('confignr zs\n')

outfile_skippedfiles = open(outfilename_skippedfiles, 'w')

for confignr in confignrs:
    print('On config number:', confignr)
    if long==True:
        infilename_all  = endlocation_in+'long/'+'all_confignr'+str(confignr)+'_long.lammpstrj'
        infilename_free = endlocation_in+'long/'+'freeatom_confignr'+str(confignr)+'_long.lammpstrj'
        plotname_dirs   = endlocation+'dxdydzR2_seed'+str(confignr)+'_long.png'
        plotname_testsect = endlocation+'testsectioned_seed'+str(confignr)+'_long.png'
    else:
        infilename_all  = endlocation_in+'all_confignr'+str(confignr)+'.lammpstrj'
        infilename_free = endlocation_in+'freeatom_confignr'+str(confignr)+'.lammpstrj'
        plotname_dirs   = endlocation+'dxdydzR2_seed'+str(confignr)+'.png'
        plotname_testsect = endlocation+'testsectioned_seed'+str(confignr)+'.png'

    
    
    #print('infilename_all:',infilename_all)
    
    # Read in:
    #### Automatic part
    ## Find the extent of the polymers: Max z-coord of beads in the chains
    try:
        infile_all = open(infilename_all, "r")
    except:
        try:
            infilename_all = endlocation_in+'long/'+'all_confignr'+str(confignr)+'.lammpstrj'
            infilename_free = endlocation_in+'long/'+'freeatom_confignr'+str(confignr)+'.lammpstrj'
            infile_all = open(infilename_all, "r")
        except:
            print('Oh, lammpstrj-file! Where art thou?')
            outfile_skippedfiles.write('%i\n' % confignr)
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
    #print('dt:', dt)
    #print('times[1]-times[0]:',times[1]-times[0])
    times_single = np.arange(Nsteps)#*dt
    times_single_real = np.arange(Nsteps)*dt
    
    time_end = time.process_time()
    
    if max(zs_fortesting)>zhigh:
        continue
    if min(zs_fortesting)<zlow:
        continue
    
    Nread += 1 # Counting the number of files we analyze
    pos_inpolymer = []
    
    #before_in
    
    # I do not take into account that the free bead can enter the polymer grid and then exit again. If that happens, there might be discontinuities or weird kinks in the data (if you look at the graphs...) # I do not have this in my test-dataset. Maybe make the bead lighter and see what happens?
    
    Nin   = 0
    maxzpol = 0
    brokenthis = 0
    for i in range(counter):
        thesepos = positions[i]
        z        = thesepos[2]
        if z>extent_polymers: # If the polymer is in bulk # We don't want it to go back and forth between brush and bulk # That will cause discontinuities in our data
            brokenthis = 1
            outfile_cuts.write('%i %i ' % (confignr, i))
            for iex in range(Nin): # Ugh, stupid loop structure... 
                posnow = positions[iex] 
                outfile_cuts.write('%.16e ' % posnow[2])
            outfile_cuts.write('\n')
            # Saving info on exit times:
            exittimes.append(i*dt) # Change to time step
            exitzs.append(z)
            break
        else:
            pos_inpolymer.append(thesepos)
            if z>maxzpol:
                maxzpol = z
            Nin+=1
    if brokenthis==0:
        outfile_noexit.write('%i ' % confignr) 
        for inoex in range(Nin):
            posnow = positions[inoex] 
            outfile_noexit.write('%.16e ' % posnow[2])
        outfile_noexit.write('\n')
    
    # testh_times:
    for i in range(counter):
        thesepos = positions[i]
        z        = thesepos[2]
        if z>testh: # If the polymer is in bulk # Measuring the time when is first reaches height h.
            testh_times.append(i*dt)
            exiths.append(z)
            break
    
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
    gamma_temp = []
    
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
        dx   = dist[0]         # Signed
        dy   = dist[1]
        dz   = dist[2]
        R2   = np.dot(dist,dist) # Squared
        dx2  = dx*dx
        dy2  = dy*dy
        dz2  = dz*dz
        gamma = (R2-dz2)/R2
        # Velocity:
        vxi  = vx[i]
        vyi  = vy[i]
        vzi  = vz[i]
        vx2i = vxi*vxi # squared
        vy2i = vyi*vyi
        vz2i = vzi*vzi
        v2i  = vx2i + vy2i + vz2i
        # All together:
        allRs.append(R2)
        alldxs.append(dx2)
        alldys.append(dy2)
        alldzs.append(dz2)
        alltimes.append(i)
        # Averages:
        averageR2s[i] +=R2   # Distance
        averagedx2s[i]+=dx2
        averagedy2s[i]+=dy2
        averagedz2s[i]+=dz2
        averagedxs[i] +=dx   # Distance, signed
        averagedys[i] +=dy
        averagedzs[i] +=dz
        averagevs[i]  += np.sqrt(vx2i+vy2i+vz2i)  # Velocity
        averagevxs[i] += vxi
        averagevys[i] += vyi
        averagevzs[i] += vzi
        averagev2s[i] += vx2i+vy2i+vz2i            # Velocity, squared
        averagevx2s[i]+= vx2i
        averagevy2s[i]+= vy2i
        averagevz2s[i]+= vz2i
        averagevparallel[i]  += np.sqrt(vx2i+vy2i)
        averagevparallel2[i] += vx2i+vy2i
        averagedparallel2[i] += dx2+dy2 # Distance
        average_counter[i] +=1
        # Separated by seed:
        R_temp.append(R2)
        dx_temp.append(dx2)
        dy_temp.append(dy2)
        dz_temp.append(dz2)
        gamma_temp.append(gamma)
        step_temp.append(i)
    Rs_byseed.append(R_temp)   # Don't need these as of yet. Might use them later.
    dxs_byseed.append(dx_temp)
    dys_byseed.append(dy_temp)
    dzs_byseed.append(dz_temp)
    gamma_byseed.append(gamma_temp)
    times_byseed.append(step_temp)
    gamma_avg = np.mean(gamma_temp)
    gamma_avgs.append(gamma_avg)
    
    coeffs = np.polyfit(step_temp,R_temp,1)
    a = coeffs[0]
    b = coeffs[1]
    
    linefit = a*np.array(step_temp)+b
    
    rmsd_with_linefit = rmsd(R_temp,linefit)
    
    single_slopes.append(a)
    rmsds.append(rmsd_with_linefit)
    Nins.append(Nin)
    
    if plotdirs==True:
        plt.figure(figsize=(6,5))
        plt.plot(step_temp, R_temp, label=r'$<R^2>$')
        plt.plot(step_temp, dx_temp, label=r'$<dx^2>$')
        plt.plot(step_temp, dy_temp, label=r'$<dy^2>$')
        plt.plot(step_temp, dz_temp, label=r'$<dz^2>$')
        plt.xlabel(r'Step number')
        plt.ylabel(r'Distance$^2$ [in unit length]')
        plt.title('Random walk in polymer brush, config %i, d = %i nm' % (confignr,spacing))
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.savefig(plotname_dirs)
    
    # Partitioned:
    # Npartitions
    part_walks    = []
    part_steps    = []
    part_start    = 0
    #smallest_part_walk = int(math.floor(Nin/Npartitions))
    for i in range(Npartitions):
        this_part    = []
        these_steps  = []
        part_start   = startpoints[i]
        len_thiswalk = lengths[i]
        part_end     = part_start+len_thiswalk
        if part_start>(Nin-1):              # Do not start this section if the bead has left the brush.
            break
        rstart = pos_inpolymer[part_start]
        for j in range(len_thiswalk):
            #print('i:',i,', j:',j)
            if part_start+j>(Nin-1):
                break
            rthis =  pos_inpolymer[part_start+j]
            drvec = rthis - rstart
            dr2   = np.dot(drvec,drvec)
            dz2   = drvec[2]*drvec[2]
            dpar2 = drvec[0]*drvec[0]+drvec[1]*drvec[1]
            '''
            print('part_start:',part_start, '; part_start+j:', part_start+j, ';   j:',j)
            print('rthis:',rthis)
            print('rstart:',rstart)
            print('drvec:', drvec)
            print('dr2:', dr2)
            '''
            average_walks[i][j]    +=dr2
            average_walks_z2[i][j] +=dz2
            average_walks_p2[i][j] +=dpar2
            average_counters[i][j] +=1
            this_part.append(dr2)
            these_steps.append(j)
        part_walks.append(this_part)
        part_steps.append(these_steps)
    sections_walk.append(part_walks)
    sections_steps.append(part_steps)
    
    if test_sectioned==True and (config==1 or seed==5 or seed==11 or seed==15 or seed==17 or seed==21): # This does not make much sense now... I wanted to inspect some plots.
        plt.figure(figsize=(6,5))
        for i in range(Npartitions):
            plt.plot(part_steps[i], part_walks[i], label='Section %i' % i)
        plt.xlabel(r'Step number')
        plt.ylabel(r'Distance$^2$ [in unit length]')
        plt.title('Random walk in polymer brush, config %i, d = %i nm' % (confignr,spacing))
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.savefig(plotname_testsect)
    # Do something else than this? Use the data I already have?
    #time_beforepartition = time.process_time()
    #partition_walks = datr.partition_dist_averaged(positions, partition_walks, Nsteps, minlength) # partition_walks is being updated inside of the function
    #time_afterpartition = time.process_time()
    #print('Time, partitioning:',time_afterpartition-time_beforepartition)

outfile_cuts.close()
outfile_noexit.close()
outfile_skippedfiles.close()

print('filescounter:', filescounter)
maxz_av /= filescounter

print('maxz_av:', maxz_av)

newlen = len(gamma_avgs)

allRs      = np.array(allRs)
alltimes   = np.array(alltimes)
gamma_avgs = np.array(gamma_avgs)
single_slopes = np.array(single_slopes)
rmsds         = np.array(rmsds)
Nins          = np.array(Nins)

# Sorted arrays:
sorted_gamma_avgs = np.sort(gamma_avgs)
plotgammaagainst  = np.arange(newlen)
index_sorted      = np.argsort(gamma_avgs)

# Making arrays
rmsds_sorted         = np.zeros(newlen)
seeds_sorted         = np.zeros(newlen)
Nins_sorted          = np.zeros(newlen)
single_slopes_sorted = np.zeros(newlen)

# Calculating averages:
#Know the velocity at t=0, same for all
vxi = vx[0]
vyi = vy[0]
vzi = vz[0]
averagevs[0]  = np.sqrt(vxi*vxi + vyi*vyi + vzi*vzi)
averagevxs[0] = vxi
averagevys[0] = vyi
averagevzs[0] = vzi
averagevparallel[0] = np.sqrt(vxi*vxi+vyi*vyi)

Ninbrush = Nsteps # Default, in case it does not exit the brush
for i in range(1,Nsteps):
    counter = average_counter[i]
    if counter!=0:
        # Distance squared
        averageR2s[i]/=counter
        averagedx2s[i]/=counter
        averagedy2s[i]/=counter
        averagedz2s[i]/=counter
        averagedparallel2[i]/= counter
        # Distance, signed:
        averagedxs[i]/=counter
        averagedys[i]/=counter
        averagedzs[i]/=counter
        # Velocity
        averagevs[i]/=counter
        averagevxs[i]/=counter
        averagevys[i]/=counter
        averagevzs[i]/=counter
        averagevparallel[i]/=counter
        # Velocity squared
        averagev2s[i]/=counter
        averagevx2s[i]/=counter
        averagevy2s[i]/=counter
        averagevz2s[i]/=counter
        averagevparallel2[i]/=counter
    else:
       Ninbrush = i-1
       break

times_single      = times_single[0:Ninbrush]
times_single_real = times_single_real[0:Ninbrush]
# Distance squared
averageR2s  = averageR2s[0:Ninbrush]
averagedx2s = averagedx2s[0:Ninbrush]
averagedy2s = averagedy2s[0:Ninbrush]
averagedz2s = averagedz2s[0:Ninbrush]
averagedparallel2 = averagedparallel2[0:Ninbrush]
# Distance to the first power
averagedxs = averagedxs[0:Ninbrush]
averagedys = averagedys[0:Ninbrush]
averagedzs = averagedzs[0:Ninbrush]
# Velocity to the first power
averagevs  = averagevs[0:Ninbrush]
averagevxs = averagevxs[0:Ninbrush]
averagevys = averagevys[0:Ninbrush]
averagevzs = averagevzs[0:Ninbrush]
averagevparallel = averagevparallel[0:Ninbrush]
# Velocity squared
averagev2s  = averagev2s[0:Ninbrush]
averagevx2s = averagevx2s[0:Ninbrush]
averagevy2s = averagevy2s[0:Ninbrush]
averagevz2s = averagevz2s[0:Ninbrush]
averagevparallel2 = averagevparallel2[0:Ninbrush]

print('len(averageR2s):',len(averageR2s))
print('Nsteps:', Nsteps)

# Sectioned walks
for i in range(Npartitions):
    for j in range(lengths[i]):
        if average_counters[i][j]!=0:
            #print('Before division: Average_walks=', average_walks[i][j], '     average_counters=',average_counters[i][j])
            average_walks[i][j]   /=average_counters[i][j]
            average_walks_z2[i][j]/=average_counters[i][j]
            average_walks_p2[i][j]/=average_counters[i][j]
            #print('After division: Average_walks=', average_walks[i][j], '     average_counters=',average_counters[i][j])

if Nseeds<12:
    print('gamma_avgs:',gamma_avgs)
    print('index_sorted:',index_sorted)
# Opening file
outfile_gamma = open(outfilename_gamma, 'w')
outfile_gamma.write('This is an attempt to find the relative contribution of dz to R. The smaller the value in the second column, the more important dz is.\nOrder: <(R^2(n)-dz^2(n))/R^2(n)>_n, slope a of line fit, rmsd R^2(n) line fit, Nin, confignr\n')
for i in range(newlen):
    index = index_sorted[i]
    rmsds_sorted[i] = rmsds[index]
    seeds_sorted[i] = confignrs[index]
    Nins_sorted[i]  = Nins[index]
    single_slopes_sorted[i] = single_slopes[index]
    outfile_gamma.write('%.16f %.16f %.16f %i %i\n' % (sorted_gamma_avgs[i], single_slopes_sorted[i], rmsds_sorted[i], Nins_sorted[i], seeds_sorted[i]))
outfile_gamma.close()


## Finding the diffusion coefficient (in polymer)
# Finding the last 25% of the data:               # Do this later! In loop! # But do now just for testing purposes.
# Performing the line fit:
coeffs_poly, covs = polyfit(alltimes, allRs, 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a_poly = coeffs_poly[0]
b_poly = coeffs_poly[1]
rms_D_poly = np.sqrt(covs[0,0])/6.
rms_b_poly = np.sqrt(covs[1,1])
D_poly = a_poly/6.

'''
coeffs, covs = np.polyfit(insteps, indata, 1, full=False, cov=True) # Using all the data in the fit
    av_slope[j]  = coeffs[0]
    rms_slope[j] = np.sqrt(covs[0,0])   
'''
 
fit_poly = a_poly*alltimes+b_poly

# Line fit to averaged data:
coeffs_poly, covs = polyfit(times_single, averageR2s, 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a_poly_av = coeffs_poly[0]
b_poly_av = coeffs_poly[1]
rms_D_poly_av = np.sqrt(covs[0,0])/6.
rms_b_poly_av = np.sqrt(covs[1,1])
D_poly_av = a_poly_av/6.
 
fit_poly_av = a_poly_av*times_single+b_poly_av

'''
print('Fit, all data:')
print('b_poly (should be small):', b_poly)
print('D_poly:', D_poly)

print('Fit, average:')
print('b_poly (should be small):', b_poly_av)
print('D_poly:', D_poly_av)
'''

## Average, SI units:
# Distance:
averageR2s_SI  = averageR2s*unitlength**2
averagedx2s_SI = averagedx2s*unitlength**2
averagedy2s_SI = averagedy2s*unitlength**2
averagedz2s_SI = averagedz2s*unitlength**2
averagedparallel2_SI = averagedparallel2*unitlength**2
# Velocity:
averagevs_SI  = averagevs*unitlength/timestepsize
averagevxs_SI = averagevxs*unitlength/timestepsize
averagevys_SI = averagevys*unitlength/timestepsize
averagevzs_SI = averagevzs*unitlength/timestepsize
averagevparallel_SI = averagevparallel*unitlength/timestepsize

#print('----------------------------------')
#print(' averageR2s_SI[2]:', averageR2s_SI[2])

coeffs_poly, covs = polyfit(times_single_real, averageR2s_SI, 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a_poly_SI = coeffs_poly[0]
b_poly_SI = coeffs_poly[1]
rms_D_poly_SI = np.sqrt(covs[0,0])/6.
rms_b_poly_SI = np.sqrt(covs[1,1])
D_poly_SI = a_poly_SI/6.

#print(' averageR2s_SI[2]:', averageR2s_SI[2])

#print('----------------------------------')
fit_poly_SI = a_poly_SI*times_single_real+b_poly_SI

#print('Fit, average, SI:')
#print('b_poly (should be small):', b_poly_SI)
#print('D_poly:', D_poly_SI)

# Do not perform fit here:
'''
outfile = open(outfilename, 'w')
outfile.write('D_poly: %.16f %.16f\n' % (D_poly,rms_D_poly))
outfile.write('b_poly: %.16f %.16f\n' % (b_poly,rms_b_poly))
outfile.close()
'''

outfile_ds = open(outfilename_ds,'w')
outfile_ds.write('Time step; Time; <R^2>; <dx^2>; <dy^2>; <dz^2>; <dx^2+dy^2>\n')
for i in range(len(averageR2s)):
    outfile_ds.write('%i %.16e %.16e %.16e %.16e %.16e %.16e\n' % (times_single[i], times_single_real[i], averageR2s_SI[i], averagedx2s_SI[i], averagedy2s_SI[i], averagedz2s_SI[i], averagedparallel2_SI[i]))
outfile_ds.close()

xpos = []
ypos = []
zpos = []

outfilename_xoft     = endlocation+'xoft_config%i.txt' % confignr
outfilename_yoft     = endlocation+'yoft_config%i.txt' % confignr
outfilename_zoft     = endlocation+'zoft_config%i.txt' % confignr


outfile_xoft = open(outfilename_xoft,'w')
outfile_yoft = open(outfilename_yoft,'w')
outfile_zoft = open(outfilename_zoft,'w')

print('len(positions):',len(positions))
print('len(times_single):',len(times_single))
'''# I still have some trouble when I don't have enough files (maybe when the last one is not made yet?)
Nloop = len(positions)
if len(times_single)<Nloop:
    Nloop = len(times_single)
'''
for i in range(len(times_single)):
    pos = positions[i]
    xthis = pos[0]
    ythis = pos[1]
    zthis = pos[2]
    xpos.append(xthis)
    ypos.append(ythis)
    zpos.append(zthis)
    outfile_xoft.write('%i %.16f\n' % (times_single[i],xthis))
    outfile_yoft.write('%i %.16f\n' % (times_single[i],ythis))
    outfile_zoft.write('%i %.16f\n' % (times_single[i],zthis))
outfile_xoft.close()
outfile_yoft.close()
outfile_zoft.close()



# Trajectory in plane
plt.figure(figsize=(6,5))
plt.plot(xpos, ypos, '.')
plt.xlabel(r'x')
plt.ylabel(r'y')
plt.title('Trajectory in xy-plane')
plt.tight_layout()
plt.savefig(plotname_traj_xy)

plt.figure(figsize=(6,5))
plt.plot(xpos, zpos, '.')
plt.xlabel(r'x')
plt.ylabel(r'z')
plt.title('Trajectory in xz-plane')
plt.tight_layout()
plt.savefig(plotname_traj_xz)

plt.figure(figsize=(6,5))
plt.plot(ypos, zpos, '.')
plt.xlabel(r'y')
plt.ylabel(r'z')
plt.title('Trajectory in yz-plane')
plt.tight_layout()
plt.savefig(plotname_traj_yz)

### x, y, z vs t
plt.figure(figsize=(6,5))
plt.plot(times_single, xpos, '.')
plt.xlabel(r't')
plt.ylabel(r'x')
plt.title('x(t)')
plt.tight_layout()
plt.savefig(plotname_traj_xt)

plt.figure(figsize=(6,5))
plt.plot(times_single, ypos, '.')
plt.xlabel(r't')
plt.ylabel(r'y')
plt.title('y(t)')
plt.tight_layout()
plt.savefig(plotname_traj_yt)

plt.figure(figsize=(6,5))
plt.plot(times_single, zpos, '.')
plt.xlabel(r't')
plt.ylabel(r'z')
plt.title('z(t)')
plt.tight_layout()
plt.savefig(plotname_traj_zt)


outfile_sections = open(outfilename_sections, 'w')
for i in range(Npartitions):
    outfile_sections.write('Section %i\n' % i) # When reading from the file afterwards, I can extract the section number by if words[0]=='Section' and use that to divide the data into sections
    for j in range(lengths[i]):
        average_walks_SI[i][j]    = average_walks[i][j]*unitlength 
        average_walks_z2_SI[i][j] = average_walks_z2[i][j]*unitlength 
        average_walks_p2_SI[i][j] = average_walks_p2[i][j]*unitlength 
        time_walks_SI[i][j]       = steps[i][j]*timestepsize
        #print('time_walks_SI[i][j]:', time_walks_SI[i][j], ': steps[i][j]:',steps[i][j], ';   timestepsize:',timestepsize, '; steps[i][j]*timestepsize:', steps[i][j]*timestepsize)
        outfile_sections.write('%.16f %.16f %.16f %.16f\n' % (time_walks_SI[i][j],average_walks_SI[i][j],average_walks_p2_SI[i][j],average_walks_z2_SI[i][j]))
outfile_sections.close()


outfile_maxz = open(outfilename_maxz, 'w')
outfile_maxz.write('%.16f' % maxz_av)         # Should I have some rms value here too?
outfile_maxz.close()

outfile_dt = open(outfilename_dt, 'w')
outfile_dt.write('%.16e' % dt)
outfile_dt.close()

if popup_plots==True:
    maxind = 100#00
    plt.figure(figsize=(6,5))
    plt.plot(times_single_real, averagevs_SI, label=r'$<R^2>$')
    plt.plot(times_single_real, averagevxs_SI, label=r'$<dx^2>$')
    plt.plot(times_single_real, averagevys_SI, label=r'$<dy^2>$')
    plt.plot(times_single_real, averagevzs_SI, label=r'$<dz^2>$')
    plt.plot(times_single_real, averagevparallel_SI, label=r'$<dx^2+dy^2>$')
    plt.xlabel(r'Index (s)')
    plt.ylabel(r'Distance$^2$ [in unit length]')
    plt.title(r'RMSD in brush, d = %i nm, $\sigma_b=%.2f$, SI' % (spacing,psigma))
    plt.tight_layout()
    plt.legend(loc='upper left')
    plt.axis([0, times_single_real[maxind], min(averagevxs_SI[0:maxind]), max(averagevs_SI[0:maxind])])
    plt.show()

'''
plt.figure()
for i in range(Npartitions):
    plt.plot(steps[:][i],average_walks[:][i], label='Start at step %i' % startpoints[i])
plt.xlabel(r'Time step')
plt.ylabel(r'Distance$^2$ [m]')
plt.title('RMSD in brush, d = %i nm, SI, sectioning' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.show()


maxind = 10000
plt.figure(figsize=(6,5))
plt.plot(times_single_real, averageR2s_SI, label=r'$<R^2>$')
plt.plot(times_single_real, averagedz2s_SI, label=r'$<dz^2>$')
plt.plot(times_single_real, averagedparallel2_SI, label=r'$<dx^2+dy^2>$')
plt.xlabel(r'Index (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %i nm, $\sigma_b=%.2f$, SI' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.axis([0, times_single_real[maxind], 0, max(averageR2s_SI[0:maxind])])
plt.show()
'''
# To determine the range
xmax_plot = 100
plt.figure(figsize=(6,5))
plt.plot(times_single, averageR2s, label=r'$<R^2>$')
plt.plot(times_single, averagedx2s, label=r'$<dx^2>$')
plt.plot(times_single, averagedy2s, label=r'$<dy^2>$')
plt.plot(times_single, averagedz2s, label=r'$<dz^2>$')
plt.plot(times_single, averagedparallel2, label=r'$<dx^2+dy^2>$')
plt.xlabel(r'Index (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.axis([0, xmax_plot, 0, max(averageR2s[0:xmax_plot])])
plt.savefig(plotname_short_all)

plt.figure(figsize=(6,5))
plt.plot(alltimes, allRs, ',', label='Data, brush')
plt.plot(times_single, fit_poly_av, '--', label='Fit, brush')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.plot(alltimes, allRs, ',', label='Data, brush')
plt.plot(times_single, averageR2s, label='Average, brush')
plt.plot(alltimes, fit_poly, '--', label='Fit, data, brush')
plt.plot(times_single, fit_poly_av, '--', label='Fit, average, brush')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_all)

plt.figure(figsize=(6,5))
plt.plot(times_single, averageR2s, label=r'$<R^2>$')
plt.plot(times_single, averagedx2s, label=r'$<dx^2>$')
plt.plot(times_single, averagedy2s, label=r'$<dy^2>$')
plt.plot(times_single, averagedz2s, label=r'$<dz^2>$')
plt.plot(times_single, averagedparallel2, label=r'$<dx^2+dy^2>$')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'Averaged RMSD in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_parallel_orthogonal)

plt.figure(figsize=(6,5))
plt.plot(times_single, averagedzs, ',')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'$<dz^2>$ in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.savefig(plotname_orthogonal)

plt.figure(figsize=(6,5))
plt.plot(times_single, averagedxs, label=r'$<dx>$')
plt.plot(times_single, averagedys, label=r'$<dy>$')
plt.plot(times_single, averagedzs, label=r'$<dz>$')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance [in unit length]')
plt.title(r'Averaged RMSD in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_dx_dy_dz)

plt.figure(figsize=(6,5))
plt.plot(times_single, averagedparallel2, ',')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'$<dx^2+dy^2>$ in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.savefig(plotname_parallel)

plt.figure(figsize=(6,5))
plt.plot(plotgammaagainst, sorted_gamma_avgs)
plt.xlabel(r'Different runs')
plt.ylabel(r'$\gamma = <\frac{R^2-dz^2}{R^2}>$')
plt.title(r'Random walk in polymer brush, $\gamma$; d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.savefig(plotname_gamma)

plt.figure(figsize=(6,5))
plt.plot(times_single_real, averageR2s_SI, ',', label='Average, brush')
plt.plot(times_single_real, fit_poly_SI, '--', label='Fit, average, brush')
plt.xlabel(r'Time (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.legend(loc='upper left')
plt.savefig(plotname_SI)

plt.figure()
plt.plot(times_single, averagevs, label=r'$<v>$')
plt.plot(times_single, averagevxs, label=r'$<v_x>$')
plt.plot(times_single, averagevys, label=r'$<v_y>$')
plt.plot(times_single, averagevzs, label=r'$<v_z>$')
plt.plot(times_single, averagevparallel, label=r'$<v_\parallel>$')
plt.xlabel(r'Step number')
plt.ylabel(r'Velocity [in unit length/time]')
plt.title('Averaged velocity in brush, d = %i nm' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_velocity)

plt.figure()
plt.plot(times_single, averagev2s, label=r'$<v^2>$')
plt.plot(times_single, averagevx2s, label=r'$<v_x^2>$')
plt.plot(times_single, averagevy2s, label=r'$<v_y^2>$')
plt.plot(times_single, averagevz2s, label=r'$<v_z^2>$')
plt.plot(times_single, averagevparallel2, label=r'$<v_x^2+vy^2>$')
plt.xlabel(r'Step number')
plt.ylabel(r'Velocity [(length/time)$^2$]')
plt.title('Averaged velocity in brush, d = %i nm' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_velocity_sq)

plt.figure()
plt.plot(times_single_real, averagevs_SI, label=r'v')
plt.plot(times_single_real, averagevxs_SI, label=r'vx')
plt.plot(times_single_real, averagevys_SI, label=r'vy')
plt.plot(times_single_real, averagevzs_SI, label=r'vz')
plt.xlabel(r'Time [s]')
plt.ylabel(r'Velocity [m/s]')
plt.title('Averaged velocity in brush, d = %i nm, SI' % spacing)
plt.tight_layout()
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.legend(loc='upper left')
plt.savefig(plotname_velocity_SI)




'''
print('Original:')
print('steps[0][:]:',steps[0][:])
print('steps[:][0]:',steps[:][0])
print('steps[1][:]:',steps[1][:])
print('steps[:][1]:',steps[:][1])
print('steps[2][:]:',steps[2][:])
print('steps[:][2]:',steps[:][2])


print('SI!!!!:')
print('time_walks_SI[0][:]:',time_walks_SI[0][:])
print('time_walks_SI[:][0]:',time_walks_SI[:][0])
print('time_walks_SI[1][:]:',time_walks_SI[1][:])
print('time_walks_SI[:][1]:',time_walks_SI[:][1])
print('time_walks_SI[2][:]:',time_walks_SI[2][:])
print('time_walks_SI[:][2]:',time_walks_SI[:][2])


print('average_walks_SI[0][:]:',average_walks_SI[0][:])
print('average_walks_SI[:][0]:',average_walks_SI[:][0])
print('average_walks_SI[1][:]:',average_walks_SI[1][:])
print('average_walks_SI[:][1]:',average_walks_SI[:][1])
print('average_walks_SI[2][:]:',average_walks_SI[2][:])
print('average_walks_SI[:][2]:',average_walks_SI[:][2])
'''
print('average_walks:',average_walks)
print('average_walks_SI:',average_walks_SI)
#'''

plt.figure()
for i in range(Npartitions):
    plt.plot(time_walks_SI[:][i],average_walks_SI[:][i], label='Start at step %i' % startpoints[i])
plt.xlabel(r'Time [s]')
plt.ylabel(r'Distance$^2$ [m]')
plt.title('RMSD in brush, d = %i nm, SI, sectioning' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.savefig(plotname_sectioned_average)

plt.figure()
for i in range(Npartitions):
    plt.plot(steps[:][i],average_walks_SI[:][i], label='Start at step %i' % startpoints[i])
plt.xlabel(r'Time step')
plt.ylabel(r'Distance$^2$ [m]')
plt.title('RMSD in brush, d = %i nm, SI, sectioning' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.savefig(plotname_sectioned_average_vs_steps)

#plotname_sectioned_average_vs_steps

#print('counters:',average_counter)

print('spacing:', spacing)
print('psigma:', psigma)

outfile_exittimes = open(outfilename_exittimes,'w')
outfile_exittimes.write('Exit time | Value of z at exit time\n')
for i in range(len(exittimes)):
    outfile_exittimes.write('%.16e %.16e\n' % (exittimes[i],exitzs[i]))
outfile_exittimes.close()

# Exit times
plt.figure()
plt.plot(np.sort(exittimes), 'o') # I'll see how this works, if not I can always replot
plt.ylabel(r'Exit time [s]')
plt.title('Exit time for dynamic brush, d = %i nm' % spacing)
plt.tight_layout()
#plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.savefig(plotname_exittimes)

Nexits = len(exittimes)
Nbins  = int(Nexits/10)      # Don't know if this is the best solution...
hist, bin_edges = np.histogram(exittimes, bins=20)


# Creating histogram 
fig, axs = plt.subplots(1, 1, 
                        figsize =(10, 7),  
                        tight_layout = True) 
  
axs.hist(exittimes, bins = 20) 
plt.xlabel(r'Exit time [s]',fontsize=15)
plt.ylabel('Number of exits',fontsize=15)
plt.title('Exit time for dynamic brush, d = %i nm, histogram' % spacing,fontsize=15)
plt.savefig(plotname_exittimes_binned)



## Distribution of times when walker reaches height testh:
Nth = len(testh_times)

outfile_th = open(outfilename_th,'w')
outfile_th.write('Files read: %i\n' % Nread)
for i in range(Nth):
    outfile_th.write('%.16e %.16e\n' % (testh_times[i],exiths[i]))
outfile_th.close()

hist_th, bin_edges_th = np.histogram(testh_times, bins=20)

# Creating histogram II
# (https://stackoverflow.com/questions/22241240/how-to-normalize-a-histogram-in-python)
#define width of each column
width = bin_edges_th[1]-bin_edges_th[0]
#standardize each column by dividing with the maximum height
hist_th = hist_th/float(Nread)
#plot
plt.bar(bin_edges_th[:-1],hist_th,width = width)
plt.xlabel(r'$t_h$ [s]')
plt.ylabel(r'No. of exits/$N_{sims}$') 
plt.title('$t_h$ in system size by d = %i nm, h = %.1f' % (spacing,testh)) 
plt.savefig(plotname_th_hist)


print('hist, exitzs:',hist)
print('hist, th:',hist_th)

print('bin_edges_exitzs:',bin_edges)
print('bin_edges_th:',bin_edges_th)

