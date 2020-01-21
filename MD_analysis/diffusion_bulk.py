from skimage import morphology, measure, io, util   # For performing morphology operations (should maybe import more from skimage...)
from mpl_toolkits.mplot3d import Axes3D             # Plotting in 3D
import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
from scipy.ndimage import measurements, convolve    # Should have used the command below, but changing now will lead to confusion
from scipy import ndimage                           # For Euclidean distance measurement
import maintools_percolation as perctools
import data_treatment as datr
import numpy as np
import random
import math
import time
import os
import glob
import copy

# Get these from datr:
'''
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

def partition_holders(Nthr,Nsteps,minlength):   # Perhaps I should have this in it's own module. Right now it is just a copy of the version I have in percolation_from_MD.py 
    interval     = minlength                           # Spacing between start points
    end_interval = Nsteps-interval+1             # +1 because of the indexing and how arange works (but will check and see later)
    startpoints  = np.arange(0,end_interval, interval) # Points to start the calculation of the rmsd.
    numberofsamples = len(startpoints)                 # Number of such walks. Will be the number of lists
    len_all = 0
    partitioned_walks_holder = []
    steps_holder             = []
    lengths                  = []
    for i in range(numberofsamples):
        length = Nsteps-i*interval
        len_all += length
        lengths.append(length)
        partitioned_walks_holder.append(np.zeros((Nthr,length))) # Should hold in general # Should I use Nseeds here? Or just store the average?
        steps_holder.append(np.arange(0,length))
    return steps_holder, partitioned_walks_holder, numberofsamples, len_all, lengths


def partition_rw(thesepositions,partition_walks, Nsteps, minlength, thisthr, Nthr): # Or should I use something like Nsections instead?
    interval     = minlength                           # Spacing between start points
    end_interval = Nsteps-interval+1             # +1 because of the indexing and how arange works (but will check and see later)
    startpoints  = np.arange(0,end_interval, interval) # Points to start the calculation of the rmsd.
    numberofsamples = len(startpoints)                 # Number of such walks. Will be the number of lists
    for i in range(numberofsamples): # Looping over the number of partitioned walks
        startindex = startpoints[i]
        startpoint = thesepositions[startindex]
        counter    = 1
        length     = Nsteps-i*interval
        these_rmsd = np.zeros((Nthr,length))
        this_index = startindex+counter
        while this_index<Nsteps:
            this_point = thesepositions[this_index]
            distvec    = this_point-startpoint
            rmsd       = np.dot(distvec,distvec)
            #print('thisthr:', thisthr, '; counter:', counter)
            these_rmsd[thisthr,counter] += rmsd 
            this_index+=1
            counter   +=1 
        partition_walks[i] += these_rmsd
        #steps.append(thesesteps)
        #partitioned_walks.append(these_rmsd)
    return partition_walks

def partition_holders_averaged(Nsteps,minlength):   # Perhaps I should have this in it's own module. Right now it is just a copy of the version I have in percolation_from_MD.py 
    interval     = minlength                           # Spacing between start points
    end_interval = Nsteps-interval+1             # +1 because of the indexing and how arange works (but will check and see later)
    startpoints  = np.arange(0,end_interval, interval) # Points to start the calculation of the rmsd.
    numberofsamples = len(startpoints)                 # Number of such walks. Will be the number of lists
    len_all = 0
    partitioned_walks_holder = []
    steps_holder             = []
    lengths                  = []
    for i in range(numberofsamples):
        length = Nsteps-i*interval
        len_all += length
        lengths.append(length)
        partitioned_walks_holder.append(np.zeros(length)) # Smaller than it's original counterpart
        steps_holder.append(np.arange(0,length))
    return steps_holder, partitioned_walks_holder, numberofsamples, len_all, lengths
'''

# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
spacing = 10 #100
psigma  = 1
# I need to set the file name in an easier way, but for now I just use this:  ## Might want to add a loop too, if I have more files...

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
plotdirs = False
startpart = '_'#'_withsubstrate_'#'_'#
if startpart=='_':
    parentfolder = 'Pure_bulk/'
else:
    parentfolder = 'Bulk_substrate/'
test_sectioned = True
seeds  = np.arange(1,1001)#[23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]
Nseeds = len(seeds)
Nsteps = 20001
minlength = 5 # For sectioning the data
unitlength   = 1e-9
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime
print('timestepsize:', timestepsize)
Npartitions  = 5 # For extracting more walks from one file (but is it really such a random walk here...?)
nextpart = 'ljunits_'
namemid  = 'Langevin_scaled_T3'#_pmass1.5'
nameend  = '_sect_placeexact_ljepsilon0.730372054992096_ljcut1p122'
namebase = 'bulkdiffusion'+startpart+'ljunits_spacing%i_' % spacing + namemid +'_psigma'+str(psigma)+ nameend
folderbase = 'Bulkdiffusion'+startpart+nextpart+namemid+nameend
namebase_short = namebase # In case I ever need to shorten it
endlocation          = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/'+parentfolder+folderbase+'/Spacing'+str(spacing)+'/Sigma_bead_' + str(psigma)+'/'
outfilename          = endlocation+'lammpsdiffusion_'+namebase_short+'.txt'
outfilename_ds       = endlocation+'lammpsdiffusion_'+namebase+'_av_ds.txt'
outfilename_gamma    = endlocation+'lammpsdiffusion_'+namebase_short+'_zimportance.txt'
outfilename_sections = endlocation+'lammpsdiffusion_'+namebase_short+'_sections.txt'
plotname             = endlocation+'lammpsdiffusion_'+namebase_short+'.png'
plotname_gamma       = endlocation+'lammpsdiffusion_'+namebase_short+'_zimportance.png'
plotname_SI          = endlocation+'lammpsdiffusion_'+namebase_short+'_SI.png'
plotname_velocity    = endlocation+'lammpsdiffusion_'+namebase_short+'_velocity.png'
plotname_velocity_SI = endlocation+'lammpsdiffusion_'+namebase_short+'_velocity_SI.png'
plotname_parallel_orthogonal = endlocation+'lammpsdiffusion_'+namebase+'_par_ort.png'
plotname_sectioned_average = endlocation+'lammpsdiffusion_'+namebase+'_sections.png'
# Setting arrays
# Prepare for sectioning distance data:
steps, partition_walks, numberofsamples, len_all, lengths = datr.partition_holders_averaged(Nsteps,minlength)
# These are all squared:
# All together:
allRs     = [] # Distance
alldxs    = []
alldys    = []
alldzs    = []
'''
allvs     = [] # Velocity
allvxs    = []
allvys    = []
allvzs    = []
'''
# Averages:
averageRs = np.zeros(Nsteps)   # Distances
averagedxs = np.zeros(Nsteps)
averagedys = np.zeros(Nsteps)
averagedzs = np.zeros(Nsteps)
averagevs  = np.zeros(Nsteps)  # Velocities
averagevxs = np.zeros(Nsteps)
averagevys = np.zeros(Nsteps)
averagevzs = np.zeros(Nsteps)
averagedparallel = np.zeros(Nsteps)
average_counter  = np.zeros(Nsteps)
average_walks    = copy.copy(partition_walks)
average_counters = copy.copy(partition_walks)

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
alltimes  = []
#times_single = np.zeros(Nsteps) # I set this in the loop, after having extracted dt


for seed in seeds:
    print('On seed', seed)
    seedstr = str(seed)
    #infilename_free = endlocation+'pmass'+str(pmass)+'_seed'+seedstr+'.lammpstrj' # Old file name
    infilename_free   = endlocation+'seed'+seedstr+'.lammpstrj'
    plotname_dirs     = endlocation+'lammpsdiffusion_'+namebase_short+'_dxdydzR2_seed'+seedstr+'.png'
    plotname_testsect = endlocation+'lammpsdiffusion_'+namebase_short+'_testsectioned_seed'+seedstr+'.png'
    #print('infilename:',infilename_free)
    # Read in:
    #### Automatic part
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
    vx = [] # x-component of the velocity. 
    vy = []
    vz = []
    positions = [] # Fill with xu,yu,zu-coordinates. One position for each time step. # We have the same time step for the whole of the simulation. Can bother with time later
    times     = [] # Not useful anymore?
    
    # For now: Only find the largest z-position among the beads in the chain. # Otherwise I must read off the atom number (will do that with the free bead anyways.)
    maxz = -1000 # No z-position is this small.
    # id type xu yu zu vx vy vz 
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
            #print('words:', words)
            # Find properties
            # Order:  id  type xu  yu   zu  vx  vy  vz # I don't have mol anymore!!! #### IS THIS A PURE DIFFUSION THING, OR DOES IT HAPPEN FOR SUBSTRATE_BEAD TOO?????
            #         [0] [1]  [2] [3]  [4] [5] [6] [7]
            # Order:  id  type mol ux  uy  uz  vx  vy   vz
            #         [0] [1]  [2] [3] [4] [5] [6] [7]  [8]
            ind      = int(words[0])-1 # Atom ids go from zero to N-1.
            #atomtype = int(words[1]) 
            #molID    = int(words[2])
            x        = float(words[2]) #float(words[3])
            y        = float(words[3]) #float(words[4])
            z        = float(words[4]) #float(words[5])
            positions.append(np.array([x,y,z]))
            vx.append(float(words[5]))
            vy.append(float(words[6]))
            vz.append(float(words[7]))
            #velocity.append(np.sqrt(vx*vx + vy*vy + vz*vz)) # Shouldn't I only perform this when I need it?
            counter+=1
            i+=1
    infile_free.close()
    dt = (times[1]-times[0])*timestepsize # This might be handy
    #print('dt:', dt)
    #print('times[1]-times[0]:',times[1]-times[0])
    times_single = np.arange(Nsteps)#*dt
    times_single_real = np.arange(Nsteps)*dt
    
    time_end = time.process_time()
    
    
    
    # I do not take into account that the free bead can enter the polymer grid and then exit again. If that happens, there might be discontinuities or weird kinks in the data (if you look at the graphs...) # I do not have this in my test-dataset. Maybe make the bead lighter and see what happens?
    
    Nin         = counter
    startpos_in = positions[0]
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
        # Distance
        this_in = positions[i]
        dist = this_in-startpos_in
        R2   = np.dot(dist,dist)
        dx2  = dist[0]*dist[0]
        dy2  = dist[1]*dist[1]
        dz2  = dist[2]*dist[2]
        gamma = (R2-dz2)/R2
        # Velocity:
        vxi = vx[i]
        vyi = vy[i]
        vzi = vz[i]
        # All together:
        allRs.append(R2)
        alldxs.append(dx2)
        alldys.append(dy2)
        alldzs.append(dz2)
        alltimes.append(i)
        # Averages:
        averageRs[i] += R2   # Distance
        averagedxs[i]+= dx2
        averagedys[i]+= dy2
        averagedzs[i]+= dz2
        averagevs[i] += np.sqrt(vxi*vxi + vyi*vyi + vzi*vzi)  # Velocity
        averagevxs[i]+= vxi
        averagevys[i]+= vyi
        averagevzs[i]+= vzi
        averagedparallel[i] += dx2+dy2 # Distance
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
    
    rmsd_with_linefit = datr.rmsd(R_temp,linefit)
    
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
        plt.title('RMSD in bulk, seed %s, d = %i nm' % (seedstr,spacing))
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.savefig(plotname_dirs)
    
    # Partitioned:
    # Npartitions
    part_walks    = []
    part_steps    = []
    part_start    = 0
    smallest_part_walk = int(math.floor(Nin/Npartitions))
    for i in range(Npartitions):
        this_part    = []
        these_steps  = []
        part_start   = i*smallest_part_walk
        len_thiswalk = Nin-part_start       # Do I need this?
        #print('part_start:', part_start)
        #print('len(positions):', len(positions))
        rstart = positions[part_start]
        for j in range(len_thiswalk):
            rthis = positions[part_start+j]
            drvec = rthis - rstart
            dr2   = np.dot(drvec,drvec)
            average_walks[i][j] += dr2 # Notation? [i,j] or [i][j] # Ugh, and should I have one of these for R2, dx2, dy2 and dz2?
            average_counters[i][j] += 1
            this_part.append(dr2)
            these_steps.append(j)
        part_walks.append(this_part)
        part_steps.append(these_steps)
    sections_walk.append(part_walks)
    sections_steps.append(part_steps)
    
    if test_sectioned==True and (seed==29 or seed==47 or seed==53 or seed==59 or seed==83 or seed==103):
        plt.figure(figsize=(6,5))
        for i in range(Npartitions):
            plt.plot(part_steps[i], part_walks[i], label='Section %i' % i)
        plt.xlabel(r'Step number')
        plt.ylabel(r'Distance$^2$ [in unit length]')
        plt.title('RMSD in bulk, seed %s, d = %i nm' % (seedstr,spacing))
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.savefig(plotname_testsect)
    # I tried calling the function before, but that was slow...
    #print('Hubba hubba zoot zoot')
    #time_beforepartition = time.process_time()
    #partition_walks = datr.partition_dist_averaged(positions, partition_walks, Nsteps, minlength) # partition_walks is being updated inside of the function
    #time_afterpartition = time.process_time()
    #print('Deba uba zat zat')
    #print('Time, partitioning:',time_afterpartition-time_beforepartition)

allRs      = np.array(allRs)
alltimes   = np.array(alltimes)
gamma_avgs = np.array(gamma_avgs)
single_slopes = np.array(single_slopes)
rmsds         = np.array(rmsds)
Nins          = np.array(Nins)

# Sorted arrays:
sorted_gamma_avgs = np.sort(gamma_avgs)
plotgammaagainst  = np.arange(Nseeds)
index_sorted      = np.argsort(gamma_avgs)

# Making arrays
rmsds_sorted         = np.zeros(Nseeds)
seeds_sorted         = np.zeros(Nseeds)
Nins_sorted          = np.zeros(Nseeds)
single_slopes_sorted = np.zeros(Nseeds)

## Calculating averages:
# All data
for i in range(Nsteps):
    counter = average_counter[i]
    if counter!=0:
        # Distance
        averageRs[i]/=counter
        averagedxs[i]/=counter
        averagedys[i]/=counter
        averagedzs[i]/=counter
        averagedparallel[i]/= counter
        # Velocity
        averagevs[i]/=counter
        averagevxs[i]/=counter
        averagevys[i]/=counter
        averagevzs[i]/=counter

# Sectioned walks
for i in range(Npartitions):
    for j in range(lengths[i]):
        if average_counters[i][j]!=0:
            average_walks[i][j]/=average_counters[i][j]

## Opening file
outfile_gamma = open(outfilename_gamma, 'w')
outfile_gamma.write('This is an attempt to find the relative contribution of dz to R. The smaller the value in the second column, the more important dz is.\nOrder: <(R^2(n)-dz^2(n))/R^2(n)>_n, slope a of line fit, rmsd R^2(n) line fit, Nin, seed\n')
for i in range(Nseeds):
    index = index_sorted[i]
    rmsds_sorted[i] = rmsds[index]
    seeds_sorted[i] = seeds[index]
    Nins_sorted[i]  = Nins[index]
    single_slopes_sorted[i] = single_slopes[index]
    outfile_gamma.write('%.16f %.16f %.16f %i %i\n' % (sorted_gamma_avgs[i], single_slopes_sorted[i], rmsds_sorted[i], Nins_sorted[i], seeds_sorted[i]))
outfile_gamma.close()


## Finding the diffusion coefficient (in polymer)
# Nonsense as of now
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
coeffs_poly, covs = polyfit(times_single, averageRs, 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a_poly_av = coeffs_poly[0]
b_poly_av = coeffs_poly[1]
rms_D_poly_av = np.sqrt(covs[0,0])/6.
rms_b_poly_av = np.sqrt(covs[1,1])
D_poly_av = a_poly_av/6.
 
fit_poly_av = a_poly_av*times_single+b_poly_av


print('Fit, all data:')
print('b_poly (should be small):', b_poly)
print('D_poly:', D_poly)

print('Fit, average:')
print('b_poly (should be small):', b_poly_av)
print('D_poly:', D_poly_av)

## Average, SI units:

averageRs_SI  = averageRs*unitlength**2
averagedxs_SI = averagedxs*unitlength**2
averagedys_SI = averagedys*unitlength**2
averagedzs_SI = averagedzs*unitlength**2
averagedparallel_SI = averagedparallel*unitlength**2
# Velocity:
averagevs_SI  = averageRs*unitlength/timestepsize
averagevxs_SI = averagedxs*unitlength/timestepsize
averagevys_SI = averagedys*unitlength/timestepsize
averagevzs_SI = averagedzs*unitlength/timestepsize


coeffs_poly, covs = polyfit(times_single_real, averageRs_SI, 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a_poly_SI = coeffs_poly[0]
b_poly_SI = coeffs_poly[1]
rms_D_poly_SI = np.sqrt(covs[0,0])/6.
rms_b_poly_SI = np.sqrt(covs[1,1])
D_poly_SI = a_poly_SI/6.
 
fit_poly_SI = a_poly_SI*times_single_real+b_poly_SI

print('Fit, average, SI:')
print('b_poly (should be small):', b_poly_SI)
print('D_poly:', D_poly_SI)

outfile = open(outfilename, 'w')
outfile.write('D_poly: %.16f %.16f\n' % (D_poly,rms_D_poly))
outfile.write('b_poly: %.16f %.16f\n' % (b_poly,rms_b_poly))
outfile.close()

outfile_ds = open(outfilename_ds,'w')
outfile_ds.write('Time step; Time; <R^2>; <dx^2>; <dy^2>; <dz^2>; <dx^2+dy^2>\n')
for i in range(len(averageRs)):
    outfile_ds.write('%i %.2e %.5e %.5e %.5e %.5e %.5e\n' % (times_single[i], times_single_real[i], averageRs_SI[i], averagedxs_SI[i], averagedys_SI[i], averagedzs_SI[i], averagedparallel_SI[i]))
outfile_ds.close()

outfile_sections = open(outfilename_sections, 'w')
for i in range(Npartitions):
    outfile_sections.write('Section %i\n' % i) # When reading from the file afterwards, I can extract the section number by if words[0]=='Section' and use that to divide the data into sections
    for j in range(average_walks[i]):
        time_this = steps[i][j]*timestepsize
        dist_this = average_walks[i][j]*unitlength
        outfile_sections.write('%.16f %16.f\n' % (time_this,dist_this))
outfile_sections.close()

## Making figures.
plt.figure(figsize=(6,5))
plt.plot(alltimes, allRs, ',', label='Data, brush')
plt.plot(times_single, averageRs, ',', label='Average, brush')
plt.plot(alltimes, fit_poly, '--', label='Fit, data, brush')
plt.plot(times_single, fit_poly_av, '--', label='Fit, average, brush')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title('RMSD in bulk, d = %i nm' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.plot(times_single, averagedzs, ',', label=r'dz$^2$')
plt.plot(times_single, averagedparallel, ',', label=r'dx$^2$+dy$^2$')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title('Averaged RMSD in bulk, d = %i nm' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_parallel_orthogonal)

plt.figure(figsize=(6,5))
plt.plot(plotgammaagainst, sorted_gamma_avgs)
plt.xlabel(r'Different runs')
plt.ylabel(r'$\gamma = <\frac{R^2-dz^2}{R^2}>$')
plt.title(r'RMSD in bulk, $\gamma$; d = %i nm' % spacing)
plt.tight_layout()
plt.savefig(plotname_gamma)

plt.figure(figsize=(6,5))
plt.plot(times_single_real, averageRs_SI, ',', label='Average, brush')
plt.plot(times_single_real, fit_poly_SI, '--', label='Fit, average, brush')
plt.xlabel(r'Time [s]')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title('RMSD in bulk, d = %i nm, SI' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_SI)

plt.figure()
plt.plot(times_single, averagevs, label=r'v')
plt.plot(times_single, averagevxs, label=r'vx')
plt.plot(times_single, averagevys, label=r'vy')
plt.plot(times_single, averagevzs, label=r'vz')
plt.xlabel(r'Step number')
plt.ylabel(r'Velocity [in unit length/time]')
plt.title('Averaged velocity in bulk, system size by d = %i nm' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_velocity)

plt.figure()
plt.plot(times_single_real, averagevs_SI, label=r'v')
plt.plot(times_single_real, averagexvs_SI, label=r'vx')
plt.plot(times_single_real, averagevys_SI, label=r'vy')
plt.plot(times_single_real, averagevzs_SI, label=r'vz')
plt.xlabel(r'Time [s]')
plt.ylabel(r'Velocity [m/s]')
plt.title('Averaged velocity in bulk, system size by d = %i nm, SI' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_velocity_SI)

plt.figure()
for i in range(Npartitions):
    plt.plot()
plt.savefig(plotname_sectioned_average)

