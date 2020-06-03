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

# Limiting the walk: Not interested in what happens when the bead has moved too far in the z-direction.
maxh  = 55.940983199999756
testh = 30 # An 'arbitrary' value for now ## Print this part to file (i.e. testh_times)?
Nbins = 20 # Bins for plotting distribution of testh_times

# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
spacing = 10 #100
psigma  = 1
damp    = 10
# I need to set the file name in an easier way, but for now I just use this:  ## Might want to add a loop too, if I have more files...

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
plotdirs = False
startpart = '_'
parentfolder = 'Pure_bulk/'
test_sectioned = True
seeds  = np.arange(1,1001)#(1,1001)#(52,54)#(1,1001)#[23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]
Nseeds = len(seeds)
Nsteps = 2001#20001
testh_times  = np.zeros(Nseeds)
Npartitions  = 5 # For extracting more walks from one file (but is it really such a random walk here...?)
minlength    = int(floor(Nsteps/Npartitions)) # For sectioning the data
#writeevery   = 10                  # I'm writing to file every this many time steps ### THIS IS AUTOMATED!!!
unitlength   = 1e-9
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime#*writeevery # LAMMPS writes the actual time step to file, but I recreate the time array using np.arange. Therefore I need to multiply with writeevery

print('timestepsize:', timestepsize)
#nextpart = 'ljunits_'
#namemid  = 'Langevin_scaled_T3'#_pmass1.5'
#nameend  = '_sect_placeexact_ljepsilon0.730372054992096_ljcut1p122'
#namebase = 'bulkdiffusion'+startpart+'ljunits_spacing%i_' % spacing + namemid +'_psigma'+str(psigma)+ nameend
#folderbase = 'Bulkdiffusion'+startpart+nextpart+namemid+nameend
#namebase_short = namebase # In case I ever need to shorten it
#namebase_out   = namebase_short+'_cut'
filestext      = '_seed'+str(seeds[0])+'to'+str(seeds[-1])
#endlocation          = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/'+parentfolder+folderbase+'/Spacing'+str(spacing)+'/Sigma_bead_' + str(psigma)+'/'
inlocation           = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing%i/damp%i_diffseedLgv/Pure_bulk/Sigma_bead_' % (spacing,damp)+str(psigma) + '/'
endlocation          = inlocation + 'maxh'+str(maxh)+'/'
outfilename          = endlocation+'all_'+filestext+'.txt'  #'lammpsdiffusion_'+namebase_out+filestext+'.txt'
outfilename_ds       = endlocation+'av_ds'+filestext+'.txt'  #'lammpsdiffusion_'+namebase_out+filestext+'_av_ds.txt'
outfilename_gamma    = endlocation+'zimportance'+filestext+'.txt' #'lammpsdiffusion_'+namebase_out+filestext+'_zimportance.txt'
outfilename_sections = endlocation+'sections'+filestext+'.txt' #'lammpsdiffusion_'+namebase_out+filestext+'_sections.txt'
plotname_dx_dy_dz    = endlocation+'dx_dy_dz'+filestext+'.png' #'lammpsdiffusion_'+namebase_out+filestext+'_dx_dy_dz.png' 
outfilename_th       = endlocation+'th'+filestext+'.txt' #'lammpsdiffusion_'+namebase_out+filestext+'_th.txt'
plotname             = endlocation+'all_'+filestext+'.png' #'lammpsdiffusion_'+namebase_out+filestext+'.png'
plotname_gamma       = endlocation+'zimportance'+filestext+'.png' #'lammpsdiffusion_'+namebase_out+filestext+'_zimportance.png'
plotname_SI          = endlocation+'SI'+filestext+'.png' #'lammpsdiffusion_'+namebase_out+filestext+'_SI.png'
plotname_velocity    = endlocation+'velocity'+filestext+'.png' #'lammpsdiffusion_'+namebase_out+filestext+'_velocity.png'
plotname_velocity_SI = endlocation+'velocity_SI'+filestext+'.png' #'lammpsdiffusion_'+namebase_out+filestext+'_velocity_SI.png'
plotname_velocity_sq = endlocation+'velocity_sq'+filestext+'.png' #'lammpsdiffusion_'+namebase_out+filestext+'_velocity_sq.png'
plotname_parallel_orthogonal = endlocation+'par_ort'+filestext+'.png' #'lammpsdiffusion_'+namebase_out+filestext+'_par_ort.png'
plotname_sectioned_average = endlocation+'sections'+filestext+'.png' #'lammpsdiffusion_'+namebase_out+filestext+'_sections.png'
plotname_traj_xy = endlocation+'traj_xy'+filestext+'.png' #'lammpsdiffusion_'+namebase_out+filestext+'_traj_xy.png'
plotname_traj_xz = endlocation+'traj_xz'+filestext+'.png' #'lammpsdiffusion_'+namebase_out+filestext+'_traj_xz.png'
plotname_traj_yz = endlocation+'traj_yz'+filestext+'.png' #'lammpsdiffusion_'+namebase_out+filestext+'_traj_yz.png'
plotname_traj_xt = endlocation+'traj_xt'+filestext+'.png' #'lammpsdiffusion_'+namebase_out+filestext+'_traj_xy.png'
plotname_traj_yt = endlocation+'traj_yt'+filestext+'.png' #'lammpsdiffusion_'+namebase_out+filestext+'_traj_xz.png'
plotname_traj_zt = endlocation+'traj_zt'+filestext+'.png' #'lammpsdiffusion_'+namebase_out+filestext+'_traj_yz.png'
plotname_th_hist = endlocation+'th_hist'+filestext+'.png' #'lammpsdiffusion_'+namebase_out+filestext+'_th_hist.png'
# Setting arrays
# Prepare for sectioning distance data:
time_walks_SI, steps, partition_walks, numberofsamples, len_all, lengths, startpoints = datr.partition_holders_averaged(Nsteps,minlength)
# These are all squared:
# All together:
allRs     = [] # Distance
alldxs    = []
alldys    = []
alldzs    = []
allvs     = [] # Velocity
allvxs    = []
allvys    = []
allvzs    = []

# Averages:
averageR2s  = np.zeros(Nsteps)   # Distances, squared
averagedx2s = np.zeros(Nsteps)
averagedy2s = np.zeros(Nsteps)
averagedz2s = np.zeros(Nsteps)
averageRs  = np.zeros(Nsteps)   # Distances
averagedxs = np.zeros(Nsteps)
averagedys = np.zeros(Nsteps)
averagedzs = np.zeros(Nsteps)
averagevs  = np.zeros(Nsteps)  # Velocities
averagevxs = np.zeros(Nsteps)
averagevys = np.zeros(Nsteps)
averagevzs = np.zeros(Nsteps)
averagev2s  = np.zeros(Nsteps)  # Velocities
averagevx2s = np.zeros(Nsteps)
averagevy2s = np.zeros(Nsteps)
averagevz2s = np.zeros(Nsteps)
averagevparallel2 = np.zeros(Nsteps) # Velocities_squared
averagevparallel  = np.zeros(Nsteps)
averagedparallel  = np.zeros(Nsteps)
averagedparallel2 = np.zeros(Nsteps)
average_counter   = np.zeros(Nsteps)
average_walks     = copy.deepcopy(partition_walks)
average_walks_SI  = copy.deepcopy(partition_walks)
average_counters  = copy.deepcopy(partition_walks)

# For testing:
# Negative elements:
vx_neg = np.zeros(Nsteps)
vy_neg = np.zeros(Nsteps)
vz_neg = np.zeros(Nsteps)
v_all = np.zeros(Nsteps)
vx_neg_fraction = np.zeros(Nsteps)
vy_neg_fraction = np.zeros(Nsteps)
vz_neg_fraction = np.zeros(Nsteps)
# Plot appearance:
averagevparallel_fake = np.zeros(Nsteps)

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
    infilename_free   = inlocation+'seed'+seedstr+'.lammpstrj'
    plotname_dirs     = inlocation+'dxdydzR2_seed'+seedstr+'.png'
    plotname_testsect = inlocation+'testsectioned_seed'+seedstr+'.png'
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
    
    Nin   = 0
    maxzbulk = 0
    pos_inbulk = []
    for i in range(counter):
        thesepos = positions[i]
        z        = thesepos[2]
        absz     = abs(z)
        if absz>maxh: # To match the bulk diffusion with the brush one
            #print('seed',seed,', z>maxh, breaking. z =', z,'. Nin=',Nin)
            break
        else:
            pos_inbulk.append(thesepos)
            if absz>maxzbulk:
                maxzbulk = absz
            Nin+=1
    
    # testh_times:
    for i in range(counter):
        thesepos = positions[i]
        z        = thesepos[2]
        if abs(z)>testh: # If the polymer is in bulk # Measuring the time when is first reaches height h.
            testh_times[seed-1] = i*dt
            break
    
    # I do not take into account that the free bead can enter the polymer grid and then exit again. If that happens, there might be discontinuities or weird kinks in the data (if you look at the graphs...) # I do not have this in my test-dataset. Maybe make the bead lighter and see what happens?
    
    startpos_in = positions[0]
    #print('test:::', startpos_in)
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
        dx   = dist[0] # dx
        dy   = dist[1]
        dz   = dist[2]
        R2   = np.dot(dist,dist)
        dx2  = dx*dx
        dy2  = dy*dy
        dz2  = dz*dz
        gamma = (R2-dz2)/R2
        # Velocity:
        vxi  = vx[i]
        vyi  = vy[i]
        vzi  = vz[i]
        vx2i = vxi*vxi
        vy2i = vyi*vyi
        vz2i = vzi*vzi
        # All together:
        allRs.append(R2)
        alldxs.append(dx2)
        alldys.append(dy2)
        alldzs.append(dz2)
        alltimes.append(i)
        # Averages:
        averageR2s[i] += R2   # Distance squared
        averagedx2s[i]+= dx2
        averagedy2s[i]+= dy2
        averagedz2s[i]+= dz2
        averagedxs[i] += dx   # dx
        averagedys[i] += dy
        averagedzs[i] += dz
        averagevs[i]  += np.sqrt(vx2i + vy2i + vz2i)  # Velocity
        averagevxs[i] += vxi
        averagevys[i] += vyi
        averagevzs[i] += vzi
        averagev2s[i] += vx2i + vy2i + vz2i           # Velocity squared
        averagevx2s[i]+= vx2i
        averagevy2s[i]+= vy2i
        averagevz2s[i]+= vz2i
        averagevparallel2[i] += vx2i+vy2i
        averagevparallel[i]  += np.sqrt(vx2i+vy2i)
        averagedparallel2[i] += dx2+dy2 # Distance
        averagedparallel[i]  += np.sqrt(dx2+dy2) # Distance
        averagevparallel_fake[i] += vxi+vyi
        average_counter[i]  +=1
        # For testing purposes:
        if vxi<0:
            vx_neg[i] += 1
        if vyi<0:
            vy_neg[i] += 1
        if vzi<0:
            vz_neg[i] += 1
        v_all[i] += 1
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
    for i in range(Npartitions):
        this_part    = []
        these_steps  = []
        part_start   = startpoints[i]
        len_thiswalk = lengths[i]
        part_end     = part_start+len_thiswalk
        if part_start>(Nin-1):              # Do not start this section if the bead has left the brush.
            break
        rstart = pos_inbulk[part_start]
        for j in range(len_thiswalk):
            #print('i:',i,', j:',j)
            if part_start+j>(Nin-1):
                break
            rthis =  pos_inbulk[part_start+j] # Could have just used positions, but...
            drvec = rthis - rstart
            dr2   = np.dot(drvec,drvec)
            '''
            print('part_start:',part_start, '; part_start+j:', part_start+j, ';   j:',j)
            print('rthis:',rthis)
            print('rstart:',rstart)
            print('drvec:', drvec)
            print('dr2:', dr2)
            '''
            average_walks[i][j]    +=dr2#= average_walks[i][j] + dr2 # Notation? [i,j] or [i][j] # Ugh, and should I have one of these for R2, dx2, dy2 and dz2?
            average_counters[i][j] +=1#= average_counters[i][j] + 1
            this_part.append(dr2)
            these_steps.append(j)
        part_walks.append(this_part)
        part_steps.append(these_steps)
    sections_walk.append(part_walks)
    sections_steps.append(part_steps)
    
    #print('Nin:', Nin, '; shape(part_walks):', shape(part_walks), len(part_walks))
    if test_sectioned==True and len(part_walks)==Npartitions and (seed==29 or seed==47 or seed==53 or seed==59 or seed==83 or seed==103):
        plt.figure(figsize=(6,5))
        for i in range(Npartitions):
            print('i:',i)
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

for i in range(1,Nsteps):    
    #print('Step',i,', neg vz:', vz_neg[i], ' out of', v_all[i])  
    vx_neg_fraction[i] = vx_neg[i]/v_all[i]
    vy_neg_fraction[i] = vy_neg[i]/v_all[i]
    vz_neg_fraction[i] = vz_neg[i]/v_all[i]

'''
plt.figure()
plt.plot(times_single,vx_neg, label=r'Neg $v_x$')
plt.plot(times_single,vy_neg, label=r'Neg $v_y$')
plt.plot(times_single,vz_neg, label=r'Neg $v_z$')
plt.plot(times_single,v_all, label=r'No of $v$s')
plt.xlabel('Time step')
plt.ylabel('Occurences')
plt.title(r'How many negative $v$-components?')
plt.legend(loc='upper left')
plt.tight_layout()
plt.show()


plt.figure()
plt.plot(times_single,vx_neg_fraction, label=r'Neg $v_x$')
plt.plot(times_single,vy_neg_fraction, label=r'Neg $v_y$')
plt.plot(times_single,vz_neg_fraction, label=r'Neg $v_z$')
plt.xlabel('Time step')
plt.ylabel('Fraction')
plt.title(r'Fraction of negative $v$-components')
plt.legend(loc='upper left')
plt.tight_layout()
plt.show()
'''

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
# Know the velocity at t=0:
vxi = vx[0]
vyi = vy[0]
vzi = vz[0]
averagevs[0] = np.sqrt(vxi*vxi + vyi*vyi + vzi*vzi)
averagevxs[0]= vxi
averagevys[0]= vyi
averagevzs[0]= vzi
averagevparallel[0] = np.sqrt(vxi*vxi+vyi*vyi)
for i in range(Nsteps):
    counter = average_counter[i]
    if counter!=0:
        # Distance
        averageRs[i]/=counter
        averagedxs[i]/=counter
        averagedys[i]/=counter
        averagedzs[i]/=counter
        averagedparallel[i]/= counter
        # Distance squared
        averageR2s[i]/=counter
        averagedx2s[i]/=counter
        averagedy2s[i]/=counter
        averagedz2s[i]/=counter
        averagedparallel2[i]/= counter
        # Velocity
        averagevs[i]/=counter
        averagevxs[i]/=counter
        averagevys[i]/=counter
        averagevzs[i]/=counter
        averagevparallel[i]/=counter
        averagevparallel_fake[i]/=counter
        # Velocity squared
        averagev2s[i]/=counter
        averagevx2s[i]/=counter
        averagevy2s[i]/=counter
        averagevz2s[i]/=counter
        averagevparallel2[i]/=counter

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
coeffs_poly, covs = polyfit(times_single, averageR2s, 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
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


coeffs_poly, covs = polyfit(times_single_real, averageR2s_SI, 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
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
    outfile_ds.write('%i %.16e %.16e %.16e %.16e %.16e %.16e\n' % (times_single[i], times_single_real[i], averageR2s_SI[i], averagedx2s_SI[i], averagedy2s_SI[i], averagedz2s_SI[i], averagedparallel2_SI[i]))
outfile_ds.close()

times_last = []
xpos = []
ypos = []
zpos = []

for i in range(Nin):
    pos= positions[i]
    xpos.append(pos[0])
    ypos.append(pos[1])
    zpos.append(pos[2])
    times_last.append(i)

'''
plt.figure(figsize=(6,5))
plt.plot(alltimes, allRs, ',', label='Data, bulk')
plt.plot(times_single, averageR2s, ',', label='Average, bulk')
plt.plot(alltimes, fit_poly, '--', label='Fit, data, bulk')
plt.plot(times_single, fit_poly_av, '--', label='Fit, average, bulk')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title('RMSD in bulk, d = %i nm' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.show()
'''

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
plt.title('Trajectory in zy-plane')
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
plt.plot(times_last, xpos, '.')
plt.xlabel(r't')
plt.ylabel(r'x')
plt.title('x(t)')
plt.tight_layout()
plt.savefig(plotname_traj_xt)

plt.figure(figsize=(6,5))
plt.plot(times_last, ypos, '.')
plt.xlabel(r't')
plt.ylabel(r'y')
plt.title('y(t)')
plt.tight_layout()
plt.savefig(plotname_traj_yt)

plt.figure(figsize=(6,5))
plt.plot(times_last, zpos, '.')
plt.xlabel(r't')
plt.ylabel(r'z')
plt.title('z(t)')
plt.tight_layout()
plt.savefig(plotname_traj_zt)



'''
average_walks_SI = copy.copy(partition_walks)
average_counters = copy.copy(partition_walks)
time_walks_SI    = copy.copy(steps)
'''

outfile_sections = open(outfilename_sections, 'w')
for i in range(Npartitions):
    outfile_sections.write('Section %i\n' % i) # When reading from the file afterwards, I can extract the section number by if words[0]=='Section' and use that to divide the data into sections
    for j in range(lengths[i]):
        average_walks_SI[i][j] = average_walks[i][j]*unitlength
        time_walks_SI[i][j]    = steps[i][j]*timestepsize
        outfile_sections.write('%.16f %16.f\n' % (time_walks_SI[i][j],average_walks_SI[i][j]))
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
outfile_sections.close()

## Making figures.

plt.figure(figsize=(6,5))
plt.plot(alltimes, allRs, ',', label='Data, bulk')
plt.plot(times_single, averageR2s, ',', label='Average, bulk')
plt.plot(alltimes, fit_poly, '--', label='Fit, data, bulk')
plt.plot(times_single, fit_poly_av, '--', label='Fit, average, bulk')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title('RMSD in bulk, d = %i nm' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.plot(times_single, averageR2s, label=r'$<dR^2>$') #label=r'$\braket{dR^2}$')
plt.plot(times_single, averagedx2s, label=r'$<dx^2>$') #label=r'$\braket{dx^2}$')
plt.plot(times_single, averagedy2s, label=r'$<dy^2>$')#label=r'$\braket{dy^2}$')
plt.plot(times_single, averagedz2s, label=r'$<dz^2>$')#label=r'$\braket{dz^2}$')
plt.plot(times_single, averagedparallel2,label=r'$<dx^2+dy^2>$') #label=r'$\braket{dx^2+dy^2}')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title('Averaged RMSD in bulk, d = %i nm' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_parallel_orthogonal)

plt.figure(figsize=(6,5))
plt.plot(times_single, averagedxs, label=r'$<dx>$')
plt.plot(times_single, averagedys, label=r'$<dy>$')
plt.plot(times_single, averagedzs, label=r'$<dz>$')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance [in unit length]')
plt.title(r'Averaged RMSD in bulk, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_dx_dy_dz)

plt.figure(figsize=(6,5))
plt.plot(plotgammaagainst, sorted_gamma_avgs)
plt.xlabel(r'Different runs')
plt.ylabel(r'$\gamma = <\frac{R^2-dz^2}{R^2}>$')
plt.title(r'RMSD in bulk, $\gamma$; d = %i nm' % spacing)
plt.tight_layout()
plt.savefig(plotname_gamma)

plt.figure(figsize=(6,5))
plt.plot(times_single_real, averageR2s_SI, ',', label='Average, bulk')
plt.plot(times_single_real, fit_poly_SI, '--', label='Fit, average, bulk')
plt.xlabel(r'Time [s]')
plt.ylabel(r'Distance$^2$ [m]')
plt.title('RMSD in bulk, d = %i nm, SI' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.savefig(plotname_SI)

plt.figure()
plt.plot(times_single, averagevs, label=r'$<v>$')
plt.plot(times_single, averagevxs, label=r'$<v_x>$')
plt.plot(times_single, averagevys, label=r'$<v_y>$')
plt.plot(times_single, averagevzs, label=r'$<v_z>$')
plt.plot(times_single, averagevparallel, label=r'$<v_\parallel>$')
plt.plot(times_single, averagevparallel_fake, label=r'$<v_x+v_y>$')
plt.xlabel(r'Step number')
plt.ylabel(r'Velocity [in unit length/time]')
plt.title('Averaged velocity in bulk, system size by d = %i nm' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_velocity)

plt.figure()
plt.plot(times_single_real, averagevs_SI, label=r'$<v>$')
plt.plot(times_single_real, averagevxs_SI, label=r'$<v_x>$')
plt.plot(times_single_real, averagevys_SI, label=r'$<v_y>$')
plt.plot(times_single_real, averagevzs_SI, label=r'$<v_z>$')
plt.plot(times_single_real, averagevparallel_SI, label=r'$<v_\parallel>$')
plt.xlabel(r'Time [s]')
plt.ylabel(r'Velocity [m/s]')
plt.title('Averaged velocity in bulk, system size by d = %i nm, SI' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.savefig(plotname_velocity_SI)

plt.figure()
plt.plot(times_single, averagev2s, label=r'$<v^2>$')
plt.plot(times_single, averagevx2s, label=r'$<v_x^2>$')
plt.plot(times_single, averagevy2s, label=r'$<v_y^2>$')
plt.plot(times_single, averagevz2s, label=r'$<v_z^2>$')
plt.plot(times_single, averagevparallel2, label=r'$<v_x^2+vy^2>$')
plt.xlabel(r'Step number')
plt.ylabel(r'Velocity [(length/time)$^2$]')
plt.title('Averaged velocity in bulk, d = %i nm' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_velocity_sq)


plt.figure()
for i in range(Npartitions):
    plt.plot(time_walks_SI[i][:],average_walks_SI[i][:], label='Start at step %i' % startpoints[i])
plt.xlabel(r'Time [s]')
plt.ylabel(r'Distance$^2$ [m]')
plt.title('RMSD in bulk, system size by d = %i nm, SI, sectioning' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.savefig(plotname_sectioned_average)


print('averagedxs[-1]:',averagedxs[-1])
print('averagedys[-1]:',averagedys[-1])
print('averagedzs[-1]:',averagedzs[-1])

print('mean(averagedxs):',np.mean(averagedxs))
print('mean(averagedys):',np.mean(averagedys))
print('mean(averagedzs):',np.mean(averagedzs))

## Distribution of times when it reaches height testh:
#testh_times[seed]
# Binning:
maxtime_testh = max(testh_times)
mintime_testh = min(testh_times)
binlength     = (maxtime_testh-mintime_testh)/float(Nbins)
bins_ht       = np.zeros(Nbins)
bincenters    = np.zeros(Nbins)
centerstart   = mintime_testh+0.5*binlength


# Finding bin values + writing to file
outfile_th = open(outfilename_th, 'w')
for i in range(Nseeds):
    val = testh_times[i]
    outfile_th.write('%.16f\n' % val)
    for j in range(Nbins):
        if val>(mintime_testh+j*binlength) and val<(mintime_testh+(j+1)*binlength):
            bins_ht[j] += 1
            continue

Nbinned = np.sum(bins_ht) # Use this or use Nseeds?
for i in range(Nbins):
    bins_ht[i] /= Nseeds#Nbinned
    bincenters[i] = centerstart +i*binlength

plt.figure()
plt.plot(bincenters, bins_ht)
plt.xlabel(r'$t_h$ [s]')
plt.ylabel(r'Hits/$N_{sims}$')
plt.title('$t_h$ in system size by d = %i nm' % spacing)
plt.tight_layout()
plt.savefig(plotname_th_hist)


