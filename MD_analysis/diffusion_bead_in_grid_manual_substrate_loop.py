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

# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
spacing = 10
psigma  = 2.5
density = 0.238732414637843 # Yields mass 1 for bead of radius 1 nm
#pmass   = 1.5
print('spacing:', spacing)
print('psigma:', psigma)
# I need to set the file name in an easier way, but for now I just use this:  ## Might want to add a loop too, if I have more files...

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them
plotseed = 0
plotdirs = False
test_sectioned = False
#seeds  = [23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]
confignrs    = np.arange(1,22) 
Nseeds       = len(confignrs)      # So that I don't have to change that much
Nsteps       = 20001
writeevery   = 10                  # I'm writing to file every this many time steps
unitlength   = 1e-9
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime*writeevery
Npartitions  = 5 # For extracting more walks from one file (but is it really such a random walk here...?)
minlength    = int(floor(Nsteps/Npartitions)) # For sectioning the data
print('timestepsize:', timestepsize)
Npartitions = 5 # For extracting more walks from one file (but is it really such a random walk here...?)
## Weird cutoff (bead):
#namebase    = '_quadr_M9N101_ljunits_spacing%i_Langevin_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_theta0is180_pmass1.5_sect_placeexact_ljcut1p122' %spacing
#folderbase  = 'Part_in_chgr_subst_all_quadr_M9N101_ljunits_Langevin_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_theta0is180_pmass1.5_sect_placeexact_ljcut1p122'
# Usual cutoff (bead):
foldername  = 'Quadr_M9N101_ljunits_Langevin_scaled_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_ljcut1p122'
endlocation = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Brush/'+foldername+'/Spacing'+str(spacing)+'/Sigma_bead_'+str(psigma)+'/'
namebase_start = '_quadr_M9N101_ljunits_'
folderbase_mid = 'Langevin_scaled_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_pmass1.5'
#folderbase_end = '_ljcut1p122'
namebase = namebase_start+'spacing%i_' % spacing + folderbase_mid + '_psigma' +str(psigma)+'_sect_placeexact_ljcut1p122'
#namebase = 'spacing%i_Langevin_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_pmass1.5_psigma' % spacing +str(psigma)+'_pstdcutoff_sect_placeexact_ljcut1p122'
#folderbase = 'Part_in_chgr_subst'+namebase_start+folderbase_mid+folderbase_end

#endlocation       = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Brush/'+foldername+'/Spacing'+str(spacing)+'/Sigma_bead_'+str(psigma)
outfilename          = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'.txt'
outfilename_ds       = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'_av_ds.txt'
outfilename_gamma    = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'_zimportance'+'.txt'
outfilename_sections = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'_sections'+'.txt'
plotname             = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'.png'
plotname_all         = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'_all.png'
plotname_gamma       = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'_zimportance'+'.png'
plotname_SI          = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'_SI.png'
plotname_parallel_orthogonal = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'_par_ort.png' 
plotname_parallel    = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'_par.png' 
plotname_orthogonal  = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'_ort.png'
plotname_short_all   = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'_short_all.png'
plotname_velocity    = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'_velocity.png'
plotname_velocity_SI = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'_velocity_SI.png'
plotname_sectioned_average = endlocation+'lammpsdiffusion_'+namebase+'_sections.png'


## Setting arrays
# Prepare for sectioning distance data:
steps, partition_walks, numberofsamples, len_all, lengths, startpoints = datr.partition_holders_averaged(Nsteps,minlength)
# These are all squared:
# All together:
allRs     = []
alldxs    = []
alldys    = []
alldzs    = []
# Averages:
averageRs = np.zeros(Nsteps)   # Distances
averagedxs = np.zeros(Nsteps)
averagedys = np.zeros(Nsteps)
averagedzs = np.zeros(Nsteps)
averagevs  = np.zeros(Nsteps)  # Velocities
averagevxs = np.zeros(Nsteps)
averagevys = np.zeros(Nsteps)
averagevzs = np.zeros(Nsteps)
averagedparallel = np.zeros(Nsteps) # Distances
average_counter  = np.zeros(Nsteps)
average_walks    = copy.copy(partition_walks)
average_walks_SI = copy.copy(partition_walks)
average_counters = copy.copy(partition_walks) # Okay, this name is confusing.
time_walks_SI    = copy.copy(steps)

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


for confignr in confignrs:
    print('On config number:', confignr)
    infilename_all  = endlocation+'all_density'+str(density)+'_confignr'+str(confignr)+'.lammpstrj'
    infilename_free = endlocation+'freeatom_density'+str(density)+'_confignr'+str(confignr)+'.lammpstrj'
    plotname_dirs   = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'_dxdydzR2_seed'+str(confignr)+'.png'
    plotname_testsect = endlocation+'lammpsdiffusion_qdrgr_'+namebase+'_testsectioned_seed'+str(confignr)+'.png'
    
    #print('infilename_all:',infilename_all)
    
    # Read in:
    #### Automatic part
    ## Find the extent of the polymers: Max z-coord of beads in the chains
    infile_all = open(infilename_all, "r")
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
                if z>maxz:
                    maxz = z
            i+=1
    extent_polymers = maxz
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
    freeatom_positions = [] # Fill with xu,yu,zu-coordinates. One position for each time step. # We have the same time step for the whole of the simulation. Can bother with time later
    times              = [] # Not useful anymore?
    
    # For now: Only find the largest z-position among the beads in the chain. # Otherwise I must read off the atom number (will do that with the free bead anyways.)
    maxz = -1000 # No z-position is this small.
    
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
            # Find properties
            # Order:  id  type mol ux  uy  uz  vx  vy   vz
            #         [0] [1]  [2] [3] [4] [5] [6] [7]  [8]
            ind      = int(words[0])-1 # Atom ids go from zero to N-1.
            #atomtype = int(words[1]) 
            #molID    = int(words[2])
            x        = float(words[3])
            y        = float(words[4])
            z        = float(words[5])
            freeatom_positions.append(np.array([x,y,z]))
            vx.append(float(words[6]))                       # Check the output format!!!!
            vy.append(float(words[7]))
            vz.append(float(words[8]))
            counter+=1
            i+=1
    infile_free.close()
    dt = (times[1]-times[0])*timestepsize # This might be handy
    #print('dt:', dt)
    #print('times[1]-times[0]:',times[1]-times[0])
    times_single = np.arange(Nsteps)#*dt
    times_single_real = np.arange(Nsteps)*dt
    
    time_end = time.process_time()
    
    pos_inpolymer = []
    
    #before_in
    
    # I do not take into account that the free bead can enter the polymer grid and then exit again. If that happens, there might be discontinuities or weird kinks in the data (if you look at the graphs...) # I do not have this in my test-dataset. Maybe make the bead lighter and see what happens?
    
    Nin   = 0
    maxzpol = 0
    for i in range(counter):
        thesepos = freeatom_positions[i]
        z        = thesepos[2]
        if z>extent_polymers: # If the polymer is in bulk # We don't want it to go back and forth between brush and bulk # That will cause discontinuities in our data
            break
        else:
            pos_inpolymer.append(thesepos)
            if z>maxzpol:
                maxzpol = z
            Nin+=1
    
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
        averageRs[i] +=R2   # Distance
        averagedxs[i]+=dx2
        averagedys[i]+=dy2
        averagedzs[i]+=dz2
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
    smallest_part_walk = int(math.floor(Nin/Npartitions))
    for i in range(Npartitions):
        this_part    = []
        these_steps  = []
        part_start   = i*smallest_part_walk
        len_thiswalk = Nin-part_start       # Do I need this?
        #print('part_start:', part_start)
        #print('len(pos_inpolymer):', len(pos_inpolymer))
        rstart = pos_inpolymer[part_start]
        for j in range(len_thiswalk):
            rthis =  pos_inpolymer[part_start+j]
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
    time_beforepartition = time.process_time()
    partition_walks = datr.partition_dist_averaged(positions, partition_walks, Nsteps, minlength) # partition_walks is being updated inside of the function
    time_afterpartition = time.process_time()
    print('Time, partitioning:',time_afterpartition-time_beforepartition)

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

# Calculating averages:
Ninbrush = Nsteps # Default, in case it does not exit the brush
for i in range(1,Nsteps):
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
    else:
       Ninbrush = i-1
       break

times_single = times_single[0:Ninbrush]
times_single_real = times_single_real[0:Ninbrush]
averageRs  = averageRs[0:Ninbrush]
averagedxs = averagedxs[0:Ninbrush]
averagedys = averagedys[0:Ninbrush]
averagedzs = averagedzs[0:Ninbrush]
averagedparallel = averagedparallel[0:Ninbrush]

print('len(averageRs):',len(averageRs))
print('Ninbrush:',Ninbrush)
print('Nsteps:', Nsteps)

# Opening file
outfile_gamma = open(outfilename_gamma, 'w')
outfile_gamma.write('This is an attempt to find the relative contribution of dz to R. The smaller the value in the second column, the more important dz is.\nOrder: <(R^2(n)-dz^2(n))/R^2(n)>_n, slope a of line fit, rmsd R^2(n) line fit, Nin, confignr\n')
for i in range(Nseeds):
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


print('----------------------------------')
print(' averageRs_SI[2]:', averageRs_SI[2])

coeffs_poly, covs = polyfit(times_single_real, averageRs_SI, 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a_poly_SI = coeffs_poly[0]
b_poly_SI = coeffs_poly[1]
rms_D_poly_SI = np.sqrt(covs[0,0])/6.
rms_b_poly_SI = np.sqrt(covs[1,1])
D_poly_SI = a_poly_SI/6.

print(' averageRs_SI[2]:', averageRs_SI[2])

print('----------------------------------')
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
    outfile_ds.write('%i %.16e %.16e %.16e %.16e %.16e %.16e\n' % (times_single[i], times_single_real[i], averageRs_SI[i], averagedxs_SI[i], averagedys_SI[i], averagedzs_SI[i], averagedparallel_SI[i]))
outfile_ds.close()


outfile_sections = open(outfilename_sections, 'w')
for i in range(Npartitions):
    outfile_sections.write('Section %i\n' % i) # When reading from the file afterwards, I can extract the section number by if words[0]=='Section' and use that to divide the data into sections
    for j in range(lengths[i]):
        average_walks_SI[i][j] = steps[i][j]*timestepsize
        time_walks_SI[i][j] = average_walks[i][j]*unitlength
        outfile_sections.write('%.16f %16.f\n' % (time_walks_SI[i][j],average_walks_SI[i][j]))
outfile_sections.close()


'''
maxind = 10000
plt.figure(figsize=(6,5))
plt.plot(times_single_real, averageRs_SI, label=r'$<R^2>$')
plt.plot(times_single_real, averagedzs_SI, label=r'$<dz^2>$')
plt.plot(times_single_real, averagedparallel_SI, label=r'$<dx^2+dy^2>$')
plt.xlabel(r'Index (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in bulk, d = %i nm, $\sigma_b=%.2f$, SI' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.axis([0, times_single_real[maxind], 0, max(averageRs_SI[0:maxind])])
plt.show()
'''
# To determine the range
xmax_plot = 10000
plt.figure(figsize=(6,5))
plt.plot(times_single, averageRs, label=r'$<R^2>$')
plt.plot(times_single, averagedzs, label=r'$<dz^2>$')
plt.plot(times_single, averagedparallel, label=r'$<dx^2+dy^2>$')
plt.xlabel(r'Index (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in bulk, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.axis([0, xmax_plot, 0, max(averageRs[0:xmax_plot])])
plt.savefig(plotname_short_all)

plt.figure(figsize=(6,5))
plt.plot(alltimes, allRs, ',', label='Data, brush')
plt.plot(times_single, fit_poly_av, '--', label='Fit, brush')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in bulk, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.plot(alltimes, allRs, ',', label='Data, brush')
plt.plot(times_single, averageRs, ',', label='Average, brush')
plt.plot(alltimes, fit_poly, '--', label='Fit, data, brush')
plt.plot(times_single, fit_poly_av, '--', label='Fit, average, brush')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in bulk, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_all)

plt.figure(figsize=(6,5))
plt.plot(times_single, averagedzs, ',', label=r'dz$^2$')
plt.plot(times_single, averagedparallel, ',', label=r'dx$^2$+dy$^2$')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'Averaged RMSD in bulk, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_parallel_orthogonal)

plt.figure(figsize=(6,5))
plt.plot(times_single, averagedzs, ',')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'$<dz^2>$ in bulk, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.savefig(plotname_orthogonal)

plt.figure(figsize=(6,5))
plt.plot(times_single, averagedparallel, ',')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'$<dx^2+dy^2>$ in bulk, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
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
plt.plot(times_single_real, averageRs_SI, ',', label='Average, brush')
plt.plot(times_single_real, fit_poly_SI, '--', label='Fit, average, brush')
plt.xlabel(r'Time (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in bulk, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
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
plt.plot(times_single_real, averagevxs_SI, label=r'vx')
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
    plt.plot(time_walks_SI[i][:],average_walks_SI[i][:], label='Start at step %i' % startpoints[i])
plt.xlabel(r'Time [s]')
plt.ylabel(r'Distance$^2$ [m]')
plt.title('RMSD in bulk, system size by d = %i nm, SI, sectioning' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_sectioned_average)



print('counters:',average_counter)

print('spacing:', spacing)
print('psigma:', psigma)
