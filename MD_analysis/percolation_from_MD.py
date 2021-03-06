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

# Function for curve_fit
'''
def costheta_exponential(s,P):
    return np.exp(-s/P)
'''

# Checking when a (decreasing) quantity dips below a given value
def crude_pcapproach_readoffgraph(p, Pi, threshold, confidence_range):
    Np = len(p)
    passed_start = 0
    passed_thr   = 0
    deltaPi      = 0.5*confidence_range
    endrangeset  = False
    startrangeset = False
    for i in range(Np-1):
        if Pi[i]<=(threshold+deltaPi) and passed_start==0: # Should maybe use 0.5*1/e or something instead of 0.2 for confidence_range?
            startrange    = p[i]
            passed_start  = 1
            startrangeset = True
        if Pi[i]<=threshold and passed_thr==0:   # The moment Pi dips below the threshold, take the corresponding p to be pc
            pc = p[i]
            slope = (Pi[i+1]-Pi[i-1])/(p[i+1]-p[i-1])
        if Pi[i]<=(threshold-deltaPi):
            endrange = p[i]
            endrangeset = True
            break
    # These are here because we might not have chosen an appropriate range
    if endrangeset==False:
        endrange = p[-1]
    if startrangeset==False:
        startrange = p[0]
    rangewidth = endrange-startrange
    return pc, slope, startrange, endrange, rangewidth, startrangeset, endrangeset # I don't need to print this much...

def partition_holders(Nthr,Nstepsrandom,minlength):
    interval     = minlength                           # Spacing between start points
    end_interval = Nstepsrandom-interval+1             # +1 because of the indexing and how arange works (but will check and see later)
    startpoints  = np.arange(0,end_interval, interval) # Points to start the calculation of the rmsd.
    numberofsamples = len(startpoints)                 # Number of such walks. Will be the number of lists
    len_all = 0
    partitioned_walks_holder = []
    steps_holder             = []
    lengths                  = []
    for i in range(numberofsamples):
        length = Nstepsrandom-i*interval
        len_all += length
        lengths.append(length)
        partitioned_walks_holder.append(np.zeros((Nthr,length))) # Should hold in general
        steps_holder.append(np.arange(0,length)) # Rewrite this slightly... Use arange
    return steps_holder, partitioned_walks_holder, numberofsamples, len_all, lengths

def partition_rw(thesepositions,partition_walks, Nstepsrandom, minlength, thisthr, Nthr): # Or should I use something like Nsections instead?
    interval     = minlength                           # Spacing between start points
    end_interval = Nstepsrandom-interval+1             # +1 because of the indexing and how arange works (but will check and see later)
    startpoints  = np.arange(0,end_interval, interval) # Points to start the calculation of the rmsd.
    numberofsamples = len(startpoints)                 # Number of such walks. Will be the number of lists
    for i in range(numberofsamples): # Looping over the number of partitioned walks
        startindex = startpoints[i]
        startpoint = thesepositions[startindex]
        counter    = 1
        length     = Nstepsrandom-i*interval
        these_rmsd = np.zeros((Nthr,length))
        this_index = startindex+counter
        while this_index<Nstepsrandom:
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
    
# Differentiation
def diff_by_middlepoint(h,f):
    # Differentiating using the
    # symmetric difference quotient:
    #   f' = (f(x+2h)-f(x))/2h
    # The error is proportional with h
    # The advantage is that we can use this blindly even when
    # there is an uneven spacing
    
    N = len(f)-2
    df = zeros(N)
    hout = h[1:N]
    
    for i in range(1,N):
        df[i]= (f[i+1]-f[i-1])/(h[i+1]-h[i-1])
    return hout, df


def diff_by_secant(h,f):
    # Differentiating using the secant method
    #   f' = (f(x+h)-f(x))/h
    # The error is proportional with h
    # The advantage is that we can use this blindly even when
    # there is an uneven spacing

    N = len(f)-1
    df = zeros(N)
    hout = h[0:N]

    for i in range(N):
        df[i]= (f[i+1]-f[i])/(h[i+1]-h[i])
    return hout, df

def element_difference(h,f):
    # The difference between adjacent elements in the input array

    N = len(f)-1
    df = zeros(N)
    hout = h[1:N+1]

    for i in range(N):
        df[i]= f[i+1]-f[i]
    return hout, df

def av_and_rms(f):
    av  = mean(f)
    N   = len(f)
    rms = 0
    for i in range(len(f)):
        rms += (av-f[i])**2
    rms = np.sqrt(rms/N)
    return av, rms

start_time = time.process_time()

### To make data into binary:
# Threshold:
thr   = 0.4 # Set any number here
Nthr  = 101
thrs  = np.linspace(20,40,Nthr) # Or something
#thrs  = np.array([2.1, 15.0])
Nthr  = len(thrs)                  # Safeguard

# 
minlength = 200  # Determines how small walks we should section the walk into

# Pi
# Pi_av, percolation probability vs threshold:
Pi_av   = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
Pi_av_x = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
Pi_av_y = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
Pi_av_z = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
Pi_av_one = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds

# Pi_rms, rms of Pi_av:
Pi_rms   = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
Pi_rms_x = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
Pi_rms_y = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
Pi_rms_z = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
Pi_rms_one = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds

# P
# P_av, percolation probability vs threshold:
P_av   = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
P_av_x = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
P_av_y = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
P_av_z = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
P_av_one = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds

# P_rms, rms of P_av:
P_rms   = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
P_rms_x = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
P_rms_y = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
P_rms_z = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
P_rms_one = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds

counter = np.zeros(Nthr)

# Bools for plotting
issphere        = False # For plotting if we have a spherical pore
plotconfig      = True
plotbool        = False


### Grid setting
zbuffer     = 10 # Don't know if the 'buffer' is big enough
Nvoxels     = 50

### Input file parameters
Kangle      = 20
Kbond       = 200
T           = 310
M           = 9
N           = 101 ###
ljdebye     = 1.042
epsilon     = ljdebye
sigma       = 1
ljcutoff    = 1.12246204830937
debyecutoff = 3
factor      = 0.05#250
#Kbond       = 2000#Kangle*factor
#Kangle      = Kbond*factor
charge      = -1
spacing     = 40
gridspacing = spacing
K           = Kangle # Because we used this notation earlier, but using Kangle is less confusing
wallenergy  = 1.042
dielectric  = 50 #1 2 10 100

### Input and output file names

''' # Testing; <ree2> is correct for completely straight chains.
N  = 100 # Number of bond vectors = number of units - 1 # Per chain
M  = 9   # Number of chains
L2 = 3   # Determining the shape of the L1xL2 polymer grid on the surface
gridspacing = 40 # The spacing in our grid
Nelem   = M*(N+1)
N       = 101
infilename   = 'chaingrids_totallystraight_N%i_Nchains%i_Ly%i_twofixed_test.lammpstrj'  % (Nelem,M,L2)    # Actual run
outfilename  = 'voxelmap_chaingrid_quadratic_M%iN%i_totallystraight_test.txt' % (M,N)
#'''

#	     #  chaingrid_quadratic_M9N101_Langevin_Kangle100_Kbond2000_factor0.05_T310_theta0is180_twofirst_are_fixed
#gridspacing  = 40
'''   # This one is for varying the factor.
infilename   = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,Kangle,Kbond,factor,T)
outfilename  = 'voxelmap_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,Kangle,Kbond,factor,T)
#'''

'''   # Changing the dielectric constant
#	     # 'chaingrid_quadratic_M9N101_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_dielectric10_T310_theta0is180_twofirst_are_fixed.lammpstrj'
#            # 'chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj'
#            # 'chaingrid_quadratic_M9N101_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_dielectric100_T310_theta0is180_twofirst_are_fixed.lammpstrj'
infilename   = 'chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
outfilename  = 'voxelmap_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
#'''


''' # With wall potential:
infilename   = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
% (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
outfilename  = 'voxelmap_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
#'''

# Output names for code testing:
#'''
#infilename   = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,spacing,Kangle,Kbond,charge,T)
#infilename      = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)

infilename_base = 'voxelmap_test_short'#_TESTII'#'spherepore_three'#'equal_distances_four'#'halving_distances'#'voxelmap_test_short'
#'''
# Varying the grid spacing # THIS IS NOW THE STANDARD FILE NAMES.
'''
infilename   = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,spacing,Kangle,Kbond,charge,T)
outfilename  = 'voxelmap_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,Kangle,Kbond,charge,T)
#'''

# Varying the charge:
#chaingrid_quadratic_M9N101_Langevin_Kangle${Kangle}_Kbond${Kbond}_charge${charge}_T$T_theta0is180_twofirst_are_fixed.lammpstrj
''' # String appending method
#               chaingrid_quadratic_M9N101_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge0_T310_theta0is180_twofirst_are_fixed
infilename   = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % T     # Actual run
outfilename  = 'voxelmap_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
#'''

''' # String appending method
infilename   = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % T     # Actual run
outfilename  = 'voxelmap_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
#'''

''' # Sprintf method
infilename   = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,Kangle,Kbond,factor,T)     # Actual run
outfilename  = 'voxelmap_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,Kangle,Kbond,factor,T)
#'''

'''
infilename   = 'chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.lammpstrj' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)     # Actual run
outfilename  = 'voxelmap_chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.txt' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)
'''
### More on names:
name_end                 = '_thr%ito%i' % (thrs[0],thrs[-1])
foldername               = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Voxelmatrices/Frames_'+infilename_base+'/'
plotfoldername           = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Voxelmatrices/Frames_'+infilename_base+'/Plots/'
distfoldername           = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Voxelmatrices/Frames_'+infilename_base+'/Distmatrices/'
percfoldername           = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Voxelmatrices/Frames_'+infilename_base+'/Percolation/'
percfoldername_totalbase = percfoldername + infilename_base # Do I actually need this?
infilename_totalbase     = foldername + infilename_base
plotname_totalbase       = plotfoldername + infilename_base
distfilename_totalbase   = distfoldername + infilename_base
name_end                 = name_end+'_TESTING_RANDOM_WALK'
outfilename_percdata_Pi  = percfoldername_totalbase + '_percolation_data'+name_end+'_Pi.txt'
outfilename_percdata_P   = percfoldername_totalbase + '_percolation_data'+name_end+'_P.txt'
outfilename_percdata_pc  = percfoldername_totalbase + '_percolation_data'+name_end+'_pc.txt'
outfilename_randomwalk   = percfoldername_totalbase + '_voxellated'+name_end+'_random_walk.txt'
outfilename_randomwalk_partiions = percfoldername_totalbase + '_voxellated'+name_end+'_random_walk_partitions.txt'
plotname_percolation_alldirs = percfoldername + infilename_base + '_percolation_total'+name_end
plotname_percolation_xdir    = percfoldername + infilename_base + '_percolation_xdir'+name_end
plotname_percolation_ydir    = percfoldername + infilename_base + '_percolation_ydir'+name_end
plotname_percolation_zdir    = percfoldername + infilename_base + '_percolation_zdir'+name_end
plotname_percolation_everyonetogether = percfoldername + infilename_base + '_percolation_everythingtogether'+name_end
plotname_randomwalk_entire   = percfoldername + infilename_base + '_rw_entire'+name_end
plotname_randomwalk_test     = percfoldername + infilename_base + '_test'+name_end
plotname_randomwalk_walks    = percfoldername + infilename_base + '_rws'+name_end
plotname_randomwalk_slope    = percfoldername + infilename_base + '_rwsslope'+name_end
#
plotname_percolation_P_alldirs = percfoldername + infilename_base + '_P_total'+name_end
plotname_percolation_P_xdir    = percfoldername + infilename_base + '_P_xdir'+name_end
plotname_percolation_P_ydir    = percfoldername + infilename_base + '_P_ydir'+name_end
plotname_percolation_P_zdir    = percfoldername + infilename_base + '_P_zdir'+name_end
plotname_percolation_P_everyonetogether = percfoldername + infilename_base + '_P_everythingtogether'+name_end
### For each frame:

### Loading data
# Extracting the numbers I'm after to make the file names (because otherwise the order might be off...)
allthetimesteps = []
regex = re.compile(r'\d+')
# List all files in a directory using os.listdir
basepath = foldername#'/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Voxelmatrices/Frames_voxelmap_test_short/'
for entry in os.listdir(basepath):
    if os.path.isfile(os.path.join(basepath, entry)):
        filename = entry#str(entry)
        thesenumbers = regex.findall(filename)
        if len(thesenumbers)!=0:
            allthetimesteps.append(thesenumbers[-1]) # I get the numbers out this way :D Now, I only need to make sure that the timestep is the last number in the filename. Should be easy
                                                     # Oooooh, I also need to remove repetitions of the numbers so that the program only performs the operation ONCE for each number
#print('allthetimesteps:',allthetimesteps)
allthetimesteps = np.array(allthetimesteps)  # This is highly redundant
timestepnumbers = np.unique(allthetimesteps) # And I want the unique time step numbers
#timestepnumbers = np.array(['2340000','3340000'])        # FOR TESTING ONLY!!! UNCOMMENT AS SOON AS THE PLOTS LOOK THE SAME!
Nsteps          = len(timestepnumbers)

# To compute time averages and rms:
# Pi
Pi_av_store   = np.zeros((Nsteps,Nthr)) # Needs to be of the same length as the number of thresholds
Pi_av_x_store = np.zeros((Nsteps,Nthr)) # Needs to be of the same length as the number of thresholds
Pi_av_y_store = np.zeros((Nsteps,Nthr)) # Needs to be of the same length as the number of thresholds
Pi_av_z_store = np.zeros((Nsteps,Nthr)) # Needs to be of the same length as the number of thresholds
Pi_av_one_store = np.zeros((Nsteps,Nthr)) # Needs to be of the same length as the number of thresholds

# P
P_av_store   = np.zeros((Nsteps,Nthr)) # Needs to be of the same length as the number of thresholds
P_av_x_store = np.zeros((Nsteps,Nthr)) # Needs to be of the same length as the number of thresholds
P_av_y_store = np.zeros((Nsteps,Nthr)) # Needs to be of the same length as the number of thresholds
P_av_z_store = np.zeros((Nsteps,Nthr)) # Needs to be of the same length as the number of thresholds
P_av_one_store = np.zeros((Nsteps,Nthr)) # Needs to be of the same length as the number of thresholds


## Adding random walk on the binary matrix:
Nstepsrandom = 1000#10001 # Don't know how many steps are appropriate
Nstarray     = np.arange(0,Nstepsrandom+1)
Rrandom      = np.zeros((Nstepsrandom+1,Nthr))
steps, partition_walks, numberofsamples, len_all, lengths = partition_holders(Nthr,Nstepsrandom,minlength)

# Testing part:
Nsteps = 6



#print('!!! timestepnumbers:',timestepnumbers)
# Can loop from here:
for timeind in range(Nsteps):
    timestep = timestepnumbers[timeind]
    infilename_vmat = infilename_totalbase +'_vox_matrix_timestep'+timestep+'.npy' # _vox_matrix_timestep
    infilename_x    = infilename_totalbase +'_x_timestep'+timestep+'.npy'
    #print('TIMESTEP:', timestep)
    
    #print('infilename_vmat:', infilename_vmat)
    #print('infilename_x:', infilename_x)
    
    #outfilename_txt  = outfilename_npy+'.txt'
    infilename_x    = infilename_totalbase+'_x.npy' # What to do with these for NPT... (in make_voxmap_several.py)
    infilename_y    = infilename_totalbase+'_y.npy'
    infilename_z    = infilename_totalbase+'_z.npy'
    #infilename_vox  = infilename_base+'_vox.txt'
    #infilename_vmat = infilename_totalbase+'_vox_matrix.npy'
    
    #outfile_txt  = open(outfilename_txt,'r')
    x_vals    = np.load(infilename_x) # Do I really need all these? Or just the spacing to translate the distances?
    y_vals    = np.load(infilename_y)
    z_vals    = np.load(infilename_z)
    vmat      = np.load(infilename_vmat) # I guess I should use this one...
    
    #print('shape, vmat:', np.shape(vmat))
    
    matsize   = np.size(vmat)
    dx_map    = x_vals[1] - x_vals[0] # Grid size
    
    dims      = np.shape(vmat) 
    Lx        = dims[0]
    Ly        = dims[1]
    Lz        = dims[2]
    #print('Lx:', Lx, 'Ly:', Ly, 'Lz:', Lz)
    
    #print('dx_map',dx_map)
    
    # Agh, double loop! :(
    # This makes True and False. Would rather have 0 and 1, I guess:
    # Should I just make a bunch of binary matrices?
    for thrind in range(Nthr):
        thr = thrs[thrind]
        #print('thrind:', thrind, '; thr:', thr)
        #print('infilename_vmat:', infilename_vmat)
        
        #print('vmat:',vmat)
        #print('min(vmat):', np.amin(vmat))
        
        vmat_binary      = vmat > thr # This yields 'solid' 0 and 'pore' 1. perctools finds clusters of value 1 (True) and checks if these are connected
        #vmat_binary      = vmat > thr \subsubsection*{make\_voxelmap\_severalframes\_jit_lean.py}# Testing.
        #print(vmat_binary)
        thesepositions, thiswalk, doesitcount = perctools.randomwalk_on_3D_clusters(Nstepsrandom,vmat_binary)
        partition_walks = partition_rw(thesepositions, partition_walks, Nstepsrandom, minlength, thrind, Nthr)
        counter[thrind] += doesitcount
        #print('thiswalk:', thiswalk)
        Rrandom[:,thrind] = Rrandom[:,thrind] + thiswalk # Casting problem here, that is why this is 'ugly'
        
        if plotbool==True and timeind==(Nsteps-1): # Save the last time step for testing
            thrindstring     = '_binary_thr' + str(thr) + '_timestep'+str(timestep)
            plotname_this    = percfoldername + thrindstring + '.png'
            
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            pos = np.where(vmat_binary==False)         # Will only plot the solid. Can infer the percolation from that 
            if issphere:
                pos = np.where(vmat_binary==True)
            ax.scatter(pos[0], pos[1], pos[2], c='black')
            plt.savefig(plotname_this)
        
        # Pi
        perc_all, perc_one, perc_x, perc_y, perc_z = perctools.is_the_given_matrix_percolating(vmat_binary)
        Pi_av_store[timeind,thrind]      = perc_all
        Pi_av_one_store[timeind,thrind]  = perc_one
        Pi_av_x_store[timeind,thrind]    = perc_x
        Pi_av_y_store[timeind,thrind]    = perc_y
        Pi_av_z_store[timeind,thrind]    = perc_z
        Pi_av[thrind]                   += perc_all
        Pi_av_one[thrind]               += perc_one
        Pi_av_x[thrind]                 += perc_x
        Pi_av_y[thrind]                 += perc_y
        Pi_av_z[thrind]                 += perc_z
        
        P_all, P_one, P_x, P_y, P_z = perctools.size_of_percolating_cluster_given_matrix(vmat_binary)
        P_av_store[timeind,thrind]      = P_all
        P_av_one_store[timeind,thrind]  = P_one
        P_av_x_store[timeind,thrind]    = P_x
        P_av_y_store[timeind,thrind]    = P_y
        P_av_z_store[timeind,thrind]    = P_z
        P_av[thrind]                   += P_all
        P_av_one[thrind]               += P_one
        P_av_x[thrind]                 += P_x
        P_av_y[thrind]                 += P_y
        P_av_z[thrind]                 += P_z
    # Loop over thresholds finished
    
    new_time = time.process_time()
    if timeind/Nsteps*100==0 or timeind/Nsteps*100==10  or timeind/Nsteps*100==25 or timeind/Nsteps*100==50 or timeind/Nsteps*100==75 or timeind/Nsteps*100==90:
        print('Time step no.', timeind, ' of', Nsteps, ' done. Time:', new_time-start_time, " s")
    
# Loop over time steps finished
# We use that to get the averages, so I don't think I will do much analysis before this point.

# Finding averages and rms values: 
# Pi
Pi_av     /= float(Nsteps)
Pi_av_one /= float(Nsteps)
Pi_av_x   /= float(Nsteps)
Pi_av_y   /= float(Nsteps)
Pi_av_z   /= float(Nsteps)

# P
P_av     /= float(Nsteps)
P_av_one /= float(Nsteps)
P_av_x   /= float(Nsteps)
P_av_y   /= float(Nsteps)
P_av_z   /= float(Nsteps)

#Rrandom  /= Nsteps # Do this?
print('Rrandom[:,-1]:',Rrandom[:,-1])

for i in range(Nthr):
    if counter[i]>0:
        Rrandom[:,i] /= counter[i] # Or this?

Npartitions = math.floor(Nstepsrandom/minlength)

# Finding the slope
av_slope    = np.zeros(Nthr)
rms_slope   = np.zeros(Nthr)
store_slope = np.zeros((Npartitions,Nthr))
fitsteps    = np.arange(0,minlength)
slope_variance_fit = np.zeros(Nthr)


for i in range(Npartitions):
    partition_walks[i] /= Nsteps

#thr
# Make lists
for j in range(Nthr):
    indata = np.zeros(len_all) # Store for plotting?
    insteps = np.zeros(len_all)
    start_ind = 0
    for i in range(Npartitions):
        end_ind = start_ind + lengths[i]
        this_walk   = partition_walks[i]
        these_steps = steps[i]
        indata[start_ind:end_ind] = this_walk[j,:]
        insteps[start_ind:end_ind] = these_steps
        start_ind = end_ind
    if j==2:
        plt.figure()
        plt.plot(insteps,indata, ',')
        plt.xlabel('Step')
        plt.ylabel(r'Distance$^2$ [voxel length]')
        plt.title('Sectioned walks, thr=%i' % thrs[2])
        plt.savefig(plotname_randomwalk_test)
    #### This part seems wrong:
    #fitpart = this_walk[j,0:minlength] # Works #... Should I find for the first 300 steps too? ... Then I can't use the last part, 'cause that only has 200 elements.
    #coeffs, covs = np.polyfit(fitsteps, fitpart, 1, full=False, cov=True)
    ####
    coeffs, covs = np.polyfit(insteps, indata, 1, full=False, cov=True) # Using all the data in the fit
    av_slope[j]  = coeffs[0]
    rms_slope[j] = np.sqrt(covs[0,0])
''' # Deprecated:
for i in range(Npartitions):
    # Averaging
    partition_walks[i] /= Nsteps
    this_walk = partition_walks[i] #
    # Finding slope
    for j in range(Nthr):
        fitpart = this_walk[j,0:minlength] # Works #... Should I find for the first 300 steps too? ... Then I can't use the last part, 'cause that only has 200 elements.
        coeffs, covs = np.polyfit(fitsteps, fitpart, 1, full=False, cov=True)
        slope_this   = coeffs[0]
        av_slope[j] += slope_this
        store_slope[i,j] = slope_this
        slope_variance_fit[j] += covs[0,0] 
        
av_slope /= Npartitions
slope_variance_fit /= Npartitions

for j in range(Nthr):
    for i in range(Npartitions):
        rms_slope[j] += (av_slope[j] - store_slope[i,j])**2 # Performing the sum
    rms_slope[j] = np.sqrt(rms_slope[j]/(Npartitions-1))    # Finding the rms by dividing and then taking the root
    slope_variance_fit[j] = np.sqrt(slope_variance_fit[j])  # I'm not sure exactly how correct this is... Doing this assumes the means are the same... But how do you find the mean of a parameter estimate...? Is that the coefficient it returns?
'''

print('Rrandom[:,-1]:',Rrandom[:,-1])

new_time = time.process_time()
print('Time:', new_time-start_time, " s")

outfile_Pi = open(outfilename_percdata_Pi,'w') # Need to print P too. Do that in a separate file?
outfile_P  = open(outfilename_percdata_P, 'w')
outfile_randomwalk = open(outfilename_randomwalk,'w')
outfile_randomwalk_partitions = open(outfilename_randomwalk_partiions,'w')

# Making a header
outfile_Pi.write('Threshold; Pi_av, Pi_rms, Pi_av, Pi_one_rms, Pi_av_x, Pi_rms_x, Pi_av_y, Pi_rms_y, Pi_av_z, Pi_rms_z\n')
outfile_P.write('Threshold; P_av, P_rms, P_av, P_one_rms, P_av_x, P_rms_x, P_av_y, P_rms_y, P_av_z, P_rms_z\n')
outfile_randomwalk.write('Number of time steps (for averaging); Number of steps in random walk:\n%i %i\n' % (Nsteps,Nstepsrandom))
outfile_randomwalk_partitions.write('Time steps: %i ; Number of steps in random walk: %i ; minmum length of section: %i\nThreshold; average slope of fit, rms of slope from fit (uncertainty from time steps)\n' % (Nsteps,Nstepsrandom, minlength)) #, variance of slope from the fit (uncertainty in fitting)\n' 

for thrind in range(Nthr):
    for timeind in range(Nsteps):
        # Pi
        Pi_rms[thrind]     += (Pi_av_store[timeind,thrind] - Pi_av[thrind])**2
        Pi_rms_one[thrind] += (Pi_av_one_store[timeind,thrind] - Pi_av_one[thrind])**2
        Pi_rms_x[thrind]   += (Pi_av_x_store[timeind,thrind] - Pi_av_x[thrind])**2
        Pi_rms_y[thrind]   += (Pi_av_y_store[timeind,thrind] - Pi_av_y[thrind])**2
        Pi_rms_z[thrind]   += (Pi_av_z_store[timeind,thrind] - Pi_av_z[thrind])**2
        
        # P
        P_rms[thrind]     += (P_av_store[timeind,thrind] - P_av[thrind])**2
        P_rms_one[thrind] += (P_av_one_store[timeind,thrind] - P_av_one[thrind])**2
        P_rms_x[thrind]   += (P_av_x_store[timeind,thrind] - P_av_x[thrind])**2
        P_rms_y[thrind]   += (P_av_y_store[timeind,thrind] - P_av_y[thrind])**2
        P_rms_z[thrind]   += (P_av_z_store[timeind,thrind] - P_av_z[thrind])**2
    # Pi    
    Pi_rms[thrind]     = np.sqrt(Pi_rms[thrind]/(Nsteps-1))
    Pi_rms_one[thrind] = np.sqrt(Pi_rms_one[thrind]/(Nsteps-1))
    Pi_rms_x[thrind]   = np.sqrt(Pi_rms_x[thrind]/(Nsteps-1))
    Pi_rms_y[thrind]   = np.sqrt(Pi_rms_y[thrind]/(Nsteps-1))
    Pi_rms_z[thrind]   = np.sqrt(Pi_rms_z[thrind]/(Nsteps-1))
    
    # P
    P_rms[thrind]     = np.sqrt(P_rms[thrind]/(Nsteps-1))
    P_rms_one[thrind] = np.sqrt(P_rms_one[thrind]/(Nsteps-1))
    P_rms_x[thrind]   = np.sqrt(P_rms_x[thrind]/(Nsteps-1))
    P_rms_y[thrind]   = np.sqrt(P_rms_y[thrind]/(Nsteps-1))
    P_rms_z[thrind]   = np.sqrt(P_rms_z[thrind]/(Nsteps-1))
    # Writing Pi
    outfile_Pi.write('%.5f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n' % (thrs[thrind], Pi_av[thrind], Pi_rms[thrind], Pi_av_one[thrind], Pi_rms_one[thrind], Pi_av_x[thrind], Pi_rms_x[thrind], Pi_av_y[thrind], Pi_rms_y[thrind], Pi_av_z[thrind], Pi_rms_z[thrind]))
    # Writing P
    outfile_P.write('%.5f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n' % (thrs[thrind], P_av[thrind], P_rms[thrind], P_av_one[thrind], P_rms_one[thrind], P_av_x[thrind], P_rms_x[thrind], P_av_y[thrind], P_rms_y[thrind], P_av_z[thrind], P_rms_z[thrind]))
    # Writing random walk
    outfile_randomwalk.write('%.5f' % thrs[thrind])
    for i in range(N):
        outfile_randomwalk.write(' %.16f' % Rrandom[i,thrind]) # Could use some rms value
    outfile_randomwalk.write('\n')
    outfile_randomwalk_partitions.write('%.5f %.16f %.16f\n' % (thrs[thrind], av_slope[thrind], rms_slope[thrind])) # Deprecated --> #, slope_variance_fit[thrind]))

outfile_Pi.close()
outfile_P.close()
outfile_randomwalk.close()
outfile_randomwalk_partitions.close()

# Finding pc (crude approach):
threshold_for_pc = 0.5
confidence_range = np.exp(-1) # 0.1
pc, slope, startrange, endrange, rangewidth, isthestartrangereal, istheendrangereal = crude_pcapproach_readoffgraph(thrs, Pi_av_z, threshold_for_pc, confidence_range) # Should this give more information? Should perhaps give the span of p for which Pi is +/-0.1 off this value # Or slope in this point? # Yeah, slope would be good. # Should have both, just in case. I do not write much to this file anyways, and will not have many files of this kind.

outfile_pc = open(outfilename_percdata_pc,'w')
outfile_pc.write('%.16f %.2e %.5f %.5f %.5f' % (pc, slope, startrange, endrange, rangewidth))
outfile_pc.close()

# Plotting:
'''
# Pi
plt.figure(figsize=(6,5))
plt.errorbar(thrs, Pi_av, yerr=Pi_rms, capsize=2) #fmt="none", capsize=2)#, label='Values')
plt.plot(thrs, Pi_av, 'o')
plt.xlabel('Threshold [in voxel lengths]')
plt.ylabel(r'Percolation probability, $\Pi$')
plt.title('Percolation probability vs threshold value, all directions')
plt.tight_layout()
plt.savefig(plotname_percolation_alldirs)

plt.figure(figsize=(6,5))
plt.errorbar(thrs, Pi_av_x, yerr=Pi_rms_x, capsize=2) #fmt="none", capsize=2)#, label='Values')
plt.plot(thrs, Pi_av_x, 'o')
plt.xlabel('Threshold [in voxel lengths]')
plt.ylabel(r'Percolation probability, $\Pi$')
plt.title('Percolation probability vs threshold value, x-direction')
plt.tight_layout()
plt.savefig(plotname_percolation_xdir)


plt.figure(figsize=(6,5))
plt.errorbar(thrs, Pi_av_y, yerr=Pi_rms_y, capsize=2) #fmt="none", capsize=2)#, label='Values')
plt.plot(thrs, Pi_av_y, 'o')
plt.xlabel('Threshold [in voxel lengths]')
plt.ylabel(r'Percolation probability, $\Pi$')
plt.title('Percolation probability vs threshold value, y-direction')
plt.tight_layout()
plt.savefig(plotname_percolation_ydir)


plt.figure(figsize=(6,5))
plt.errorbar(thrs, Pi_av_z, yerr=Pi_rms_z, capsize=2) #fmt="none", capsize=2)#, label='Values')
plt.plot(thrs, Pi_av_z, 'o')
plt.xlabel('Threshold [in voxel lengths]')
plt.ylabel(r'Percolation probability, $\Pi$')
plt.title('Percolation probability vs threshold value, z-direction')
plt.tight_layout()
plt.savefig(plotname_percolation_zdir)

plt.figure(figsize=(6,5))
plt.plot(thrs, Pi_av, label='all dirs')
plt.plot(thrs, Pi_av_one, label='one dir')
plt.plot(thrs, Pi_av_x, label='x-dir')
plt.plot(thrs, Pi_av_y, label='y-dir')
plt.plot(thrs, Pi_av_z, label='z-dir')
plt.xlabel('Threshold [in voxel lengths]')
plt.ylabel(r'Size of spanning cluster, $P$')
plt.title('$P$ vs threshold value, all directions')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(plotname_percolation_everyonetogether)


# P

plt.figure(figsize=(6,5))
plt.errorbar(thrs, P_av, yerr=P_rms, capsize=2) #fmt="none", capsize=2)#, label='Values')
plt.plot(thrs, P_av, 'o')
plt.xlabel('Threshold [in voxel lengths]')
plt.ylabel(r'Size of spanning cluster, $P$')
plt.title('$P$ vs threshold value, all directions')
plt.tight_layout()
plt.savefig(plotname_percolation_P_alldirs)

plt.figure(figsize=(6,5))
plt.errorbar(thrs, P_av_x, yerr=P_rms_x, capsize=2) #fmt="none", capsize=2)#, label='Values')
plt.plot(thrs, P_av_x, 'o')
plt.xlabel('Threshold [in voxel lengths]')
plt.ylabel(r'Percolation probability, $\Pi$')
plt.title('Percolation probability vs threshold value, x-direction')
plt.tight_layout()
plt.savefig(plotname_percolation_P_xdir)


plt.figure(figsize=(6,5))
plt.errorbar(thrs, P_av_y, yerr=P_rms_y, capsize=2) #fmt="none", capsize=2)#, label='Values')
plt.plot(thrs, P_av_y, 'o')
plt.xlabel('Threshold [in voxel lengths]')
plt.ylabel(r'Size of spanning cluster, $P$')
plt.title('$P$ vs threshold value, y-direction')
plt.tight_layout()
plt.savefig(plotname_percolation_P_ydir)


plt.figure(figsize=(6,5))
plt.errorbar(thrs, P_av_z, yerr=P_rms_z, capsize=2) #fmt="none", capsize=2)#, label='Values')
plt.plot(thrs, P_av_z, 'o')
plt.xlabel('Threshold [in voxel lengths]')
plt.ylabel(r'Size of spanning cluster, $P$')
plt.title('$P$ vs threshold value, z-direction')
plt.tight_layout()
plt.savefig(plotname_percolation_P_zdir)

plt.figure(figsize=(6,5))
plt.plot(thrs, P_av, label='all dirs')
plt.plot(thrs, P_av_one, label='one dir')
plt.plot(thrs, P_av_x, label='x-dir')
plt.plot(thrs, P_av_y, label='y-dir')
plt.plot(thrs, P_av_z, label='z-dir')
plt.xlabel('Threshold [in voxel lengths]')
plt.ylabel(r'Size of spanning cluster, $P$')
plt.title('$P$ vs threshold value')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(plotname_percolation_P_everyonetogether)
'''

plt.figure(figsize=(6,5))
plt.plot(Nstarray, Rrandom[:,0], label='thr=%i' % thrs[0])
plt.plot(Nstarray, Rrandom[:,10], label='thr=%i' % thrs[10])
plt.plot(Nstarray, Rrandom[:,-1], label='thr=%i' % thrs[-1])
#plt.plot(Nstarray, Rrandom[:,1000], label='Step 1000')
#plt.plot(Nstarray, Rrandom[:,10000], label='Step 10000')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in voxel lengths]')
plt.title('Random walk in pore space')
plt.legend(loc='upper left')
plt.tight_layout()
#plt.show()
plt.savefig(plotname_randomwalk_entire)


rw_longest = partition_walks[0]
rw_3 = partition_walks[1]
rw_2 = partition_walks[2]
rw_1 = partition_walks[3]
rw_shortest = partition_walks[4]

steps_longest = steps[0]
steps_3 = steps[1]
steps_2 = steps[2]
steps_1 = steps[3]
steps_shortest = steps[4]

plt.figure(figsize=(6,5))
plt.plot(steps_longest, rw_longest[2,:])
plt.plot(steps_1, rw_1[2,:])
plt.plot(steps_2, rw_2[2,:])
plt.plot(steps_3, rw_3[2,:])
plt.plot(steps_shortest, rw_shortest[2,:])
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in voxel lengths]')
plt.title('Random walk in pore space, thr=%i' % thrs[2])
plt.tight_layout()
#plt.savefig(plotname_randomwalk_walks)
plt.show()

plt.figure(figsize=(6,5))
plt.errorbar(thrs, av_slope, yerr=rms_slope, capsize=2)
plt.xlabel(r'Threshold')
plt.ylabel(r'Slope, Distance$^2$ [in voxel lengths]')
plt.title('Random walk in pore space: Slope vs threshold')
plt.tight_layout()
plt.savefig(plotname_randomwalk_slope)
#plt.show()
 
# No need for BoundingBox or stuff like that. Should I label clusters, or just let the morhpology functions do their work?
# Should use a sphere for structuring element, I guess. Can see if ndimage provides that or if I should come up with it myself.

# Look at the morphology section of scipy.ndimage!:
# https://docs.scipy.org/doc/scipy/reference/ndimage.html
# scipy.ndimage.binary_opening
# scipy.ndimage.binary_closing
# Can choose different structuring elements to get different information.
# Other commands (from Svenn-Arne's page https://dragly.org/2013/03/25/working-with-percolation-clusters-in-python/):
# Furthermore we want to label each cluster. This is performed by the scipy.ndimage.measurements.label function:
# lw, num = measurements.label(z)
# After this we may want to extract some properties of the clusters, such as the area. This is possible using the scipy.ndimage.measurements.sum function:
# area = measurements.sum(z, lw, index=arange(lw.max() + 1))

