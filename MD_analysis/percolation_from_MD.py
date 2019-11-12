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
Nthr  = 5
thrs  = np.linspace(2.1,100.0,Nthr) # Or something

# Pi_av, percolation probability vs threshold:
Pi_av   = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
Pi_av_x = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
Pi_av_y = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
Pi_av_z = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds

# Pi_rms, rms of Pi_av:
Pi_rms   = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
Pi_rms_x = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
Pi_rms_y = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds
Pi_rms_z = np.zeros(Nthr) # Needs to be of the same length as the number of thresholds

# Bools for plotting
issphere        = False # For plotting if we have a spherical pore
plotconfig      = True
plotbool        = True


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
foldername             = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Voxelmatrices/Frames_'+infilename_base+'/'
plotfoldername         = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Voxelmatrices/Frames_'+infilename_base+'/Plots/'
distfoldername         = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Voxelmatrices/Frames_'+infilename_base+'/Distmatrices/'
percfoldername         = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Voxelmatrices/Frames_'+infilename_base+'/Percolation/'
infilename_totalbase   = foldername + infilename_base
plotname_totalbase     = plotfoldername + infilename_base
distfilename_totalbase = distfoldername + infilename_base
outfilename_percdata   = percfoldername + '_percolation_data.txt'
plotname_percolation_alldirs = percfoldername + infilename_base + '_percolation_total'
plotname_percolation_xdir    = percfoldername + infilename_base + '_percolation_xdir'
plotname_percolation_ydir    = percfoldername + infilename_base + '_percolation_ydir'
plotname_percolation_zdir    = percfoldername + infilename_base + '_percolation_zdir'
plotname_percolation_everyonetogether = percfoldername + infilename_base + '_percolation_everythingtogether'
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
print('allthetimesteps:',allthetimesteps)
allthetimesteps = np.array(allthetimesteps)  # This is highly redundant
timestepnumbers = np.unique(allthetimesteps) # And I want the unique time step numbers
#timestepnumbers = np.array(['2340000','3340000'])        # FOR TESTING ONLY!!! UNCOMMENT AS SOON AS THE PLOTS LOOK THE SAME!
Nsteps          = len(timestepnumbers)

# To compute time averages and rms:
Pi_av_store   = np.zeros((Nsteps,Nthr)) # Needs to be of the same length as the number of thresholds
Pi_av_x_store = np.zeros((Nsteps,Nthr)) # Needs to be of the same length as the number of thresholds
Pi_av_y_store = np.zeros((Nsteps,Nthr)) # Needs to be of the same length as the number of thresholds
Pi_av_z_store = np.zeros((Nsteps,Nthr)) # Needs to be of the same length as the number of thresholds



print('!!! timestepnumbers:',timestepnumbers)
# Can loop from here:
for timeind in range(Nsteps):
    timestep = timestepnumbers[timeind]
    infilename_vmat = infilename_totalbase +'_vox_matrix_timestep'+timestep+'.npy' # _vox_matrix_timestep
    infilename_x    = infilename_totalbase +'_x_timestep'+timestep+'.npy'
    print('TIMESTEP:', timestep)
    
    print('infilename_vmat:', infilename_vmat)
    print('infilename_x:', infilename_x)
    
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
    print('Lx:', Lx, 'Ly:', Ly, 'Lz:', Lz)
    
    #print('dx_map',dx_map)
    
    # Agh, double loop! :(
    # This makes True and False. Would rather have 0 and 1, I guess:
    # Should I just make a bunch of binary matrices?
    for thrind in range(Nthr):
        thr = thrs[thrind]
        #print('infilename_vmat:', infilename_vmat)
        
        #print('vmat:',vmat)
        #print('min(vmat):', np.amin(vmat))
        
        vmat_binary      = vmat > thr # This yields 'solid' 0 and 'pore' 1. perctools finds clusters of value 1 (True) and checks if these are connected
        #vmat_binary      = vmat > thr # Testing.
        #print(vmat_binary)
        
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
        
        
        perc_all, perc_x, perc_y, perc_z = perctools.is_the_given_matrix_percolating(vmat_binary)
        Pi_av_store[timeind,thrind]      = perc_all
        Pi_av_x_store[timeind,thrind]    = perc_x
        Pi_av_y_store[timeind,thrind]    = perc_y
        Pi_av_z_store[timeind,thrind]    = perc_z
        Pi_av[thrind]                   += perc_all
        Pi_av_x[thrind]                 += perc_x
        Pi_av_y[thrind]                 += perc_y
        Pi_av_z[thrind]                 += perc_z
    # Loop over thresholds finished
# Loop over time steps finished
# We use that to get the averages, so I don't think I will do much analysis before this point.

# Finding Pi_av and Pi_rms: 

Pi_av   /= float(Nsteps)
Pi_av_x /= float(Nsteps)
Pi_av_y /= float(Nsteps)
Pi_av_z /= float(Nsteps)

outfile = open(outfilename_percdata,'w')

# Making a header
outfile.write('Threshold; Pi_av, Pi_rms, Pi_av_x, Pi_rms_x, Pi_av_y, Pi_rms_y, Pi_av_z, Pi_rms_z\n')

for thrind in range(Nthr):
    for timeind in range(Nsteps):
        Pi_rms[thrind]   += (Pi_av_store[timeind,thrind] - Pi_av[thrind])**2
        Pi_rms_x[thrind] += (Pi_av_x_store[timeind,thrind] - Pi_av_x[thrind])**2
        Pi_rms_y[thrind] += (Pi_av_y_store[timeind,thrind] - Pi_av_y[thrind])**2
        Pi_rms_z[thrind] += (Pi_av_z_store[timeind,thrind] - Pi_av_z[thrind])**2
    Pi_rms[thrind]   = np.sqrt(Pi_rms[thrind]/(Nsteps-1))
    Pi_rms_x[thrind] = np.sqrt(Pi_rms_x[thrind]/(Nsteps-1))
    Pi_rms_y[thrind] = np.sqrt(Pi_rms_y[thrind]/(Nsteps-1))
    Pi_rms_z[thrind] = np.sqrt(Pi_rms_z[thrind]/(Nsteps-1))
    outfile.write('%.5f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n' % (thrs[thrind], Pi_av[thrind], Pi_rms[thrind], Pi_av_x[thrind], Pi_rms_x[thrind], Pi_av_y[thrind], Pi_rms_y[thrind], Pi_av_z[thrind], Pi_rms_z[thrind]))

 
# Set of plot commands:
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
plt.errorbar(thrs, Pi_av, yerr=Pi_rms, capsize=2, label='all dirs')
plt.errorbar(thrs, Pi_av_x, yerr=Pi_rms_x, capsize=2, label='x-dir')
plt.errorbar(thrs, Pi_av_y, yerr=Pi_rms_y, capsize=2, label='y-dir')
plt.errorbar(thrs, Pi_av_z, yerr=Pi_rms_z, capsize=2, label='z-dir')
plt.plot(thrs, Pi_av_z, 'o')
plt.xlabel('Threshold [in voxel lengths]')
plt.ylabel(r'Percolation probability, $\Pi$')
plt.title('Percolation probability vs threshold value')
plt.legend(loc='center right')
plt.tight_layout()
plt.savefig(plotname_percolation_everyonetogether)

print('thrs:', thrs)
print('Pi_av_x_store[:,0]:',Pi_av_x_store[:,0])
print('Pi_av_x_store[:,1]:',Pi_av_x_store[:,1])
print('Pi_av_x_store[:,2]:',Pi_av_x_store[:,2])
print('Pi_av_x_store[:,3]:',Pi_av_x_store[:,3])
print('Pi_av_x_store[:,4]:',Pi_av_x_store[:,4], '\n')


print('Pi_av_y_store[:,0]:',Pi_av_y_store[:,0])
print('Pi_av_y_store[:,1]:',Pi_av_y_store[:,1])
print('Pi_av_y_store[:,2]:',Pi_av_y_store[:,2])
print('Pi_av_y_store[:,3]:',Pi_av_y_store[:,3])
print('Pi_av_y_store[:,4]:',Pi_av_y_store[:,4], '\n')

print('Pi_av_z_store[:,0]:',Pi_av_z_store[:,0])
print('Pi_av_z_store[:,1]:',Pi_av_z_store[:,1])
print('Pi_av_z_store[:,2]:',Pi_av_z_store[:,2])
print('Pi_av_z_store[:,3]:',Pi_av_z_store[:,3])
print('Pi_av_z_store[:,4]:',Pi_av_z_store[:,4])

 
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

