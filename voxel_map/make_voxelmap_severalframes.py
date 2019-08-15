from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt               # To plot
from scipy.optimize import curve_fit
from pathlib import Path
import numpy as np
import random
import math
import time
import os
import glob

# Function for curve_fit
def costheta_exponential(s,P):
    return np.exp(-s/P)

start_time = time.process_time()

### There is a possible issue here if I want to go from NVT to NPT. That will determine how much I will put in the loop... 

### Number of time frames
Nframes     = 3                                 # Number of time steps we want to include here
totframes   = 10000000/10000                    # Total number of time steps in the file # This should be 1000
startframe  = 100#250
endframe    = 200#750# totframes-1
indices     = np.linspace(0,totframes-1,Nframes)  # Might need to take floor() of the values here, just in case
for i in range(Nframes):
    indices[i] = math.floor(indices[i])

### Grid setting
zbuffer     = 10 # Don't know if the 'buffer' is big enough
Nvoxels     = 55

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

### For selecting time steps:
dt                  = 0.00045       # The time step of our simulation. 0.00045 ns default for nano
skiplines           = 9             # If we hit 'ITEM:', skip this many steps...

### Input and output file names

''' # Testing; <ree2> is correct for completely straight chains.
N  = 100 # Number of bond vectors = number of units - 1 # Per chain
M  = 9   # Number of chains
L2 = 3   # Determining the shape of the L1xL2 polymer grid on the surface
gridspacing = 40 # The spacing in our grid
Nelem   = M*(N+1)
N       = 101
infilename        = 'chaingrids_totallystraight_N%i_Nchains%i_Ly%i_twofixed_test.lammpstrj'  % (Nelem,M,L2)    # Actual run
outfilename_base  = 'voxelmap_chaingrid_quadratic_M%iN%i_totallystraight_test' % (M,N)
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
infilename        = 'chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
outfilename_base  = 'voxelmap_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
#'''


''' # With wall potential:
infilename        = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
% (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
outfilename_base  = 'voxelmap_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
#'''

# Output names for code testing:
#'''
#infilename   = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,spacing,Kangle,Kbond,charge,T)
infilename      = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
# chaingrid_quadratic_M9N101_gridspacing40_Langevin_wall1.042_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T310_theta0is180_twofirst_are_fixed.lammpstrj
outfilename_base  = 'voxelmap_test_short'

#'''
# Varying the grid spacing # THIS IS NOW THE STANDARD FILE NAMES.
'''
infilename        = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,spacing,Kangle,Kbond,charge,T)
outfilename_base  = 'voxelmap_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,spacing,Kangle,Kbond,charge,T)
#'''

# Varying the charge:
#chaingrid_quadratic_M9N101_Langevin_Kangle${Kangle}_Kbond${Kbond}_charge${charge}_T$T_theta0is180_twofirst_are_fixed.lammpstrj
''' # String appending method
#               chaingrid_quadratic_M9N101_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge0_T310_theta0is180_twofirst_are_fixed
infilename        = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % T     # Actual run
outfilename_base  = 'voxelmap_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed' % T
#'''

''' # String appending method
infilename        = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % T     # Actual run
outfilename_base  = 'voxelmap_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed' % T
#'''

''' # Sprintf method
infilename        = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,Kangle,Kbond,factor,T)     # Actual run
outfilename_base  = 'voxelmap_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,Kangle,Kbond,factor,T)
#'''

'''
infilename        = 'chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.lammpstrj' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)     # Actual run
outfilename_base  = 'voxelmap_chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)
'''

# Foldername
foldername = 'Voxelmatrices/Frames_'+outfilename_base+'/'

# The log file in LAMMPS contains a lot of lines detailing the initiation of the system. We want to access the information that comes after that.
startat             = 50            # To equilibrate. The number of ns we want to skip
dt                  = 0.00045       # The time step of our simulation. 0.00045 ns default for nano
skiplines           = 9             # If we hit 'ITEM:', skip this many steps...

#### Automatic part
print("infilename:", infilename)
infile = open(infilename, "r")
lines = infile.readlines()
# Getting the number of lines, etc.
totlines = len(lines)         # Total number of lines
lineend = totlines-1          # Index of last element

# Extracting the number of atoms:
words = lines[3].split()
#print("words:", words)
#print("words[0]:", words[0])
Nall = int(words[0])
print('Nall:',Nall)
lineindices = indices*(Nall+skiplines) # Should work...
#print('indices[1]:', indices[1])
#print('lineindices[1]:', lineindices[1])
#print('indices*(Nall+skiplines):',indices*(Nall+skiplines))

## Making an array of indices:
# Not sure if this is really the most efficient way of doing things
# Probably is if I have a lot of time steps
index1 = np.zeros(M*N) # The chain index
index2 = np.zeros(M*N) # The atom (in chain) index
counter = 0
for i in range(M):
    for j in range(N):
        index1[counter] = i
        index2[counter] = j
        counter += 1

# These I need.
xes     = np.zeros((M,N)) # Is storing the position after chain AND number in chain really neccessary? And should I perhaps store as a vector?
ys      = np.zeros((M,N))
zs      = np.zeros((M,N))
posvecs = np.zeros((Nall,3))
    
tempzs = np.zeros(M)  # Do I need this?
counter = 1 

# Loop here
for lindex in lineindices:
    lindex = int(math.ceil(lindex))#math.floor(lindex)
    ### Extracting time and volume data from file: ###
    ## Time
    #lindex should yield ['ITEM:', 'TIMESTEP']
    tline     = lines[lindex+1].split()
    timestep  = float(tline[0])         # This should give us the time
    
    ## Volume
    words = lines[lindex+5].split()
    xmin  = float(words[0])
    xmax  = float(words[1])
    
    words = lines[lindex+6].split()
    ymin  = float(words[0])
    ymax  = float(words[1])
    
    words = lines[lindex+7].split()
    zmin  = float(words[0])
    zmax  = float(words[1])
    
    Lx = xmax-xmin
    Ly = ymax-ymin
    Lz = zmax-zmin

    # Maybe I should do something with xmax, ymax, zmax?
    
    halfLx = 0.5*Lx
    halfLy = 0.5*Ly
    halfLz = 0.5*Lz
    
    #xes = np.zeros((N, ?))
    #ys  = np.zeros((N, ?))
    #zs  = np.zeros((N, ?))
    print('lindex:',lindex)
    print('lines[lindex]:', lines[lindex])
    print('Going into the coord.-gathering loop:')
    for i in range(Nall): # This should work, but check!
        words   = lines[i+lindex+skiplines].split()
        #print('words:',words)
        #print('i:', i)
        ind   = int(words[0])-1 # Atom ids go from zero to N-1.
        x     = float(words[3])
        y     = float(words[4])
        z     = float(words[5])
        # Making array indices:
        ind1  = int(index1[ind])
        ind2  = int(index2[ind])
        # Add values to lists
        xes[ind1,ind2] = x
        ys[ind1,ind2]  = y
        zs[ind1,ind2]  = z
        posvecs[ind,0] = x
        posvecs[ind,1] = y
        posvecs[ind,2] = z
    
    '''
    for i in range(M*N):
        print("Atomid-1 =", i, ": chain:", index1[i], ", atom nr.", index2[i])
    '''
    
    #print("xes:", xes)
    
    #print("xes[-1]:",xes[-1])
    #print("ys[-1]:",ys[-1])
    #print("zs[-1]:",zs[-1])
    #print("len, xes:", Nt)
    
    #print("xes[0]:", xes[0])
    #print("ys[0]:", ys[0])
    #print("zs[0]:", zs[0])
    #print("xes:", xes)
    
    ### Make voxel grid ###
    ####  It's enough to make this grid once for NVT, but not for NPT ####
    print('Make voxel grid')
    
    print('time step:', timestep)
    
    '''
    # Print this to a plot for comparison:    # DO THIS IN OVITO INSTEAD!
    plotname = outfilename_npy + '_atomplot.png'
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    size_plot = len(xes_makegrid)
    m = np.zeros(shape = (N*M, N*M, N*M))
    for i in range(size_plot):
        location = (xes[i],ys[i],zs[i])    
        m[location] = 1
    pos = np.where(m==1)
    ax.scatter(pos[0], pos[1], pos[2], c='black')
    plt.savefig(plotname)
    '''
    
    '''
    print('zs_makegrid:',zs_makegrid)
    print('len(zs):',len(zs))
    print('size(zs):',np.size(zs))#.size())
    '''
    ## I want to cut the box a bit: I want to reduce the amount of empty space above my chains
    # Find max value of zs_makegrid:
    zmaxes = []
    for i in range(M):
        zmax_temp = max(zs[i])
        zmaxes.append(zmax_temp)
    zmax_makegrid = max(zmaxes)
    
    # Set edges of grid # I should probably print these, along with information about the number of voxels
    minx_grid          = xmin
    maxx_grid          = xmax
    miny_grid          = ymin
    maxy_grid          = ymax
    minz_grid          = zmin
    maxz_grid          = zmax_makegrid+zbuffer # Don't know if the 'buffer' is big enough
    if abs(maxz_grid-Lx)<zbuffer:              # I don't want to make the voxel map bigger than the simulation box
        maxz_grid          = zmax_makegrid
    
    # Make voxel grid, store information about positions
    # Is the information about the positions stored as the corners or the centres?
    # It is probably easiest to work with the centres when finding the distance to the closest atom.
    # But what if the atom is WITHIN the voxel? Should I just set the value to 0?
    # Should be quite easy to check this? But the shape of the voxel causes some trouble. It is not a sphere, so the distance from the centre of the voxel to a point on one of its edges is not the same for all points. So I can't just compare lengths.
    # Maybe it is easier to just store the distance to the centre even if the atom is inside the voxel.
    # Should I just give a voxel size in the beginning? But that might not fit with the simulation box, and I want to set the number of voxels.
    # But it is difficult to fit a set number of voxels... The simulation box does not have the same length in all direction.
    # OR should I set number of voxels in x-direction and then add voxels (cubes) of that size until I've covered all the simulation box (up until maxz_grid)?
    
    len_voxel = Lx/float(Nvoxels) # This is the length of all sides.
    ns        = np.arange(Nvoxels)
    x_centres = 0.5*(2*ns+1)*len_voxel+minx_grid # Do I really need this one, though?...
    
    #'''
    #print('x_centres:',x_centres) 
    print('Lx:', Lx)
    print('Ly:', Ly)
    print('Lz:', Lz)
    #'''
    # Figure out how many voxels in the other directions
    # In y-direction (guess I could just set this automatically, at least as long as I have constant volume) -- (But if I don't have constant volume, then I'm gonna be in trouble anyways)
    edgehit   = 0
    y_centres = []
    i         = 0
    while edgehit==0:
        centrevalue = 0.5*(2*i+1)*len_voxel+miny_grid # Possibly add starting y-value
        edgevalue   = len_voxel*(i+1)+miny_grid
        y_centres.append(centrevalue)
        i += 1.0
        #print('edgevalue:', edgevalue, 'maxy_grid:', maxy_grid)
        if (maxy_grid-centrevalue)<0.5*len_voxel or abs(maxy_grid-edgevalue)<0.1*len_voxel: # Could have abs() as well, but this works fine since we break it here
            edgehit = 1
    
    y_centres = np.array(y_centres)
    #print('y_centres:', y_centres)
    
    # In z-direction
    edgehit   = 0
    z_centres = []
    i         = 0
    while edgehit==0:
        centrevalue = 0.5*(2*i+1)*len_voxel+minz_grid
        edgevalue   = len_voxel*(i+1)+minz_grid
        z_centres.append(centrevalue)
        i += 1.0
        # Test this:
        if (maxz_grid-centrevalue)<0.5*len_voxel or abs(maxz_grid-edgevalue)<1e-15: # Could have abs() as well, but this works fine since we break it here
            edgehit = 1
    
    z_centres = np.array(z_centres)
    #print('z_centres:', z_centres)
    
    Nx = len(x_centres)
    Ny = len(y_centres)
    Nz = len(z_centres)
    voxN = Nx*Ny*Nz
    
    #'''
    print('xmin:',xmin)
    print('ymin:',ymin)
    print('zmin:',zmin)
    print('len_voxel:',len_voxel)
    #'''
    
    # Regarding indices: Should I have one running index, or one index for each coordinate?
    # Need to work at the coordinate level at some point. I can start with that and tweak the indices later
    # BUT I probably want to do this in a vectorized fashion.
    # I should make a map of the running index to the voxel centre coordinates.
    # First: Make array of positions. # Or don't I need that?
    
    print('voxN:',voxN)

    ### Fill voxel grid ###
    # Should I make an array with xpos, ypos, zpos, voxelvalue? # Or will that take up too much space?
    outarray = np.zeros((voxN,4))
    voxmat   = np.zeros((Nx,Ny,Nz))
    # This is probably the bottleneck
    # Print to file here, I guess...
    for tindex in range(1): # Only look at the first frame. Can generalize.
        print('Find voxel values')
        l = 0 # Can't even remember what this was for...
        n = 0
        voxelvalues = np.zeros(voxN)
        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    smallest_dist = 1e10 # No distance is this big.
                    xc            = x_centres[i]
                    yc            = y_centres[j]
                    zc            = z_centres[k]
                    centrevec     = np.array([xc,yc,zc])
                    # Loop over all the atoms. Oh vey.
                    for m in range(Nall):
                        vecthis = posvecs[m,:]
                        distvec = centrevec-vecthis
                        dotprod = np.dot(distvec,distvec)
                        if dotprod<smallest_dist:
                            smallest_dist = dotprod
                    '''
                    for m in range(M):       # Loop over the chains 
                        x_thischain = xes[m] # Change these if I use multiple timesteps
                        y_thischain = ys[m]  # Change these if I use multiple timesteps
                        z_thischain = zs[m]  # Change these if I use multiple timesteps
                        for nc in range(N):  # Loop over the atoms in a chain
                            vecthis = np.array([x_thischain[nc],y_thischain[nc],z_thischain[nc]]) # Is this really cost-efficient? I should probably just import the data into position arrays...
                            distvec = centrevec-vecthis
                            dotprod = np.dot(distvec,distvec)
                            if dotprod<smallest_dist:
                                smallest_dist = dotprod
                    '''
                    #print('n:',n)
                    #print('smallest_dist:', smallest_dist)
                    voxval         = np.sqrt(smallest_dist)
                    voxelvalues[n] = voxval
                    voxmat[i,j,k]  = voxval
                    outarray[n,0]  = xc
                    outarray[n,1]  = yc
                    outarray[n,2]  = zc
                    outarray[n,3]  = voxval
                    n+=1
    # It's running, so that's a plus. The distances seem way to big, though...
    print('The absolutely smallest distance:', min(voxelvalues))
    print('max(voxelvalues):', max(voxelvalues))
    
    ### Make the folder for the file if it does not already exist
    p = Path(foldername)
    p.mkdir(exist_ok=True)
    
    ### Opening files
    outfilename_npy  = foldername + outfilename_base+'_timestep%i' % timestep
    outfilename_x    = foldername + outfilename_base+'_x_timestep%i' % timestep
    outfilename_y    = foldername + outfilename_base+'_y_timestep%i' % timestep
    outfilename_z    = foldername + outfilename_base+'_z_timestep%i' % timestep
    outfilename_vox  = foldername + outfilename_base+'_vox_timestep%i' % timestep
    outfilename_vmat = foldername + outfilename_base+'_vox_matrix_timestep%i' % timestep
    
    #outfile_txt  = open(outfilename_txt,'w')
    outfile_x    = open(outfilename_x,'w')
    outfile_y    = open(outfilename_y,'w')
    outfile_z    = open(outfilename_z,'w')
    outfile_vox  = open(outfilename_vox,'w')
    outfile_vmat = open(outfilename_vmat,'w')
    
    
    np.save(outfilename_npy,outarray)
    np.save(outfilename_x,x_centres)
    np.save(outfilename_y,y_centres)
    np.save(outfilename_z,z_centres)
    np.save(outfilename_vox,voxelvalues)
    np.save(outfilename_vmat,voxmat)
    
    new_time = time.process_time()
    print('Time step no.', counter, ' of', Nframes, ' done. Time:', new_time-start_time, " s")
    #for i in range(voxN): # Why this no work?
    #   outfile_txt.write('%.16f %.16f %.16f %.16f\n' % (outarray[n,0],outarray[n,1],outarray[n,2],outarray[n,3]))
    #outfile_txt.close()

# Plotting and printing
'''
plt.figure(figsize=(6,5))
plt.errorbar(separation, costheta[-1,:], yerr=costheta_rms[-1,:], fmt="none", capsize=2, label='Values')
plt.plot(separation, costheta[-1,:], 'o')
plt.plot(separation, fittedgraph, label='Fit')
#plt.plot(x_persistencelength, y_persistencelength, '--')
#plt.plot(x_eline, y_eline, '--')
plt.xlabel(r'Bond distance', fontsize=16)
plt.ylabel(r'<cos($\theta$)>', fontsize=16)
plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
plt.legend(loc="upper right")
plt.title(r'<cos($\theta$)> vs separation (last chain)', fontsize=16)
plt.savefig(plot1name)

plt.figure(figsize=(6,5))
plt.plot(costheta_allvalues[5], '.')
plt.xlabel(r'Index', fontsize=16)
plt.ylabel(r'cos($\theta$)', fontsize=16)
plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
plt.title(r'Values of cos($\theta$) for separation 5', fontsize=16)
plt.savefig(plot2name)

outfile8 = open(outfilename8,'w')
for i in range(M):
    outfile8.write('%.16e %.16e\n' % (bondlength_chain_vals[i],bondlength_chain_vals_rms[i]))
outfile8.write('System wide: %.16e %.16e' % (bondlength_av_system, bondlength_rms_system))
'''
infile.close()
print("Script finished.")
