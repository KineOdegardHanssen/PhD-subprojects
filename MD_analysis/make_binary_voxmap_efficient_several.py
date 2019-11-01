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

### CODE DON'T WORK FOR NPT!!!

def make_3Dlist(a, b, c): # Borrowed this from GeeksforGeeks (and renamed it)
    lst = [[ [[] for col in range(c)] for col in range(b)] for row in range(a)] 
    return lst 

# Function taking a box and calculating the distance between its atoms and a box centre
def find_atomdists_givenbox(smallest_dist, natoms_box, searchmore, thisbox, centrevec):
    for l in range(natoms_box):
        vecthis = posvecs[thisbox[l]] # Getting the position vector of an atom in the box
        distvec = centrevec-vecthis
        dotprod = np.dot(distvec,distvec)
        if dotprod<smallest_dist:
            smallest_dist = dotprod
            if dotprod < qrtdistance:
                searchmore = False     # If the atom is sufficiently close to the centre of the box, I don't need to search the neighbouring boxes
    return smallest_dist, searchmore

def order_distances(Nx,Ny,Nz): # Is this too big now?
    # This is actually complicated (because we have cut the system so that generally Nx=Ny!=Nz). Do the easy part first.
    ## Making lists
    distance_indices    = []
    voxcentres_distance = []
    ## Getting elements
    for deltax in range(-Nx+1,Nx): # Should this go from -Nx to Nx? Then this is huuuuge. # Maybe half will do? # At least in the x- and y- directions...
        for deltay in range(-Ny+1,Ny):
            for deltaz in range(-Nz+1,Nz):
                vcd = np.sqrt(deltax**2+deltay**2+deltaz**2) # Distance between centres of voxels
                voxcentres_distance.append(vcd)
                distance_indices.append([deltax,deltay,deltaz])
    ## Co-sort the arrays. Order voxcentres_distance by increasing values without messing with the nice correspondence between that array and distance_indices
    voxcentres_distance, distance_indices = (list(t) for t in zip(*sorted(zip(voxcentres_distance, distance_indices))))
    ## Make an array of the different distances
    ## AND count the number of different distances
    voxcentres_distance_array = np.array(voxcentres_distance)
    voxcentres_distance_slim, number_each_dist = np.unique(voxcentres_distance_array, return_counts=True)
    print('voxcentres_distance_slim',voxcentres_distance_slim)
    print('voxcentres_distance_array',voxcentres_distance_array)
    ## Gather tha distance indices into one neat array as well:
    gathered_indices = []
    start_counter    = 0
    #print('distance_indices:',distance_indices)
    #print('shape(distance_indices):',np.shape(distance_indices))
    for i in range(len(voxcentres_distance_slim)):
        nthis = number_each_dist[i]-1 # I need the -1 for some reason... Must have messed up at some point.
        end_counter = start_counter + nthis
        #gathered_indices.append([distance_indices[start_counter:end_counter]]) # Unfortunately, I am working with a list here...
        temp_list = []
        for j in range(nthis+1):
            temp_list.append(distance_indices[start_counter+j])
            if i==0 or i==1 or i==2 or i==3:
                print('i =', i, ', j =', j, '; index: ', start_counter+j, '; dist =', voxcentres_distance[start_counter+j])
        gathered_indices.append(temp_list)
        # Test of code:
        '''
        if i==0 or i==1 or i==2 or i==3:
            print('i =', i)
            print('start_counter:', start_counter)
            print('end_counter:', end_counter)
            print('nthis:', nthis)
            print('voxcentres_distance[',start_counter,']:',voxcentres_distance[start_counter])
            print('; voxcentres_distance[',start_counter+nthis,']:',voxcentres_distance[start_counter+nthis])
            print('; voxcentres_distance[',start_counter+nthis+1,']:',voxcentres_distance[start_counter+nthis+1])
        '''
        start_counter = end_counter+1
    ## return
    return voxcentres_distance_slim, gathered_indices, number_each_dist
    
    
                


# Function for curve_fit
def costheta_exponential(s,P):
    return np.exp(-s/P)

start_time = time.process_time()

### There is a possible issue here if I want to go from NVT to NPT. That will determine how much I will put in the loop... 

### Number of time frames
Nframes     = 3                                 # Number of time steps we want to include here
totframes   = 10000000/10000                    # Total number of time steps in the file # This should be 1000
startframe  = 300#0#100#250
endframe    = 450#totframes-1#200#750# 
indices     = np.linspace(startframe,endframe,Nframes)#np.linspace(0,totframes-1,Nframes)  # Might need to take floor() of the values here, just in case
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
outfilename_base  = 'voxelmap_test_short_TESTII'

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
xes     = np.zeros(Nall) # Is storing the position after chain AND number in chain really neccessary? And should I perhaps store as a vector?
ys      = np.zeros(Nall)
zs      = np.zeros(Nall)
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
    print('TIMESTEP:', timestep)
    
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
    #print('lines[lindex]:', lines[lindex])
    print('Going into the coord.-gathering loop:')
    for i in range(Nall): # This should work, but check!
        words   = lines[i+lindex+skiplines].split()
        #print('words:',words)
        #print('i:', i)
        ind   = int(words[0])-1 # Atom ids go from zero to N-1.
        x     = float(words[3])-xmin # Shift so that box corner is at (0,0,0) and all elements of any other position in the box is positive
        y     = float(words[4])-ymin
        z     = float(words[5])-zmin
        # Add values to lists
        xes[ind] = x
        ys[ind]  = y
        zs[ind]  = z
        posvecs[ind,0] = x
        posvecs[ind,1] = y
        posvecs[ind,2] = z
    
    Natoms = len(xes)
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
    print('zs_makegrid:',zs_makegrid)
    print('len(zs):',len(zs))
    print('size(zs):',np.size(zs))#.size())
    '''
    ## I want to cut the box a bit: I want to reduce the amount of empty space above my chains
    # Find max value of zs_makegrid:
    zmax_makegrid = max(zs)
    
    # Set edges of grid # I should probably print these, along with information about the number of voxels
    minx_grid          = 0
    maxx_grid          = xmax-xmin
    miny_grid          = 0
    maxy_grid          = ymax-ymin
    minz_grid          = 0
    maxz_grid          = zmax_makegrid-zmin+zbuffer # Don't know if the 'buffer' is big enough
    if abs(maxz_grid-Lx)<zbuffer:                   # I don't want to make the voxel map bigger than the simulation box
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
    
    len_voxel = Lx/float(Nvoxels) # This is the length of all sides. # The voxel length should be the same for all time frames
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
    
    #print('voxN:',voxN)
    
    ### Place atoms in boxes ###
    
    ## Ready the distances:
    print('About to enter order_distances')
    voxcentres_distance_slim, gathered_indices, number_each_dist = order_distances(Nx,Ny,Nz)
    Ndistances = len(voxcentres_distance_slim)
    
    box        = make_3Dlist(Nx, Ny, Nz) # Make box to use
    
    # Point
    
    
    ## Manual:
    #list.append(elem) -- adds a single element to the end of the list. ...
    #list.insert(index, elem) -- inserts the element at the given index, shifting elements to the right.
    #list.extend(list2) adds the elements in list2 to the end of the list.
    
    # Loop over all atoms and place them in box # Index the atoms from 1 to N+1?
    for i in range(Natoms):
        # Finding x,y,z-coordinates
        nx = int(xes[i]/len_voxel)
        ny = int(ys[i]/len_voxel)
        nz = int(zs[i]/len_voxel)
        #print('Nx:', Nx, 'Ny:', Ny, 'Nz:', Nz)
        #print('nx:', nx, '; ny:', ny, '; nz:', nz)
        #print('shape(box):', np.shape(box))
        box[nx][ny][nz].append(i)     # I hope this works...
    
    
    # What Anders wrote down. I don't really understand all of it...:
    # It doesn't seem like linked lists are that useful in Python, more in C++. 
    '''
    particle_counter = 1
    for i in range(Natoms):
        # Finding x,y,z-coordinates
        nx = int(xes[i]/len_voxel)
        ny = int(ys[i]/len_voxel)
        nz = int(zs[i]/len_voxel)
        p1 = box[nx,ny,nz]          # Hmmm...  
        box[nx,ny,nz] = particle_counter
        point[particle_counter] = p1
        particle_counter +=1
    '''
    
    # Moving on:
    ### Fill voxel grid ###
    # Should I make an array with xpos, ypos, zpos, voxelvalue? # Or will that take up too much space?
    outarray = np.zeros((voxN,4))
    voxmat   = np.zeros((Nx,Ny,Nz))
    # This is probably the bottleneck
    # Print to file here, I guess...
    print('Find voxel values')
    l = 0 # Can't even remember what this was for...
    n = 0
    print('posvecs[25,:]:',posvecs[25,:])
    voxelvalues = np.zeros(voxN)
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                smallest_dist  = 1e20 # No distance is this big.
                searchmore     = True
                atomsfound     = False
                xc             = x_centres[i]
                yc             = y_centres[j]
                zc             = z_centres[k]
                centrevec      = np.array([xc,yc,zc]) # Should probably just have stored this at once...
                qrtdistance    = len_voxel/4.         # If the atom is closer to the centre than half the box length, we don't need to search another box.
                thisbox        = box[i][j][k]
                natoms_box     = len(thisbox)
                atomsfound_box = [False,False,False]
                atomsfound_nbbox_index = False
                if natoms_box==0:
                    atomsfound = False
                else:                   # Only inspect box if there are atoms in it.
                    atomsfound = True
                    atomsfound_box = [i,j,k]
                    atomsfound_nbbox_index = 0
                    smallest_dist, searchmore = find_atomdists_givenbox(smallest_dist,natoms_box,searchmore,thisbox,centrevec)
                
                nbn = 0 # neighbournumber 
                while searchmore==True:
                    # Adding +/-1 to all the indices should yield the closest atom. 
                    # But we might hit a boundary at some point. 
                    nbn += 1 # Searching the boxes that are just a little bit further out (i.e. going from distance 0 to 1 or 1 to sqrt(2))
                    
                    if nbn==Ndistances:
                        print('WARNING: Failed to end atom search. Might not have found any atom.')
                        print('smallest_dist:', smallest_dist)
                        break
                    
                    if atomsfound_nbbox_index!=False:
                        if nbn-atomsfound_nbbox_index==2:
                            break # We have found an atom, and only need to search in the boxes closest to that one # Could chech nbn-atomsfound_nbbox+1 and set searchmore to False...
                    
                    #voxcentres_distance_slim, gathered_indices, number_each_dist
                    
                    distance_to_voxel  = voxcentres_distance_slim[nbn]
                    indices_nbn_voxels = gathered_indices[nbn]
                    N_nbn_voxels       = number_each_dist[nbn]
                    
                    for indices in indices_nbn_voxels:
                        # Extract indices
                        deltai = indices[0]
                        deltaj = indices[1]
                        deltak = indices[2]
                        # Check if this neighbour exists:
                        newi = i+deltai
                        newj = j+deltaj
                        newk = k+deltak
                        if newi<Nx and newi>=0 and newj<Ny and newj>=0 and newk<Nz and newk>=0:      
                            thisbox    = box[newi][newj][newk]
                            natoms_box = len(thisbox) 
                            if natoms_box>0:
                                smallest_dist, searchmore = find_atomdists_givenbox(smallest_dist,natoms_box,searchmore,thisbox,centrevec)
                                atomsfound_box = [newi,newj,newk]
                                atomsfound_nbbox_index = nbn
                                atomsfound = True # This here?
                    #atomsfound_now = False
                        
                    # Need another test to turn searchmore off.
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
    outfilename_txt  = foldername + outfilename_base+'_timestep%i' % timestep+'_voxarray.txt'
    
    print('!!! outfilename_npy:', outfilename_npy)
    # Plotting
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
    #'''
    
    outfile_txt  = open(outfilename_txt,'w')
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
    for i in range(voxN): # Why this no work?
       outfile_txt.write('%i %.16f\n' % (i,outarray[i,3]))
    outfile_txt.close()
    print('DONE WITH ONE, ONTO THE NEXT!\n')

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
