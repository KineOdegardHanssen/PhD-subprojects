import matplotlib.pyplot as plt               # To plot
import numpy as np
import random
import math

# 21.03.19: Making a completely random chain for input to LAMMPS.
# I probably don't need a lot of these variables anymore, but this is not a bottleneck
# so I don't bother cleaning it up

# Set variables
N             = 100 # Number of bond vectors = number of units - 1 # Per chain
M             = 9   # Number of chains
L2            = 3   # Determining the shape of the L1xL2 polymer grid on the surface
charge        = -1
gridspacing   = 4 # The spacing in our grid # Default is 40
mass          = 0.00081
Nelem         = M*(N+1)
Nbonds        = M*N
Nangles       = M*(N-1)
xpos          = np.zeros(Nelem)
ypos          = np.zeros(Nelem)
zpos          = np.zeros(Nelem)
atomtypes_all = np.zeros(Nelem)
molID         = np.zeros(Nelem)
qs            = np.zeros(Nelem)
vxs           = np.zeros(Nelem)
vys           = np.zeros(Nelem)
vzs           = np.zeros(Nelem)
vectorlengths = np.zeros(N*M)
vectors       = []
positions     = []

k = 0
l = 0

k = 0
l = 0

for j in range(M):
    # Starting the chain
    x_accm = (j//L2)*gridspacing
    y_accm = j%L2*gridspacing
    z_accm = 0
    positions.append([x_accm,y_accm,0])
    xpos[k]  = x_accm
    ypos[k]  = y_accm
    zpos[k]  = 0
    vxs[k]   = np.random.normal(0,1)
    vys[k]   = np.random.normal(0,1)
    vzs[k]   = np.random.normal(0,1)
    qs[k]    = charge
    molID[k] = j+1
    atomtypes_all[k] = 1
    
    k += 1
    
    # We want another fixed atom per chain:
    chain_twofirst = 0
    # Set up system
    # (Generate random configuration to start with)
    # Should add a constraint to avoid z = 0 plane (or x = 0, or y = 0...)
    for i in range(N):
        # I could maybe flip the sign of z0, but that does not seem very ergodic... 
        # I probably want a Monte Carlo-style test to set the random chain in a realistic way
        z_acc_test = -1
        '''# There is an intialization problem. Random vectors: Some atoms get too close to the wall... I should start from straight chains just to be sure.
        while z_acc_test<0:
            x0 = np.random.normal(0,1)
            y0 = np.random.normal(0,1)
            z0  = np.random.normal(0,1)+2
            den = np.sqrt(x0**2+y0**2+z0**2)
            x = x0/den
            y = y0/den
            z = z0/den
            z_acc_test = z_accm+z
        '''
        x = 0
        y = 0
        z = 1
        x_accm += x
        y_accm += y
        z_accm += z
        qs[k]   = charge
        xpos[k] = x_accm
        ypos[k] = y_accm
        zpos[k] = z_accm
        vectorlengths[l] = np.sqrt(x**2+y**2+z**2)
        # All atoms are of type 2 (=moving):
        atomtypes_all[k] = 2
        # Except the first two:
        if chain_twofirst==0:
            atomtypes_all[k] = 1
            chain_twofirst   = 1 # This test is only true once per chain
        molID[k]         = j+1
        # Set velocity...
        # Should set the velocity according to some temperature, but fix it later
        # Or maybe I don't need to set it here...
        vxs[k] = np.random.normal(0,1)
        vys[k] = np.random.normal(0,1)
        vzs[k] = np.random.normal(0,1)
        qs[k] = charge
        k += 1
        l += 1


### Data on system
'''
xmax = max(xpos)*3.0
ymax = max(ypos)*3.0
zmax = max(zpos)*3.0 
xmin = min(xpos)*3.0
ymin = min(ypos)*3.0
zmin = min(zpos)*3.0 # Change something here...
'''

'''
xmax = max(xpos)*3.0+20.0
ymax = max(ypos)*3.0+20.0
zmax = max(zpos)*3.0
zmin = min(zpos)*3.0 # Change something here...
xmin = min(xpos)*3.0-20.0
ymin = min(ypos)*3.0-20.0
'''

'''
xmax = max(xpos)+2.0*gridspacing
ymax = max(ypos)+2.0*gridspacing
zmax = max(zpos)*3.0
zmin = min(zpos)*3.0 # Change something here...
xmin = min(xpos)-2.0*gridspacing
ymin = min(ypos)-2.0*gridspacing
'''

# Kan evt. ha +- 0.5*gridspacing slik at det ser ut som om systemet er uniformt og uendelig... Best dersom systemet er ganske stort
xmax = max(xpos)+0.5*gridspacing
ymax = max(ypos)+0.5*gridspacing
zmax = max(zpos)*3.0
zmin = min(zpos)*3.0 # Change something here...
xmin = min(xpos)-0.5*gridspacing
ymin = min(ypos)-0.5*gridspacing

# Change names of atom-style, etc

### Print to file
#outfilename = 'data.randomchain_N%i' % N+1
outfilename = 'testsystem_chaingrid_straight.lammpstrj'#'data.chaingrids_N%i_Nchains%i_Ly%i_gridspacing%i_twofixed_charge%i_mass%.5f' % (Nelem,M,L2,gridspacing,charge,mass)
outfile = open(outfilename,'w')

for j in range(100):
    outfile.write('ITEM: TIMESTEP\n%i\nITEM: NUMBER OF ATOMS\n%i\nITEM: BOX BOUNDS pp pp pp\n' % (j,M*(N+1)))
    outfile.write('%.16e %.16e xlo xhi\n' % (xmin, xmax))
    outfile.write('%.16e %.16e ylo yhi\n' % (ymin, ymax))
    outfile.write('%.16e %.16e zlo zhi\n' % (zmin, zmax))
    outfile.write('ITEM: ATOMS id type mol x y z vx vy vz \n')
    for i in range(M*(N+1)):
        # ITEM: ATOMS id type mol x y z vx vy vz
        outfile.write('%i %i %i %.16e %.16e %.16e %.16e %.16e %.16e\n' % ((i+1), atomtypes_all[i], molID[i], xpos[i], ypos[i], zpos[i], vxs[i], vys[i], vzs[i]))

# And the last command:
# Close file
outfile.close()
