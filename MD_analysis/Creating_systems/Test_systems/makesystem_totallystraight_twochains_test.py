import matplotlib.pyplot as plt               # To plot
import numpy as np
import random
import math

# 21.03.19: Making a completely random chain for input to LAMMPS.
# I probably don't need a lot of these variables anymore, but this is not a bottleneck
# so I don't bother cleaning it up

# Set variables
N  = 100 # Number of bond vectors = number of units - 1 # Per chain
M  = 4   # Number of chains
L2 = 1   # Determining the shape of the L1xL2 polymer grid on the surface
gridspacing = 40 # The spacing in our grid
Nelem   = M*(N+1)
Nbonds  = M*N
Nangles = M*(N-1)
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
    qs[k]    = -1
    molID[k] = j+1
    atomtypes_all[k] = 1
    
    k += 1
    
    # We want another fixed atom per chain:
    chain_twofirst = 0
    # Set up system
    # (Generate random configuration to start with)
    # Should add a constraint to avoid z = 0 plane (or x = 0, or y = 0...)
    x = 0
    y = np.sin(j*np.pi/6.)
    z = np.cos(j*np.pi/6.)
    for i in range(N):
        # I could maybe flip the sign of z0, but that does not seem very ergodic... 
        # I probably want a Monte Carlo-style test to set the random chain in a realistic way    
        x_accm += x
        y_accm += y
        z_accm += z
        qs[k]   = -1
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
        qs[k] = -1
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
k = 0
#outfilename = 'data.randomchain_N%i' % N+1
outfilename = 'chaingrids_totallystraight_N%i_Nchains%i_Ly%i_twofixed_test.lammpstrj' % (Nelem,M,L2)
outfile = open(outfilename,'w')
outfile.write('ITEM: TIMESTEP\n%i\nITEM: NUMBER OF ATOMS\n%i\nITEM: BOX BOUNDS pp pp pp\n' % (k,Nelem))
outfile.write('%.16e %.16e xlo xhi\n' % (xmin, xmax))
outfile.write('%.16e %.16e ylo yhi\n' % (ymin, ymax))
outfile.write('%.16e %.16e zlo zhi\n' % (zmin, zmax))
outfile.write('ITEM: ATOMS id type x y z vx vy vz \n')
for i in range(Nelem):
    # ITEM: ATOMS id type x y z vx vy vz
    outfile.write('%i %i %i %.16e %.16e %.16e %.16e %.16e %.16e\n' % ((i+1), atomtypes_all[i], molID[i], xpos[i], ypos[i], zpos[i], vxs[i], vys[i], vzs[i]))

# Write that step again too
for j in range(2000):
    k+=1
    outfile.write('ITEM: TIMESTEP\n%i\nITEM: NUMBER OF ATOMS\n%i\nITEM: BOX BOUNDS pp pp pp\n' % (k,Nelem))
    outfile.write('%.16e %.16e xlo xhi\n' % (xmin, xmax))
    outfile.write('%.16e %.16e ylo yhi\n' % (ymin, ymax))
    outfile.write('%.16e %.16e zlo zhi\n' % (zmin, zmax))
    outfile.write('ITEM: ATOMS id type x y z vx vy vz \n')
    for i in range(Nelem):
        # ITEM: ATOMS id type x y z vx vy vz
        outfile.write('%i %i %i %.16e %.16e %.16e %.16e %.16e %.16e\n' % ((i+1), atomtypes_all[i], molID[i], xpos[i], ypos[i], zpos[i], vxs[i], vys[i], vzs[i]))


# And the last command:
# Close file
outfile.close()



'''
outfile.write('LAMMPS data file via python scripting, version 11 Aug 2017, timestep = 0\n\n%i atoms \n%i bonds \n2 atom types\n1 bond types \n%i angles \n1 angle types\n\n' % (Nelem,Nbonds, Nangles))
outfile.write('%.16e %.16e xlo xhi\n' % (xmin, xmax))
outfile.write('%.16e %.16e ylo yhi\n' % (ymin, ymax))
outfile.write('%.16e %.16e zlo zhi\n\n' % (zmin, zmax))
outfile.write('Masses\n\n')
outfile.write('1 1\n2 1\n') # Should maybe automate this, but...

#outfile.write('Atoms # angle\n\n') # Original line
outfile.write('\nAtoms # full\n\n')

for i in range(Nelem):
    # For atom_type full:
    # atom-ID molecule-ID atom-type q x y z
    outfile.write('%i %i %i %i %.16e %.16e %.16e\n' % ((i+1), molID[i], atomtypes_all[i], qs[i], xpos[i], ypos[i], zpos[i]))

# Have the velocities here?
# Argh, the velocities
outfile.write('\nVelocities\n\n')
for i in range(Nelem):
    outfile.write('%i %.16e %.16e %.16e\n' % ((i+1), vxs[i], vys[i], vzs[i]))

# Have bonds here?

outfile.write('\nBonds\n\n')

counter = 0
for j in range(M):
    for i in range(N):
        counter += 1
        k = j*(N+1)+i+1
        outfile.write('%i 1 %i %i\n' % (counter, k, (k+1)))

outfile.write('\nAngles\n\n')

counter = 0
for j in range(M):
    for i in range(N-1):
        counter += 1
        k = j*(N+1)+i+1
        outfile.write('%i 1 %i %i %i\n' % (counter, k, (k+1), (k+2)))
'''

# And the last command:
# Close file
outfile.close()

