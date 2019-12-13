import matplotlib.pyplot as plt               # To plot
import numpy as np
import random
import math

# 21.03.19: Making a completely random chain for input to LAMMPS.
# I probably don't need a lot of these variables anymore, but this is not a bottleneck
# so I don't bother cleaning it up

#
onechain = False

# Set variables
N             = 100 # Number of bond vectors = number of units - 1 # Per chain
if onechain==True:
    M             = 1   # Number of chains
    L2            = 1   # Determining the shape of the L1xL2 polymer grid on the surface
    gridspacing   = 300 # Huuuge spacing in the grid in order to only study ONE chain (this will affect the box size)
else:
    M             = 9   # Number of chains
    L2            = 3   # Determining the shape of the L1xL2 polymer grid on the surface
    gridspacing   = 5   # The spacing in our grid # Default is 40
substratecharge = 0
charge          = 0
mass            = 1.5
Nelem           = M*(N+1)
Nall            = 1
xpos            = []
ypos            = []
zpos            = []

substrate_set = np.zeros((M,2))


for j in range(M):
    # Starting the chain
    x_accm = (j//L2)*gridspacing
    y_accm = j%L2*gridspacing
    z_accm = 0
    xpos.append(x_accm)
    ypos.append(y_accm)
    zpos.append(0)
    
### Data on system

# Kan evt. ha +- 0.5*gridspacing slik at det ser ut som om systemet er uniformt og uendelig... Best dersom systemet er ganske stort
xmax = max(xpos)+0.5*gridspacing
ymax = max(ypos)+0.5*gridspacing
zmax = N*3.0
zmin = 0
xmin = min(xpos)-0.5*gridspacing
ymin = min(ypos)-0.5*gridspacing


# Change names of atom-style, etc

### Print to file
outfilename = 'data.bulkdiffusion_gridspacing%i_mass%.1f' % (gridspacing,mass)
outfile = open(outfilename,'w')
outfile.write('LAMMPS data file via python scripting, version 11 Aug 2017, timestep = 0\n\n%i atoms\n1 atom types\n\n' % (Nall))
outfile.write('%.16e %.16e xlo xhi\n' % (xmin, xmax))
outfile.write('%.16e %.16e ylo yhi\n' % (ymin, ymax))
outfile.write('%.16e %.16e zlo zhi\n\n' % (zmin, zmax))
outfile.write('Masses\n\n')
outfile.write('1 %.5f\n' % (mass))

outfile.write('\nAtoms # atomic\n\n')

# atom-ID atom-type x y z
outfile.write('1 1 1.0 1.0 1.0\n')

outfile.write('\nVelocities\n\n')
outfile.write('1 1e-5 1e-5 0.01\n')

# Have bonds here?
'''
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

