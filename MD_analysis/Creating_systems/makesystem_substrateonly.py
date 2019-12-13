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
    gridspacing   = 100   # The spacing in our grid # Default is 40
mass            = 1.5
smass           = 1.0
xpos            = []
ypos            = []
zpos            = []
atomtypes_all   = []
vxs             = []
vys             = []
vzs             = []

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
xmin = min(xpos)-0.5*gridspacing
ymin = min(ypos)-0.5*gridspacing
zmax = N*3.0
zmin = 0

xpos            = []
ypos            = []
zpos            = []

atomtypes_all.append(1)
xpos.append(1)
ypos.append(1)
zpos.append(1)
vxs.append(1e-5)
vys.append(1e-5)
vzs.append(0.01)

# Now I know how big the system is. Next: Setting the substrate.
# Which atom type? 4, since I use 3 for the unbound bead?
d    = 1 # One nm (=sigma? in between the beads? Or less? I guess an atom could get pretty close to the substrate if it went in at the right angle/orientation...)
Lx   = xmax-xmin
Ly   = ymax-ymin
Lz   = zmax-zmin # Don't need this now...
Nx   = int(math.floor(Lx/d))
Ny   = int(math.floor(Ly/d))
Nall = 1
for i in range(Nx):
    for j in range(Ny):
        x = i*d+xmin
        y = j*d+ymin
        atomtypes_all.append(2)
        xpos.append(x)
        ypos.append(y)
        zpos.append(0)
        vxs.append(np.random.normal(0,1))
        vys.append(np.random.normal(0,1))
        vzs.append(np.random.normal(0,1))
        Nall += 1
        

# Change names of atom-style, etc

### Print to file

outfilename = 'data.bulkdiffusion_withsubstrate_gridspacing%i_mass%.1f' % (gridspacing,mass)
outfile = open(outfilename,'w')
outfile.write('LAMMPS data file via python scripting, version 11 Aug 2017, timestep = 0\n\n%i atoms \n2 atom types\n\n' % Nall)
outfile.write('%.16e %.16e xlo xhi\n' % (xmin, xmax))
outfile.write('%.16e %.16e ylo yhi\n' % (ymin, ymax))
outfile.write('%.16e %.16e zlo zhi\n\n' % (zmin, zmax))
outfile.write('Masses\n\n')
outfile.write('1 %.5f\n2 %.5f\n' % (mass,smass))

#outfile.write('Atoms # angle\n\n') # Original line
outfile.write('\nAtoms # atomic\n\n')

for i in range(Nall):
    # For atom_type full:
    # atom-ID molecule-ID atom-type q x y z
    outfile.write('%i %i %.16e %.16e %.16e\n' % ((i+1), atomtypes_all[i], xpos[i], ypos[i], zpos[i]))

# Have the velocities here?
# Argh, the velocities
outfile.write('\nVelocities\n\n')
for i in range(Nall):
    outfile.write('%i %.16e %.16e %.16e\n' % ((i+1), vxs[i], vys[i], vzs[i]))

'''
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

