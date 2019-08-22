import matplotlib.pyplot as plt               # To plot
import numpy as np
import random
import math

#outfilename = 'data.hexgrid_'  # Actually, I might move this a bit further down...

mass        = 0.00081
charge      = -1.9 
gridspacing = 10
### Get tethering points:
infilename = 'data.hexgrid_box2x2_gridspacing%i_mass0.00081' % gridspacing
infile     = open(infilename, 'r')

lines = infile.readlines()

xes = []
ys  = []

start   = 0
M       = 0    # Number of chains
N       = 101  # Number of atoms in chain
for line in lines:
    words = line.split()
    #print(words)
    if len(words)>0:
        if words[0]=='Velocities':
            break
        if start>0 and len(words)>0:
            # Update number of chains:
            M+=1
            # Append coordinates:
            xes.append(float(words[2]))
            ys.append(float(words[3]))
        if words[0]=='Atoms':
            start = 1

# Setting the number of elements, bonds and angles
Nelem   = M*N
Nbonds  = M*(N-1)
Nangles = M*(N-2)

# Make a straight chain of 101 atoms for each tethering point:

n = 0 # Atom number
# Should make some neighbour list as well...     # I can easily deal with the neighbour list later...
# I should make a matrix for the positions...    # Or should I unpack them? I guess that is just as easy... Have a running counter
xpos          = np.zeros(Nelem)
ypos          = np.zeros(Nelem)
zpos          = np.zeros(Nelem)
molID         = np.zeros(Nelem)
atomtypes_all = np.zeros(Nelem)                  # The fixed monomers are type 1. The moving ones are of type 2.s
molID         = np.zeros(Nelem)
qs            = np.zeros(Nelem)
vxs           = np.zeros(Nelem)
vys           = np.zeros(Nelem)
vzs           = np.zeros(Nelem)

counter = 0
for i in range(M):
    # Tethering point
    zcoord = 0
    x = xes[i]
    y = ys[i]
    xpos[counter]  = x
    ypos[counter]  = y
    zpos[counter]  = 0
    molID[counter] = i+1
    qs[counter]    = charge
    atomtypes_all[counter] = 1
    counter+=1
    zcoord +=1
    for j in range(N-1):
        xpos[counter]  = x
        ypos[counter]  = y
        zpos[counter]  = zcoord
        molID[counter] = i+1
        qs[counter]    = charge
        if j==0:
            atomtypes_all[counter] = 1
        else:
            atomtypes_all[counter] = 2
        counter+=1
        zcoord +=1


######################################################################## Dette er lagt til, men ikke fikset ##################################################################################
# Kan evt. ha +- 0.5*gridspacing slik at det ser ut som om systemet er uniformt og uendelig... Best dersom systemet er ganske stort
xmax = max(xpos)+1.0*gridspacing # Var 0.5*gridspacing
ymax = max(ypos)+1.0*gridspacing
zmax = max(zpos)*3.0
zmin = min(zpos)*3.0 # Change something here...
xmin = min(xpos)-1.0*gridspacing
ymin = min(ypos)-1.0*gridspacing

# Change names of atom-style, etc

### Print to file
#outfilename = 'data.randomchain_N%i' % N+1
outfilename = 'data.hexatethered_N%i_Nchains%i_gridspacing%i_twofixed_charge%.1f_mass%.5f' % (Nelem,M,gridspacing,charge,mass)
outfile = open(outfilename,'w')
outfile.write('LAMMPS data file via python scripting, version 11 Aug 2017, timestep = 0\n\n%i atoms \n%i bonds \n2 atom types\n1 bond types \n%i angles \n1 angle types\n\n' % (Nelem,Nbonds, Nangles))
outfile.write('%.16e %.16e xlo xhi\n' % (xmin, xmax))
outfile.write('%.16e %.16e ylo yhi\n' % (ymin, ymax))
outfile.write('%.16e %.16e zlo zhi\n\n' % (zmin, zmax))
outfile.write('Masses\n\n')
outfile.write('1 %.5f\n2 %.5f\n' % (mass,mass)) # Should maybe automate this, but...

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
    for i in range(N-1):
        counter += 1
        k = j*N+i+1
        outfile.write('%i 1 %i %i\n' % (counter, k, (k+1)))

outfile.write('\nAngles\n\n')

counter = 0
for j in range(M):
    for i in range(N-2):
        counter += 1
        k = j*N+i+1
        outfile.write('%i 1 %i %i %i\n' % (counter, k, (k+1), (k+2)))


# And the last command:
# Close file
outfile.close()



 
