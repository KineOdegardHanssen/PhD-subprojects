from numpy import *
import random
import sys
import math

## Placing the bead in the brush
# I should make a lot of files at once here.

# System info
spacing = 5
charge  = -1

# Number of files out, etc.
Nfiles = 21 # The number of files I had before
seed   = 23 # Or somethin' 

# Could I just write chunks of text to file without touching it? Use line in lines to write? I only really need the atom positions now.

brush_beads = [] # We know where the substrate atoms are (at z=0), but we should know where the beads are so that we can avoid overlap

'''
# Line 0: 
LAMMPS data file via write_data, version 11 Aug 2017, timestep = 200000 

# Line 2 to 7 (enumeration starting at 0)
1134 atoms
4 atom types
900 bonds
1 bond types
891 angles
1 angle types
'''

infilename = 'data.chaingrids_substrate_N909_Nchains9_Ly3_gridspacing%i_twofixed_charge%i_mass1_equilibrated' % (spacing, charge)
infile     = open(infilename, 'r')
lines      = infile.readlines()

N_atoms      = float(lines[2].split()[0]) # Line 2
N_atomtypes  = float(lines[3].split()[0]) # Line 3
N_bonds      = float(lines[4].split()[0]) # Line 4
N_bondtypes  = float(lines[5].split()[0]) # Line 5
N_angles     = float(lines[6].split()[0]) # Line 6
N_angletypes = float(lines[7].split()[0]) # Line 7


#atom_nr   = []
molID     = np.zeros(N_atoms)
atom_type = np.zeros(N_atoms)
qs        = np.zeros(N_atoms)
xpos      = np.zeros(N_atoms)
ypos      = np.zeros(N_atoms)
zpos      = np.zeros(N_atoms)
bond_atom1 = [] # Could I perhaps just write lines directly to file?
bond_atom2 = []


# Just rewrite the first lines to file?


# ((i+1), molID[i], atomtypes_all[i], qs[i], xpos[i], ypos[i], zpos[i]))
# Have switches instead?
keyword = 'None'
for line in lines:
    words = line.split()
    '''
    if words[0]=='Atoms': # Do I need these?
        keyword = 'Atoms'
    '''
    if words[0]=='Velocities':
        break # I'll do this if I can write read lines directly to file.
        #keyword = 'Velocities'
    '''
    if words[0]=='Bonds':
        keyword = 'Bonds'
    if words[0]=='Angles':
        keyword = 'Angles'
    '''
    if len(words)==10:
        i = int(words[0])
        this_atom_type = int(words[2])
        molID[i-1]     = int(words[1])
        atom_type[i-1] = this_atom_type
        qs[i-1]        = int(words[3])
        xpos[i-1]      = int(words[4])
        ypos[i-1]      = int(words[5])
        zpos[i-1]      = int(words[6])
        # Check if this is a brush atom, then append it to the brush list if it is.
        if this_atom_type==1 or this_atom_type==2:
            brush_beads.append(i-1)
    '''
    if len(words)==4:
        if keyword=='Velocities':
        if keyword=='Bonds':
    if len(words)==5:
        if keyword=='Angles':
    '''
# ten elements in the lines that we want. Not that number of elements in any other line.

# atom style full: atom-ID molecule-ID atom-type q x y z

# I want to keep the information about the velocities and the bonds and angles too

### Part where I randomly place the diffusing bead.

for i in range(Nfiles): # I need to get the box edges here
    xran = random.uniform(xmin,xmax)
    yran = random.uniform(ymin,ymax)
    zran = random.uniform(1,3) # Random placement not too far from the substrate.


