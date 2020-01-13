import numpy as np
import random
import sys
import math

## Placing the bead in the brush
# I should make a lot of files at once here.

# System info
spacing = 5
charge  = -1
T       = 3

# Number of files out, etc.
Nfiles = 21 # The number of files I had before
seed   = 23 # Or somethin' 

# Technicalities: Placement of bead
tol        = 0.2 # The closest we'll allow the bead to be to a brush 'atom'
substr_tol = tol # the closest we'll allow the bead to be to the substrate
distrange  = 1   # The spread of z-values of the bead. E.g. from 0.5 to 1.5. Will be from substr_tol to substr_tol+distrange

# 
Ncopylines = 36 # The number of lines we are going to read from the data.-infile and write to the data.-outfile # I will update this with the number of atoms later, so it should be sort of automatic

# Could I just write chunks of text to file without touching it? Use line in lines to write? I only really need the atom positions now.

vx = np.sqrt(T)
vy = np.sqrt(T)
vz = np.sqrt(T)

xran = np.zeros(Nfiles)
yran = np.zeros(Nfiles)
zran = np.zeros(Nfiles)

brush_beads = [] # We know where the substrate atoms are (at z=0), but we should know where the beads are so that we can avoid overlap

'''
# Line 0: 
LAMMPS data file via write_data, version 11 Aug 2017, timestep = 200000 

# Line 2 to 7 (enumeration starting at 0)
1134 atoms # And I'll add one... Can't just copy this line...
4 atom types
900 bonds
1 bond types
891 angles
1 angle types

# Line 9 to 11 (in 0-enumeration)
-2.5000000000000000e+00 1.2500000000000000e+01 xlo xhi
-2.5000000000000000e+00 1.2500000000000000e+01 ylo yhi
0.0000000000000000e+00 3.0000000000000000e+02 zlo zhi

'''

foldername = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Brush/Quadr_M9N101_ljunits_Langevin_scaled_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T%i_ljcut1p122/Spacing%i/Initial_configurations/' % (T,spacing)
#foldername = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Brush/Quadr_M9N101_ljunits_Langevin_scaled_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T%i_ljcut1p122/' % (T)
infilename = 'data.chaingrids_substrate_N909_Nchains9_Ly3_gridspacing%i_twofixed_charge%i_mass1_equilibrated' % (spacing, charge)
infile     = open(infilename, 'r')
lines      = infile.readlines()
Nlines     = len(lines)            # I will probably need this to write as many lines to file as necessary
infile.close()

N_atoms      = int(lines[2].split()[0]) # Line 2
N_atomtypes  = int(lines[3].split()[0]) # Line 3
N_bonds      = int(lines[4].split()[0]) # Line 4
N_bondtypes  = int(lines[5].split()[0]) # Line 5
N_angles     = int(lines[6].split()[0]) # Line 6
N_angletypes = int(lines[7].split()[0]) # Line 7


xmin = float(lines[9].split()[0]) 
xmax = float(lines[9].split()[1])
ymin = float(lines[10].split()[0])
ymax = float(lines[10].split()[1])
zmin = float(lines[9].split()[0])
zmax = float(lines[9].split()[1])

Ncopylines += N_atoms

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
    if len(words)>0 and words[0]=='Velocities':
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
        molID[i-1]     = int(words[1])
        this_atom_type = int(words[2])
        xposition      = float(words[4])
        yposition      = float(words[5])
        zposition      = float(words[6])
        atom_type[i-1] = this_atom_type
        qs[i-1]        = float(words[3])
        xpos[i-1]      = xposition
        ypos[i-1]      = yposition
        zpos[i-1]      = zposition
        # Check if this is a brush atom, then append it to the brush list if it is.
        if this_atom_type==1 or this_atom_type==2:
            brush_beads.append(np.array([xposition,yposition,zposition]))
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

beadpos = []

i = 0
while i<Nfiles: # I need to get the box edges here # USE CONTINUE??????
    breakit = False
    xran[i] = random.uniform(xmin,xmax)
    yran[i] = random.uniform(ymin,ymax)
    zran[i] = random.uniform(substr_tol,substr_tol+distrange) # Random placement not too far from the substrate.
    rran = np.array([xran[i],yran[i],zran[i]])
    for rbr in brush_beads:
        distvec = rbr-rran
        dist2   = np.dot(distvec,distvec)
        if dist2<tol:
            breakit = True
            break
    if breakit==True:
        continue
    else:
        beadpos.append(rran)
        i+=1

print('beadpos:',beadpos)

freebead_number = N_atoms+1
############# Write to file. Multiple times. #############
for i in range(Nfiles):
    filename = foldername + infilename + '_bead_file%i' % (i+1)
    outfile = open(filename, 'w')
    # Write lines from infile to the outfile.
    for j in range(Ncopylines+1):
        # if-test to update number of atoms?
        if len(lines[j].split())==2 and lines[j].split()[1]=='atoms':
            outfile.write('%i atoms\n' % (N_atoms+2))
        else:
            outfile.write(lines[j])
    outfile.write('%i %i 3 0 %.16e %.16e %.16e 0 0 0\n' % (freebead_number, max(molID)+1,xran[i],yran[i],zran[i])) # atom-ID molecule-ID atom-type q x y z # And the three last numbers I don't know. Flags? Counters of how many times they have crossed a border?
    for j in range(N_atoms+3):                        # Writing velocities to file.
        outfile.write(lines[j+Ncopylines+1]) # -1 ?
    # How to pick the velocity? Just choose one?
    outfile.write('%i %.16e %.16e %.16e\n' % (freebead_number,vx,vy,vz))
    for j in range(N_atoms+4+Ncopylines,Nlines):      # Write the rest of the lines to file. Everything from here on can be copy-pasted
        outfile.write(lines[j])
    outfile.close()                                   # This should be important since I write to multiple files
    # Write to meta file too?
