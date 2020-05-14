import numpy as np
import random
import sys
import math

## Placing the bead in the brush
# I should make a lot of files at once here.

# System info
radius  = 2
spacing = 3
charge  = -1
T       = 3

# Number of files out, etc.
Nconfigs    = 100
Nplacements = 10
Nfiles = Nplacements # 21 # The number of files I had before
seed   = 23 # Or somethin' 

# Technicalities: Placement of bead
tol        = 0.2*radius # The closest we'll allow the bead to be to a brush 'atom'
substr_tol = 1.1*radius # the closest we'll allow the bead to be to the substrate # Previously: 0.7
distrange  = 0.5*radius  # The spread of z-values of the bead. E.g. from 0.5 to 1.5. Will be from substr_tol to substr_tol+distrange. # Previously: 0.5

#          # Is this not as set as I thought? It was 36 before:
startlines = 23 # The number of lines we are going to read from the data.-infile and write to the data.-outfile # I will update this with the number of atoms later, so it should be sort of automatic

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
infoldername  = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/Spacing%i/Initial_configs/Before_bead/' % spacing #'/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Initial_configurations/Spacing%i/' % spacing
outfoldername = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/Spacing%i/Radius' % spacing +str(radius)+'/Initial_configs/' 


for i in range(Nconfigs):
    outfilename = 'data.config%i_beadplacement' % (i+1)
    infilename  = infoldername+'data.config%i' % (i+1)
    infile      = open(infilename, 'r')
    lines       = infile.readlines()
    Nlines      = len(lines)            # I will probably need this to write as many lines to file as necessary
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
    
    #Ncopylines += N_atoms
    
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
        if len(words)>0 and words[0]=='Velocities':
            break # I'll do this if I can write read lines directly to file.
        
        if len(words)==10:
            j = int(words[0])
            molID[j-1]     = int(words[1])
            this_atom_type = int(words[2])
            xposition      = float(words[4])
            yposition      = float(words[5])
            zposition      = float(words[6])
            atom_type[j-1] = this_atom_type
            qs[j-1]        = float(words[3])
            xpos[j-1]      = xposition
            ypos[j-1]      = yposition
            zpos[j-1]      = zposition
            # Check if this is a brush atom, then append it to the brush list if it is.
            if this_atom_type==1 or this_atom_type==2:
                brush_beads.append(np.array([xposition,yposition,zposition]))
    # ten elements in the lines that we want. Not that number of elements in any other line.
    
    # atom style full: atom-ID molecule-ID atom-type q x y z
    
    # I want to keep the information about the velocities and the bonds and angles too
    
    ### Part where I randomly place the diffusing bead.
    
    beadpos = []
    
    j = 0
    while j<Nfiles: # I need to get the box edges here # USE CONTINUE??????
        breakit = False
        xran[j] = random.uniform(xmin,xmax)
        yran[j] = random.uniform(ymin,ymax)
        zran[j] = random.uniform(substr_tol,substr_tol+distrange) # Random placement not too far from the substrate.
        rran = np.array([xran[j],yran[j],zran[j]])
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
            j+=1
    
    print('beadpos:',beadpos)
    
    freebead_number = N_atoms+1
    ############# Write to file. Multiple times. #############
    for j in range(Nfiles):
        filename = outfoldername + outfilename + '%i' % (j+1)
        outfile = open(filename, 'w')
        # Write lines from infile to the outfile.
        for k in range(startlines-1):
            # if-test to update number of atoms?
            if len(lines[k].split())==2 and lines[k].split()[1]=='atoms':
                outfile.write('%i atoms\n' % (N_atoms+1))
            else:
                outfile.write(lines[k])
        for k in range(startlines-1,startlines+N_atoms-1):
            words = lines[k].split()
            #print('len(words):',len(words))
            if len(words)>0:
                outfile.write('%i %i %i %.16e %.16e %.16e %.16e\n' % (int(words[0]), int(words[1]), int(words[2]), float(words[3]), float(words[4]), float(words[5]), float(words[6])))
            else:    # I don't want an empty line between this one and the next
                #print('I broke.')
                break
        outfile.write('%i %i 3 0 %.16e %.16e %.16e\n' % (freebead_number, max(molID)+1,xran[j],yran[j],zran[j])) # atom-ID molecule-ID atom-type q x y z # And the three last numbers I don't know. Flags? Counters of how many times they have crossed a border?
        for k in range(N_atoms+3):                        # Writing velocities to file.
            outfile.write(lines[k+startlines+N_atoms-1])
        # How to pick the velocity? Just choose one?
        outfile.write('%i %.16e %.16e %.16e\n' % (freebead_number,vx,vy,vz))
        for k in range(2*N_atoms+4+startlines-2,Nlines):      # Write the rest of the lines to file. Everything from here on can be copy-pasted
            outfile.write(lines[k])
        outfile.close()                                   # This should be important since I write to multiple files
        # Write to meta file too?
