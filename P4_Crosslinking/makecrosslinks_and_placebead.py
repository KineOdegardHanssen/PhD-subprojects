import numpy as np
import random
import sys
import math


bondcutoff = 1.7 # Too big?
bondcutoff2 = bondcutoff**2
bondfraction = 0.2

rng = np.random.default_rng()

## Placing the bead in the brush
# I should make a lot of files at once here.

# System info
spacing = 7
charge  = -1
T       = 3

print('d =',spacing)

confignrs = np.arange(1,3)
# Number of files out, etc.
Nfiles = 1000 # 21 # The number of files I had before

# Technicalities: Placement of bead
tol        = 0.2 # The closest we'll allow the bead to be to a brush 'atom'
substr_tol = 1.1 # the closest we'll allow the bead to be to the substrate # Previously: 0.7
distrange  = 0.5  # The spread of z-values of the bead. E.g. from 0.5 to 1.5. Will be from substr_tol to substr_tol+distrange. # Previously: 0.5

#          # Is this not as set as I thought? It was 36 before:
startlines = 23 # The number of lines we are going to read from the data.-infile and write to the data.-outfile # I will update this with the number of atoms later, so it should be sort of automatic

# Could I just write chunks of text to file without touching it? Use line in lines to write? I only really need the atom positions now.

vx = np.sqrt(T)
vy = np.sqrt(T)
vz = np.sqrt(T)

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
infoldername  = '/home/kine/Documents/Backup2_P2_PolymerMD/P3_Crosslinking/Spacing'+str(spacing)+'/Sigma_free_1/Initial_configurations/Before_bead/'
outfoldername = '/home/kine/Documents/Backup2_P2_PolymerMD/P3_Crosslinking/Spacing'+str(spacing)+'/Sigma_free_1/Initial_configurations/'
for confignr in confignrs:
    outfilename = 'data.crl_fb_confignr%i' % confignr
    infilename  = infoldername+'data.eq_nobead_config%i' % confignr
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
    
    Lx = xmax-xmin
    Ly = ymax-ymin
    Lz = zmax-zmin
    Lxd2 = 0.5*Lx
    
    pLx = np.array([Lx,0,0])
    mLx = np.array([-Lx,0,0])
    pLy = np.array([0,Ly,0])
    mLy = np.array([0,-Ly,0])
    pLxpLy = np.array([0,Lx,Lz])
    pLxmLy = np.array([0,Lx,-Lz])
    mLxpLy = np.array([0,-Lx,Lz])
    mLxmLy = np.array([0,-Lx,-Lz])
    
    molID     = np.zeros(N_atoms)
    atom_type = np.zeros(N_atoms)
    qs        = np.zeros(N_atoms)
    xpos      = np.zeros(N_atoms)
    ypos      = np.zeros(N_atoms)
    zpos      = np.zeros(N_atoms)
    bond_atom1 = [] # Could I perhaps just write lines directly to file?
    bond_atom2 = []
    
    brush_beads = np.zeros((909,3))
    brush_beads_matrix = np.zeros((9,99,3))
    linked_beads = np.zeros((9,99))
    atomnr_beads = np.zeros((9,99))
    chaincounter = np.zeros(9)
    brushbeadcounter = 0
    
    keyword = 'None'
    for line in lines:
        words = line.split()
        if len(words)>0 and words[0]=='Velocities':
            break # I'll do this if I can write read lines directly to file.
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
                brush_beads[brushbeadcounter,0]= xposition
                brush_beads[brushbeadcounter,1]= yposition
                brush_beads[brushbeadcounter,2]= zposition
                brushbeadcounter+=1
                if this_atom_type==2:
                    thismolID = int(molID[i-1])-1
                    #print('thismolID:',thismolID)
                    thiscc = int(chaincounter[thismolID])
                    #print('thiscc:',thiscc)
                    brush_beads_matrix[thismolID,thiscc,0] = xposition
                    brush_beads_matrix[thismolID,thiscc,1] = yposition
                    brush_beads_matrix[thismolID,thiscc,2] = zposition
                    atomnr_beads[thismolID,thiscc] = i-1
                    chaincounter[thismolID] += 1

    ### Part where I crosslink
    # Pick random chain
    # Pick random bead in that chain
    # See if it can crosslink to some bead on another chain
    #      # Cannot crosslink if they are too far apart
    #      # Cannot crosslink if one of them already has a crosslink <--- Need some extra array for this. Value 0 or 1 probably easiest. 0: No bond; 1: A bond.
    #      #                                                         <--- Check if linked_beads[i,j]=0 before doing ANYTHING.
    #      # Need to store info on extra bonds somewhere. Prob need a list. Append bondnr,bondtype,bead1,bead2(strictly speaking only the last, but these are things that'll enter the data.-file)
    # Will I stop at a certain number of crosslinks?
    # Will I stop when we have drawn X (non-crosslinked) beads? <---Probably this.
    #     # X should be quite large compared to the number Nmcb of mobile chain beads. Prob quite a bit smaller than Nmcb. 75%? Will fail to set a bond sometime
    # Should I shuffle chain numbers a bit? No bias whatsoever in which chain gets bonded? Probably.
    # bondcutoff
    
    # WHAT ABOUT THE PERIODIC IMAGES????
    newbonds = []
    Nmcb = 99*9 # Number of moving chain beads
    Nca  = int(math.ceil(bondfraction*Nmcb)) # Number of crosslinking attempts
    beadchains = np.arange(0,9)
    #print('beadchains:',beadchains)
    for i in range(Nca):
        rchain = random.randint(0,8)# Random chain
        rbead  = random.randint(0,98)
        if linked_beads[rchain,rbead]==1:
            i-=1
            continue
        rng.shuffle(beadchains) # Shuffling
        success = False
        for j in range(9):
            beadchainj =beadchains[j]
            if beadchainj!=rchain: # Will not crosslink a chain to itself:
                for k in range(99):
                    if linked_beads[beadchainj,k]==0 and success==False:
                        diff = brush_beads_matrix[rchain,rbead,:]-brush_beads_matrix[beadchainj,k,:]
                        dist2 = np.dot(diff,diff)
                        if dist2<bondcutoff2:
                            success = True
                        elif dist2>Lxd2:
                            # Periodic images
                            diff_mx1 = brush_beads_matrix[rchain,rbead,:]-brush_beads_matrix[beadchainj,k,:]-pLx    
                            diff_mx2 = brush_beads_matrix[rchain,rbead,:]-brush_beads_matrix[beadchainj,k,:]-mLx
                            diff_my1 = brush_beads_matrix[rchain,rbead,:]-brush_beads_matrix[beadchainj,k,:]-pLy
                            diff_my2 = brush_beads_matrix[rchain,rbead,:]-brush_beads_matrix[beadchainj,k,:]-mLy    
                            diff_mxy1 = brush_beads_matrix[rchain,rbead,:]-brush_beads_matrix[beadchainj,k,:]-pLxpLy
                            diff_mxy2 = brush_beads_matrix[rchain,rbead,:]-brush_beads_matrix[beadchainj,k,:]-mLxpLy
                            diff_mxy3 = brush_beads_matrix[rchain,rbead,:]-brush_beads_matrix[beadchainj,k,:]-pLxmLy
                            diff_mxy4 = brush_beads_matrix[rchain,rbead,:]-brush_beads_matrix[beadchainj,k,:]-mLxmLy
                            dist2_mx1 = np.dot(diff_mx1,diff_mx1)    
                            dist2_mx2 = np.dot(diff_mx1,diff_mx2)    
                            dist2_my1 = np.dot(diff_my1,diff_my1)
                            dist2_my2 = np.dot(diff_my1,diff_my2)
                            dist2_mxy1 = np.dot(diff_mxy1,diff_mxy1)
                            dist2_mxy2 = np.dot(diff_mxy2,diff_mxy2)
                            dist2_mxy3 = np.dot(diff_mxy3,diff_mxy3)    
                            dist2_mxy4 = np.dot(diff_mxy4,diff_mxy4)
                            if dist2_mx1<bondcutoff2:
                                success = True
                            elif dist2_mx2<bondcutoff2:
                                success = True
                            elif dist2_my1<bondcutoff2:
                                success = True
                            elif dist2_my2<bondcutoff2:
                                success = True    
                            elif dist2_mxy1<bondcutoff2:
                                success = True    
                            elif dist2_mxy2<bondcutoff2:
                                success = True
                            elif dist2_mxy3<bondcutoff2:
                                success = True
                            elif dist2_mxy4<bondcutoff2:
                                success = True    
                        if success==True:
                            linked_beads[rchain,rbead]=1
                            linked_beads[beadchainj,k]=1
                            bead1 = atomnr_beads[rchain,rbead]
                            bead2 = atomnr_beads[beadchainj,k]
                            thebeads = [bead1,bead2]
                            firstbead = min(thebeads)
                            secondbead = max(thebeads)
                            newbonds.append([firstbead, secondbead])
                            break
    
    # Print number of crosslinking bonds created
    print('Number of crosslinking bonds:',len(newbonds))
    
    ### Part where I randomly place the diffusing bead.
    
    i = 0
    beadok = False
    while beadok==False: # I need to get the box edges here # USE CONTINUE??????
        breakit = False
        xran = random.uniform(xmin,xmax)
        yran = random.uniform(ymin,ymax)
        zran = random.uniform(substr_tol,substr_tol+distrange) # Random placement not too far from the substrate.
        rran = np.array([xran,yran,zran])
        for rbr in brush_beads:
            distvec = rbr-rran   
            dist2   = np.dot(distvec,distvec)
            if dist2<tol:
                breakit = True
                break   
        if breakit==True:
            continue
        else:
            beadpos = rran
            beadok = True
    print('beadpos:',beadpos)
    
    freebead_number = N_atoms+1
    ############# Write to file. #############
    filename = outfoldername + outfilename
    outfile = open(filename, 'w')
    # Write lines from infile to the outfile.
    for j in range(startlines-1):
        # if-test to update number of atoms?    
        if len(lines[j].split())==2 and lines[j].split()[1]=='atoms':
            outfile.write('%i atoms\n' % (N_atoms+1))
        elif len(lines[j].split())==2 and lines[j].split()[1]=='bonds':
            outfile.write('%i bonds\n' % (N_bonds+len(newbonds)))
        else:
            outfile.write(lines[j])
    for j in range(startlines-1,startlines+N_atoms-1):
        words = lines[j].split()    
        #print('len(words):',len(words))
        if len(words)>0:
            outfile.write('%i %i %i %.16e %.16e %.16e %.16e\n' % (int(words[0]), int(words[1]), int(words[2]), float(words[3]), float(words[4]), float(words[5]), float(words[6])))
        else:    # I don't want an empty line between this one and the next    
            #print('I broke.')
            break
    outfile.write('%i %i 3 0 %.16e %.16e %.16e\n' % (freebead_number, max(molID)+1,beadpos[0],beadpos[1],beadpos[2])) # atom-ID molecule-ID atom-type q x y z # And the three last numbers I don't know. Flags? Counters of how many times they have crossed a border?
    for j in range(N_atoms+3):                        # Writing velocities to file.
        outfile.write(lines[j+startlines+N_atoms-1])
    # How to pick the velocity? Just choose one?
    outfile.write('%i %.16e %.16e %.16e\n' % (freebead_number,vx,vy,vz))
    for j in range(2*N_atoms+4+startlines-2,2*N_atoms+4+startlines+1+N_bonds):      # Write original bonds to file
        outfile.write(lines[j])
    for j in range(len(newbonds)):                                                 # Write new bonds to file
        currbond = newbonds[j]
        outfile.write('%i 2 %i %i\n' % (N_bonds+j+1,currbond[0],currbond[1]))
    for j in range(2*N_atoms+4+startlines+1+N_bonds,Nlines): # Assuming no angle constraint on new bonds, write the rest of the lines to file.
        outfile.write(lines[j])
    outfile.close()                                   # This should be important since I write to multiple files
