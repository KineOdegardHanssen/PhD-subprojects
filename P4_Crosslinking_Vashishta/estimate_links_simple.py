import numpy as np
import random
import sys
import math


# System info
spacing = 3
radius  = 1
charge  = -1
T       = 3

# Number of files out, etc.
Nfiles = 1000
Ncr    = 5
pad    = 10
Ncrfr  = 50 # Or scale to system size

cutoff2 = 1.0 # Possibly change

confignrs = [1,2,3]

for confignr in confignrs:
    infilename = '/home/kine/Documents/P4_Crosslinking_Vashishta/Spacing'+str(spacing)+'/Sigma_free_'+str(radius)+'/Charge%i/all_confignr%i_TEST.lammpstrj' % (charge,confignr)
    infile = open(infilename,'r')
    
    lines = infile.readlines() # This takes some time
    # Getting the number of lines, etc.
    totlines = len(lines)         # Total number of lines
    lineend = totlines-1          # Index of last element
    
    # Extracting the number of atoms:
    words = lines[3].split()
    Nall = int(words[0])
    N    = Nall
    
    skiplines   = 9             # If we hit 'ITEM:', skip this many steps...
    skipelem    = 0
    sampleevery = 0
    i           = int(math.ceil(skipelem*(Nall+9)))
    skiplines  += (Nall+skiplines)*sampleevery # Check!
    
    # Setting arrays for treatment:
    
    # Setting arrays for treatment:
    positions_crlinker = [] # Fill with xu,yu,zu-coordinates. One position for each time step. # We have the same time step for the whole of the simulation. Can bother with time later
    positions_linksite = []
    positions_crlinker_this = []
    positions_linksite_this = []
    times     = [] # Not useful anymore?
        
    # For now: Only find the largest z-position among the beads in the chain. # Otherwise I must read off the atom number (will do that with the free bead anyways.)
    maxz = -1000 # No z-position is this small.
    
    crlinkcounter = 0
    linksitecounter = 0
    while i<totlines:
        words = lines[i].split()
        if (words[0]=='ITEM:'):
            if words[1]=='TIMESTEP':
                words2 = lines[i+1].split() # The time step is on the next line
                t = float(words2[0])
                times.append(t)
                positions_crlinker.append(positions_crlinker_this) 
                positions_linksite.append(positions_linksite_this)
                positions_crlinker_this = []
                positions_linksite_this = []
                i+=skiplines
            elif words[1]=='NUMBER': # These will never kick in. 
                i+=1
            elif words[1]=='BOX':
                i+=1
            elif words[1]=='ATOMS':
                i+=1
        elif len(words)<8:
            i+=1
        else:
            # Find properties
            # Order:  id  type mol ux  uy  uz  vx  vy   vz
            #         [0] [1]  [2] [3] [4] [5] [6] [7]  [8]
            atomtype = int(words[1]) 
            #molID    = int(words[2])
            x        = float(words[3])
            y        = float(words[4])
            z        = float(words[5])
            if atomtype==5:
                positions_linksite_this.append(np.array([x,y,z]))
                linksitecounter+=1
            elif atomtype==6:
                positions_crlinker_this.append(np.array([x,y,z]))
                crlinkcounter+=1
            i+=1
    infile.close()
    
    print('len(times):',len(times))
    print('len(positions_crlinker):',len(positions_crlinker))
    print('len(positions_linksite):',len(positions_linksite))
    print('len(positions_crlinker[0]):',len(positions_crlinker[0]))
    print('len(positions_linksite[0]):',len(positions_linksite[0]))
    
    Nt = len(times)
    Nlinks = np.zeros(Nt-1)
    dotprods = []
    for i in range(1,Nt):
        thistime = times[i]
        positions_linksite_this = positions_linksite[i]
        positions_crlinker_this = positions_crlinker[i]
        for j in range(Ncrfr):
            thisr = positions_crlinker_this[j]
            team = 0
            for k in range(Ncr*9):
                thislinkr = positions_linksite_this[k]
                diff = thisr-thislinkr
                dotprod = np.dot(diff,diff)
                dotprods.append(dotprod)
                if dotprod<cutoff2:
                    team+=1
            if team>=2:
                Nlinks[i]+=1

    avgNlinks = np.mean(Nlinks)
    dotprods = np.array(dotprods)
    print('avg(dotprods):',np.mean(dotprods), '; min(dotprods):',min(dotprods), '; max(dotprods):',max(dotprods))
    print('File ', confignr, '; avg number of links per time frame:', avgNlinks)
    print('Nlinks:', Nlinks)
