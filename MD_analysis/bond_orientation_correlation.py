import matplotlib.pyplot as plt               # To plot
import numpy as np
import random
import math
import time

def find_rms(inarray):
    rms     = 0
    N       = len(inarray)
    #print("inarray:", inarray)
    average = np.mean(inarray)
    for k in range(N):
        rms += (inarray[k]-average)**2 
    rms = np.sqrt(rms/float(N))
    return average, rms

## Set parameters:
Lx = 3
Ly = 3
d  = 2 # Dim.

gridspacing = 2
spacing     = gridspacing
T           = 310
M           = Lx*Ly
N           = 101
Nb          = N-1
ljdebye     = 1.042
epsilon     = ljdebye
sigma       = 1
ljcutoff    = 1.12246204830937
debyecutoff = 3
factor      = 0.1#250
Kbond       = 200#Kangle*factor
Kangle      = Kbond*factor
charge      = -1
K           = Kangle # Because we used this notation earlier, but using Kangle is less confusing
oneframe    = 2
wallenergy  = 1.042
dielectric  = 1000

# The LAMMPS output files contains a some lines with details the system size, number of atom types etc.. We want to access the information that comes after that.
linestart           = 6 #17            # The line the data starts at.
printeverynthstep   = 100
dt                  = 0.00045       # The time step of our simulation. 0.00045 ns default for nano
skiplines           = 9             # If we hit 'ITEM:', skip this many steps...
skipelem            = 0#10#1000#10000#10000#90000 # The number of elements we skip in order to equilibrate (set to 0 if the .lammpstrj file should be equilibrated)
sampleevery         = 0#1 # Sample every 10th of the values written to file # The program gets WAY too slow if this is too small.
timefac             = dt*printeverynthstep*1e-9*sampleevery

## Extracting neighbours:
infilename = 'chainneighbourinfo_%ix%i.txt' % (Lx,Ly)
infile     = open(infilename,'r')

# Extracting length of distance vector:
header = infile.readline()
hwords = header.split()
length = int(hwords[0])

# Preparing correlationvector:
bond_corr     = np.zeros(length)
bond_corr_rms = np.zeros(length)
distances     = np.zeros(length)

# Making array for pairs:
if d==2:
    maxelements = 2*Lx*Ly-(Lx+Ly)
if d==1:
    maxelements = Lx-1
neighbours  = np.zeros((length,maxelements,2))

# Reading data:
lines   = infile.readlines()
ind     = 0
counter = 0

for line in lines:
    words = line.split()
    print(words)
    if len(words)>0:
        if words[0]=='PAIRS:':
            # Extract index of neighbours:
            ind = int(words[2])
            distances[ind] = np.sqrt(int(words[4])**2+int(words[5])**2)
            counter = 0
        else:
            # Add into arrays:
            if words[0]!='#':
                neighbours[ind][counter][0] = int(words[0])
                neighbours[ind][counter][1] = int(words[1])
                counter += 1
infile.close()            
print(neighbours[0])     # All the nearest-neighbour pairs
print(neighbours[0,:,0]) # All the lowest indices of the nearest neighbour pairs
print(neighbours[0,1])   # The second nearest neighbour pair in the list
print(neighbours[0,-1])  # The last nearest neighbour pair in the list

# Include some test to chech if the pair indices are 0 and 0. Maybe perform a break if they are.

#################################################################################################################################################################################################

# Function for curve_fit
def decaying_exponential(s,P):
    return np.exp(-s/P)

start_time = time.process_time()

### Input and output file names
#	     #  chaingrid_quadratic_M9N101_Langevin_Kangle100_Kbond2000_factor0.05_T310_theta0is180_twofirst_are_fixed
'''
infilename   = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge-1_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,spacing,Kangle,Kbond,T)
plotname     = 'cos_bond_correlation_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge-1_T%i_theta0is180_twofirst_are_fixed_LONG.png' % (M,N,spacing,Kangle,Kbond,T)
outfilename  = 'cos_bond_correlation_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge-1_T%i_theta0is180_twofirst_are_fixed_LONG.txt' % (M,N,spacing,Kangle,Kbond,T)
outfilename2 = 'cos_bond_correlation_oneframe%i_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge-1_T%i_theta0is180_twofirst_are_fixed_LONG.txt' % (oneframe,M,N,spacing,Kangle,Kbond,T)
outfilename3 ='cos_bond_correlation_nn_bybondnr_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge-1_T%i_theta0is180_twofirst_are_fixed_LONG.txt' % (M,N,spacing,Kangle,Kbond,T)
#'''


# Tester:
testtype     = 'rotatingarc'#'rotatingarc'#'rotatingarc_withtime'#'quartercirc'#'straight'
infilename   = 'testsystem_chaingrid_'+testtype+'.lammpstrj'
plotname     = 'cos_bond_correlation_testsystem_chaingrid_'+testtype+'.png'
outfilename  = 'cos_bond_correlation_testsystem_chaingrid_'+testtype+'.txt'
outfilename2 = 'cos_bond_correlation_oneframe%i_testsystem_chaingrid_' % oneframe +testtype+'.txt'
outfilename3 ='cos_bond_correlation_nn_bybondnr_testsystem_chaingrid_'+testtype+'.txt'

''' # Varying the dielectric constant:
infilename   = 'chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
plotname     = 'cos_bond_correlation_chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
outfilename  = 'cos_bond_correlation_chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
outfilename2 = 'cos_bond_correlation_oneframe%i_chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (oneframe,M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
outfilename3 ='cos_bond_correlation_nn_bybondnr_chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
#'''

''' # With a wall potential:
infilename   = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
plotname     = 'cos_bond_correlation_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
outfilename  = 'cos_bond_correlation_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
outfilename2 = 'cos_bond_correlation_oneframe%i_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (oneframe,M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
outfilename3 ='cos_bond_correlation_nn_bybondnr_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)

#'''


'''
infilename   = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,Kangle,Kbond,factor,T)
plotname     = 'cos_bond_correlation_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,spacing,Kangle,Kbond,factor,T)
(M,N,spacing,Kangle,Kbond,factor,T)
outfilename  = 'cos_bond_correlation_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,Kangle,Kbond,factor,T)
#'''
# No problem with bond lengths here...

# Varying the grid spacing # THIS IS NOW THE STANDARD FILE NAMES.
'''
infilename   = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,spacing,Kangle,Kbond,charge,T)
plotname     = 'bond_correlation_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,spacing,Kangle,Kbond,charge,T)
outfilename  = 'bond_correlation_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,Kangle,Kbond,charge,T)'''

# Testing:
'''
#infilename  = 'chaingrids_totallystraight_N909_Nchains9_Ly3_twofixed_test.lammpstrj'
infilename  = 'chaingrids_totallystraight_N404_Nchains4_Ly1_twofixed_test.lammpstrj'
plotname    = 'test.png'
outfilename = 'test.txt'
#'''


# Varying the charge:
#chaingrid_quadratic_M9N101_Langevin_Kangle${Kangle}_Kbond${Kbond}_charge${charge}_T$T_theta0is180_twofirst_are_fixed.lammpstrj
''' # String appending method
#               chaingrid_quadratic_M9N101_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge0_T310_theta0is180_twofirst_are_fixed
infilename   = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % T     # Actual run
plotname     = 'bond_correlation_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.png' % T
outfilename  = 'bond_correlation_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
#'''

''' # String appending method
infilename   = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % T     # Actual run
plotname     = 'bond_correlation_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.png' % T
outfilename  = 'bond_lengths_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
#'''

''' # Sprintf method
infilename   = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,Kangle,Kbond,factor,T)     # Actual run
plotname     = 'bond_correlation_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,Kangle,Kbond,factor,T)
outfilename  = 'bond_correlation_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,Kangle,Kbond,factor,T)
#'''

'''
infilename   = 'chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.lammpstrj' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)     # Actual run
plotname     = 'bond_correlation_chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.png' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)
outfilename  = 'bond_correlation_chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.txt' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)
#'''

'''
infilename   = 'chaingrid_quadratic_M'+str(M)+'N'+str(N)+'_ljdebye_epsilon'+str(epsilon)+'_sigma'+str(sigma)+'_ljcutoff'+str(ljcutoff)+'_kappa1_debyecutoff'+str(debyecutoff)+'_Langevin_wall'+str(ljdebye)+'_Kangle'+str(Kangle)+'_Kbond'+str(Kbond)+'_T'+str(T)+'_theta0is180_twofirst_are_fixed.lammpstrj'     # Actual run
plotname     = 'bond_correlation_chaingrid_quadratic_M'+str(M)+'N'+str(N)+'_ljdebye_epsilon'+str(epsilon)+'_sigma'+str(sigma)+'_ljcutoff'+str(ljcutoff)+'_kappa1_debyecutoff'+str(debyecutoff)+'_Langevin_wall'+str(ljdebye)+'_Kangle'+str(Kangle)+'_Kbond'+str(Kbond)+'_T'+str(T)+'_theta0is180_twofirst_are_fixed.png'
outfilename  = 'bond_correlation_chaingrid_quadratic_M'+str(M)+'N'+str(N)+'_ljdebye_epsilon'+str(epsilon)+'_sigma'+str(sigma)+'_ljcutoff'+str(ljcutoff)+'_kappa1_debyecutoff'+str(debyecutoff)+'_Langevin_wall'+str(ljdebye)+'_Kangle'+str(Kangle)+'_Kbond'+str(Kbond)+'_T'+str(T)+'_theta0is180_twofirst_are_fixed.txt'
'''


### Opening files
outfile  = open(outfilename,'w')
outfile2 = open(outfilename2,'w')
outfile3 = open(outfilename3,'w')

#### Automatic part
infile = open(infilename, "r")
lines = infile.readlines()
# Getting the number of lines, etc.
totlines = len(lines)         # Total number of lines
lineend = totlines-1          # Index of last element

# Extracting the number of atoms:
words = lines[3].split()
print("words:", words)
print("words[0]:", words[0])
Nall = int(words[0])

words = lines[5].split()
xmin = float(words[0])
xmax = float(words[1])

words = lines[6].split()
ymin = float(words[0])
ymax = float(words[1])

words = lines[7].split()
zmin = float(words[0])
zmax = float(words[1])

Lx = xmax-xmin
Ly = ymax-ymin
Lz = zmax-zmin

halfLx = 0.5*Lx
halfLy = 0.5*Ly
halfLz = 0.5*Lz


## Making an array of indices:
# Not sure if this is really the most efficient way of doing things
# Probably is if I have a lot of time steps
index1 = np.zeros(M*N) # The chain index
index2 = np.zeros(M*N) # The atom (in chain) index
counter = 0
for i in range(M):
    for j in range(N):
        index1[counter] = i
        index2[counter] = j
        counter += 1

'''
for i in range(M*N):
    print("Atomid-1 =", i, ": chain:", index1[i], ", atom nr.", index2[i])
'''

xes = []
ys  = []
zs  = []

x_curr = np.zeros((M,N))
y_curr = np.zeros((M,N))
z_curr = np.zeros((M,N))

times = []

print("infilename:", infilename)

#skipelem = 0#10000#10000#10000#90000
#totlines = N+9+9
# Extracting data
i = int(math.ceil(skipelem*(Nall+9)))
j = 0
testcounter = 0
#totlines = 3*(N+skiplines)
skiplines += (Nall+skiplines)*sampleevery # Check!
while i<totlines:
    words = lines[i].split()
    if i==int(math.ceil(skipelem*(N+9))):
        print("In loop, words[0]=", words[0])
        print("First line I read:", words)
    if (words[0]=='ITEM:' and words[1]=='TIMESTEP'):
        if words[1]=='TIMESTEP':
            words2 = lines[i+1].split() # The time step is on the next line
            t = float(words2[0])
            times.append(t)
            #print("New time step")
            #print("Loop time step:", j)
            if t!=0:
                #print("x_curr:",x_curr)
                #print("No. of atoms processed in time step:", counter)
                counter=0
                xes.append(x_curr)
                ys.append(y_curr)
                zs.append(z_curr)
                # I don't think I ACTUALLY need to reset these, but it's probably no problem:
                x_curr    = np.zeros((M,N))
                y_curr    = np.zeros((M,N))
                z_curr    = np.zeros((M,N))
                j+=1
            i+=skiplines
            #testcounter += 1
            #if(testcounter==10):
            #    i=totlines
            #print("Next i:", i)
        elif words[1]=='NUMBER':
            i+=7
        elif words[1]=='BOX':
            i+=5
        elif words[1]=='ATOMS':
            i+=1
    elif len(words)<8:
        i+=1
    else:
        # Find properties
        # Order:  id  type mol  x   y   z  vx   vy   vz
        #         [0] [1]  [2] [3] [4] [5] [6] [7]  [8]
        #print(words)
        # Extracting:
        #print("words:",words)
        #print("words[0]:",words[0])
        ind   = int(words[0])-1 # Atom ids go from zero to N-1.
        x     = float(words[3])
        y     = float(words[4])
        z     = float(words[5])
        # Making array indices:
        #print("ind:", ind)
        ind1  = int(index1[ind])
        ind2  = int(index2[ind])
        # Add to lists
        x_curr[ind1,ind2] = x # Does this work? Check...
        y_curr[ind1,ind2] = y
        z_curr[ind1,ind2] = z
        i+=1 #Should this be in the else-block?
        counter+=1
        '''
        if(ind==749):
            print("ind=749")
        if(ind==2):
            print("ind=2")
        if(ind==5):
            print("ind=5")
        if(ind==108):
            print("ind=108")
        if(ind==75):
            print("ind=75")
        if(ind==302):
            print("ind=302")
        '''
        #if z>200:
        #    print("z:", z, "; index:", ind, "frame:", t/100.)
        #    break
        #print("In extraction part")

# Only if I have ONE short test file:
#xes.append(x_curr)
#ys.append(y_curr)
#zs.append(z_curr)

infile.close()

#print("xes:",xes)
#print("xes[0]:", xes[0])
#print("max(xes):", max(xes))

# SEEMS to work fine up until here.

# Finding the bond vectors for each time step.
# This will be an awfully blocky code, it seems
# Removing the first redundant element of position lists (uncomment if short test-file)
xes.pop(0)
ys.pop(0)
zs.pop(0)

Nt = len(xes) # Could have used j (Nt = j)

#print("xes:", xes)

#print("xes[-1]:",xes[-1])
#print("ys[-1]:",ys[-1])
#print("zs[-1]:",zs[-1])
#print("len, xes:", Nt)

#print("xes[0]:", xes[0])
#print("ys[0]:", ys[0])
print("zs[0]:", zs[0])
#print("xes:", xes)

# Oh, vey.
bondvecs        = np.zeros((Nt, M, N-1, 3))
length_bondvecs = np.zeros((Nt, M, N-1))

# Crap. I need another loop and even more indices. Yikes.
bindex1 = 0
bindex2 = 0
for i in range(Nt):    # Loop over time steps
    for k in range(M): # Loop over chains
        xes_thistime = xes[i] # [i,k] # I think I need to extract this twice, because it is a list of arrays (a double one at that...)
        ys_thistime  = ys[i]
        zs_thistime  = zs[i]
        #print("len(xes_this):",len(xes_thistime))
        #print(xes_thistime)
        #print("shape(xes_this):",np.size(xes_this))
        xes_this     = xes_thistime[k,:]
        ys_this      = ys_thistime[k,:]
        zs_this      = zs_thistime[k,:]
        #print("len(xes_this):",len(xes_this))
        x1 = xes_this[0]
        y1 = ys_this[0]
        z1 = zs_this[0]
        bindex1 = 0
        for j in range(N-1): # Looping over atoms in the chain to get bonds
            bindex2 = j+1
            x2 = xes_this[j+1]
            y2 = ys_this[j+1]
            z2 = zs_this[j+1]
            dx = x2-x1
            dy = y2-y1
            dz = z2-z1
            if abs(dx)>halfLx:
                if dx>0:
                    dx-=Lx
                else:
                    dx+=Lx
            if abs(dy)>halfLy:
                if dy>0:
                    dy-=Ly
                else:
                    dy+=Ly
            if abs(dz)>halfLz:
                if dz>0:
                    dz-=Lz
                else:
                    dz+=Lz
            bondvecs[i,k,j,0] = dx # x-component of bond vector j for time i
            bondvecs[i,k,j,1] = dy # y-component of bond vector j for time i
            bondvecs[i,k,j,2] = dz # z-component of bond vector j for time i
            veclen2                = dx**2+dy**2+dz**2
            length_bondvecs[i,k,j] = np.sqrt(veclen2)
            
            if(length_bondvecs[i,k,j]>100):
                print("Loooong bond")
                print("Index1:", bindex1, "; index2:", bindex2)
                if(abs(dx)>100):
                    print("x1:", x1, "; x2:", x2, "dx:", dx)
                elif(abs(dy)>100):
                    print("y1:", y1, "; y2:", y2, "dy:", dy)
                else:
                    print("z1:", z1, "; z2:", z2, "dz:", dz)
                print("i:",i)
            
            # Test:
            #print("Index1:", bindex1, "; index2:", bindex2)
            chainno1 = index1[bindex1]
            chainno2 = index1[bindex2]
            if abs(chainno1-chainno2)>0:
                print("Dangah dangah! Made bond between atoms in different chains! Difference:", chainno1-chainno2)
            # Next bond: The end particle is now the start particle 
            x1 = x2
            y1 = y2
            z1 = z2
            bindex1 = bindex2

print("Made bonds")
end_time = time.process_time()

print("Time:", end_time-start_time, " s")
############################################################
# Have the bonds. Can perform the dot products now.
# Ignore the first bonds in each chain.

# I'll just write some stuff, then I know what I need to include:
# bond_corr
# To calculate the standard deviation:
bond_corr_all = []
counters      = np.zeros(length)
# Other measures: correlation for ONE frame and correlation for each bond.
# One frame:
bond_corr_oneframe     = np.zeros(length)
bond_corr_oneframe_rms = np.zeros(length)
bond_corr_oneframe_all = []
# All bonds:
bond_corr_allbonds_nn     = np.zeros(Nb)
bond_corr_allbonds_nn_rms = np.zeros(Nb)
bond_corr_allbonds_nn_all = np.zeros((Nb,Nt*len(neighbours[0])))
kcounters                 = np.zeros(Nb)

#dpcounter     = 0
for l in range(length):                  # Looping over nearest neighb., next nearest neighb., etc.
    nchains          = neighbours[l]
    temparr          = []
    temparr_oneframe = []
    for i in range(Nt):                  # Looping over all time steps. This is not that neccessary...
        for j in range(maxelements):     # Looping over all elements in the pairs.
            chainpair = nchains[j]
            #if l==0:
            #    print("chain pair, l= 0", chainpair)
            chain1    = int(chainpair[0])
            chain2    = int(chainpair[1])
            if chain1==0 and chain2==0:  # 'We're done here'
                break
            for k in range(Nb):#(1,Nb):        # Looping over all bonds in the chains (except the first one, which is fixed)
                # perform a sum. Maybe find a way to take the average.
                bond1         = bondvecs[i,chain1,k]
                bond2         = bondvecs[i,chain2,k]
                len1          = length_bondvecs[i,chain1,k]
                len2          = length_bondvecs[i,chain2,k]
                dotprod       = np.dot(bond1,bond2)/(len1*len2)  
                bond_corr[l] += dotprod # Index this with k for l=0, and I should be ok.
                temparr.append(dotprod)
                # Other measurements:
                # Looking at ONE frame:
                if i==oneframe:
                    bond_corr_oneframe[l] += dotprod
                    temparr_oneframe.append(dotprod)
                # Looking at all bonds:
                if l==0:
                    bond_corr_allbonds_nn[k] += dotprod
                    #print('kcounters[k]:',int(kcounters[k]))
                    bond_corr_allbonds_nn_all[k,int(kcounters[k])] = dotprod
                    kcounters[k] +=1 
                '''
                if dotprod>1:
                    #print("dot product>1!:", dotprod)
                    dpcounter += 1
                
                if l==0:
                    #print("l:", l, "; dotprod:", dotprod)
                    print("l=",l, "; chain1:", chain1, "; chain 2:", chain2, "; dotprod:", dotprod)
                '''
                counters[l] += 1
    bond_corr_all.append(temparr)
    bond_corr_oneframe_all.append(temparr_oneframe)
#print("dpcounter:", dpcounter)
#print("sum(counters):", sum(counters))

### Testing ###
print("bond_corr_allbonds_nn:", bond_corr_allbonds_nn)
print("bond_corr_allbonds_nn_all[2,:]",bond_corr_allbonds_nn_all[2,:])

### Find average and rms ###
for i in range(length):
    #print(bond_corr_all)
    if len(bond_corr_all[i])==0:
        break
    bond_corr[i], bond_corr_rms[i] = find_rms(bond_corr_all[i])

for i in range(length):
    #print(bond_corr_all)
    if len(bond_corr_oneframe_all[i])==0:
        break
    bond_corr_oneframe[i], bond_corr_oneframe_rms[i] = find_rms(bond_corr_oneframe_all[i])

for k in range(Nb):
    if len(bond_corr_allbonds_nn_all[k,:])==0:
        break
    bond_corr_allbonds_nn[k], bond_corr_allbonds_nn_rms[k] = find_rms(bond_corr_allbonds_nn_all[k,:])

#print("bond_corr_allbonds_nn:", bond_corr_allbonds_nn)

### Write to file ###

for i in range(length):
    outfile.write('%.16e %.16e %.16e\n' % (distances[i], bond_corr[i], bond_corr_rms[i]))
outfile.close()

for i in range(length):
    outfile2.write('%.16e %.16e %.16e\n' % (distances[i], bond_corr_oneframe[i], bond_corr_oneframe_rms[i]))
outfile2.close()

for i in range(Nb):
    outfile3.write('%i %.16e %.16e\n' % (i, bond_corr_allbonds_nn[i], bond_corr_allbonds_nn_rms[i]))
outfile3.close()


#divby = (Nb-1)*Nt* # Check that this is equal to counter. # This is difficult since we have a different number of pairs.


end_time = time.process_time()
print("Time:", end_time-start_time, " s")

#print("costheta_rms:", costheta_rms)

print(np.mean(length_bondvecs))


'''# e line:
x_eline = np.zeros(2)
y_eline = np.zeros(2)+ 1./np.exp(1)
x_eline[0] = 0
x_eline[1] = 100
'''

# Maybe plotting, printing and fitting?
plt.figure(figsize=(6,5))
plt.errorbar(distances, bond_corr, bond_corr_rms, fmt="none", capsize=2, label='Values')
plt.plot(distances, bond_corr, 'o')
#plt.plot(separation, fittedgraph, label='Fit')
#plt.plot(x_persistencelength, y_persistencelength, '--')
#plt.plot(x_eline, y_eline, '--')
plt.xlabel(r'Chain distance $\Delta$', fontsize=16)
plt.ylabel(r'<$\cos\theta_{\Delta}$>', fontsize=16)
plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
plt.legend(loc="upper right")
plt.title(r'<$\cos\theta_{\Delta}$> vs chain distance $\Delta$', fontsize=16)
plt.savefig(plotname)

'''
plt.figure(figsize=(6,5))
#plt.plot(distances_unpacked, costheta_unpacked, 'o', label='Values')
plt.plot(separation, fittedgraph, label='Fit, unit b')
plt.plot(separation, fittedgraph2, label='Fit, exact b')
#plt.plot(x_persistencelength, y_persistencelength, '--')
#plt.plot(x_eline, y_eline, '--')
plt.xlabel(r'Bond distance', fontsize=16)
plt.ylabel(r'<cos($\theta$)>', fontsize=16)
plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
plt.legend(loc="upper right")
plt.title(r'<cos($\theta$)> vs separation, actual length', fontsize=16)
plt.savefig(plot3name)
'''
