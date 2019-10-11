import matplotlib.pyplot as plt               # To plot
from scipy.optimize import curve_fit
import numpy as np
import random
import math
import time

# Function for curve_fit
def costheta_exponential(s,P):
    return np.exp(-s/P)

start_time = time.process_time()

arctest     = False
### Input file parameters
Kangle      = 20#125000
Kbond       = 200#0
factor      = Kangle/float(Kbond)
T           = 310
M           = 9
N           = 101 ###
ljdebye     = 1.042
epsilon     = ljdebye
sigma       = 1
kappa       = 1
ljcutoff    = 1.12246204830937
debyecutoff = 3
#factor      = 0.05#250
#Kbond       = 2000#Kangle*factor
#Kangle      = Kbond*factor
charge      = -1
mass        = 1
spacing     = 300#4#0
gridspacing = spacing
K           = Kangle # Because we used this notation earlier, but using Kangle is less confusing
wallenergy  = 1.042
dielectric  = 1#50 #1 2 10 100
systemtest  = False
lotsofchains = False

### Input and output file names

spacings = [1,2,3,4,5,6,7,8,10,15,40,100]
outfilename = 'crossings_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed_varygridspacing.txt' % (M,N,Kangle,Kbond,charge,T)
outfile     = open(outfilename,'w')
outfile.write('begin{table}\n\centering\n\caption{The amount of the chains that are under the z=0-plane per time}\n begin{tabular}{r|c|c|c|c|c|c}\nSpacing ')
for spacing in spacings:
    outfile.write('& %i ' % spacing)
outfile.write('\ \ \n\hline\n')
outfile.write('$f_{under}$')


for spacing in spacings:
    #basename     = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
    basename     = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,spacing,Kangle,Kbond,charge,T)
    #basename     = 'chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)
    #basename     = 'chaingrid_quadratic_M'+str(M)+'N'+str(N)+'_ljdebye_epsilon'+str(epsilon)+'_sigma'+str(sigma)+'_ljcutoff'+str(ljcutoff)+'_kappa1_debyecutoff'+str(debyecutoff)+'_Langevin_wall'+str(ljdebye)+'_Kangle'+str(Kangle)+'_Kbond'+str(Kbond)+'_T'+str(T)+'_theta0is180_twofirst_are_fixed'     # Actual run
    #basename     = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed' % (M,N,Kangle,Kbond,factor,T)
    #basename2    = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed' % (M,N,gridspacing,Kangle,Kbond,factor,T)
    #basename      = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa%i_debyecutoff%i_charge%i_mass%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,gridspacing,Kangle,Kbond,kappa,debyecutoff,charge,mass,T)
    basename2     = basename
    
    infilename   = basename+'.lammpstrj'
    
    # 
    skiplines           = 9             # If we hit 'ITEM:', skip this many steps...
    skipelem            = 0#10#1000#10000#10000#90000 # The number of elements we skip in order to equilibrate (set to 0 if the .lammpstrj file should be equilibrated)
    sampleevery         = 0#1 # Sample every 10th of the values written to file # The program gets WAY too slow if this is too small.
    
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
            
    xes = []
    ys  = []
    zs  = []
    
    x_curr = np.zeros((M,N))
    y_curr = np.zeros((M,N))
    z_curr = np.zeros((M,N))
    tempzs = np.zeros(M)
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
    
    if systemtest==False:
        # Removing the first redundant element of position lists (neccessary if t=0 is not the first time step in the file/that we read)
        xes.pop(0)
        ys.pop(0)
        zs.pop(0)
    if sampleevery==0: # Because sampleevery!=0 can mess with the program, and leave an empty x_curr after the read-in loop
        xes.append(x_curr)
        ys.append(y_curr)
        zs.append(z_curr)
    
    infile.close()
    print('Done with read-in loop')
    
    # Finding the bond vectors for each time step.
    # This will be an awfully blocky code, it seems
    
    Nt = len(xes) # Could have used j (Nt = j)
    
    # A test:
    infilename_zendtest = 'test_zend_from_find_persistencelength.txt'
    infile_zendtest     = open(infilename_zendtest,'w')
    # Crap. I need another loop and even more indices. Yikes.
    bindex1 = 0
    bindex2 = 0
    counter_all  = 0
    chainisunder_counter = 0
    for i in range(Nt):    # Loop over time steps
        for k in range(M): # Loop over chains
            xes_thistime = xes[i] # [i,k] # I think I need to extract this twice, because it is a list of arrays (a double one at that...)
            ys_thistime  = ys[i]
            zs_thistime  = zs[i]
            xes_this     = xes_thistime[k,:]
            ys_this      = ys_thistime[k,:]
            zs_this      = zs_thistime[k,:]
            x1 = xes_this[0]
            y1 = ys_this[0]
            z1 = zs_this[0]
            bindex1 = 0
            chainisunder = False # Bool for knowing if the current part of the chain is under or over the z=0-plane.
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
                        chainisunder = True
                    else:
                        dz+=Lz
                        chainisunder = False
                        
                if chainisunder==True:
                    chainisunder_counter += 1
                # Test:
                #print("Index1:", bindex1, "; index2:", bindex2)
                chainno1 = index1[bindex1]
                chainno2 = index1[bindex2]
                counter_all +=1
                if abs(chainno1-chainno2)>0:
                    print("Dangah dangah! Made bond between atoms in different chains! Difference:", chainno1-chainno2)
                # Next bond: The end particle is now the start particle 
                x1 = x2
                y1 = y2
                z1 = z2
                bindex1   = bindex2
    chainisunder_counter /= counter_all#(M*(N-1)*Nt)
    outfile.write(' & %.2e' % chainisunder_counter)
outfile.write('\ \ \n')

outfile.write('\end{tabular}}\n\label{table:endz_chains_something}\n\end{table}') 
outfile.close()
end_time = time.process_time()

print("Time:", end_time-start_time, " s")


