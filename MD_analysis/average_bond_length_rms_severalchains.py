import matplotlib.pyplot as plt               # To plot
from scipy.optimize import curve_fit
import numpy as np
import random
import math
import time


start_time = time.process_time()

M        = 9
N        = 101
kangle   = 20
kbond    = 200
factors  = [0.1,1,10,100, 250]
charges  = [0,-1,-5,-10]
#spacings = [1,2,3,4,5,10,40,100]
N        = 101
spacing  = 40
spacings = [1,2,3,4,5,6,7,8,10,15,40,100]#[1,2,3,4,5,8,10,15,40,100]
lengths  = [21,61,101,141]
lensp    = [1,3,5,7]
krfacs   = [0.01, 0.05, 0.10, 0.50, 1.00]
kangles  = [20, 100, 200, 1000, 2000]
charge   = -1
T        = 310

chargesims  = False
factorsims  = False
spacesims   = True
lengthsims  = False # Varying the length and the spacing
lensfxsims  = False  # Varying the length, fixing the spacing
krfxfacsims = False

if chargesims==True:
    Narray  = len(charges)
if factorsims==True:
    Narray  = len(factors)
if spacesims==True:
    Narray = len(spacings)
if lengthsims==True:
    Narray = len(lengths)
if lensfxsims==True:
    Narray = len(lengths)
if krfxfacsims==True:
    kbond = 2000
    Narray = len(krfacs)


bondlengths = np.zeros(Narray)
bondrms     = np.zeros(Narray)

if chargesims==True:
    outfilename = 'bond_lengths_and_rms_Kangle%i_Kbond%i_T%i_adebye_kappa1_debyecutoff3_varycharge.txt' % (kangle,kbond,T)
    outfile     = open(outfilename, 'w')
    outfile.write('Order of elements: Charge,   Average bond length,   Root mean square bond length\n') ####### Rewrite output if you change this!!!

if factorsims==True:
    outfilename = 'bond_lengths_and_rms_Kangle%i_T%i_angle_stretching_only_varyfactor.txt' % (kangle,T)
    outfile     = open(outfilename, 'w')
    outfile.write('Order of elements: Factor,   Average bond length,   Root mean square bond length\n')

if spacesims==True:
    outfilename  = 'bond_lengths_and_rms_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varyspacing_N%i.txt' % (kangle,kbond,T,N)
    outfile      = open(outfilename, 'w')
    outfile.write('Order of elements: Spacing,   Average bond length,   Root mean square bond length\n')

if lengthsims==True:
    outfilename  = 'bond_lengths_and_rms_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varylength_contourdsep20.txt' % (kangle,kbond,T)
    outfile      = open(outfilename, 'w')
    outfile.write('Order of elements: N,   Average bond length,   Root mean square bond length\n')

if lensfxsims==True:
    outfilename  = 'bond_lengths_and_rms_quadratic_gridspacing%i_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varylength.txt' % (spacing,kangle,kbond,T)
    outfile      = open(outfilename, 'w')
    outfile.write('Order of elements: N,   Average bond length,   Root mean square bond length\n')

if krfxfacsims==True:
    outfilename  = 'bond_lengths_and_rms_spacing%i_Kbond%i_T%i_angle_stretching_only_varyfactor_Krfixed.txt' % (spacing,kbond,T)
    outfile      = open(outfilename, 'w')
    outfile.write('Order of elements: Factor,   Average bond length,   Root mean square bond length\n')

k = 0
#for kangle in kangles:
for i in range(Narray):
    print((k+1)/float(Narray))
    #kbond  = kangle*factor[i]
    if chargesims==True:
        charge = charges[i]
        infilename   = 'bond_lengths_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,kangle,kbond) +str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        outfilename2 = 'bond_lengths_average_and_rms_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,kangle,kbond) +str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
    if factorsims==True:
        factor       = factors[i]
        infilename   = 'bond_lengths_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,kangle,kbond) +str(factor)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        outfilename2 = 'bond_lengths_average_and_rms_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,kangle,kbond) +str(factor)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
    
    if spacesims==True:
        spacing    = spacings[i]
        infilename = 'bond_lengths_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        outfilename2 = 'bond_lengths_average_and_rms_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
    if lengthsims==True:
        N          = lengths[i]
        spacing    = lensp[i]
        infilename = 'bond_lengths_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        outfilename2 = 'bond_lengths_average_and_rms_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
    
    if lensfxsims==True:
        N          = lengths[i]
        infilename = 'bond_lengths_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        outfilename2 = 'bond_lengths_average_and_rms_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
    
    if krfxfacsims==True:
        factor      = krfacs[i]
        kangle      = kangles[i]
        #	       costhetas_neighbourbonds_chaingrid_quadratic_M9N101_gridspacing40_Langevin_Kangle2000_Kbond2000_factor1.00_T310_theta0is180_twofirst_are_fixed.txt
        infilename  = 'bond_lengths_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
        outfilename2 = 'bond_lengths_average_and_rms_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
    
    infile = open(infilename,'r')
    
    header = infile.readline()
    words  = header.split()
    Nt     = int(words[1])
    Nbonds = int(words[3])
    
    lines  = infile.readlines()
    Nlines = len(lines)
    
    i = 0
    
    #totalav      = 0
    allbonds     = []
    thistimestep = []
    alltimesteps = []
    while i<Nlines:
        words = lines[i].split()
        if words[0]=='Time:':
            i+=1
            alltimesteps.append(thistimestep)
            thistimestep = []
        else:
            bondlength = float(words[2])
            #totalav += bondlength
            allbonds.append(bondlength)
            thistimestep.append(bondlength)
            i+=1
    infile.close()
    allbonds = np.array(allbonds)
    totalav = np.mean(allbonds)
    #totalav /= len(totalav)
    
    alltimesteps.append(thistimestep)
    #print(alltimesteps[0])
    alltimesteps.pop(0)
    
    #print("len(alltimesteps):", len(alltimesteps))
    
    # rms, the whole simulation:
    total_rms = 0
    for i in range(Nbonds):
        total_rms += (allbonds[i]-totalav)**2
    total_rms =  np.sqrt(total_rms/Nbonds)
    
    Nb = Nbonds # Had some trouble here, not going to rewrite it all #len(alltimesteps[2])
    #alltimesteps = np.array(alltimesteps)
    each_average = np.zeros(Nt) #np.mean(alltimesteps, axis=1)
    timestep_rms = np.zeros(Nt)
    for i in range(Nt):
        thistimestep = alltimesteps[i]
        for j in range(Nb):
            each_average[i] += thistimestep[j]
        each_average[i] /=  Nb
    timestep_rms = np.zeros(Nt)
    Nrms = len(each_average) # This should be Nt
    # rms, per time step
    for i in range(Nt):
        thistimestep = alltimesteps[i]
        for j in range(Nb):
            timestep_rms[i] += (thistimestep[j]-each_average[i])**2
        timestep_rms[i] =  np.sqrt(timestep_rms[i]/Nb)
    
    
    # Print to file: Properties of the entire simulation # I should probably exclude some of the first time steps...
    if factorsims==True or krfxfacsims==True:
        outfile.write('%.16f %.16f %.16f\n' % (factor, totalav, total_rms))
    if chargesims==True:
        outfile.write('%.16f %.16f %.16f\n' % (charge, totalav, total_rms))
    if spacesims==True:
        outfile.write('%.16f %.16f %.16f\n' % (spacing, totalav, total_rms))
    if lengthsims==True or lensfxsims==True:
        outfile.write('%.16f %.16f %.16f\n' % (N, totalav, total_rms))
    
    # Print to file: Properties of this time step
    outfile2 = open(outfilename2, 'w')
    outfile2.write('Order of elements: Timestep number,   Average bond length,   Root mean square bond length\n')
    for i in range(Nt):
        outfile2.write('%i %.16f %.16f\n' % (i, each_average[i], timestep_rms[i]))
    outfile2.close()
    k += 1
    
outfile.close()
