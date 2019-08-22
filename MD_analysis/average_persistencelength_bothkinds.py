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
factors  = [0.1,1,10,100,250]
charges  = [0,-1,-5,-10]
#spacings = [1,2,3,4,5,10,40,100]
#N        = 21
spacing  = 40
spacings = [2,3,4,5,8,10,15,40,100]#[1,2,3,4,5,8,10,15,40,100]
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
    outfilename  = 'persistence_length_varycharge_chaingrid_quadratic_M9N101_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_T%i_theta0is180_twofirst_are_fixed.txt' % (kangle,kbond,T)
    outfile      = open(outfilename, 'w')
    outfile.write('Order of elements: Charge,   Average persistence length,   Root mean square persistence length\n')

if factorsims==True:
    outfilename  = 'persistence_length_varycharge_neighbourbonds_Kangle%i_T%i_angle_stretching_only_varyfactor.txt' % (kangle,T)
    outfile      = open(outfilename, 'w')
    outfile.write('Order of elements: Factor,   Average persistence length,   Root mean square persistence length\n')

if spacesims==True:
    outfilename  = 'persistence_length_neighbourbonds_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_N%i_varyspacing_II.txt' % (kangle,kbond,T,N)
    outfile      = open(outfilename, 'w')
    outfile.write('Order of elements: Spacing,   Average persistence length,   Root mean square persistence length\n')

if lengthsims==True:
    outfilename  = 'persistence_length_varycharge_neighbourbonds_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varylength_contourdsep20.txt' % (kangle,kbond,T)
    outfile      = open(outfilename, 'w')
    outfile.write('Order of elements: N,   Average persistence length,   Root mean square persistence length\n')

if lensfxsims==True:
    outfilename  = 'persistence_length_quadratic_gridspacing%i_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varylength.txt' % (spacing,kangle,kbond,T)
    outfile      = open(outfilename, 'w')
    outfile.write('Order of elements: N,   Average bond length,   Root mean square bond length\n')

if krfxfacsims==True:
    outfilename  = 'persistence_length_neighbourbonds_spacing%i_Kbond%i_T%i_angle_stretching_only_varyfactor_Krfixed.txt' % (spacing,kbond,T)
    outfile      = open(outfilename, 'w')
    outfile.write('Order of elements: Factor,   Average bond length,   Root mean square bond length\n')

k = 0
#for kangle in kangles:
for i in range(Narray):
    print((k+1)/float(Narray))
    
    if chargesims==True:
        charge = charges[i]
        infilename   = 'persistence_length_chaingrid_quadratic_M9N101_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge%i_T310_theta0is180_twofirst_are_fixed.txt' % charge 
        infilename2  = 'persistence_length_actualdistance_chaingrid_quadratic_M9N101_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge%i_T310_theta0is180_twofirst_are_fixed.txt' % charge    
    
    if factorsims==True:
        factor = factors[i]
        kbond  = kangle*factors[i]
        infilename   = 'persistence_length_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,kangle,kbond) +str(factor)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename2  = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,kangle,kbond) +str(factor)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
    
    if spacesims==True:
        spacing      = spacings[i]
        infilename   = 'persistence_length_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename2  = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
    
    if lengthsims==True:
        N            = lengths[i]
        spacing      = lensp[i]
        infilename   = 'persistence_length_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename2  = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
    
    if lensfxsims==True:
        N            = lengths[i]
        infilename   = 'persistence_length_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename2  = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
    
    if krfxfacsims==True:
        factor       = krfacs[i]
        kangle       = kangles[i]
        #	       costhetas_neighbourbonds_chaingrid_quadratic_M9N101_gridspacing40_Langevin_Kangle2000_Kbond2000_factor1.00_T310_theta0is180_twofirst_are_fixed.txt
        infilename   = 'persistence_length_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
        infilename2  = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
    
    #print('infilename:', infilename)
    # First file:
    infile = open(infilename,'r')
    
    line   = infile.readline()
    words  = line.split()
    pl_av  = float(words[8])
    pl_rms = float(words[9])
    
    infile.close()
    
    infile2 = open(infilename,'r')
    
    line   = infile2.readline()
    words  = line.split()
    pl_av  = float(words[8])
    pl_rms = float(words[9])
    
    infile2.close()
    
    if chargesims==True:
        outfile.write('%.16f %.16f %.16f %.16f\n' % (charge, pl_av, pl_rms, pl_rms2))
    if factorsims==True or krfxfacsims==True:
        outfile.write('%.16f %.16f %.16f %.16f\n' % (factor, pl_av, pl_rms, pl_rms2))
    if spacesims==True:
        outfile.write('%.16f %.16f %.16f %.16f\n' % (spacing, pl_av, pl_rms, pl_rms2))
    if lengthsims==True or lensfxsims==True:
        outfile.write('%.16f %.16f %.16f %.16f\n' % (N, pl_av, pl_rms, pl_rms2))
    k+=1
    
outfile.close()
