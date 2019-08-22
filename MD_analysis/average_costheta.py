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
    outfilename  = 'average_costhetas_neighbourbonds_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varycharge.txt' % (kangle,kbond,T)
    outfilename2 = 'ree_average_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varycharge.txt' % (kangle,kbond,T)
    outfile      = open(outfilename, 'w')
    outfile2     = open(outfilename2, 'w')
    outfile.write('Order of elements: Charge,   Average bond length,   Root mean square bond length\n')

if factorsims==True:
    outfilename  = 'average_costhetas_neighbourbonds_Kangle%i_T%i_angle_stretching_only_varyfactor.txt' % (kangle,T)
    outfilename2 = 'ree_average_Kangle%i_T%i_angle_stretching_only_varyfactor.txt' % (kangle,T)
    outfile      = open(outfilename, 'w')
    outfile2     = open(outfilename2, 'w')
    outfile.write('Order of elements: Factor,   Average bond length,   Root mean square bond length\n')

if spacesims==True:
    outfilename  = 'average_costhetas_neighbourbonds_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varyspacing_N%i.txt' % (kangle,kbond,T,N)
    outfilename2 = 'ree_average_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varyspacing_N%i.txt' % (kangle,kbond,T,N)
    outfile      = open(outfilename, 'w')
    outfile2     = open(outfilename2, 'w')
    outfile.write('Order of elements: Spacing,   Average cos(theta_1),   Root mean square cos(theta_1)\n')
    outfile2.write('Order of elements: Spacing,   Average ree2,   Root mean square ree2\n')

if lengthsims==True:
    outfilename  = 'average_costhetas_neighbourbonds_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varylength_contourdsep20.txt' % (kangle,kbond,T)
    outfilename2 = 'ree_average_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varylength_contourdsep20.txt' % (kangle,kbond,T)
    outfile      = open(outfilename, 'w')
    outfile2     = open(outfilename2, 'w')
    outfile.write('Order of elements: N,   Average cos(theta_1),   Root mean square cos(theta_1)\n')
    outfile2.write('Order of elements: N,   Average ree2,   Root mean square ree2\n')

if lensfxsims==True:
    outfilename  = 'average_costhetas_neighbourbonds_quadratic_gridspacing%i_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varylength.txt' % (spacing,kangle,kbond,T)
    outfilename2 = 'ree_average_quadratic_gridspacing%i_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varylength.txt' % (spacing,kangle,kbond,T)
    outfile      = open(outfilename, 'w')
    outfile2     = open(outfilename2, 'w')
    outfile.write('Order of elements: N,   Average bond length,   Root mean square bond length\n')
    outfile2.write('Order of elements: N,   Average ree2,   Root mean square ree2\n')

if krfxfacsims==True:
    outfilename  = 'average_costhetas_neighbourbonds_spacing%i_Kbond%i_T%i_angle_stretching_only_varyfactor_Krfixed.txt' % (spacing,kbond,T)
    outfilename2 = 'ree_average_spacing%i_Kbond%i_T%i_angle_stretching_only_varyfactor_Krfixed.txt' % (spacing,kbond,T)
    outfile      = open(outfilename, 'w')
    outfile2     = open(outfilename2, 'w')
    outfile.write('Order of elements: Factor,   Average bond length,   Root mean square bond length\n')

k = 0
#for kangle in kangles:
for i in range(Narray):
    print((k+1)/float(Narray))
    
    if chargesims==True:
        charge      = charges[i]
        infilename  = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,kangle,kbond) +str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename2 = 'ree_average_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,kangle,kbond) +str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
    
    if factorsims==True:
        kbond       = kangle*factors[i]
        infilename  = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,kangle,kbond) +str(factor)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename2 = 'ree_average_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,kangle,kbond) +str(factor)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
    
    if spacesims==True:
        spacing     = spacings[i]
        infilename  = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename2 = 'ree_average_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        
    if lengthsims==True:
        N           = lengths[i]
        spacing     = lensp[i]
        infilename  = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename2 = 'ree_average_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
    
    if lensfxsims==True:
        N           = lengths[i]
        infilename  = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename2 = 'ree_average_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
    
    if krfxfacsims==True:
        factor      = krfacs[i]
        kangle      = kangles[i]
        #	       costhetas_neighbourbonds_chaingrid_quadratic_M9N101_gridspacing40_Langevin_Kangle2000_Kbond2000_factor1.00_T310_theta0is180_twofirst_are_fixed.txt
        infilename  = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
        infilename2 = 'ree_average_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
    print('infilename:', infilename)
    # First file:
    infile = open(infilename,'r')
    
    lines  = infile.readlines()
    
    #print(infilename)
    counter      = 0
    costheta_av  = 0
    costheta_rms = 0
    for line in lines:
        words = line.split()
        if len(words)!=0:
            #print("words:", words)
            counter      += 1
            costheta_av  +=float(words[0])
            costheta_rms +=float(words[1])
    
    infile.close()
    costheta_av  /= counter
    costheta_rms /= counter
    
    # Second file:
    infile2 = open(infilename2,'r')
    
    lines  = infile2.readlines()
    
    counter      = 0
    ree2_av      = 0
    ree2_rms     = 0
    for line in lines:
        words = line.split()
        if len(words)!=0:
            counter   += 1
            ree2_av   +=float(words[4])
            ree2_rms  +=float(words[5])
    
    infile2.close()
    ree2_av  /= counter
    ree2_rms /= counter
    
    if chargesims==True:
        outfile.write('%.16f %.16f %.16f\n' % (charge, costheta_av, costheta_rms))
        outfile2.write('%.16f %.16f %.16f\n' % (charge, ree2_av, ree2_rms))
    if factorsims==True or krfxfacsims==True:
        outfile.write('%.16f %.16f %.16f\n' % (factor, costheta_av, costheta_rms))
        outfile2.write('%.16f %.16f %.16f\n' % (factor, ree2_av, ree2_rms))
    if spacesims==True:
        outfile.write('%.16f %.16f %.16f\n' % (spacing, costheta_av, costheta_rms))
        outfile2.write('%.16f %.16f %.16f\n' % (spacing, ree2_av, ree2_rms))
    if lengthsims==True or lensfxsims==True:
        outfile.write('%.16f %.16f %.16f\n' % (N, costheta_av, costheta_rms))
        outfile2.write('%.16f %.16f %.16f\n' % (N, ree2_av, ree2_rms))
    k+=1
    
outfile.close()
outfile2.close()
