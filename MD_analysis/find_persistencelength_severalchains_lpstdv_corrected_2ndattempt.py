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
Kangle      = 25
Kbond       = 2000
T           = 310
M           = 9
N           = 101 ###
ljdebye     = 1.042
epsilon     = ljdebye
sigma       = 1
ljcutoff    = 1.12246204830937
debyecutoff = 3
factor      = 0.05#250
#Kbond       = 2000#Kangle*factor
#Kangle      = Kbond*factor
charge      = -1
spacing     = 40
gridspacing = spacing
K           = Kangle # Because we used this notation earlier, but using Kangle is less confusing
wallenergy  = 1.042
dielectric  = 1#50 #1 2 10 100

### Input and output file names

''' # Testing; <ree2> is correct for completely straight chains.
N  = 100 # Number of bond vectors = number of units - 1 # Per chain
M  = 9   # Number of chains
L2 = 3   # Determining the shape of the L1xL2 polymer grid on the surface
gridspacing = 40 # The spacing in our grid
Nelem   = M*(N+1)
N       = 101
infilename   = 'chaingrids_totallystraight_N%i_Nchains%i_Ly%i_twofixed_test.lammpstrj'  % (Nelem,M,L2)    # Actual run
plot1name    = 'costheta_vs_separation_chaingrid_quadratic_M%iN%i_totallystraight_test.png' % (M,N)
plot2name    = 'angledistr_sep5_chaingrid_quadratic_M%iN%i_totallystraight_test.png' % (M,N)
plot3name    = 'costheta_vs_separation_actualdistance_chaingrid_quadratic_M%iN%i_totallystraight_test.png' % (M,N)
outfilename  = 'bond_lengths_chaingrid_quadratic_M%iN%i_totallystraight_test.txt' % (M,N)
outfilename2 = 'persistence_length_chaingrid_quadratic_M%iN%i_totallystraight_test.txt' % (M,N)
outfilename3 = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_totallystraight_test.txt' % (M,N)
outfilename4 = 'ree_last_chaingrid_quadratic_M%iN%i_totallystraight_test.txt' % (M,N)
outfilename5 = 'ree_average_chaingrid_quadratic_M%iN%i_totallystraight_test.txt' % (M,N)
outfilename6 = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_totallystraight_test.txt' % (M,N)
outfilename7 = 'endzs_chaingrid_quadratic_M%iN%i_totallystraight_test.txt' % (M,N)
outfilename8 = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_totallystraight_test.txt' % (M,N)
#'''

#	     #  chaingrid_quadratic_M9N101_Langevin_Kangle100_Kbond2000_factor0.05_T310_theta0is180_twofirst_are_fixed
#gridspacing  = 40
'''   # This one is for varying the factor.
infilename   = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,Kangle,Kbond,factor,T)
plot1name    = 'costheta_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,spacing,Kangle,Kbond,factor,T)
plot2name    = 'angledistr_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,spacing,Kangle,Kbond,factor,T)
plot3name    = 'costheta_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,spacing,Kangle,Kbond,factor,T)
outfilename  = 'bond_lengths_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,Kangle,Kbond,factor,T)

outfilename2 = 'persistence_length_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,Kangle,Kbond,factor,T)
outfilename3 = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,Kangle,Kbond,factor,T)
outfilename4 = 'ree_last_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,Kangle,Kbond,factor,T)
outfilename5 = 'ree_average_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,Kangle,Kbond,factor,T)
outfilename6 = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,Kangle,Kbond,factor,T)
outfilename7 = 'endzs_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,Kangle,Kbond,factor,T)
outfilename8 = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,Kangle,Kbond,factor,T)
#'''

'''   # Changing the dielectric constant
#	     # 'chaingrid_quadratic_M9N101_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_dielectric10_T310_theta0is180_twofirst_are_fixed.lammpstrj'
#            # 'chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj'
#            # 'chaingrid_quadratic_M9N101_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_dielectric100_T310_theta0is180_twofirst_are_fixed.lammpstrj'
infilename   = 'chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
plot1name    = 'costheta_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
plot2name    = 'angledistr_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
plot3name    = 'costheta_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)

outfilename  = 'bond_lengths_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
outfilename2 = 'persistence_length_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
outfilename3 = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
outfilename4 = 'ree_last_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
outfilename5 = 'ree_average_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
outfilename6 = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
outfilename7 = 'endzs_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
outfilename8 = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
#'''


# Tester:
arctest      = True
testtype     = 'straight'#'rotatingarc_withtime'#'rotatingarc_withtime'#'rotatingarc'#'quartercirc'#'straight'
infilename   = 'testsystem_chaingrid_'+testtype+'.lammpstrj'
plot1name    = 'costheta_testsystem_chaingrid_'+testtype+'.png'
plot2name    = 'angledistr_testsystem_chaingrid_'+testtype+'.png'
plot3name    = 'costheta_actualdistance_testsystem_chaingrid_'+testtype+'.png'  
outfilename  = 'bond_lengths_testsystem_chaingrid_'+testtype+'.txt'

outfilename2 = 'persistence_length_testsystem_chaingrid_'+testtype+'.txt'
outfilename3 = 'costhetas_neighbourbonds_testsystem_chaingrid_'+testtype+'.txt'
outfilename4 = 'ree_last_testsystem_chaingrid_'+testtype+'.txt'
outfilename5 = 'ree_average_testsystem_chaingrid_'+testtype+'.txt'
outfilename6 = 'persistence_length_actualdistance_testsystem_chaingrid_'+testtype+'.txt'
outfilename7 = 'endzs_testsystem_chaingrid_'+testtype+'.txt'
outfilename8 = 'bondlengths_systemwide_testsystem_chaingrid_'+testtype+'.txt'  



''' # With wall potential: #### This is what I've been using recently (good for voxellation)
infilename   = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
plot1name    = 'costheta_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.png'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
plot2name    = 'angledistr_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.png'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
plot3name    = 'costheta_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.png'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
outfilename  = 'bond_lengths_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)

outfilename2 = 'persistence_length_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
outfilename3 = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
outfilename4 = 'ree_last_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
outfilename5 = 'ree_average_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
outfilename6 = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
outfilename7 = 'endzs_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
outfilename8 = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
#'''

# Output names for code testing:
'''
infilename   = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,spacing,Kangle,Kbond,charge,T)
plot1name    = 'costheta_test_short.png'
plot2name    = 'angledistr_test_short.png'
plot3name    = 'costheta_actualdistance_test_short.png'
outfilename  = 'bond_lengths_test_short.txt'

outfilename2 = 'persistence_length_test_short.txt'
outfilename3 = 'costhetas_neighbourbonds_test_short.txt'
outfilename4 = 'ree_last_test_short.txt'
outfilename5 = 'ree_average_test_short.txt'
outfilename6 = 'persistence_length_actualdistance_test_short.txt' 
outfilename7 = 'endzs_test_short.txt'
outfilename8 = 'bondlengths_systemwide_short.txt'
#'''
# Varying the grid spacing # THIS IS NOW THE STANDARD FILE NAMES.
'''
infilename   = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,spacing,Kangle,Kbond,charge,T)
plot1name    = 'costheta_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,spacing,Kangle,Kbond,charge,T)
plot2name    = 'angledistr_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,spacing,Kangle,Kbond,charge,T)
plot3name    = 'costheta_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,spacing,Kangle,Kbond,charge,T)
outfilename  = 'bond_lengths_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,Kangle,Kbond,charge,T)

outfilename2 = 'persistence_length_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,Kangle,Kbond,charge,T)
outfilename3 = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,Kangle,Kbond,charge,T)
outfilename4 = 'ree_last_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,Kangle,Kbond,charge,T)
outfilename5 = 'ree_average_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,Kangle,Kbond,charge,T)
outfilename6 = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed_varygridspacing.txt' % (M,N,spacing,Kangle,Kbond,charge,T)
outfilename7 = 'endzs_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed_varygridspacing.txt' % (M,N,spacing,Kangle,Kbond,charge,T)
outfilename8 = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed_varygridspacing.txt' % (M,N,spacing,Kangle,Kbond,charge,T)
#'''

# Varying the charge:
#chaingrid_quadratic_M9N101_Langevin_Kangle${Kangle}_Kbond${Kbond}_charge${charge}_T$T_theta0is180_twofirst_are_fixed.lammpstrj
''' # String appending method
#               chaingrid_quadratic_M9N101_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge0_T310_theta0is180_twofirst_are_fixed
infilename   = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % T     # Actual run
plot1name    = 'costheta_vs_separation_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.png' % T
plot2name    = 'angledistr_sep5_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.png' % T
plot3name    = 'costheta_vs_separation_actualdistance_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.png' % T
outfilename  = 'bond_lengths_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
outfilename2 = 'persistence_length_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
outfilename3 = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
outfilename4 = 'ree_last_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
outfilename5 = 'ree_average_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
outfilename6 = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
outfilename7 = 'endzs_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
outfilename8 = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
#'''


''' # String appending method
infilename   = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % T     # Actual run
plot1name    = 'costheta_vs_separation_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.png' % T
plot2name    = 'angledistr_sep5_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.png' % T
plot3name    = 'costheta_vs_separation_actualdistance_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.png' % T
outfilename  = 'bond_lengths_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
outfilename2 = 'persistence_length_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
outfilename3 = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
outfilename4 = 'ree_last_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
outfilename5 = 'ree_average_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
outfilename6 = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
outfilename7 = 'endzs_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
outfilename8 = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed.txt' % T
#'''


''' # Sprintf method
infilename   = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,Kangle,Kbond,factor,T)     # Actual run
plot1name    = 'costheta_vs_separation_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,Kangle,Kbond,factor,T)
plot2name    = 'angledistr_sep5_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,Kangle,Kbond,factor,T)
plot3name    = 'costheta_vs_separation_actualdistance_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,Kangle,Kbond,factor,T)
outfilename  = 'bond_lengths_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,Kangle,Kbond,factor,T)
outfilename2 = 'persistence_length_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,Kangle,Kbond,factor,T)
outfilename3 = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,Kangle,Kbond,factor,T)
outfilename4 = 'ree_last_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,Kangle,Kbond,factor,T)
outfilename5 = 'ree_average_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,Kangle,Kbond,factor,T)
outfilename6 = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,Kangle,Kbond,factor,T)
outfilename7 = 'endzs_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,Kangle,Kbond,factor,T)
outfilename8 = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,Kangle,Kbond,factor,T)
#'''


'''
infilename   = 'chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.lammpstrj' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)     # Actual run
plot1name    = 'costheta_vs_separation_chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.png' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)
plot2name    = 'angledistr_sep5_chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.png' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)
plot3name    = 'costheta_vs_separation_actualdistance_chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.png' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)
outfilename  = 'bond_lengths_chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.txt' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)
outfilename2 = 'persistence_length_chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.txt' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)
outfilename3 = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.txt' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)
outfilename4 = 'ree_last_chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.txt' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)
outfilename5 = 'ree_average_chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.txt' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)
outfilename6 = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.txt' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)
outfilename7 = 'endzs_chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.txt' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)
outfilename8 = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180.txt' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)
#'''

'''
infilename   = 'chaingrid_quadratic_M'+str(M)+'N'+str(N)+'_ljdebye_epsilon'+str(epsilon)+'_sigma'+str(sigma)+'_ljcutoff'+str(ljcutoff)+'_kappa1_debyecutoff'+str(debyecutoff)+'_Langevin_wall'+str(ljdebye)+'_Kangle'+str(Kangle)+'_Kbond'+str(Kbond)+'_T'+str(T)+'_theta0is180_twofirst_are_fixed.lammpstrj'     # Actual run
plot1name    = 'costheta_vs_separation_chaingrid_quadratic_M'+str(M)+'N'+str(N)+'_ljdebye_epsilon'+str(epsilon)+'_sigma'+str(sigma)+'_ljcutoff'+str(ljcutoff)+'_kappa1_debyecutoff'+str(debyecutoff)+'_Langevin_wall'+str(ljdebye)+'_Kangle'+str(Kangle)+'_Kbond'+str(Kbond)+'_T'+str(T)+'_theta0is180_twofirst_are_fixed.png'
plot2name    = 'angledistr_sep5_chaingrid_quadratic_M'+str(M)+'N'+str(N)+'_ljdebye_epsilon'+str(epsilon)+'_sigma'+str(sigma)+'_ljcutoff'+str(ljcutoff)+'_kappa1_debyecutoff'+str(debyecutoff)+'_Langevin_wall'+str(ljdebye)+'_Kangle'+str(Kangle)+'_Kbond'+str(Kbond)+'_T'+str(T)+'_theta0is180_twofirst_are_fixed.png'
plot3name    = 'costheta_vs_separation_actualdistance_chaingrid_quadratic_M'+str(M)+'N'+str(N)+'_ljdebye_epsilon'+str(epsilon)+'_sigma'+str(sigma)+'_ljcutoff'+str(ljcutoff)+'_kappa1_debyecutoff'+str(debyecutoff)+'_Langevin_wall'+str(ljdebye)+'_Kangle'+str(Kangle)+'_Kbond'+str(Kbond)+'_T'+str(T)+'_theta0is180_twofirst_are_fixed.png'
outfilename  = 'bond_lengths_chaingrid_quadratic_M'+str(M)+'N'+str(N)+'_ljdebye_epsilon'+str(epsilon)+'_sigma'+str(sigma)+'_ljcutoff'+str(ljcutoff)+'_kappa1_debyecutoff'+str(debyecutoff)+'_Langevin_wall'+str(ljdebye)+'_Kangle'+str(Kangle)+'_Kbond'+str(Kbond)+'_T'+str(T)+'_theta0is180_twofirst_are_fixed.txt'
outfilename2 = 'persistence_length_chaingrid_quadratic_'+str(M)+'N'+str(N)+'_ljdebye_epsilon'+str(epsilon)+'_sigma'+str(sigma)+'_ljcutoff'+str(ljcutoff)+'_kappa1_debyecutoff'+str(debyecutoff)+'_Langevin_wall'+str(ljdebye)+'_Kangle'+str(Kangle)+'_Kbond'+str(Kbond)+'_T'+str(T)+'_theta0is180_twofirst_are_fixed.txt'
outfilename3 = 'costhetas_neighbourbonds_chaingrid_quadratic_'+str(M)+'N'+str(N)+'_ljdebye_epsilon'+str(epsilon)+'_sigma'+str(sigma)+'_ljcutoff'+str(ljcutoff)+'_kappa1_debyecutoff'+str(debyecutoff)+'_Langevin_wall'+str(ljdebye)+'_Kangle'+str(Kangle)+'_Kbond'+str(Kbond)+'_T'+str(T)+'_theta0is180_twofirst_are_fixed.txt'
outfilename4 = 'ree_last_chaingrid_quadratic_'+str(M)+'N'+str(N)+'_ljdebye_epsilon'+str(epsilon)+'_sigma'+str(sigma)+'_ljcutoff'+str(ljcutoff)+'_kappa1_debyecutoff'+str(debyecutoff)+'_Langevin_wall'+str(ljdebye)+'_Kangle'+str(Kangle)+'_Kbond'+str(Kbond)+'_T'+str(T)+'_theta0is180_twofirst_are_fixed.txt'
outfilename5 = 'ree_average_chaingrid_quadratic_M'+str(M)+'N'+str(N)+'_ljdebye_epsilon'+str(epsilon)+'_sigma'+str(sigma)+'_ljcutoff'+str(ljcutoff)+'_kappa1_debyecutoff'+str(debyecutoff)+'_Langevin_wall'+str(ljdebye)+'_Kangle'+str(Kangle)+'_Kbond'+str(Kbond)+'_T'+str(T)+'_theta0is180_twofirst_are_fixed.txt'
outfilename6 = 'persistence_length_actualdistance_chaingrid_quadratic_M'+str(M)+'N'+str(N)+'_ljdebye_epsilon'+str(epsilon)+'_sigma'+str(sigma)+'_ljcutoff'+str(ljcutoff)+'_kappa1_debyecutoff'+str(debyecutoff)+'_Langevin_wall'+str(ljdebye)+'_Kangle'+str(Kangle)+'_Kbond'+str(Kbond)+'_T'+str(T)+'_theta0is180_twofirst_are_fixed.txt'
outfilename7 = 'endzs_chaingrid_quadratic_M'+str(M)+'N'+str(N)+'_ljdebye_epsilon'+str(epsilon)+'_sigma'+str(sigma)+'_ljcutoff'+str(ljcutoff)+'_kappa1_debyecutoff'+str(debyecutoff)+'_Langevin_wall'+str(ljdebye)+'_Kangle'+str(Kangle)+'_Kbond'+str(Kbond)+'_T'+str(T)+'_theta0is180_twofirst_are_fixed.txt'
outfilename8 = 'bondlengths_systemwide_chaingrid_quadratic_M'+str(M)+'N'+str(N)+'_ljdebye_epsilon'+str(epsilon)+'_sigma'+str(sigma)+'_ljcutoff'+str(ljcutoff)+'_kappa1_debyecutoff'+str(debyecutoff)+'_Langevin_wall'+str(ljdebye)+'_Kangle'+str(Kangle)+'_Kbond'+str(Kbond)+'_T'+str(T)+'_theta0is180_twofirst_are_fixed.txt'
'''


### Opening files
outfile3 = open(outfilename3,'w')
outfile4 = open(outfilename4,'w')
outfile5 = open(outfilename5,'w')
outfile6 = open(outfilename6,'w')

# The log file in LAMMPS contains a lot of lines detailing the initiation of the system. We want to access the information that comes after that.
linestart           = 6 #17            # The line the data starts at.
printeverynthstep   = 100
L                   = 10            # The size of our simulation box. LxLxL
startat             = 50            # To equilibrate. The number of ns we want to skip
dt                  = 0.00045       # The time step of our simulation. 0.00045 ns default for nano
skiplines           = 9             # If we hit 'ITEM:', skip this many steps...
skipelem            = 0#10#1000#10000#10000#90000 # The number of elements we skip in order to equilibrate (set to 0 if the .lammpstrj file should be equilibrated)
sampleevery         = 0#1 # Sample every 10th of the values written to file # The program gets WAY too slow if this is too small.
timefac             = dt*printeverynthstep*1e-9*sampleevery

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

#xes = np.zeros((N, ?))
#ys  = np.zeros((N, ?))
#zs  = np.zeros((N, ?))

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
xes.append(x_curr)
ys.append(y_curr)
zs.append(z_curr)

infile.close()

#print("xes:",xes)
#print("xes[0]:", xes[0])
#print("max(xes):", max(xes))

# SEEMS to work fine up until here.

# Finding the bond vectors for each time step.
# This will be an awfully blocky code, it seems
# Removing the first redundant element of position lists (uncomment if short test-file)

'''
xes.pop(0)
ys.pop(0)
zs.pop(0)
'''

Nt = len(xes) # Could have used j (Nt = j)

#print("xes:", xes)

#print("xes[-1]:",xes[-1])
#print("ys[-1]:",ys[-1])
#print("zs[-1]:",zs[-1])
#print("len, xes:", Nt)

#print("xes[0]:", xes[0])
#print("ys[0]:", ys[0])
#print("zs[0]:", zs[0])
#print("xes:", xes)

# Oh, vey.
bondvecs        = np.zeros((Nt, M, N-1, 3))
length_bondvecs = np.zeros((Nt, M, N-1))
ree_vec         = np.zeros((Nt, M))
ree2_vec        = np.zeros((Nt, M))
endzs		= np.zeros((Nt, M))   # I can use ree_dz to find this.

# A test:
infilename_zendtest = 'test_zend_from_find_persistencelength.txt'
infile_zendtest     = open(infilename_zendtest,'w')
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
        ree_dx  = 0
        ree_dy  = 0
        ree_dz  = 0
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
            ree_dx         += dx
            ree_dy         += dy
            ree_dz         += dz
            #tempvec = np.array([dx,dy,dz])
            #length_bondvecs[i,j] = np.linalg.norm(tempvec)
            veclen2              = dx**2+dy**2+dz**2
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
            bindex1   = bindex2
        ree_len2      = ree_dx**2+ree_dy**2+ree_dz**2
        ree_vec[i,k]  = np.sqrt(ree_len2)
        ree2_vec[i,k] = ree_len2
        endzs[i,k]    = ree_dz
        infile_zendtest.write('%.16f\n' % ree_dz)
    #print(ree_dz)
infile_zendtest.close()

print("Made bonds")
end_time = time.process_time()

print("Time:", end_time-start_time, " s")



# Finding the persistence length
# Want <cos(theta)> of sites that are at distance i apart
# We want this to be a time average
# Distance 0: <cos(theta)> = 1 # Can just set that and not calculate it
Nb = N-1
costheta = np.zeros((M,Nb))
counters = np.zeros((M,Nb))
costheta[:,0] = 1            # Distance 0: <cos(theta)> = 1 # Can just set that and not calculate it
costheta_allvalues   = []    # Sorted after distance between bonds
costheta_unpacked    = []    
distances_unpacked   = []
costheta_chainsorted = []    # costheta_allvalues for each chain

# This pfart is the bottleneck: # Aaaaand I'm adding another loop to the bottleneck. Great! # It also messes up the data structure. I think I need to add another list...
# OLD:
'''
for l in range(M):           # Looping over the chains
    costheta_allvalues = []
    for j in range(1,Nb):    # The distance between the bonds to be measured
        costheta_j = []
        for i in range(Nt):  # Looping over the time steps
            for k in range(Nb-j): # The position of the first bond
                bond1        = bondvecs[i,l,k]
                bond2        = bondvecs[i,l,k+j]
                distance     = np.sum(length_bondvecs[i,l,(k+1):(k+j+1)]) # Getting the distance between the two vectors
                length1      = length_bondvecs[i,l,k]
                length2      = length_bondvecs[i,l,k+j] 
                dotprod      = np.dot(bond1,bond2)         # I only find the angle between bond vectors in the same time step
                ct           = dotprod/(length1*length2)
                if(length1*length2==0):
                    a = 1
                    #print("ct:", ct, "; i =", i)
                costheta[l,j] += ct
                counters[l,j] += 1
                costheta_j.append(ct)
                costheta_unpacked.append(ct)
                distances_unpacked.append(distance)
                # Should I save all the contributions to costheta so that I can find the rms?
        costheta_allvalues.append(costheta_j) # Should I not have indented this?
    costheta_chainsorted.append(costheta_allvalues)
'''

costheta_allvalues_separation    = []
costheta_allvalues_chain         = []
costheta_chainsorted_separations = []
costheta_chainsorted_distances   = []
costheta_chainmaster             = []
costheta_timemaster              = []

# This part is the bottleneck: # Aaaaand I'm adding another loop to the bottleneck. Great! # It also messes up the data structure. I think I need to add another list...
for i in range(Nt):  # Looping over the time steps
    for l in range(M):           # Looping over the chains
        costheta_allvalues                     = []
        costheta_allvalues_unzipped            = []
        costheta_allvalues_separation_unzipped = []
        costheta_allvalues_distance_unzipped   = []
        for j in range(1,Nb):    # The distance between the bonds to be measured
            costheta_j = []
            for k in range(Nb-j): # The position of the first bond
                bond1        = bondvecs[i,l,k]
                bond2        = bondvecs[i,l,k+j]
                distance     = np.sum(length_bondvecs[i,l,(k+1):(k+j+1)]) # Getting the distance between the two vectors
                length1      = length_bondvecs[i,l,k]
                length2      = length_bondvecs[i,l,k+j] 
                dotprod      = np.dot(bond1,bond2)         # I only find the angle between bond vectors in the same time step
                ct           = dotprod/(length1*length2)
                if(length1*length2==0):
                    a = 1
                    #print("ct:", ct, "; i =", i)
                costheta[l,j] += ct
                counters[l,j] += 1
                costheta_j.append(ct)
                #costheta_unpacked.append(ct)
                #distances_unpacked.append(distance)
                costheta_allvalues_unzipped.append(ct)
                costheta_allvalues_separation_unzipped.append(j)
                costheta_allvalues_distance_unzipped.append(distance)
                # Should I save all the contributions to costheta so that I can find the rms?
            # These are not in use anymore:
            costheta_allvalues.append(costheta_j) # Should I not have indented this?
            costheta_allvalues_separation.append(j)
            costheta_allvalues_chain.append(l)
        # But these are:
        costheta_chainmaster.append(l)
        costheta_timemaster.append(i)
        costheta_chainsorted.append(costheta_allvalues_unzipped)
        costheta_chainsorted_separations.append(costheta_allvalues_separation_unzipped)
        costheta_chainsorted_distances.append(costheta_allvalues_distance_unzipped)
    # Do something similar to costheta_unpacked.
    #costheta_unpacked.append(costheta_)
#costheta_unpacked  = np.array(costheta_unpacked)
#distances_unpacked = np.array(distances_unpacked)

for j in range(M):
    for i in range(1,Nb):
        costheta[j,i] /= Nt*(Nb-i) # Could also have a counter, that's probably wise
        #print("counter[i]-(Nt*(Nb-i)):", counters[i]-(Nt*(Nb-i))) # Seems fine

costheta1_av_system  = np.mean(costheta[:,1])
costheta1_counter    = 0
# Finding the rms values:
costheta_rms         = np.zeros((M,Nb))
costheta1_rms_system = 0
#costheta_rms[0] = 0   # No variation here. # Don't actually need to set this
# Should I have binned the values? Some of the bins are really small, so I guess it is superfluous
thiscounter = 0
''' # Deprecated:
for l in range(Nt):
    for k in range(M):
        costheta_allvalues = costheta_chainsorted[k]
        thiscounter += 1
        for i in range(1,Nb):
            cai = costheta_allvalues[i-1]
            Nthis = len(cai)
            for j in range(Nthis):
                costheta_rms[k,i] += (costheta[k,i]-cai[j])**2    # Should work up until here...
            costheta_rms[k,i] =  np.sqrt(costheta_rms[k,i]/Nthis)
'''
Nthis = np.zeros((M,Nb))
for i1 in range(len(costheta_chainsorted)):
    costheta_allvalues      = costheta_chainsorted[i1]
    costheta_allseparations = costheta_chainsorted_separations[i1]
    k                       = costheta_chainmaster[i1]
    thiscounter += 1
    for i2 in range(len(costheta_allvalues)):
        cai = costheta_allvalues[i2]
        i = costheta_allseparations[i2]
        Nthis[k,i] += 1 
        costheta_rms[k,i] += (costheta[k,i]-cai)**2    # Should work up until here...
        if i==1:
            costheta1_rms_system += (costheta1_av_system-cai)**2
            costheta1_counter    += 1

costheta1_rms_system      = np.sqrt(costheta1_rms_system/(costheta1_counter-1))
for k in range(M):
    for i in range(1,Nb):
        costheta_rms[k,i] = np.sqrt(costheta_rms[k,i]/(Nthis[k,i]-1))


print("Found cos(theta)")

for i in range(M):
    outfile3.write('%.16f %.16f\n' % (costheta[i,1],costheta_rms[i,1]))
outfile3.write('System wide: %.16f %.16f' % (costheta1_av_system,costheta1_rms_system))
outfile3.close()



#### Write ree last to file
for i in range(M):
    outfile4.write('%.16e\n' % (ree_vec[-1,i]))

outfile4.close()
#### Finding ree average and stdv                   # This might take some time
ree_av         = np.zeros(M)
ree2_av        = np.zeros(M)
ree_tot_av     = 0
ree2_tot_av    = 0

for i in range(M):
    ree_av[i]  = np.mean(ree_vec[:,i])
    ree2_av[i] = np.mean(ree2_vec[:,i])

ree_tot_av     = np.mean(ree_av)
ree2_tot_av    = np.mean(ree2_av)

#print("ree_vec:", ree_vec[:,0])
#print("That was ree_vec[:,0]")

# Stdv
ree_stdv       = np.zeros(M)
ree2_stdv      = np.zeros(M)
ree_tot_stdv   = 0
ree2_tot_stdv  = 0
for j in range(M):
    for i in range(Nt):
        ree_stdv[j]   += (ree_vec[i,j]  - ree_av[j])**2
        ree2_stdv[j]  += (ree2_vec[i,j] - ree2_av[j])**2
        ree_tot_stdv  += (ree_vec[i,j]  - ree_tot_av)**2
        ree2_tot_stdv += (ree2_vec[i,j] - ree2_tot_av)**2
    ree_stdv[j]  = np.sqrt(ree_stdv[j]/float(Nt-1))
    ree2_stdv[j] = np.sqrt(ree2_stdv[j]/float(Nt-1))
ree_tot_stdv  = np.sqrt(ree_tot_stdv/(M*Nt-1))
ree2_tot_stdv = np.sqrt(ree2_tot_stdv/(M*Nt-1))

for i in range(M):
    outfile5.write('<r_ee>: %.16e %.16e <r_ee^2>: %.16e %.16e\n' % (ree_av[i],ree_stdv[i],ree2_av[i],ree2_stdv[i]))
outfile5.write('System total: %.16e %.16e %.16e %.16e' % (ree_tot_av,ree_tot_stdv,ree2_tot_av,ree2_tot_stdv))
outfile5.close()


# For plotting:
separation         = np.arange(Nb)
# Make each bond count as separation 1
persistencelength  = np.zeros(M)
pl_stdv            = np.zeros(M)
allvals_pl         = np.zeros((M,Nt)) 
pl_stdv_all        = 0
pl_all             = 0         
# Make each bond count with its length
persistencelength2 = np.zeros(M)
pl_stdv2           = np.zeros(M)
allvals_pl2        = np.zeros((M,Nt)) 
pl_stdv_all2       = 0
pl_all2            = 0 

#Ntimeav = M*Nt/100. # 100 bins
for k in range(len(costheta_chainmaster)):
    i = costheta_chainmaster[k]
    j = costheta_timemaster[k]
    costhetas_this        = costheta_chainsorted[k]
    separations_this      = costheta_chainsorted_separations[k]
    distances_this        = costheta_chainsorted_distances[k]
    # First: Separation fit:
    popt, pcov            = curve_fit(costheta_exponential, separations_this, costhetas_this)
    this_pl               = popt[0]
    uncertainty           = np.sqrt(pcov[0])
    persistencelength[i] += this_pl
    #pl_stdv[i]           = np.sqrt(pcov[0])
    allvals_pl[i,j]       = this_pl
    '''
    print("This persistence length:", this_pl, "; pcov:", uncertainty)
    if this_pl>20:
        separations_this = np.array(separations_this)
        fittedgraph = costheta_exponential(separations_this, this_pl)
        x_eline = np.zeros(2)
        y_eline = np.zeros(2)+ 1./np.exp(1)
        x_eline[0] = 0
        x_eline[1] = 100
        plt.figure(figsize=(6,5))
        plt.plot(separations_this, costhetas_this, label='Values')
        plt.plot(separations_this, fittedgraph, label='Fit')#plt.plot(x_persistencelength, y_persistencelength, '--')#plt.plot(x_eline, y_eline, '--')
        plt.plot(x_eline, y_eline, '--')
        plt.xlabel(r'Bond distance', fontsize=16)
        plt.ylabel(r'<cos($\theta$)>', fontsize=16)
        plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
        plt.legend(loc="upper right")
        plt.title(r'<cos($\theta$)> vs separation (last chain)', fontsize=16)
        plt.show()
    '''
    # Then: Distance fit:
    popt2, pcov2           = curve_fit(costheta_exponential, distances_this, costhetas_this)
    this_pl2               = popt2[0]
    persistencelength2[i] += this_pl2
    #pl_stdv[i]           = np.sqrt(pcov[0])
    allvals_pl2[i,j]       = this_pl2
persistencelength /=float(Nt)
persistencelength2/=float(Nt)

stdv_test = False
counter   = 1
runcumul  = 0
runningav = 0
pl_all  = np.mean(persistencelength)
pl_all2 = np.mean(persistencelength2)
for i in range(M):
    for j in range(Nt):
        # Separations:
        pl_stdv[i]   += (persistencelength[i]-allvals_pl[i,j])**2
        pl_stdv_all  += (pl_all-allvals_pl[i,j])**2
        if stdv_test==True:
            runcumul  += allvals_pl[i,j]
            runningav = runcumul/counter 
            if allvals_pl[i,j]>20:
                print("l:,",i,"t:",j,"stdv,",counter," values:",np.sqrt(pl_stdv_all/float(counter)), "; pl_av =", pl_all, "; this pl-value:", allvals_pl[i,j], "av. of these:", runningav )
        counter+=1
        # Distances:
        pl_stdv2[i]  += (persistencelength2[i]-allvals_pl2[i,j])**2
        pl_stdv_all2 += (pl_all2-allvals_pl2[i,j])**2
    pl_stdv[i]        = np.sqrt(pl_stdv[i]/float(Nt-1))
    pl_stdv2[i]       = np.sqrt(pl_stdv2[i]/float(Nt-1))
pl_stdv_all           = np.sqrt(pl_stdv_all/float(M*Nt-1))
pl_stdv_all2          = np.sqrt(pl_stdv_all2/float(M*Nt-1))
fittedgraph = costheta_exponential(separation, persistencelength[-1]) # I don't need THAT many graphs

# Standard deviation in the mean persistence length:
pl_stdv_mean  = 0
pl_stdv_mean2 = 0
for i in range(M):
    pl_stdv_mean  += (persistencelength[i]-pl_all)**2
    pl_stdv_mean2 += (persistencelength2[i]-pl_all2)**2
pl_stdv_mean  = np.sqrt(pl_stdv_mean/(M-1))
pl_stdv_mean2 = np.sqrt(pl_stdv_mean2/(M-1))

print("pl_stdv_mean:", pl_stdv_mean)
print("pl_stdv_mean2:", pl_stdv_mean2)


lastpls      = allvals_pl[:,-1]
mean_lastpls = np.mean(lastpls)
lastpls_rms  = 0
for i in range(M):
    lastpls_rms += (lastpls[i]-mean_lastpls)**2
lastpls_rms = np.sqrt(lastpls_rms/(M-1))

lastpls2      = allvals_pl2[:,-1]
mean_lastpls2 = np.mean(lastpls2)
lastpls_rms2  = 0
for i in range(M):
    lastpls_rms2 += (lastpls2[i]-mean_lastpls2)**2
lastpls_rms2 = np.sqrt(lastpls_rms2/(M-1))

print("rms from the last values:", lastpls_rms)
print("Nt:", Nt)
'''
# Use the actual bond length # I think maybe this is too computationally heavy. I guess I should try once to see if the results are fairly similar?
popt2, pcov2       = curve_fit(costheta_exponential, distances_unpacked, costheta_unpacked)
persistencelength2 = popt2[0]
pl_stdv2           = np.sqrt(pcov2[0])
'''
fittedgraph2 = costheta_exponential(separation, persistencelength2[-1])


print("K:", K)
print("Persistence length (from separation fit):", np.mean(pl_all), " (mean of chains)")
print("Persistence length 2 (from distance fit):", pl_all2, "\ndiff pl 1 and pl2:", (pl_all-pl_all2))
print("Percent deviation, diff pl 1 and pl2:", (pl_all-pl_all2)/pl_all)

end_time = time.process_time()
print("Time:", end_time-start_time, " s")

#print("costheta_rms:", costheta_rms)

print(np.mean(length_bondvecs))

# Persistence length line
x_persistencelength = np.zeros(2) + persistencelength[-1]
y_persistencelength = np.zeros(2)
y_persistencelength[0] = 0
y_persistencelength[1] = 1

# e line:
x_eline = np.zeros(2)
y_eline = np.zeros(2)+ 1./np.exp(1)
x_eline[0] = 0
x_eline[1] = 100

############## FIX THE PLOTS AT SOME POINT!!!
# Maybe plotting, printing and fitting?
plt.figure(figsize=(6,5))
plt.errorbar(separation, costheta[-1,:], yerr=costheta_rms[-1,:], fmt="none", capsize=2, label='Values')
plt.plot(separation, costheta[-1,:], 'o')
if arctest==False:
    plt.plot(separation, fittedgraph, label='Fit')
else:
    plt.plot(separation,np.cos(np.pi*separation/(2*Nb)), label='cos(pi*s/2N)')
#plt.plot(x_persistencelength, y_persistencelength, '--')
#plt.plot(x_eline, y_eline, '--')
plt.xlabel(r'Bond distance', fontsize=16)
plt.ylabel(r'<cos($\theta$)>', fontsize=16)
plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
plt.legend(loc="upper right")
plt.title(r'<cos($\theta$)> vs separation (last chain)', fontsize=16)
plt.savefig(plot1name)

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
plt.title(r'<cos($\theta$)> vs separation', fontsize=16)
plt.savefig(plot3name)

'''
plt.figure(figsize=(6,5))
plt.plot(costheta_allvalues[5], '.')
plt.xlabel(r'Index', fontsize=16)
plt.ylabel(r'cos($\theta$)', fontsize=16)
plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
plt.title(r'Values of cos($\theta$) for separation 5', fontsize=16)
plt.savefig(plot2name)
'''
outfile  = open(outfilename,'w')

# Have a costheta chain- and system file too.
bondlength_chain_vals     = np.zeros(M)
bondlength_chain_vals_rms = np.zeros(M)
bondlength_av_system      = np.mean(length_bondvecs)
bondlength_rms_system     = 0 


outfile.write('Nt: %i Nbonds: %i\n' % (Nt, N-1))
for i in range(Nt):
    outfile.write('Time: %.16e i: %i\n' % (times[i],i))
    for k in range(M):
        for j in range(N-1):
            #print("length_bondvecs[i,k,j]:",length_bondvecs[i,k,j])
            bondlength_chain_vals[k] += length_bondvecs[i,k,j]
            bondlength_rms_system    += (bondlength_av_system-length_bondvecs[i,k,j])**2
            outfile.write('%i %i %.16e\n' % (j+1, k,length_bondvecs[i,k,j]))
bondlength_chain_vals /= Nt*(N-1)
bondlength_rms_system  = np.sqrt(bondlength_rms_system/(M*Nt*(N-1)-1))
print('M*Nt*(N-1)-1:',M*Nt*(N-1)-1)
outfile.close()

print('SEARCHING FOR THE ERROR IN BOND LENGTH RMS:')
print('BOND LENGTH AV. FOR CHAIN:',bondlength_chain_vals)
for i in range(Nt):
    for k in range(M):
        for j in range(N-1):
            #print("length_bondvecs[i,k,j]:",length_bondvecs[i,k,j])
            bondlength_chain_vals_rms[k] += (bondlength_chain_vals[k]-length_bondvecs[i,k,j])**2 # Do I get an overflow?
            if (bondlength_chain_vals[k]-length_bondvecs[i,k,j])**2!=0.0:
                print('Contrib:',(bondlength_chain_vals[k]-length_bondvecs[i,k,j])**2)
                print('accumulated:', bondlength_chain_vals_rms[k])
for i in range(M):
    bondlength_chain_vals_rms[i] = np.sqrt(bondlength_chain_vals_rms[k]/(Nt*(N-1)-1))#(M*Nt*(N-1)-1)) # Shouldn't divide by M here...
print('BOND LENGTH AV. FOR CHAIN:',bondlength_chain_vals)
print('Nt*(N-1)-1:',Nt*(N-1)-1)

outfile8 = open(outfilename8,'w')
for i in range(M):
    outfile8.write('%.16e %.16e\n' % (bondlength_chain_vals[i],bondlength_chain_vals_rms[i]))
outfile8.write('System wide: %.16e %.16e' % (bondlength_av_system, bondlength_rms_system))

outfile2 = open(outfilename2,'w')
outfile2.write('The average persistence length for all chains: K=%i: %.16e %.16e\n' % (K, pl_all, pl_stdv_all))
for i in range(M):
    outfile2.write('%.16e %.16e %.16e\n' % (persistencelength[i], pl_stdv[i], persistencelength[i]*np.mean(length_bondvecs[:,i,:])))
outfile2.write('RMS of average persistence length of chains: %.16e (Taking the rms of the average chain persistencelengths)\n' % pl_stdv_mean) # I can use this to see if one chain is shorter than the others. Or I could just look at the chains lol

outfile2.close()

outfile6.write('The average persistence length2 for all chains: K=%i: %.16e %.16e\n' % (K, pl_all2, pl_stdv_all2))
for i in range(M):
    outfile6.write('%.16e %.16e\n' % (persistencelength2[i], pl_stdv2[i]))
outfile6.write('RMS of average persistence length2 of chains: %.16e (Taking the rms of the average chain persistencelengths)\n' % pl_stdv_mean2)
outfile6.close()

outfile7 = open(outfilename7, 'w')

print("Largest end z:",max(endzs[0]))
print("Largest end z:",max(endzs[1]))
print("Largest end z:",max(endzs[2]))
print("Largest end z:",max(endzs[3]))
print("Largest end z:",max(endzs[4]))
print("Largest end z:",max(endzs[5]))
print("Largest end z:",max(endzs[6]))
print("Largest end z:",max(endzs[7]))
print("Largest end z:",max(endzs[8]))

print("End zs:",endzs[1])

lastz_by_chain = np.zeros(M)
for i in range(Nt):
    thistimestep_lastz = endzs[i]
    for j in range(M):
        lastz_by_chain[j] += thistimestep_lastz[j] # Could take np.mean instead, maybe?
lastz_by_chain /= Nt
lastz_totalav   = np.mean(lastz_by_chain)

lastz_totalrms     = 0
lastz_by_chain_rms = np.zeros(M)
for i in range(Nt):
    thistimestep_lastz = endzs[i]
    for j in range(M):
        lastz_by_chain_rms[j] += (lastz_by_chain[j]-thistimestep_lastz[j])**2
        lastz_totalrms        += (lastz_totalav-thistimestep_lastz[j])**2
lastz_totalrms = np.sqrt(lastz_totalrms/(M*Nt-1))

for j in range(M):
    lastz_by_chain_rms[j] = np.sqrt(lastz_by_chain_rms[j]/Nt)
    outfile7.write('%.16e %.16e\n' % (lastz_by_chain[j],lastz_by_chain_rms[j]))
outfile7.write('Total average, end: %.16e %.16e' % (lastz_totalav,lastz_totalrms))

outfile7.close()
print("Script finished.")