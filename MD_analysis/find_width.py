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

# Max number of iterations in curve_fit: # (Use default unless it crashes)
thismaxfev  = 5000#400# # (default appears to be 400)

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
spacing     = 10#4#0
gridspacing = spacing
K           = Kangle # Because we used this notation earlier, but using Kangle is less confusing
wallenergy  = 1.042
dielectric  = 1#50 #1 2 10 100
systemtest  = False
lotsofchains = False
flubbedup    = False

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
#'''   # This one is for varying the factor.
# Kbond     = 2000:
#basename     = 'chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
#basename     = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed'  % (M,N,spacing,wallenergy,Kangle,Kbond,charge,T)
#basename     = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,spacing,Kangle,Kbond,charge,T)
#basename     = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,Kangle,Kbond)+ str(charge) +'_T%i_theta0is180_twofirst_are_fixed' % T
#basename     = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,Kangle,Kbond)+ str(factor) +'_T%i_theta0is180_twofirst_are_fixed' % T
#basename     = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,Kangle,Kbond,factor,T)
#basename     = 'chaingrid_quadratic_M%iN%i_ljdebye%.3f_angle_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180' % (M,N,ljdebye,ljdebye,Kangle,Kbond,T)
#basename     = 'chaingrid_quadratic_M'+str(M)+'N'+str(N)+'_ljdebye_epsilon'+str(epsilon)+'_sigma'+str(sigma)+'_ljcutoff'+str(ljcutoff)+'_kappa1_debyecutoff'+str(debyecutoff)+'_Langevin_wall'+str(ljdebye)+'_Kangle'+str(Kangle)+'_Kbond'+str(Kbond)+'_T'+str(T)+'_theta0is180_twofirst_are_fixed'     # Actual run
#basename     = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed' % (M,N,Kangle,Kbond,factor,T)
#basename2    = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed' % (M,N,gridspacing,Kangle,Kbond,factor,T)
#basename      = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa%i_debyecutoff%i_charge%i_mass%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,gridspacing,Kangle,Kbond,kappa,debyecutoff,charge,mass,T)
#basename      = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa%i_debyecutoff%i_charge%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,gridspacing,Kangle,Kbond,kappa,debyecutoff,charge,T)
#basename      = 'chaingrid_quadratic_M%iN%i_gridspacing%i_lj_epsilon%.3f_sigma%i_ljcutoff%.14f_kappa1_Langevin_Kangle%i_Kbond%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,gridspacing,epsilon,sigma,ljcutoff,Kangle,Kbond,T)
#basename      = 'chaingrid_quadratic_M%iN%i_gridspacing%i_ljdebye_epsilon%.3f_sigma%i_ljcutoff%.14f_kappa1_debyecutoff3_Langevin_Kangle%i_Kbond%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,gridspacing,epsilon,sigma,ljcutoff,Kangle,Kbond,T) # charge-1 on this one
#basename      = 'chaingrid_quadratic_M%iN%i_gridspacing%i_ljdebye_epsilon%.3f_sigma%i_ljcutoff%.14f_kappa1_debyecutoff3_Langevin_Kangle%i_Kbond%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,gridspacing,epsilon, sigma, ljcutoff,Kangle,Kbond,T)
#basename2     = 'chaingrid_quadratic_M%iN%i_gridspacing%i_lj_epsilon1p042_sigma%i_ljcutoff1p12246204830937_kappa1_Langevin_Kangle%i_Kbond%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,gridspacing,sigma,Kangle,Kbond,T)
#basename      = 'chaingrid_quadratic_M%iN%i_gridspacing%i_ljdebye_epsilon%.3f_sigma%i_ljcutoff%.14f_kappa1_debyecutoff3_charge%i_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,gridspacing,epsilon,sigma,ljcutoff,charge,epsilon,Kangle,Kbond,T) # charge-1 on this one
#basename2     = 'tester'
#basename      = 'chaingrid_quadratic_M1N101_gridspacing300_Langevin_Kbond200_wall1.042_T310_theta0is180_twofirst_are_fixed' #'onechain_M1N101_gridspacing300_Langevin_wall1.042_KfeneR0sq200_1.2_T310_theta0is180_twofirst_are_fixed'
#'onechain_M1N101_gridspacing300_Langevin_wall1.042_KfeneR0sq200_1.2_bondepsilon1.042_fenesigma0.8_T310_theta0is180_twofirst_are_fixed'


#'_onechain_M1N101_gridspacing300_Langevin_wall1.042_KfeneR0sq200_1.2_T310_theta0is180_twofirst_are_fixed.txt', '_onechain_M1N101_gridspacing300_Langevin_wall1.042_KfeneR0sq200_1.2_bondepsilon1.042_fenesigma0.8_T310_theta0is180_twofirst_are_fixed.txt']

# Just Kangle and Kbond # Spacing is 40
'''
M = 9
N = 101
T = 310
Kangle = 20
Kbond  = 2000
factor = Kangle/float(Kbond)
basename  = 'chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed' % (M,N,Kangle,Kbond,factor,T)
basename2 = basename
#chaingrid_quadratic_M9N101_Langevin_Kangle1000_Kbond2000_factor0.50_T310_theta0is180_twofirst_are_fixed.lammpstrj
#chaingrid_quadratic_M9N101_Langevin_Kangle1000_Kbond2000_factor0.5_T310_theta0is180_twofirst_are_fixed.lammpstrj
#chaingrid_quadratic_M9N101_Langevin_Kangle100_Kbond2000_factor0.05_T310_theta0is180_twofirst_are_fixed.lammpstrj
#chaingrid_quadratic_M9N101_Langevin_Kangle2000_Kbond2000_factor1.00_T310_theta0is180_twofirst_are_fixed.lammpstrj
#chaingrid_quadratic_M9N101_Langevin_Kangle200_Kbond2000_factor0.10_T310_theta0is180_twofirst_are_fixed.lammpstrj
#chaingrid_quadratic_M9N101_Langevin_Kangle20_Kbond2000_factor0.01_T310_theta0is180_twofirst_are_fixed.lammpstrj
#chaingrid_quadratic_M9N101_Langevin_Kangle20_Kbond2000_factor100_T310_theta0is180_twofirst_are_fixed.lammpstrj
'''

M = 1
basename  = 'tester_totally_straight_tiltedonly'#'tester_totally_straight'
basename2 = basename
#flubbedup = True

M = 9
basename = 'testsystem_chaingrid_straight'
basename = 'testsystem_chaingrid_straight_differentorientations'
basename = 'testsystem_chaingrid_straight_differentorientations_varying_with_time'
basename2 = basename
## FENE
'''
M             = 1
KfeneR0sq     = 200#67.5 # 200
Kangle        = 20#2000 # 200 # 2000
R0            = 1.2#1.5#1.2
sigma         = 0.8#1.0#0.8
basename      = 'onechain_M%iN%i_gridspacing%i_Langevin_nowall_Kangle%i_KfeneR0sq%i_R0%.1f_feneepsilon%.3f_fenesigma%.1f_T%i_theta0is180_twofirst_are_fixed' % (M,N,spacing,Kangle,KfeneR0sq,R0,epsilon,sigma,T)
#basename      = 'onechain_M%iN%i_gridspacing%i_Langevin_nowall_Kangle%i_KfeneR0sq%i_R0%.1f_feneepsilon%.3f_fenesigma%.1f_T%i_theta0is180_twofirst_are_fixed' % (M,N,spacing,Kangle,KfeneR0sq,R0,epsilon,sigma,T)
#basename      = 'onechain_M1N101_gridspacing300_Langevin_wall1.042_KfeneR0sq67.5_1.5_bondepsilon1.042_fenesigma1_T310_theta0is180_twofirst_are_fixed'
#basename      = 'onechain_M1N101_gridspacing300_Langevin_wall1.042_KfeneR0sq200_1.2_bondepsilon1.042_fenesigma0.8_T310_theta0is180_twofirst_are_fixed'
basename       = 'onechain_M1N101_gridspacing300_Langevin_wall1.042_KfeneR0sq67.5_1.5_bondepsilon1.042_fenesigma1_T310_theta0is180_twofirst_are_fixed' # 'onechain_M1N101_gridspacing300_Langevin_wall1.042_KfeneR0sq200_1.2_T310_theta0is180_twofirst_are_fixed'
#basename      = 'onechain_M%iN%i_gridspacing%i_Langevin_nowall_Kangle%i_KfeneR0sq%i_R0%.1f_fenesigma%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,spacing,Kangle,KfeneR0sq,R0,sigma,T)
basename2     = basename
'''

## Using Lennard-Jones units:
'''
M             = 9
spacing       = 7
gridspacing   = spacing
Kangle        = 14.0186574854529#20
Kbond         = 140.186574854529#200
K             = Kangle
T             = 3
effectivedielectric = 0.00881819074717447
##
# chaingrid_quadratic_M9N101_ljunits_gridspacing10_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_chargeelementary-1_T3_theta0is180_twofirst_are_fixed
#basename      = 'chaingrid_quadratic_M%iN%i_ljunits_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_chargeelementary%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,spacing,Kangle,Kbond,charge,T)
#basename      = 'chaingrid_quadratic_M%iN%i_ljunits_gridspacing%i_Langevin_Kangle%.13f_Kbond%.12f_debye_kappa1_debyecutoff3_chargeelementary%i_T%i_theta0is180_correctconversions_twofirst_are_fixed' % (M,N,spacing,Kangle,Kbond,charge,T)
basename      = 'chaingrid_quadratic_M%iN%i_ljunits_gridspacing%i_Langevin_Kangle%.13f_Kbond%.12f_debye_kappa1_debyecutoff3_chargeelementary%i_effectivedielectric%.17f_T%i_theta0is180' % (M,N,spacing,Kangle,Kbond,charge,effectivedielectric,T)
basename2     = basename 
#basename2     = 'chaingrid_quadratic_M%iN%i_gridspacing%i_ljdebye_epsilon%.3f_sigma%i_ljcutoff%.14f_kappa1_debyecutoff3_charge%i_Langevin_Kangle%i_Kbond%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,gridspacing,epsilon,sigma,ljcutoff,charge,Kangle,Kbond,T) # charge-1 on this one
#'''
## FENE for tethered chains ##
'''
M           = 1
spacing     = 300
R0          = 1.5
T           = 310
dt          = 0.0001
KfeneR0sq   = 200
bondepsilon = 1.042
fenesigma   = 1
gridspacing = spacing
#basename    = 'onechain_M1N101_gridspacing300_Langevin_KfeneR0sq200_R0is1.5_bondepsilon1.042_fenesigma1_T310_theta0is180_twofirst_are_fixed'
basename    = 'onechain_M%iN%i_gridspacing%i_Langevin_KfeneR0sq%i_R0is%.1f_bondepsilon%.3f_fenesigma%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,spacing,KfeneR0sq,R0,bondepsilon,fenesigma, T)
#basename    = 'onechain_M%iN%i_gridspacing%i_Langevin_KfeneR0sq%i_R0is%.1f_T%i_theta0is180_dt%.4f_twofirst_are_fixed' % (M,N,spacing,KfeneR0sq,R0, T,dt)
#basename    = 'onechain_M%iN%i_gridspacing%i_Langevin_KfeneR0sq%i_R0is%.1f_bondepsilon%.3f_fenesigma%i_T%i_theta0is180_dt%.4f_twofirst_are_fixed' % (M,N,spacing,KfeneR0sq,R0, bondepsilon, fenesigma, T,dt)
#basename    = 'onechain_M%iN%i_gridspacing%i_Langevin_KfeneR0sq%i_R0is%.1f_bondepsilon%.3f_fenesigma%i_T%i_theta0is180_dt%.4f_specialbonds_twofirst_are_fixed' % (M,N,spacing,KfeneR0sq,R0, bondepsilon, fenesigma, T,dt)
basename2   = basename
#'''
##                        ##

'''
## FENE testing ## # Yeesh... I have to do this for two bonds if I want to use this program
M = 1
N = 3
Nt = 10
basename  = 'twobonds_test_gridspacing300_Langevin_Kfene30_R0is1.5_bondepsilon1_fenesigma1_T310_dt0.0001_specialbonds011'
basename2 = basename
##              ##
'''

# LJ, Debye, bending, FENE in LJ units:
'''
KfeneR0sq   = 67.5
R0          = 1.5
Kangle      = 5#14.0186574854529
spacing     = 4
gridspacing = spacing
T           = 3
effectivedielectric = 0.00881819074717447
M           = 9
flubbedup   = False
basename  = 'chgr_quadr_M%iN%i_ljun_gridsp%i_Langevin_wall_Kangle%i_KfR0sq%.1f_R0%.1f_feps%.3f_fsgm%i_ljeps%.3f_ljsgm%i_debye_kappa%i_debcut%i_elch%i_effdiel%.17f_T%i_theta0is180' % (M,N,spacing,Kangle,KfeneR0sq,R0,epsilon,sigma,epsilon,sigma, kappa, debyecutoff,charge,effectivedielectric,T)
basename2 = basename
#basename = 'onechain_M1N101_ljunits_gridspacing300_Langevin_wall_Kangle14.0186574854529_KfeneR0sq67.5_R01.5_feneepsilon1.042_fenesigma1_debye_kappa1_debyecutoff3_chargeelementary-1_effectivedielectric0.00881819074717447_T3_theta0is180'
# chaingrid_quadratic_M9N101_ljunits_gridspacing15_Langevin_nowall_Kangle14.0186574854529_KfeneR0sq67.5_R01.5_feneepsilon1.042_fenesigma1_debye_kappa1_debyecutoff3_chargeelementary-1_effectivedielectric0.00881819074717447_T3_theta0is180.lammpstrj
# chaingrid_quadratic_M9N101_ljunits_gridspacing100_Langevin_wall_Kangle14.0186574854529_KfeneR0sq67.5_R01.5_feneepsilon1.042_fenesigma1_debye_kappa1_debyecutoff3_chargeelementary-1_effectivedielectric0.00881819074717447_T3_theta0is180
#basename    = 'chaingrid_quadratic_M%iN%i_ljunits_gridspacing%i_Langevin_wall_Kangle%.13f_KfeneR0sq%.1f_R0%.1f_feneepsilon%.3f_fenesigma%i_debye_kappa%i_debyecutoff%i_chargeelementary%i_effectivedielectric%.17f_T%i_theta0is180' % (M,N,spacing,Kangle,KfeneR0sq,R0,epsilon,sigma,kappa,debyecutoff,charge,effectivedielectric,T)
#basename2  = 'chaingrid_quadratic_M%iN%i_ljunits_gridspacing%i_Langevin_wall_Kangle%.13f_KfeneR0sq%.1f_R0%.1f_feneepsilon%.3f_fenesigma%i_debye_kappa%i_debyecutoff%i_chargeelementary%i_debconv%.17f_T%i' % (M,N,spacing,Kangle,KfeneR0sq,R0,epsilon,sigma,kappa,debyecutoff,charge,effectivedielectric,T)
#M           = 1
#basename    = 'onech_M1N101_ljun_gridsp300_Langevin_wall_Kangle5_KfR0sq67.5_R01.5_feps1.042_fsgm1_ljeps1.042_ljsgm1_debye_kappa1_debcut3_elch-1_effdiel0.00881819074717447_T3_theta0is180'
#basename    = 'chgr_quadr_M9N101_ljun_gridsp300_Langevin_wall_Kangle200_KfR0sq67.5_R01.5_feps1.042_fsgm1_ljeps1.042_ljsgm1_debye_kappa1_debcut3_elch-1_effdiel0.00881819074717447_T3_theta0is180'
#basename    = 'chaingrid_quadratic_M9N101_ljunits_gridspacing10_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_chargeelementary-1_T3_theta0is180_twofirst_are_fixed'
#basename    = 'onechain_M1N101_gridspacing300_Langevin_nowall_Kangle20_KfeneR0sq67.5_R01.5_feneepsilon1.042_fenesigma1_T310_theta0is180_twofirst_are_fixed'
#basename2   = 'test'
#'''

# Kbond     = 200:
#basename     = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.1f_T%i_theta0is180_twofirst_are_fixed' % (M,N,gridspacing,Kangle,Kbond,factor,T)
#basename2    = basename+'_test'
infilename   = basename+'.lammpstrj'
plot1name    = 'omega_vs_h_'+basename2+'.png'
plot2name    = 'omega_vs_sqrth_'+basename2+'.png'
outfilename  = 'omega_vs_h_'+basename2+'.txt'
#'''

# Tester:
'''
arctest      = True
systemtest   = True
testtype     = 'rotatingarc_withtime'#'rotatingarc_withtime'#'rotatingarc'#'quartercirc'#'straight'
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
#'''

# Tester 2:
'''
N = 100
systemtest   = True
arctest      = False
foldername   = 'MC_for_testing/'
#p = Path(foldername)
#p.mkdir(exist_ok=True)
basename     = 'M%iN%i_1000trials_freelyjointedchains' % (M,N) # /home/kine/Projects_PhD/P2_PolymerMD'
basename2    = basename
testtype     = 'MC_for_testing/'+basename
infilename   = testtype+'.lammpstrj'
plot1name    = 'omega_vs_h_'+basename+'.png'
plot2name    = 'omega_vs_sqrth_'+basename2+'.png'
outfilename  = 'omega_vs_h_'+basename+'.txt'
N = 101
#'''


# Output names for code testing: # go here!
'''
Kangle = 2000.
Kbond  = 2000.
factor = Kangle/Kbond
infilename   = 'chaingrid_quadratic_M9N101_Langevin_Kangle%i_Kbond%i_factor%.2f_T310_theta0is180_twofirst_are_fixed.lammpstrj' % (Kangle,Kbond,factor)

#infilename   = 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,spacing,Kangle,Kbond,charge,T)
#'chaingrid_quadratic_M9N101_gridspacing40_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T310_theta0is180_twofirst_are_fixed.lammpstrj' % (M,N,spacing,Kangle,Kbond,charge,T)
plot1name    = 'costheta_test_short_Kangle%i.png' % Kangle
plot2name    = 'angledistr_test_short_Kangle%i.png' % Kangle
plot3name    = 'costheta_actualdistance_test_short_Kangle%i.png' % Kangle
outfilename  = 'bond_lengths_test_short_Kangle%i.txt' % Kangle

outfilename2 = 'persistence_length_test_short_Kangle%i.txt' % Kangle
outfilename3 = 'costhetas_neighbourbonds_test_short_Kangle%i.txt' % Kangle
outfilename4 = 'ree_last_test_short_Kangle%i.txt' % Kangle
outfilename5 = 'ree_average_test_short_Kangle%i.txt' % Kangle
outfilename6 = 'persistence_length_actualdistance_test_short_Kangle%i.txt' % Kangle 
outfilename7 = 'endzs_test_short_Kangle%i.txt' % Kangle
outfilename8 = 'bondlengths_systemwide_short_Kangle%i.txt' % Kangle
#'''

# Varying the charge:
#chaingrid_quadratic_M9N101_Langevin_Kangle${Kangle}_Kbond${Kbond}_charge${charge}_T$T_theta0is180_twofirst_are_fixed.lammpstrj

print('spacing:', gridspacing)


### Opening files
outfile = open(outfilename,'w')
'''
outfile3 = open(outfilename3,'w')
outfile4 = open(outfilename4,'w')
outfile5 = open(outfilename5,'w')
outfile6 = open(outfilename6,'w')
'''

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

indexshift = 0
if M==1 or flubbedup==True:
    indexshift = -1

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
        x     = float(words[3+indexshift])
        y     = float(words[4+indexshift])
        z     = float(words[5+indexshift])
        #print('3+indexshift:', 3+indexshift)
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

#print("xes:",xes)
#print("xes[0]:", xes[0])
#print("xes[-1]:", xes[-1])
#print("max(xes):", max(xes))
print('Done with read-in loop')
end_time = time.process_time()
print("Time:", end_time-start_time, " s")

Nt = len(xes) # Could have used j (Nt = j)

# Oh, vey.  ## Not sure I need all these...
bondvecs        = np.zeros((Nt, M, N-1, 3))
length_bondvecs = np.zeros((Nt, M, N-1))
ree_vec         = np.zeros((Nt, M))
ree2_vec        = np.zeros((Nt, M))
endzs		= np.zeros((Nt, M))   # I can use ree_dz to find this.


# Should I have a main loop over time...? Maybe it saves me some time due to storing? # I guess I want things stored like this anyways... 
## Need to find end-to-end distance. Must take the (possibly) periodic boundary conditions into account.
bindex1  = 0
bindex2  = 0
rees     = np.zeros((Nt,M,3)) # End-to-end distance (I don't really need this anyways...)
x_start  = np.zeros((Nt,M,3)) # I don't think I actually need these? The information is already in x_points... Can extract these later... (Keep them for now and remove later)
x_end    = np.zeros((Nt,M,3))
x_points = np.zeros((Nt,M,N,3)) # Huuuuuge array...
array_accumul = np.zeros(3)
for i in range(Nt):    # Loop over time steps
    for k in range(M): # Loop over chains
        xes_thistime = xes[i] # [i,k] # I think I need to extract this twice, because it is a list of arrays (a double one at that...)
        ys_thistime  = ys[i]
        zs_thistime  = zs[i]
        #print("len(xes_this):",len(xes_thistime))
        #print(xes_thistime)
        #print("shape(xes_this):",np.size(xes_this))
        xes_this       = xes_thistime[k,:]
        ys_this        = ys_thistime[k,:]
        zs_this        = zs_thistime[k,:]
        #print("len(xes_this):",len(xes_this))
        x1 = xes_this[0]
        y1 = ys_this[0]
        z1 = zs_this[0]
        x_start[i,k,:] = np.array([x1,y1,z1])
        array_accumul  = x_start[i,k,:]
        x_points[i,k,0,:] = x_start[i,k,:]
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
            ree_dx           += dx
            ree_dy           += dy
            ree_dz           += dz
            array_accumul    += np.array([dx,dy,dz])
            x_points[i,k,j+1,:] = array_accumul
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
        rees[i,k,:]  = np.array([ree_dx,ree_dy,ree_dz])
        ree_len2      = ree_dx**2+ree_dy**2+ree_dz**2
        ree_vec[i,k]  = np.sqrt(ree_len2)
        ree2_vec[i,k] = ree_len2
        endzs[i,k]    = ree_dz
        x_end[i,k,:]  = rees[i,k,:]+x_start[i,k,:] # The coordinates of the end are the end-to-end vector plus starting point.

print('Done with dealing with BCs')
end_time = time.process_time()
print("Time:", end_time-start_time, " s")

### Now I have all the points after taking the BCs into account. Now I need to find the distances from each point to the end-to-end line
### I do this first and then divide into blocks h.

# Do it for each bond and each time step
distances = np.zeros((Nt,M,N))# Distance from line to point
for i in range(Nt):
    for k in range(M):
        # Extract start- and endpoints of line here 
        x_start = x_points[i,k,0,:]
        x_end   = x_points[i,k,N-1,:]
        for j in range(1,N-1): # We know that the distance is zero for the start- and the end bead.
            x_point = x_points[i,k,j,:]
            distances[i,k,j] =  np.linalg.norm(np.cross((x_point-x_start),(x_point-x_end)))/np.linalg.norm(x_end-x_start)
            print('----------------------------------------')
            print('distances[i,k,j]:',distances[i,k,j])
            print('x_start:',x_start)
            print('x_end:',x_end)
            print('x_end-x_start:',x_end-x_start)
            print('np.linalg.norm(x_end-x_start):',np.linalg.norm(x_end-x_start)) # This should not be nan, because we get a non-nan value
            print('x_point-x_start:',x_point-x_start)
            print('x_point-x_end:',x_point-x_end)
            print('np.cross((x_point-x_start),(x_point-x_end)):',np.cross((x_point-x_start),(x_point-x_end)))
            print('np.linalg.norm(np.cross((x_point-x_start),(x_point-x_end))):',np.linalg.norm(np.cross((x_point-x_start),(x_point-x_end))))
            if math.isnan(distances[i,k,j]):
                print('distances[i,k,j]:',distances[i,k,j])
                print('x_start:',x_start)
                print('x_end:',x_end)
                print('x_end-x_start:',x_end-x_start)
                print('np.linalg.norm(x_end-x_start):',np.linalg.norm(x_end-x_start)) # This should not be nan, because we get a non-nan value
                print('x_point-x_start:',x_point-x_start)
                print('x_point-x_end:',x_point-x_end)
                print('np.cross((x_point-x_start),(x_point-x_end)):',np.cross((x_point-x_start),(x_point-x_end)))
                print('np.linalg.norm(np.cross((x_point-x_start),(x_point-x_end))):',np.linalg.norm(np.cross((x_point-x_start),(x_point-x_end))))
                
print('----------------------------------------')
print('Done with finding distances')
end_time = time.process_time()
print("Time:", end_time-start_time, " s")

## partition into h:
boxes  = np.array([2,4,5,10,20,25,50]) # Number of boxes
Nboxes = len(boxes)                    # The number of different ways we divide the system into boxes
hs     = (N-1)/boxes                   # The number of atoms in one box
Nh     = len(hs)
box_average      = np.zeros((Nt, M, Nh,int(max(hs)))) # This array will have redundant entries. Need to remedy that somehow
box_sq_average   = np.zeros((Nt, M, Nh,int(max(hs)))) # This array will have redundant entries. Need to remedy that somehow
box_abs_average  = np.zeros((Nt, M, Nh,int(max(hs)))) # This array will have redundant entries. Need to remedy that somehow
box_rms          = np.zeros((Nt, M, Nh,int(max(hs)))) # This array will have redundant entries. Need to remedy that somehow
avgs_vs_time     = np.zeros((Nt,Nh))        # For binning values
avgs_sq_vs_time  = np.zeros((Nt,Nh))        # For binning values
rmses_vs_time    = np.zeros((Nt,Nh))        # For binning values
counter_t        = np.zeros((Nt,Nh))
absavgs_rms      = np.zeros(Nh)
# Should I take the average over time straight away?
for i in range(Nt):                   # Looping over time
    for k in range(M):                # Looping over chains
        for l in range(Nh):           # Looping over each choice of box size
            h = int(hs[l])
            for m in range(boxes[l]): # Looping over the boxes
                rms = 0
                #print('type k:', type(k))
                #print('type h:', type(h))
                #print('h:',h)
                #print('type m:', type(m))
                avg    = np.mean(distances[i,k,m*h:(m+1)*h])
                avg_sq = np.mean(distances[i,k,m*h:(m+1)*h]**2)
                absavg = np.mean(np.absolute(distances[i,k,m*h:(m+1)*h])) # Do I need this
                avgs_vs_time[i,l]     += np.sum(distances[i,k,m*h:(m+1)*h])
                avgs_sq_vs_time[i,l]  += np.sum(distances[i,k,m*h:(m+1)*h]**2)
                for n in range(h):    # Looping over the beads in each box
                    rms += (distances[i,k,m*h+n]-avg)**2
                    #if l==(Nh-1) and m==0:
                    #    print('n:',n, ', distance:', distances[i,k,m*h+n])
                rms = np.sqrt(rms/(h-1))         # rms for each box
                #if i==1 and l==0: # blocks in each
                #    print('rms:', rms)
                box_average[i,k,l,m]     = avg
                box_sq_average[i,k,l,m]  = avg_sq
                box_abs_average[i,k,l,m] = absavg
                box_rms[i,k,l,m]         = rms
                #if l==(Nh-1) and m==0:
                #    print('h =',h,': l==Nh-1, absavg:',absavg)
                #    print('abs(d) in m=0:',np.absolute(distances[i,k,m*h:(m+1)*h])) # Do I need this
                #print('distances:',distances[i,k,m*h:(m+1)*h])
                #print('avg distance:', avg)

for i in range(Nt):
    for l in range(Nh):
        # hs[l]: Number of distances in block, boxes[l]: Number of blocks. M: Number of chains
        # hs[l]*boxes[l]: Number of distances in one chain
        # hs[l]*boxes[l]*M: Number of distances in the system at one time step
        avgs_vs_time[i,l]    /= (hs[l]*boxes[l]*M)  
        avgs_sq_vs_time[i,l] /= (hs[l]*boxes[l]*M)
        rmses_vs_time[i,l]    = np.sqrt(avgs_sq_vs_time[i,l]-avgs_vs_time[i,l])

print('Done with first avgs and rms')
end_time = time.process_time()
print("Time:", end_time-start_time, " s")

# Taking the time average straight away.
avgs        = np.zeros((Nh,int(max(hs))))
avgs_sq     = np.zeros((Nh,int(max(hs))))
absavgs     = np.zeros((Nh,int(max(hs))))
avgs_tot    = np.zeros(Nh)
absavgs_tot = np.zeros(Nh)
avgs_sq_tot = np.zeros(Nh)
rmses       = np.zeros(Nh)
for l in range(Nh):        # Looping over each choice of box size
    h = int(hs[l])
    for i in range(Nt):    # Looping over time
        for k in range(M): # Looping over chains
            for m in range(boxes[l]): # Looping over each box
                avgs[l,m]      += np.mean(distances[i,k,m*h:(m+1)*h])
                avgs_sq[l,m]   += np.mean(distances[i,k,m*h:(m+1)*h]**2)
                absavgs[l,m]   += np.mean(np.absolute(distances[i,k,m*h:(m+1)*h]))
                avgs_tot[l]    += np.mean(distances[i,k,m*h:(m+1)*h])
                avgs_sq_tot[l] += np.mean(distances[i,k,m*h:(m+1)*h]**2)
                absavgs_tot[l] += np.mean(np.absolute(distances[i,k,m*h:(m+1)*h]))

avgs        /= Nt*M
avgs_sq     /= Nt*M
absavgs     /= Nt*M
avgs_tot    /= Nt*M
avgs_sq_tot /= Nt*M
absavgs_tot /= Nt*M

'''
for i in range(Nboxes):
    avgs_tot    /= boxes[i]
    avgs_sq_tot /= boxes[i]
    absavgs_tot /= boxes[i]
'''

rms_main_from_varianceformulation     = np.sqrt(avgs_sq_tot-(avgs_tot**2))
rms_temp_avs_from_varianceformulation = np.zeros(Nh)
rms_rms_from_varianceformulation      = np.zeros(Nh)
for i in range(Nt):
    for l in range(Nh):
        rms_temp_avs_from_varianceformulation[l] += np.sqrt(avgs_sq_vs_time[i,l]-avgs_vs_time[i,l]**2)
rms_temp_avs_from_varianceformulation /= Nt


for l in range(Nh):
    for i in range(Nt):
        rms_rms_from_varianceformulation[l] += (rmses_vs_time[i,l]-rms_temp_avs_from_varianceformulation[l])**2
    rms_rms_from_varianceformulation[l] = np.sqrt(rms_rms_from_varianceformulation[l]/(Nt*(Nt-1)))

print('Done with total averages')

end_time = time.process_time()
print("Time:", end_time-start_time, " s")

for l in range(Nh):                    # Looping over each choice of box size
    h = int(hs[l])
    for i in range(Nt):                # Looping over time
        for k in range(M):             # Looping over chains
            for m in range(boxes[l]):  # Looping over each box
                for n in range(h):     # Looping over the beads in each box
                    rmses[l] += (distances[i,k,m*h+n]-avgs_tot[l])**2 # I have used the total average # 'Cause when you calculate the rms, you should subtract with the same average for all values
    rmses[l] = np.sqrt(rmses[l]/(h*M*Nt-1))

print('Done with total rms')
end_time = time.process_time()
print("Time:", end_time-start_time, " s")

rms_av  = np.zeros(Nh)
rms_rms = np.zeros(Nh)
absavgs_correct      = np.zeros(Nh)
absavgs_correct_rms  = np.zeros(Nh)
absavgs_oneblock     = np.zeros(Nh)
absavgs_oneblock_rms = np.zeros(Nh)


for l in range(Nh):
    rms_av[l]          = np.mean(box_rms[:,:,l,0:boxes[l]])                 # Printing tests below verified that this is the right way to treat the arrays
    absavgs_correct[l] = np.mean(box_abs_average[:,:,l,0:boxes[l]])
    absavgs_oneblock[l] = np.mean(box_abs_average[:,:,l,0])                 # Just taking one block to actually get information about the scaling with h
    #if l==(Nh-1):
    #    print('box_abs_average[:,:,Nh-1,0]:',box_abs_average[:,:,l,0])
    #print('box_abs_average[2,2,l,0:boxes[l]-1]:',box_abs_average[2,2,l,0:boxes[l]-1])
    #print('box_abs_average[2,2,l,0:boxes[l]]:',box_abs_average[2,2,l,0:boxes[l]])           # Get exactly what we need
    #print('box_abs_average[2,2,l,0:boxes[l]+1]:',box_abs_average[2,2,l,0:boxes[l]+1])
    #print('box_abs_average[2,2,l,:]:',box_abs_average[2,2,l,:])
    #print('box_abs_average[0,0,l,:]:',box_abs_average[0,0,l,:])
    

for l in range(Nh):
    h = int(hs[l])
    for i in range(Nt):
        for k in range(M):
            absavgs_oneblock_rms[l] += (box_abs_average[i,k,l,0]-absavgs_oneblock[l])**2
            for m in range(boxes[l]):
                rms_rms[l] += (box_rms[i,k,l,m]-rms_av[l])**2
                absavgs_correct_rms[l] += (box_abs_average[i,k,l,m]-absavgs_correct[l])**2
    rms_rms[l]              = np.sqrt(rms_rms[l]/(h*M*Nt-1))
    absavgs_correct_rms[l]  = np.sqrt(absavgs_correct_rms[l]/(h*M*Nt-1))
    absavgs_oneblock_rms[l] = np.sqrt(absavgs_oneblock_rms[l]/(M*Nt-1))

print('Done with rms of rms')
end_time = time.process_time()
print("Time:", end_time-start_time, " s")

print('rms_av[0]:',rms_av[0])
print('rms_av:', rms_av)

sqrths = hs**(1/2.)

coeffs = np.polyfit(sqrths,rms_av,1)
a      = coeffs[0]
b      = coeffs[1]
p      = np.poly1d(coeffs)
fitarray = p(sqrths)

# Print coeffs and data points to file


'''# Don't really need this
plt.figure(figsize=(6,5))
plt.errorbar(sqrths,rms_av, yerr=rms_rms,  fmt="none", capsize=2, label='Results')
plt.plot(sqrths,rms_av, 'o')
plt.plot(sqrths,fitarray, label='Fit')
plt.xlabel(r'$h^{1/2}$', fontsize=16)
plt.ylabel(r'$\Delta\omega$', fontsize=16)
plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
plt.legend(loc="lower right")
plt.title(r'$\Delta\omega$ vs $h^{1/2}$, from average ds', fontsize=16)
plt.savefig(plot1name)
'''

'''
plt.figure(figsize=(6,5))
plt.errorbar(hs,rms_av, yerr=rms_rms,  fmt="none", capsize=2)#, label='results')
plt.plot(sqrths,rms_av, 'o')
plt.xlabel(r'Height $h$', fontsize=16)
plt.ylabel(r'$\Delta\omega$', fontsize=16)
plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
#plt.legend(loc="upper right")
plt.title(r'$\Delta\omega$ vs $h$', fontsize=16)
plt.savefig(plot1name)
'''

plt.figure(figsize=(6,5))
plt.errorbar(sqrths,rms_av, yerr=rms_rms,  fmt="none", capsize=2, label='Results')
plt.plot(sqrths,rms_av, 'o')
plt.plot(sqrths,fitarray, label='Fit')
plt.xlabel(r'$h^{1/2}$', fontsize=16)
plt.ylabel(r'$\Delta\omega$', fontsize=16)
#plt.tight_layout(pad=3.0,w_pad=0.0, h_pad=0.5)
plt.tight_layout(pad=3.0)#, w_pad=0.0, h_pad=0.5)
plt.legend(loc="lower right")
plt.title(r'$\Delta\omega$ vs $h^{1/2}$, from average ds', fontsize=16)
plt.savefig(plot2name)


'''
plt.figure(figsize=(6,5))
plt.errorbar(sqrths,rms_main_from_varianceformulation, yerr=rms_rms_from_varianceformulation,  fmt="none", capsize=2)#, label='results')
plt.xlabel(r'$h^{1/2}$', fontsize=16)
plt.ylabel(r'$\Delta\omega$', fontsize=16)
plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
#plt.legend(loc="upper right")
plt.title(r'$\Delta\omega$ vs $h^{1/2}$, from average ds', fontsize=16)
plt.savefig('fromvariance')


plt.figure(figsize=(6,5))
plt.errorbar(hs,rmses, yerr=rms_rms,  fmt="none", capsize=2)#, label='results')
plt.xlabel(r'Height $h$', fontsize=16)
plt.ylabel(r'$\Delta\omega$', fontsize=16)
plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
#plt.legend(loc="upper right")
plt.title(r'$\Delta\omega$ vs $h$', fontsize=16)
plt.savefig(plot1name)


plt.figure(figsize=(6,5))
plt.errorbar(sqrths,rmses, yerr=rms_rms,  fmt="none", capsize=2)#, label='results')
plt.xlabel(r'$h^{1/2}$', fontsize=16)
plt.ylabel(r'$\Delta\omega$', fontsize=16)
plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
#plt.legend(loc="upper right")
plt.title(r'$\Delta\omega$ vs $h$', fontsize=16)
plt.savefig(plot2name)
'''
print('Plotting done, writing to file')


outfile.write('Fit to linear function delta omega=a**x+b, x=h^(1/2.). Order: a, b:')
outfile.write('%.16f %.16f\n' % (a,b))
outfile.write('Order: h, absavgs_oneblock, absavgs_oneblock_rms, rms_av, rms_rms. Averages over time, chains and boxes\n')
#outfile.write('Order: h, absavgs_correct, absavgs_correct_rms, rms_av, rms_rms. Averages over time, chains and boxes\n')
for l in range(Nh):
    outfile.write('%i %.16f %.16f %.16f %.16f\n' % (hs[l], absavgs_oneblock[l], absavgs_oneblock_rms[l], rms_av[l], rms_rms[l]))
    #outfile.write('%i %.16f %.16f %.16f %.16f\n' % (hs[l], absavgs_correct[l], absavgs_correct_rms[l], rms_av[l], rms_rms[l]))

outfile.close()

end_time = time.process_time()
print("Time:", end_time-start_time, " s")
