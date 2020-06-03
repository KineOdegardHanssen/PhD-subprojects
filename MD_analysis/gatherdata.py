import matplotlib.pyplot as plt               # To plot
from scipy.optimize import curve_fit
import numpy as np
import random
import math
import time

start_time = time.process_time()

M           = 9
N           = 101
kangle      = 20
kbond       = 200#0
Kangle      = kangle
Kbond       = kbond
factors     = [0.1,1,10,100,250]
charges     = [0,-1,-5,-10] 
#spacings   = [1,2,3,4,5,10,40,100]
spacing     = 1
gridspacing = spacing
spacings    = [1,2,3,4,5,6,7,8,10,15,40,100,300]#[100,200,300]#[1,2,3,4,5,6,7,8,10,15,40,100] #[100,200,300]#[1,3,5,10,40]
spacings    = [3,5,7,10,40,100]
lengths     = [21,61,101,141]
lensp       = [1,3,5,7]
krfacs      = [0.01, 0.05, 0.10, 0.50, 1.00]
#kangles     = [20, 100, 200, 1000, 2000] # krfxfacsims
#kangles     = [2000,5000, 8000, 10000, 50000, 125000, 150000] # kthetasims
kangles     = [20,200,2000]
diels       = [1, 2, 10, 50, 100, 1000] # For spacing 1 nm
#diels       = [1, 2, 10, 100] # For spacing 2 nm
charge      = -1
mass        = 1
T           = 310
wall        = 1.042
epsilon     = 1.042
ljcutoff    = 1.12246204830937
sigma       = 1


chargesims  = False
factorsims  = False
spacesims   = True
lengthsims  = False # Varying the length and the spacing
lensfxsims  = False # Varying the length, fixing the spacing
krfxfacsims = False
wallsims    = False
dielsims    = False
kthetasims  = False
comparefene = False

# Variation in file names:
mass_specified = False
withwall       = False
ljsims_b_a     = False # Terms: Lennard-Jones, bond, angle
spacingfenestd = False # FENE, standard parameters
spacingfenemy  = False # FENE, my parameters
ljunits        = False
debyecorrected = False
ljdebfeneanglenowall_ljunits = False
ljdebfeneanglenowall_ljunits_corrected = True


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
if wallsims == True:
    wallenergy=1.042
    Narray = 1
if dielsims == True:
    Narray = len(diels)
if kthetasims == True:
    kbond       = 200
    gridspacing = 4
    spacing     = gridspacing
    Narray      = len(kangles)
    kangles     = np.array(kangles)
    krfacs      = kangles/float(kbond)
if comparefene == True:
    M         = 1
    filenames = ['_chaingrid_quadratic_M1N101_gridspacing300_Langevin_Kbond200_wall1.042_T310_theta0is180_twofirst_are_fixed.txt', '_onechain_M1N101_gridspacing300_Langevin_wall1.042_KfeneR0sq200_1.2_T310_theta0is180_twofirst_are_fixed.txt', '_onechain_M1N101_gridspacing300_Langevin_wall1.042_KfeneR0sq200_1.2_bondepsilon1.042_fenesigma0.8_T310_theta0is180_twofirst_are_fixed.txt', '_onechain_M1N101_gridspacing300_Langevin_wall1.042_KfeneR0sq67.5_1.5_bondepsilon1.042_fenesigma1_T310_theta0is180_twofirst_are_fixed.txt']
    types     = ['Harmonic','FENE', 'FENE+LJ', 'FENE std']
    Narray    = len(filenames)
if spacingfenestd==True:
    M         = 1
    KfeneR0sq = 67.5
    R0        = 1.5
    spacing   = 300
    fenesigma = 1
    feneepsilon = epsilon
if spacingfenemy==True:
    M         = 1
    KfeneR0sq = 200
    R0        = 1.2
    spacing   = 300
    fenesigma = 0.8
    feneepsilon = epsilon
if ljunits==True:
    kangle = 14.0186574854529
    kbond  = 140.186574854529
    T      = 3
if debyecorrected==True:
    effectivedielectric = 0.00881819074717447
if ljdebfeneanglenowall_ljunits==True:
    #(M,N,spacing,kangle,kfene,R0,epsilon,sigma,kappa,debyecutoff,charge,debconv,T)
    kangle  = 14.0186574854529
    sigma   = 1
    epsilon = 1.042
    T       = 3
    kappa   = 1
    debconv = 0.00881819074717447
    R0      = 1.5
    kfene   = 67.5
    debyecutoff = 3
    spacings    = [3,4,5,7,10,15,100,300]
    Narray      = len(spacings)
if ljdebfeneanglenowall_ljunits_corrected==True:
    print('Im in')
    kangle  = 5
    sigma   = 1
    epsilon = 1.042
    T       = 3
    kappa   = 1
    debconv = 0.00881819074717447
    R0      = 1.5
    kfene   = 67.5
    debyecutoff = 3
    spacings    = [3,5,7,10,300]#[3,4,5,7,10,15,100,300]
    Narray      = len(spacings)
    effectivedielectric = debconv

bondlengths = np.zeros(Narray)
bondrms     = np.zeros(Narray)


if chargesims==True:
    filenamebase = 'quadratic_gridspacing%i_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varycharge' % (spacing,kangle,kbond,T)
    outfilename  = 'table_systemwide_'+filenamebase+'.txt'
    outfile      = open(outfilename, 'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c}\nCharge & $tilde{l}_p$ & $l_p$ & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$ & $l_b$ & $z_{N}$ \ \ \n\hline\n')
    #outfilename2 = 'table_endzs_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,kangle,kbond) +str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varycharge.txt' % T
    outfilename2 = 'table_endzs_'+filenamebase+'txt'
    outfile2     = open(outfilename2, 'w')
    outfile2.write('resizebox{textwidth}{!}{begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')


if factorsims==True:
    filenamebase = 'quadratic_M%iN%i_Langevin_Kangle%i_T%i_angle_stretching_only_varyfactor' % (M,N,kangle,T)
    outfilename  = 'table_systemwide_'+filenamebase+'.txt'
    outfile      = open(outfilename, 'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c}\nFactor & $tilde{l}_p$ & $l_p$ & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$ & $l_b$ & $z_{N}$ \ \ \n\hline\n')
    #outfilename2 = 'table_endzs_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,kangle,kbond) +str(factor)+'_T%i_theta0is180_twofirst_are_fixed_varyspacing_N%i.txt' % T
    outfilename2 = 'table_endzs_'+filenamebase+'txt'
    outfile2     = open(outfilename2, 'w')
    outfile2.write('resizebox{textwidth}{!}{begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')


if spacesims==True:
    filenamebase = 'quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_charge%i_varyspacing' % (M,N,kangle,kbond,T,charge)
    if mass_specified==True:
        filenamebase = 'quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_charge%i_mass%i_varyspacing' % (M,N,kangle,kbond,T,charge,mass)
    if withwall==True:
        filenamebase = 'quadratic_M%iN%i_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_charge%i_varyspacing' % (M,N,wall,kangle,kbond,T,charge)
    if ljsims_b_a==True:
        filenamebase = 'quadratic_M%iN%i_lj_epsilon%.3f_sigma%i_ljcutoff%.14f_kappa1_Langevin_Kangle%i_Kbond%i_T%i_theta0is180_twofirst_are_fixed_varyspacing' % (M,N,epsilon,sigma,ljcutoff,Kangle,Kbond,T)
    if ljsims_b_a==True and withwall==True:
        filenamebase = 'quadratic_M%iN%i_ljdebye_epsilon%.3f_sigma1_ljcutoff%.14f_kappa1_debyecutoff3_charge%i_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180_twofirst_are_fixed_varyspacing' % (M,N,epsilon,ljcutoff,charge,wall,kangle,kbond,T)
    if ljunits==True: ### Add more bools here as I go along!
        filenamebase = 'quadratic_M%iN%i_ljunits_Langevin_Kangle%.13f_Kbond%.12f_debye_kappa1_debyecutoff3_chargeelementary%i_T%i_theta0is180_correctconversions_twofirst_are_fixed_varyspacing.txt' % (M,N, kangle,kbond,charge,T)
    if ljunits==True and debyecorrected==True:
        filenamebase = 'quadratic_M%iN%i_ljunits_Langevin_Kangle%.13f_Kbond%.12f_debye_kappa1_debyecutoff3_chargeelementary%i_effectivedielectric%.17f_T%i_theta0is180_varyspacing.txt' % (M,N,Kangle,Kbond,charge,effectivedielectric,T)
    if ljdebfeneanglenowall_ljunits==True: # CHECK!!!!!!!!!!!!
        filenamebase = 'quadratic_M%iN%i_ljunits_Langevin_wall_Kangle%.13f_KfeneR0sq%.1f_R0%.1f_feneepsilon%.3f_fenesigma%i_debye_kappa%i_debyecutoff%i_chargeelementary%i_debconv%.17f_T%i_varyspacing.txt' % (M,N,kangle,kfene,R0,epsilon,sigma,kappa,debyecutoff,charge,debconv,T)
    if ljdebfeneanglenowall_ljunits_corrected==True:
        print('So in')
        filenamebase = 'quadr_M%iN%i_ljun_Langevin_wall_Kangle%i_KfR0sq%.1f_R0%.1f_feps%.3f_fsgm%i_ljeps%.3f_ljsgm%i_debye_kappa%i_debcut%i_elch%i_effdiel%.17f_T%i_theta0is180_varyspacing.txt' % (M,N,kangle,kfene,R0,epsilon,sigma,epsilon,sigma, kappa, debyecutoff,charge,effectivedielectric,T)
    ## Tests done
    outfilename  = 'table_systemwide_'+filenamebase+'.txt'
    outfile      = open(outfilename, 'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c}\nSpacing & $tilde{l}_p$ & $l_p$ & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$ & $l_b$ & $z_{N}$ \ \ \n\hline\n')
    #outfilename2 = 'table_endzs_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varyspacing.txt' % T
    outfilename2 = 'table_endzs_'+filenamebase+'txt'
    outfile2     = open(outfilename2, 'w')
    outfile2.write('resizebox{textwidth}{!}{begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')


if lengthsims==True:
    filenamebase = 'quadratic_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varylength_contourdsep20' % (kangle,kbond,T)
    outfilename  = 'table_systemwide_'+filenamebase+'.txt'
    outfile      = open(outfilename, 'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c}\nN & $tilde{l}_p$ & $l_p$ & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$ & $l_b$ & $z_{N}$ \ \ \n\hline\n')
    #outfilename2 = 'table_endzs_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varylength_contourdsep20.txt' % T
    outfilename2 = 'table_endzs_'+filenamebase+'txt'
    outfile2     = open(outfilename2, 'w')
    outfile2.write('resizebox{textwidth}{!}{begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')


if lensfxsims==True:
    filenamebase = 'quadratic_gridspacing%i_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varylength' % (spacing,kangle,kbond,T)
    outfilename  = 'table_systemwide_'+filenamebase+'.txt'
    outfile      = open(outfilename, 'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c}\nN & $tilde{l}_p$ & $l_p$ & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$ & $l_b$ & $z_{N}$ \ \ \n\hline\n')
    #outfilename2 = 'table_endzs_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varylength.txt' % T
    outfilename2 = 'table_endzs_'+filenamebase+'txt'
    outfile2     = open(outfilename2, 'w')
    outfile2.write('resizebox{textwidth}{!}{begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')


if krfxfacsims==True or kthetasims==True:
    filenamebase = 'quadratic_spacing%i_Kbond%i_T%i_angle_stretching_only_varyfactor_Krfixed' % (spacing,kbond,T)
    if spacingfenestd==True:
        filenamebase = 'onechain_M%iN%i_spacing%i_Langevin_nowall_KfeneR0sq%.1f_R0%.1f_feneepsilon%.3f_fenesigma%i_T%i_theta0is180_twofirst_are_fixed_varyKangle' % (M,N,spacing,KfeneR0sq,R0,feneepsilon,fenesigma,T)
    if spacingfenemy==True:
        filenamebase = 'onechain_M%iN%i_spacing%i_Langevin_nowall_KfeneR0sq%i_R0%.1f_feneepsilon%.3f_fenesigma%.1f_T%i_theta0is180_twofirst_are_fixed_varyKangle' % (M,N,spacing,KfeneR0sq,R0,feneepsilon,fenesigma,T)
    outfilename  = 'table_systemwide_'+filenamebase+'.txt'
    outfilename  = 'table_systemwide_'+filenamebase+'.txt'
    outfile      = open(outfilename, 'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|r|c|c|c|c|c|c}\n K_theta/K_r & $K_theta$ & $tilde{l}_p$ & $l_p$ & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$ & $l_b$ & $z_{N}$ \ \ \n\hline\n')
    #outfilename2 = 'table_endzs_chaingrid_quadratic_M%iN%i_spacing%i_Langevin_Kbond%i_T%i_theta0is180_twofirst_are_fixed_varyfactor_Krfixed.txt' % (M,N,spacing,kbond,T)
    outfilename2 = 'table_endzs_'+filenamebase+'txt'
    outfile2     = open(outfilename2, 'w')
    outfile2.write('resizebox{textwidth}{!}{begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')

if wallsims==True:
    filenamebase  = 'quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed_wall' % (M,N,spacing,wallenergy,kangle,kbond,charge,T)
    outfilename  = 'table_systemwide_'+filenamebase+'.txt'
    outfile       = open(outfilename, 'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|r|c|c|c|c|c}\n $E_{wall}$ &  $tilde{l}_p$ & $l_p$ & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$ & $l_b$ & $z_{N}$ \ \ \n\hline\n')
    #outfilename2 = 'table_endzs_chaingrid_quadratic_M%iN%i_spacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed_wall.txt' % (M,N,spacing,wallenergy,kangle,kbond,charge,T)
    outfilename2 = 'table_endzs_'+filenamebase+'txt'
    outfile2     = open(outfilename2, 'w')
    outfile2.write('resizebox{textwidth}{!}{begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')

if dielsims==True:
    filenamebase = 'chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,gridspacing,Kangle,Kbond,charge,T)
    outfilename  = 'table_systemwide_'+filenamebase+'.txt'
    outfile      = open(outfilename,'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|r|c|c|c|c|c}\n $\epsilon_{diel}$ &  $tilde{l}_p$ & $l_p$ & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$ & $l_b$ & $z_{N}$ \ \ \n\hline\n')
    #outfilename2 = 'table_endz_chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,T)
    outfilename2 = 'table_endzs_'+filenamebase+'txt'
    outfile2     = open(outfilename2,'w')
    outfile2.write('resizebox{textwidth}{!}{begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')

if comparefene==True:
    filenamebase = 'compare_fene_and_harmonic'
    outfilename  = 'table_systemwide_'+filenamebase+'.txt'
    outfile      = open(outfilename, 'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|r|c|c|c|c|c|c}\n Type & $tilde{l}_p$ & $l_p$ & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$ & $l_b$ & $z_{N}$ \ \ \n\hline\n')
    #outfilename2 = 'table_endzs_chaingrid_quadratic_M%iN%i_spacing%i_Langevin_Kbond%i_T%i_theta0is180_twofirst_are_fixed_varyfactor_Krfixed.txt' % (M,N,spacing,kbond,T)
    outfilename2 = 'table_endzs_'+filenamebase+'.txt'
    outfile2     = open(outfilename2, 'w')
    outfile2.write('resizebox{textwidth}{!}{begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')


# Making lists to read in data
pl_sep_avg_array   = np.zeros(Narray)
pl_sep_rms_array   = np.zeros(Narray)
pl_dist_avg_array  = np.zeros(Narray)
pl_dist_rms_array  = np.zeros(Narray)
ree2_avg_array     = np.zeros(Narray)
ree2_rms_array     = np.zeros(Narray)
costheta_avg_array = np.zeros(Narray)  # costhetas are for neighbour bonds
costheta_rms_array = np.zeros(Narray)
endz_avg_array     = np.zeros(Narray)
endz_rms_array     = np.zeros(Narray)

k = 0
#for kangle in kangles:
for i in range(Narray):
    print((k+1)/float(Narray))
    infilename_base = 'error'
    if chargesims==True:
        charge      = charges[i]
        infilename_base  = '_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,kangle,kbond) +str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varydielectric.txt' % T
        
    if factorsims==True:
        factor = factors[i]
        kbond  = kangle*factors[i]
        infilename_base = '_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,kangle,kbond) +str(factor)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        
    if spacesims==True:
        spacing     = spacings[i]
        infilename_base = '_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        if mass_specified==True:
            infilename_base = '_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_mass%i_T%i_theta0is180_twofirst_are_fixed.txt' % (mass,T)  
        if withwall==True:
            infilename_base = '_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,wall,kangle,kbond,charge,T)
        if ljsims_b_a==True:
            infilename_base = '_chaingrid_quadratic_M%iN%i_gridspacing%i_lj_epsilon%.3f_sigma%i_ljcutoff%.14f_kappa1_Langevin_Kangle%i_Kbond%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,epsilon,sigma,ljcutoff,Kangle,Kbond,T)
        if withwall==True and ljsims_b_a==True:
            infilename_base = '_chaingrid_quadratic_M%iN%i_gridspacing%i_ljdebye_epsilon%.3f_sigma1_ljcutoff%.14f_kappa1_debyecutoff3_charge%i_Langevin_wall%.3f_Kangle%i_Kbond%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,epsilon,ljcutoff,charge,wall,kangle,kbond,T)
        if ljunits==True:
            infilename_base = '_chaingrid_quadratic_M%iN%i_ljunits_gridspacing%i_Langevin_Kangle%.13f_Kbond%.12f_debye_kappa1_debyecutoff3_chargeelementary%i_T%i_theta0is180_correctconversions_twofirst_are_fixed.txt' % (M,N, spacing,kangle,kbond,charge,T)
        if ljunits==True and debyecorrected==True:
            infilename_base = '_chaingrid_quadratic_M%iN%i_ljunits_gridspacing%i_Langevin_Kangle%.13f_Kbond%.12f_debye_kappa1_debyecutoff3_chargeelementary%i_effectivedielectric%.17f_T%i_theta0is180.txt' % (M,N,spacing,kangle,kbond,charge,effectivedielectric,T)
        if ljdebfeneanglenowall_ljunits==True:
            infilename_base = '_chaingrid_quadratic_M%iN%i_ljunits_gridspacing%i_Langevin_wall_Kangle%.13f_KfeneR0sq%.1f_R0%.1f_feneepsilon%.3f_fenesigma%i_debye_kappa%i_debyecutoff%i_chargeelementary%i_debconv%.17f_T%i.txt'  % (M,N,spacing,kangle,kfene,R0,epsilon,sigma,kappa,debyecutoff,charge,debconv,T)
        if ljdebfeneanglenowall_ljunits_corrected==True:
            print('But I just cant win')
            infilename_base = '_chgr_quadr_M%iN%i_ljun_gridsp%i_Langevin_wall_Kangle%i_KfR0sq%.1f_R0%.1f_feps%.3f_fsgm%i_ljeps%.3f_ljsgm%i_debye_kappa%i_debcut%i_elch%i_effdiel%.17f_T%i_theta0is180.txt' % (M,N,spacing,kangle,kfene,R0,epsilon,sigma,epsilon,sigma, kappa, debyecutoff,charge,effectivedielectric,T)
    if lengthsims==True:
        print('Is this were I went wrong')
        N            = lengths[i]
        spacing      = lensp[i]
        infilename_base   = '_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        
    
    if lensfxsims==True:
        print('And made the filename that dont belong')
        N            = lengths[i]
        infilename_base   = '_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
    
    if krfxfacsims==True:
        factor       = krfacs[i]
        kangle       = kangles[i]
        infilename_base   = '_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
        
    if wallsims==True:
        infilename_base   = '_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,wallenergy,kangle,kbond,charge,T)
        
    if dielsims==True:
        dielectric  = diels[i]
        infilename_base  = '_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
    
    if kthetasims==True:
        factor       = krfacs[i]
        kangle       = kangles[i]
        infilename_base   = '_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.1f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
        if spacingfenestd==True:
            infilename_base = '_onechain_M%iN%i_gridspacing%i_Langevin_nowall_Kangle%i_KfeneR0sq%.1f_R0%.1f_feneepsilon%.3f_fenesigma%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,KfeneR0sq,R0,feneepsilon,fenesigma,T)
            factor = 1.0
        if spacingfenemy==True:
            infilename_base = '_onechain_M%iN%i_gridspacing%i_Langevin_nowall_Kangle%i_KfeneR0sq%i_R0%.1f_feneepsilon%.3f_fenesigma%.1f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,KfeneR0sq,R0,feneepsilon,fenesigma,T)
            factor = 1.0
    
    if comparefene==True:
        infilename_base = filenames[i]
    
    ################################  
    
    infilename  = 'costhetas_neighbourbonds'+infilename_base
    infilename2 = 'ree_average'+infilename_base
    ###
    infilename3  = 'persistence_length'+infilename_base
    infilename4  = 'persistence_length_actualdistance'+infilename_base
    ###
    infilename5 = 'bondlengths_systemwide'+infilename_base
    infilename6 = 'endzs'+infilename_base
    
    print('infilename:', infilename)
    
    # First file: cos(theta_1)
    infile       = open(infilename,'r')
    lines        = infile.readlines()
    lastline     = lines[-1]
    words        = lastline.split()
    costheta_av  = float(words[2])
    costheta_rms = float(words[3])
    infile.close()
    
    # Second file: <ree> and <ree2>:
    infile2  = open(infilename2,'r')
    lines    = infile2.readlines()
    lastline = lines[-1]
    words    = lastline.split()
    ree2_av  = float(words[4])
    ree2_rms = float(words[5])
    infile2.close()
    
    # Third file: lp (separation):
    infile3    = open(infilename3,'r')
    line       = infile3.readline()
    words      = line.split()
    pl_sep_av  = float(words[8])
    pl_sep_rms = float(words[9])
    infile3.close()
    
    # Fourth file: lp (distance):
    infile4     = open(infilename4,'r')
    line        = infile4.readline()
    words       = line.split()
    pl_dist_av  = float(words[8])
    pl_dist_rms = float(words[9])
    infile4.close()
    
    # Fifth file: bond length
    #print('infilename5:',infilename5)
    infile5        = open(infilename5,'r')
    lines          = infile5.readlines()
    lastline       = lines[-1]
    words          = lastline.split()
    bondlength_av  = float(words[2])
    bondlength_rms = float(words[3])
    infile5.close()
    
    # Sixth file: End zs:
    infile6  = open(infilename6,'r')
    lines    = infile6.readlines()
    lastline = lines[-1]
    words    = lastline.split()
    endz_av  = float(words[3])
    endz_rms = float(words[4])
    
    
    # Read in for endzs:
    endzs_av_chain  = []
    endzs_rms_chain = [] 
    for line in lines:
        words = line.split()
        if len(words)>0 and words[0]!='Total':
            endzs_av_chain.append(float(words[0]))
            endzs_rms_chain.append(float(words[1]))
    infile6.close()
    
    # Appending some stuff to arrays:
    pl_sep_avg_array[i]   = pl_sep_av
    pl_sep_rms_array[i]   = pl_sep_rms
    pl_dist_avg_array[i]  = pl_dist_av
    pl_dist_rms_array[i]  = pl_sep_av
    ree2_avg_array[i]     = ree2_av
    ree2_rms_array[i]     = ree2_rms
    costheta_avg_array[i] = costheta_av
    costheta_rms_array[i] = costheta_rms
    endz_avg_array[i]     = endz_av
    endz_rms_array[i]     = endz_av
    
    '''
    outfile2     = open(outfilename2, 'w')
    outfile2.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')
    outfile2.write('$z_{end}$ ')
    for j in range(M):
        outfile2.write('& %.3f$\pm$%.3f ' % (endzs_av_chain[j], endzs_rms_chain[j]))
    outfile2.write('& %.3f$\pm$%.3f\n' % (endz_av, endz_rms))
    outfile2.write('\n\end{tabular}\n\label{table:quantities_vs_spacing}\n\end{table}') 
    outfile2.close()
    '''
    
    # Old printing:
    '''
    if chargesims==True:
        outfile.write('%i & %.5f$\pm$%.5f   &  %.5f$\pm$%.5f   & %.2f$\pm$%.2f & %.3f$\pm$%.3f & %.3f$\pm$%.3f & %.3f$\pm$%.3f \ \ \n' % (charge, pl_dist_av,pl_dist_rms, pl_sep_av, pl_sep_rms, ree2_av, ree2_rms, costheta_av,costheta_rms, bondlength_av, bondlength_rms, endz_av,endz_rms))
        outfile2.write('%i ' % charge)
       
    if factorsims==True or krfxfacsims==True or kthetasims==True:
        outfile.write('%.2f & %i & %.5f$\pm$%.5f   &  %.5f$\pm$%.5f   & %.2f$\pm$%.2f & %.4f$\pm$%.4f & %.3f$\pm$%.3f & %.3f$\pm$%.3f  \ \ \n' % (factor, kangle, pl_dist_av,pl_dist_rms, pl_sep_av, pl_sep_rms, ree2_av, ree2_rms, costheta_av,costheta_rms, bondlength_av, bondlength_rms, endz_av,endz_rms))
        outfile2.write('%.2f & %i ' % (factor,kangle))
    
    if spacesims==True:
        outfile.write('%i & %.5f$\pm$%.5f   &  %.5f$\pm$%.5f   & %.2f$\pm$%.2f & %.3f$\pm$%.3f & %.3f$\pm$%.3f & %.3f$\pm$%.3f  \ \ \n' % (spacing, pl_dist_av,pl_dist_rms, pl_sep_av, pl_sep_rms, ree2_av, ree2_rms, costheta_av,costheta_rms, bondlength_av, bondlength_rms, endz_av,endz_rms))
        outfile2.write('%i ' % spacing)
    
    if lengthsims==True or lensfxsims==True:
        outfile.write('%i & %.5f$\pm$%.5f   &  %.5f$\pm$%.5f   & %.2f$\pm$%.2f & %.3f$\pm$%.3f & %.3f$\pm$%.3f & %.3f$\pm$%.3f \ \ \n' % (N, pl_dist_av,pl_dist_rms, pl_sep_av, pl_sep_rms, ree2_av, ree2_rms, costheta_av,costheta_rms, bondlength_av, bondlength_rms, endz_av,endz_rms))
        outfile2.write('%i ' % N)
    
    if wallsims==True:
        outfile.write('%.3f & %.5f$\pm$%.5f   &  %.5f$\pm$%.5f   & %.2f$\pm$%.2f & %.3f$\pm$%.3f & %.3f$\pm$%.3f & %.3f$\pm$%.3f \ \ \n' % (wallenergy, pl_dist_av,pl_dist_rms, pl_sep_av, pl_sep_rms, ree2_av, ree2_rms, costheta_av,costheta_rms, bondlength_av, bondlength_rms, endz_av,endz_rms))
        outfile2.write('%.3f ' % wallenergy)
    
    if dielsims==True:
        outfile.write('%i & %.5f$\pm$%.5f   &  %.5f$\pm$%.5f   & %.2f$\pm$%.2f & %.3f$\pm$%.3f & %.3f$\pm$%.3f & %.3f$\pm$%.3f \ \ \n' % (dielectric, pl_dist_av,pl_dist_rms, pl_sep_av, pl_sep_rms, ree2_av, ree2_rms, costheta_av,costheta_rms, bondlength_av, bondlength_rms, endz_av,endz_rms))
        outfile2.write('%i ' % dielectric)
    
    if comparefene==True:
        # type    & $tilde{l}_p$             & $l_p$                   & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$     & $l_b$                          & $z_{N}$
        #types[i] & pl_dist_av+/-pl_dist_rms &  pl_sep_av+/-pl_sep_rms & ree2_av+/-ree2_rms & costheta_av+/-costheta_rms & bondlength_av+/-bondlength_rms & endz_av+/-endz_rms
        outfile.write('%s & %.5f$\pm$%.5f   &  %.5f$\pm$%.5f   & %.2f$\pm$%.2f & %.3f$\pm$%.3f & %.3f$\pm$%.3f & %.3f$\pm$%.3f \ \ \n' % (types[i], pl_dist_av,pl_dist_rms, pl_sep_av, pl_sep_rms, ree2_av, ree2_rms, costheta_av,costheta_rms, bondlength_av, bondlength_rms, endz_av, endz_rms))
        print('endz_av:',endz_av)
        outfile2.write('%s ' % types[i])
    '''
    
    # New printing (so I won't have to correct so much when putting into a table):
    
    if chargesims==True:
        outfile.write('%i & %i$\pm$%i   &  %i$\pm$%i   & %i$\pm$%i & %.2f$\pm$%.2f & %.2f$\pm$%.2f & %.1f$\pm$%.1f \ \ \n' % (charge, pl_dist_av,pl_dist_rms, pl_sep_av, pl_sep_rms, ree2_av, ree2_rms, costheta_av,costheta_rms, bondlength_av, bondlength_rms, endz_av,endz_rms))
        outfile2.write('%i ' % charge)
       
    if factorsims==True or krfxfacsims==True or kthetasims==True:
        outfile.write('%.2f & %i & %i$\pm$%i   &  %i$\pm$%i   & %i$\pm$%i & %.4f$\pm$%.4f & %.2f$\pm$%.2f & %.1f$\pm$%.1f  \ \ \n' % (factor, kangle, pl_dist_av,pl_dist_rms, pl_sep_av, pl_sep_rms, ree2_av, ree2_rms, costheta_av,costheta_rms, bondlength_av, bondlength_rms, endz_av,endz_rms))
        outfile2.write('%.2f & %i ' % (factor,kangle))
    
    if spacesims==True:
        outfile.write('%i & %.2f$\pm$%.2f   &  %.2f$\pm$%.2f   & %i$\pm$%i & %.3f$\pm$%.3f & %.2f$\pm$%.2f & %.1f$\pm$%.1f  \ \ \n' % (spacing, pl_dist_av,pl_dist_rms, pl_sep_av, pl_sep_rms, ree2_av, ree2_rms, costheta_av,costheta_rms, bondlength_av, bondlength_rms, endz_av,endz_rms))
        outfile2.write('%i ' % spacing)
    
    if lengthsims==True or lensfxsims==True:
        outfile.write('%i & %i$\pm$%i   &  %i$\pm$%i   & %i$\pm$%i & %.2f$\pm$%.2f & %.2f$\pm$%.2f & %.1f$\pm$%.1f \ \ \n' % (N, pl_dist_av,pl_dist_rms, pl_sep_av, pl_sep_rms, ree2_av, ree2_rms, costheta_av,costheta_rms, bondlength_av, bondlength_rms, endz_av,endz_rms))
        outfile2.write('%i ' % N)
    
    if wallsims==True:
        outfile.write('%.3f & %i$\pm$%i   &  %i$\pm$%i   & %i$\pm$%i & %.2f$\pm$%.2f & %.2f$\pm$%.2f & %.1f$\pm$%.1f \ \ \n' % (wallenergy, pl_dist_av,pl_dist_rms, pl_sep_av, pl_sep_rms, ree2_av, ree2_rms, costheta_av,costheta_rms, bondlength_av, bondlength_rms, endz_av,endz_rms))
        outfile2.write('%.3f ' % wallenergy)
    
    if dielsims==True:
        outfile.write('%i & %i$\pm$%i   &  %i$\pm$%i   & %i$\pm$%i & %.2f$\pm$%.2f & %.2f$\pm$%.2f & %.1f$\pm$%.1f \ \ \n' % (dielectric, pl_dist_av,pl_dist_rms, pl_sep_av, pl_sep_rms, ree2_av, ree2_rms, costheta_av,costheta_rms, bondlength_av, bondlength_rms, endz_av,endz_rms))
        outfile2.write('%i ' % dielectric)
    
    if comparefene==True:
        # type    & $tilde{l}_p$             & $l_p$                   & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$     & $l_b$                          & $z_{N}$
        #types[i] & pl_dist_av+/-pl_dist_rms &  pl_sep_av+/-pl_sep_rms & ree2_av+/-ree2_rms & costheta_av+/-costheta_rms & bondlength_av+/-bondlength_rms & endz_av+/-endz_rms
        outfile.write('%s & %.2f$\pm$%.2f   &  %.2f$\pm$%.2f   & %i$\pm$%i & %.4f$\pm$%.4f & %.2f$\pm$%.2f & %.1f$\pm$%.1f \ \ \n' % (types[i], pl_dist_av,pl_dist_rms, pl_sep_av, pl_sep_rms, ree2_av, ree2_rms, costheta_av,costheta_rms, bondlength_av, bondlength_rms, endz_av, endz_rms))
        print('endz_av:',endz_av)
        outfile2.write('%s ' % types[i])
    
    for j in range(M):
        outfile2.write('& %.3f$\pm$%.3f ' % (endzs_av_chain[j], endzs_rms_chain[j]))
    outfile2.write('& %.3f$\pm$%.3f \ \ \n' % (endz_av, endz_rms))
    
    k+=1

outfile.write('\end{tabular}\n\label{table:quantities_vs_spacing}\n\end{table}')    
outfile.close()
outfile2.write('\end{tabular}}\n\label{table:endz_chains_something}\n\end{table}') 
outfile2.close()

print('Plotting')
xarray = np.zeros(Narray)
xstrng = 'ERROR'
if chargesims==True:
    xarray = charges
    xstrng = r'Charge $q$'
if factorsims==True or krfxfacsims==True or kthetasims==True:
    xarray = kangles
    xstrng = r'$K_\theta$'
if spacesims==True:
    xarray = spacings
    xstrng = r'$l_{grid}$'
if lengthsims==True or lensfxsims==True:
    xarray = lengths
    xstrng = 'N'
if wallsims==True:
    xarray = wallenergy
if dielsims==True:
    xarray = dielsims

print('xarray:',xarray)
print('pl_sep_avg_array:',pl_sep_avg_array)
print('pl_sep_rms_array:',pl_sep_rms_array)

plt.figure(figsize=(6,5))
plt.errorbar(xarray, pl_sep_avg_array, yerr=pl_sep_rms_array, capsize=2)
plt.xlabel(xstrng)
plt.ylabel(r'$l_p$')
plt.title(r'Persistence length $l_p$')
plt.tight_layout()
plt.savefig('persistencelength_by_sep_'+filenamebase+'.png')

plt.figure(figsize=(6,5))
plt.errorbar(xarray, pl_dist_avg_array, yerr=pl_dist_rms_array, capsize=2)
plt.xlabel(xstrng)
plt.ylabel(r'$\tilde{l}_p$')
plt.title(r'Persistence length $\tilde{l}_p$')
plt.tight_layout()
plt.savefig('persistencelength_by_dist_'+filenamebase+'.png')

plt.figure(figsize=(6,5))
plt.errorbar(xarray, ree2_avg_array, yerr=ree2_rms_array, capsize=2)
plt.xlabel(xstrng)
plt.ylabel(r'$<r_{ee}^2>$')
plt.title(r'End-to-end distance $<r_{ee}^2>$')
plt.tight_layout()
plt.savefig('ree_'+filenamebase+'.png')

plt.figure(figsize=(6,5))
plt.errorbar(xarray, costheta_avg_array, yerr=costheta_rms_array, capsize=2)
plt.xlabel(xstrng)
plt.ylabel(r'$<\cos\theta(s=1)>$')
plt.title(r'<$\cos\theta$> for neighbour bonds')
plt.tight_layout()
plt.savefig('costheta_neighbonds_'+filenamebase+'.png')

plt.figure(figsize=(6,5))
plt.errorbar(xarray, endz_avg_array, yerr=endz_rms_array, capsize=2)
plt.xlabel(xstrng)
plt.ylabel(r'$<z_N>$')
plt.title(r'$z$-coordinate of last bead')
plt.tight_layout()
plt.savefig('endz_neighbonds_'+filenamebase+'.png')

print('Done')
