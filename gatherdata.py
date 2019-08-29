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
N           = 101
spacing     = 1
gridspacing = spacing
spacings    = [1,2,3,4,5,6,7,8,10,15,40,100]#[1,3,5,10,40]#[1,2,3,4,5,6,7,8,10,15,40,100]
lengths     = [21,61,101,141]
lensp       = [1,3,5,7]
krfacs      = [0.01, 0.05, 0.10, 0.50, 1.00]
#kangles     = [20, 100, 200, 1000, 2000] # krfxfacsims
kangles     = [2000,5000, 8000, 10000, 50000, 125000, 150000] # kthetasims
diels       = [1, 2, 10, 50, 100, 1000] # For spacing 1 nm
#diels       = [1, 2, 10, 100] # For spacing 2 nm
charge      = -1
T           = 310

chargesims  = False
factorsims  = False
spacesims   = False
lengthsims  = False # Varying the length and the spacing
lensfxsims  = False # Varying the length, fixing the spacing
krfxfacsims = False
wallsims    = False
dielsims    = False
kthetasims  = True

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

bondlengths = np.zeros(Narray)
bondrms     = np.zeros(Narray)


if chargesims==True:
    outfilename  = 'table_systemwide_quadratic_gridspacing%i_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varycharge.txt' % (spacing,kangle,kbond,T)
    outfile      = open(outfilename, 'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c}\nCharge & $tilde{l}_p$ & $l_p$ & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$ & $l_b$ & $z_{N}$ \ \ \n\hline\n')
    outfilename2 = 'table_endzs_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,kangle,kbond) +str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varycharge.txt' % T
    outfile2     = open(outfilename2, 'w')
    outfile2.write('resizebox{textwidth}{!}{begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')


if factorsims==True:
    outfilename  = 'table_systemwide_quadratic_M%iN%i_Langevin_Kangle%i_T%i_angle_stretching_only_varyfactor.txt' % (M,N,kangle,T)
    outfile      = open(outfilename, 'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c}\nFactor & $tilde{l}_p$ & $l_p$ & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$ & $l_b$ & $z_{N}$ \ \ \n\hline\n')
    outfilename2 = 'table_endzs_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,kangle,kbond) +str(factor)+'_T%i_theta0is180_twofirst_are_fixed_varyspacing_N%i.txt' % T
    outfile2     = open(outfilename2, 'w')
    outfile2.write('resizebox{textwidth}{!}{begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')


if spacesims==True:
    outfilename  = 'table_systemwide_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varyspacing.txt' % (M,N,kangle,kbond,T)
    outfile      = open(outfilename, 'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c}\nSpacing & $tilde{l}_p$ & $l_p$ & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$ & $l_b$ & $z_{N}$ \ \ \n\hline\n')
    outfilename2 = 'table_endzs_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varyspacing.txt' % T
    outfile2     = open(outfilename2, 'w')
    outfile2.write('resizebox{textwidth}{!}{begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')


if lengthsims==True:
    outfilename  = 'table_systemwide_quadratic_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varylength_contourdsep20.txt' % (kangle,kbond,T)
    outfile      = open(outfilename, 'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c}\nN & $tilde{l}_p$ & $l_p$ & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$ & $l_b$ & $z_{N}$ \ \ \n\hline\n')
    outfilename2 = 'table_endzs_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varylength_contourdsep20.txt' % T
    outfile2     = open(outfilename2, 'w')
    outfile2.write('resizebox{textwidth}{!}{begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')


if lensfxsims==True:
    outfilename  = 'table_systemwide_quadratic_gridspacing%i_Kangle%i_Kbond%i_T%i_debye_kappa1_debyecutoff3_varylength.txt' % (spacing,kangle,kbond,T)
    outfile      = open(outfilename, 'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c}\nN & $tilde{l}_p$ & $l_p$ & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$ & $l_b$ & $z_{N}$ \ \ \n\hline\n')
    outfilename2 = 'table_endzs_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varylength.txt' % T
    outfile2     = open(outfilename2, 'w')
    outfile2.write('resizebox{textwidth}{!}{begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')


if krfxfacsims==True or kthetasims==True:
    outfilename  = 'table_systemwide_quadratic_spacing%i_Kbond%i_T%i_angle_stretching_only_varyfactor_Krfixed.txt' % (spacing,kbond,T)
    outfile      = open(outfilename, 'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|r|c|c|c|c|c|c}\n K_theta/K_r & $K_theta$ & $tilde{l}_p$ & $l_p$ & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$ & $l_b$ & $z_{N}$ \ \ \n\hline\n')
    outfilename2 = 'table_endzs_chaingrid_quadratic_M%iN%i_spacing%i_Langevin_Kbond%i_T%i_theta0is180_twofirst_are_fixed_varyfactor_Krfixed.txt' % (M,N,spacing,kbond,T)
    outfile2     = open(outfilename2, 'w')
    outfile2.write('resizebox{textwidth}{!}{begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')

if wallsims==True:
    outfilename   = 'table_systemwide_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed_wall.txt' % (M,N,spacing,wallenergy,kangle,kbond,charge,T)
    outfile       = open(outfilename, 'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|r|c|c|c|c|c}\n $E_{wall}$ &  $tilde{l}_p$ & $l_p$ & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$ & $l_b$ & $z_{N}$ \ \ \n\hline\n')
    outfilename2 = 'table_endzs_chaingrid_quadratic_M%iN%i_spacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed_wall.txt' % (M,N,spacing,wallenergy,kangle,kbond,charge,T)
    outfile2     = open(outfilename2, 'w')
    outfile2.write('resizebox{textwidth}{!}{begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')

if dielsims==True:
    outfilename  = 'table_systemwide_chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,T)
    outfile      = open(outfilename,'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|r|c|c|c|c|c}\n $\epsilon_{diel}$ &  $tilde{l}_p$ & $l_p$ & $braket{r^2_{ee}}$ & $braket{\cos theta_1}$ & $l_b$ & $z_{N}$ \ \ \n\hline\n')
    outfilename2 = 'table_endz_chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,T)
    outfile2     = open(outfilename2,'w')
    outfile2.write('resizebox{textwidth}{!}{begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Total \ \ \n\hline\n')


k = 0
#for kangle in kangles:
for i in range(Narray):
    print((k+1)/float(Narray))
    
    if chargesims==True:
        charge      = charges[i]
        infilename  = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,kangle,kbond) +str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varydielectric.txt' % T
        infilename2 = 'ree_average_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,kangle,kbond) +str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        ###
        infilename3  = 'persistence_length_chaingrid_quadratic_M9N101_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge%i_T310_theta0is180_twofirst_are_fixed.txt' % charge 
        infilename4  = 'persistence_length_actualdistance_chaingrid_quadratic_M9N101_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge%i_T310_theta0is180_twofirst_are_fixed.txt' % charge 
        ###
        infilename5  = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,kangle,kbond) +str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varydielectric.txt' % T
        infilename6 = 'endzs_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,kangle,kbond) +str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        
    
    if factorsims==True:
        factor = factors[i]
        kbond  = kangle*factors[i]
        infilename  = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,kangle,kbond) +str(factor)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename2 = 'ree_average_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,kangle,kbond) +str(factor)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        ###
        infilename3  = 'persistence_length_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,kangle,kbond) +str(factor)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename4  = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,kangle,kbond) +str(factor)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        ###
        infilename5  = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,kangle,kbond) +str(factor)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename6 = 'endzs_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_factor' % (M,N,kangle,kbond) +str(factor)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        
        
    if spacesims==True:
        spacing     = spacings[i]
        infilename  = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename2 = 'ree_average_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        ###
        infilename3  = 'persistence_length_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename4  = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varygridspacing.txt' % T
        ###
        infilename5 = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varygridspacing.txt' % T
        infilename6 = 'endzs_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varygridspacing.txt' % T
        
        
    if lengthsims==True:
        N            = lengths[i]
        spacing      = lensp[i]
        infilename   = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename2  = 'ree_average_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        ###
        infilename3  = 'persistence_length_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename4  = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varygridspacing.txt' % T
        ###
        infilename5  = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varygridspacing.txt' % T
        infilename6  = 'endzs_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varygridspacing.txt' % T
        
    
    if lensfxsims==True:
        N            = lengths[i]
        infilename   = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename2  = 'ree_average_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        ###
        infilename3  = 'persistence_length_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed.txt' % T
        infilename4  = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varygridspacing.txt' % T
        ###
        infilename5  = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varygridspacing.txt' % T
        infilename6  = 'endzs_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+ str(charge)+'_T%i_theta0is180_twofirst_are_fixed_varygridspacing.txt' % T
        
    
    if krfxfacsims==True:
        factor       = krfacs[i]
        kangle       = kangles[i]
        infilename   = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
        infilename2  = 'ree_average_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
        infilename3  = 'persistence_length_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
        infilename4  = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
        infilename5  = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
        infilename6  = 'endzs_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.2f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
    
    if wallsims==True:
        infilename   = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,wallenergy,kangle,kbond,charge,T)
        infilename2  = 'ree_average_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,wallenergy,kangle,kbond,charge,T)
        infilename3  = 'persistence_length_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,wallenergy,kangle,kbond,charge,T)
        infilename4  = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,wallenergy,kangle,kbond,charge,T)
        infilename5  = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,wallenergy,kangle,kbond,charge,T)
        infilename6  = 'endzs_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,wallenergy,kangle,kbond,charge,T)
        
    if dielsims==True:
        dielectric  = diels[i]
        infilename  = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
        infilename2 = 'ree_average_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
        infilename3  = 'persistence_length_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
        infilename4 = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
        infilename5 = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
        infilename6  = 'endzs_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
    
    
    if kthetasims==True:
        factor       = krfacs[i]
        kangle       = kangles[i]
        infilename   = 'costhetas_neighbourbonds_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.1f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
        infilename2  = 'ree_average_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.1f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
        infilename3  = 'persistence_length_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.1f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
        infilename4  = 'persistence_length_actualdistance_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.1f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
        infilename5  = 'bondlengths_systemwide_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.1f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
        infilename6  = 'endzs_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_factor%.1f_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,spacing,kangle,kbond,factor,T)
    ################################  
    
    #print('infilename:', infilename)
    
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
    
    for j in range(M):
        outfile2.write('& %.3f$\pm$%.3f ' % (endzs_av_chain[j], endzs_rms_chain[j]))
    outfile2.write('& %.3f$\pm$%.3f \ \ \n' % (endz_av, endz_rms))
    
    k+=1


outfile.write('\end{tabular}\n\label{table:quantities_vs_spacing}\n\end{table}')    
outfile.close()
outfile2.write('\end{tabular}}\n\label{table:endz_chains_something}\n\end{table}') 
outfile2.close()
