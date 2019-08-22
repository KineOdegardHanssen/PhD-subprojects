import matplotlib.pyplot as plt               # To plot
from scipy.optimize import curve_fit
import numpy as np
import random
import math
import time

start_time = time.process_time()

M            = 9
N            = 101
kangle       = 20
kbond        = 200
Kangle       = kangle
Kbond        = kbond
#factors     = [0.1,1,10,100,250]
#charges     = [0,-1,-5,-10] 
#spacings    = [1,2,3,4,5,10,40,100]
N            = 101
spacing      = 1
gridspacing  = spacing
spacings     = [1,2,3,4,5,6,7,8,10,15,40,100]#[1,2,3,4,5,8,10,15,40,100]
#spacings     = [40]
#dielectrics  = [1,2,10,100] # For lg = 2 nm
dielectrics  = [1,2,10,50,100,1000]
#lengths     = [21,61,101,141]
#lensp       = [1,3,5,7]
#krfacs      = [0.01, 0.05, 0.10, 0.50, 1.00]
#kangles     = [20, 100, 200, 1000, 2000]
wallenergies = [1.042]
charge       = -1
T            = 310

spacesims   = False
dielsims    = True 
wallsims    = False

if spacesims == True:
    Nsp = len(spacings)
    outfilename = 'table_radgyrs_chaingrid_quadratic_M%iN%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed_varyspacing.txt' % (M,N,wallenergy,kangle,kbond,charge,T)
    outfile = open(outfilename,'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\nSpacing/Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Average \ \ \n\hline\n')

if dielsims  == True:
    Nsp = len(dielectrics)   
    outfilename   = 'table_radgyrs_chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed_varydielectric.txt' % (M,N,gridspacing,Kangle,Kbond,charge,T) 
    outfile = open(outfilename,'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\nDielectric/Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Average \ \ \n\hline\n')   

if wallsims  == True:
    Nsp = len(wallenergies)
    outfilename   = 'table_radgyrs_chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed_varywallenergy.txt' % (M,N,gridspacing,Kangle,Kbond,charge,T) 
    outfile = open(outfilename,'w')
    outfile.write('begin{table}\n\centering\n\caption{}\n begin{tabular}{r|c|c|c|c|c|c|c|c|c|c}\n $\epsilon_w$/Chain & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & Average \ \ \n\hline\n')

totalaverage = np.zeros(Nsp)
totalrms     = np.zeros(Nsp)

for i in range(Nsp):
    if spacesims == True:
        spacing = spacings[i]
        outfile.write('%i' % spacing)
        infilename = 'log.radgyrs_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge' % (M,N,spacing,kangle,kbond)+str(charge)+'_T%i_theta0is180_twofirst_are_fixed' % T
    if dielsims == True:
        dielectric = dielectrics[i]
        outfile.write('%i' % dielectric)
        infilename   = 'log.chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
    if wallsims == True:
        wallenergy = wallenergies[i]
        outfile.write('%.3f' % wallenergy)
        infilename = 'log.chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed_with_rgyr' % (M,N,spacing,wallenergy,kangle,kbond,charge,T)
    # Readying arrays:
    radgyr_average  = np.zeros(M)
    radgyr_stdv     = np.zeros(M)
    # This is really not the optimal solution:
    allradgyrs_vals = [] 
    allradgyrs_inds = []
    
    infile = open(infilename,'r')
    lines = infile.readlines()
    
    #print('infilename:', infilename)
    
    # Finding the mean and rms:
    # Finding the mean:
    starter1 = 0
    starter2 = 0
    counter  = 0
    for line in lines:
        words = line.split()
        #print('words=',words)
        #print('starter1:', starter1, '; starter2:', starter2)
        if len(words)>2:
            if words[1]=='Run' and words[2]=='and':
                # Finding the line: ####################### Run and write to file #########################################
                starter1 = 1
                #print('First mark hit')
        #if starter1==1:
        #    print(words)
        if starter1==1 and starter2==1:
            # Test if we should break:
            if len(words)>0:
                if words[0]=='WARNING:' or words[0]=='Loop':
                    break
            #print('Starting to read data')
            if len(words)==12 or len(words)==18:
                #print('I am in')
                if len(words)==12:
                    addon = 3
                else:
                    addon = 9
                for j in range(M):
                    #print(words)
                    thisvalue          = float(words[j+addon])
                    radgyr_average[j] += thisvalue
                    allradgyrs_vals.append(thisvalue)
                    allradgyrs_inds.append(j)
                counter+=1
        if starter1==1 and starter2==0:
            if len(words)>0:
                if words[0]=='Step':
                    starter2=1
                    #print('Second mark hit')
    infile.close()
    
    radgyr_average  /= counter
    totalaverage[i]  = np.mean(radgyr_average)
    # Finding the rms:
    for j in range(len(allradgyrs_vals)):
        chain               = allradgyrs_inds[j]
        val                 = allradgyrs_vals[j]
        radgyr_stdv[chain] += (radgyr_average[chain]-val)**2
        totalrms[i]         = (totalaverage[i]-val)**2
    
    totalrms[i] = np.sqrt(totalrms[i]/(counter-1))
    for j in range(M):
        radgyr_stdv[j] = np.sqrt(radgyr_stdv[j]/(counter-1))
        outfile.write(' & %.3f$\pm$%.3f' % (radgyr_average[j], radgyr_stdv[j]))
    outfile.write(' & %.4f$\pm$%.4f \ \ \n' % (totalaverage[i], totalrms[i]))
    
outfile.write('\end{tabular}\n\label{table:radgyrs_chain_and_total_something}\n\end{table}')
outfile.close()
