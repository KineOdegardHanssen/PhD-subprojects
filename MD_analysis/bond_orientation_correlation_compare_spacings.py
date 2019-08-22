import matplotlib.pyplot as plt               # To plot
import numpy as np
import random
import math
import time

Kangle      = 20
Kbond       = 200
T           = 310
M           = 9
N           = 101 ###
ljdebye     = 1.042
epsilon     = ljdebye
sigma       = 1
ljcutoff    = 1.12246204830937
debyecutoff = 3
#factor      = 0.05#250
#Kbond       = 2000#Kangle*factor
#Kangle      = Kbond*factor
charge      = -1
spacing     = 2
gridspacing = spacing
K           = Kangle # Because we used this notation earlier, but using Kangle is less confusing
wallenergy  = 1.042
dielectrics = [1,2,10,100]
spacings    = [1,4,5,7,10,100]#[1,2,3,4,5,6,7,8,10,15,40,100]

spacesims = False
dielsims  = True

if spacesims == True:
    Nsp         = len(spacings)
    outfilename = 'cos_bond_correlation_vs_spacings_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge-1_T%i_theta0is180_twofirst_are_fixed_few.txt' % (M,N,Kangle,Kbond,T)
    plotname    = 'cos_bond_correlation_vs_spacings_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge-1_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,Kangle,Kbond,T)
if dielsims  == True:
    Nsp         = len(dielectrics)
    outfilename = 'cos_bond_correlation_vs_dielectrics_chaingrid_quadratic_M%iN%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge-1_T%i_theta0is180_twofirst_are_fixed_few.txt' % (M,N,Kangle,Kbond,T)

bondcorr_avg_all = []
bondcorr_rms_all = [] 

outfile     = open(outfilename,'w')

for i in range(Nsp):
    if spacesims == True:
        spacing     = spacings[i]
        infilename  = 'cos_bond_correlation_chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge-1_T%i_theta0is180_twofirst_are_fixed_LONG.txt' % (M,N,spacing,Kangle,Kbond,T)
    if dielsims == True:
        dielectric = dielectrics[i]
        infilename = 'cos_bond_correlation_chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
    #infilename_main = 'cos_bond_correlation_chaingrid_quadratic_M9N101_gridspacing%i_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T310_theta0is180_twofirst_are_fixed_LONG.txt' % gridspacing
    '''
    print("infilename:",infilename)
    print('\n')
    print("fnamfolder:", 'cos_bond_correlation_chaingrid_quadratic_M9N101_Langevin_gridspacing1_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_dielectric1_T310_theta0is180_twofirst_are_fixed')
    print('\n\n')
    '''
    infile      = open(infilename,'r')
    lines       = infile.readlines()
    
    separations  = []
    bondcorr_av  = []
    bondcorr_rms = []
    
    for line in lines:
        words = line.split()
        if len(words)>0:
            separations.append(float(words[0]))
            bondcorr_av.append(float(words[1]))
            bondcorr_rms.append(float(words[2]))
    infile.close()
    
    separations  = np.array(separations)
    bondcorr_av  = np.array(bondcorr_av)
    bondcorr_rms = np.array(bondcorr_rms)
    Nsep = len(separations)
    
    bondcorr_avg_all.append(bondcorr_av)
    bondcorr_rms_all.append(bondcorr_rms)
    
    if i==0:
        outfile.write('begin{table}\n\centering\n\caption{}\n resizebox{textwidth}{!}{begin{tabular}{r|c|c|c|c|c|c}\n ')
        for j in range(Nsep):
            outfile.write(' & %.3f' % separations[j])
        outfile.write('\ \ \n\hline\n')
    if spacesims == True:
        outfile.write('%i' % spacing)
    if dielsims  == True:
        outfile.write('%i' % dielectric)
    for j in range(len(separations)):
        outfile.write(' & %.3f$\pm$%.3f' % (bondcorr_av[j], bondcorr_rms[j]))
    outfile.write(' \ \ \n')    

outfile.write('\end{tabular}}\n\label{table:bondcorr_vs_spacing}\n\end{table}')    
outfile.close()

# Plotting
if spacesims==True:
    plt.figure(figsize=(6,5))
    for i in range(Nsp):
        plt.errorbar(separations, bondcorr_avg_all[i], bondcorr_rms_all[i], capsize=2, label=r'$l_g=%i$' % spacings[i])
        plt.plot(separations, bondcorr_avg_all[i], 'o')
    plt.xlabel(r'Chain distance $\Delta$', fontsize=16)
    plt.ylabel(r'$<\cos\theta_{\Delta}>_{corr}$', fontsize=16)
    plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
    plt.legend(loc="lower left")
    plt.title(r'$<\cos\theta_{\Delta}>_{corr}$ vs chain distance $\Delta$', fontsize=16)
    plt.savefig(plotname)

