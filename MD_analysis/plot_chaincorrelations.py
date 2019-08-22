import matplotlib.pyplot as plt               # To plot
import numpy as np
import random
import math
import time

M           = 9
N           = 101
Kangle      = 20
Kbond       = 200
charge      = -1
dielectric  = 1
T           = 310
gridspacing = 1

'''
infilename_main = 'cos_bond_correlation_chaingrid_quadratic_M9N101_gridspacing4_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T310_theta0is180_twofirst_are_fixed_LONG.txt'
infilename_ks = 'cos_bond_correlation_nn_bybondnr_chaingrid_quadratic_M9N101_gridspacing4_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T310_theta0is180_twofirst_are_fixed_LONG.txt'
infilename_oneframe = 'cos_bond_correlation_oneframe2_chaingrid_quadratic_M9N101_gridspacing4_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T310_theta0is180_twofirst_are_fixed_LONG.txt'
plotname1 = 'cos_bond_correlation_main_vs_oneframe2_chaingrid_quadratic_M9N101_gridspacing4_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T310_theta0is180_twofirst_are_fixed.png'
plotname2 = 'cos_bond_correlation_vs_ks_chaingrid_quadratic_M9N101_gridspacing4_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T310_theta0is180_twofirst_are_fixed.png'
'''
infilename_main = 'cos_bond_correlation_chaingrid_quadratic_M9N101_gridspacing%i_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T310_theta0is180_twofirst_are_fixed_LONG.txt' % gridspacing
infilename_ks = 'cos_bond_correlation_nn_bybondnr_chaingrid_quadratic_M9N101_gridspacing%i_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T310_theta0is180_twofirst_are_fixed_LONG.txt' % gridspacing
infilename_oneframe = 'cos_bond_correlation_oneframe2_chaingrid_quadratic_M9N101_gridspacing%i_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T310_theta0is180_twofirst_are_fixed_LONG.txt' % gridspacing
plotname1 = 'cos_bond_correlation_main_vs_oneframe2_chaingrid_quadratic_M9N101_gridspacing%i_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T310_theta0is180_twofirst_are_fixed.png' % gridspacing
plotname2 = 'cos_bond_correlation_vs_ks_chaingrid_quadratic_M9N101_gridspacing%i_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T310_theta0is180_twofirst_are_fixed.png' % gridspacing


# Varying the dielectric constant:
infilename_main = 'cos_bond_correlation_chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.png' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
infilename_ks = 'cos_bond_correlation_nn_bybondnr_chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
infilename_oneframe = 'cos_bond_correlation_oneframe2_chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
plotname1 = 'cos_bond_correlation_main_vs_oneframe2_chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)
plotname2 = 'cos_bond_correlation_vs_ks_chaingrid_quadratic_M%iN%i_Langevin_gridspacing%i_Kangle%i_Kbond%i_debye_kappa1_debyecutoff3_charge%i_dielectric%i_T%i_theta0is180_twofirst_are_fixed.txt' % (M,N,gridspacing,Kangle,Kbond,charge,dielectric,T)


separations        = []
corr_main          = []
corr_main_stdv     = []
corr_oneframe      = []
corr_oneframe_stdv = []
ks                 = []
corr_ks            = []
corr_ks_stdv       = []


### Extracting data ###
## Main ##
infile_main     = open(infilename_main,'r')
lines = infile_main.readlines()

for line in lines:
    words = line.split()
    if len(words)>0:
        separations.append(float(words[0]))
        corr_main.append(float(words[1]))
        corr_main_stdv.append(float(words[2]))

separations    = np.array(separations)
corr_main      = np.array(corr_main)
corr_main_stdv = np.array(corr_main_stdv)

infile_main.close()

## One frame ##
infile_oneframe = open(infilename_oneframe,'r')
lines = infile_oneframe.readlines()

for line in lines:
    words = line.split()
    if len(words)>0:
        corr_oneframe.append(float(words[1]))
        corr_oneframe_stdv.append(float(words[2]))

corr_oneframe      = np.array(corr_oneframe)
corr_oneframe_stdv = np.array(corr_oneframe_stdv)

infile_oneframe.close()

## Numbered bonds ##
infile_ks       = open(infilename_ks,'r')
lines = infile_ks.readlines()

for line in lines:
    words = line.split()
    if len(words)>0:
        ks.append(int(words[0]))
        corr_ks.append(float(words[1]))
        corr_ks_stdv.append(float(words[2]))

corr_ks      = np.array(corr_ks)
corr_ks_stdv = np.array(corr_ks_stdv)
#corr_ks[0]   = 1 # Didn't bother setting this in bond_orientation_correlation.py

infile_ks.close()

#print('corr_main:',corr_main)

# Maybe plotting, printing and fitting?
plt.figure(figsize=(6,5))
plt.errorbar(separations, corr_main, corr_main_stdv, capsize=2, label='Main')
plt.plot(separations, corr_main, 'o')
plt.errorbar(separations, corr_oneframe, corr_oneframe_stdv, capsize=2, label='One frame')
plt.plot(separations, corr_oneframe, 'o')
#plt.plot(separation, fittedgraph, label='Fit')
#plt.plot(x_persistencelength, y_persistencelength, '--')
#plt.plot(x_eline, y_eline, '--')
plt.xlabel(r'Chain distance $\Delta$', fontsize=16)
plt.ylabel(r'<$\cos\theta_{\Delta}$>', fontsize=16)
plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
plt.legend(loc="lower left")
plt.title(r'<$\cos\theta_{\Delta}$> vs chain distance $\Delta$', fontsize=16)
plt.savefig(plotname1)

plt.figure(figsize=(6,5))
plt.errorbar(ks, corr_ks, corr_ks_stdv, fmt="none", capsize=2)
plt.plot(ks, corr_ks, 'o')
plt.xlabel(r'Bond nr. $k$', fontsize=16)
plt.ylabel(r'<$\cos\theta_{\Delta}$>', fontsize=16)
plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
plt.title(r'<$\cos\theta_{\Delta}$> vs $k$', fontsize=16)
plt.savefig(plotname2)

