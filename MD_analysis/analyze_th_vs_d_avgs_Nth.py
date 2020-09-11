import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob
import copy

# NB!!! Currently, only plots Nth and Nth/Nsim for short sims.

damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
spacings = [1,1.25,1.5,2,2.5,3,3.5,4,5,6,7,8,10,15,25,50,75,100]
longlst  = [False, False, False, False, True, False, True, False, False, False, False, False, False, False, False, False, False, False]
Nlong = 2
Nsp = len(spacings)
psigma  = 1

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T         = 3
testh     = 20     # For transport measurements (default 50, might change for tests)
confignrs = np.arange(1,1001)

## Initializing lists for the relevant values:
Nread_list  = []
Nth_list    = []
# time
th_avg_list = []
th_rms_list = []
th_max_list = []
th_min_list = []
# Distance (for verification):
zh_avg_list = []
zh_rms_list = []
zh_max_list = []
zh_min_list = []

## Readying paths and plot names:
endlocation          = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/th_vs_d/'
filestext            = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])

plot_Nth   = endlocation + 'Nth_vs_d_h%i.png' % testh
plot_ftr   = endlocation + 'fractiontoreach_vs_d_h%i.png' % testh
plot_th_av = endlocation + 'th_vs_d_av_rms_maxmin_h%i.png' % testh
plot_zh_av = endlocation + 'zh_vs_d_av_rms_maxmin_h%i.png' % testh

# Looping over the spacings:
for i in range(Nsp):
    spacing = spacings[i]
    long    = longlst[i]
    inlocation          = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp + str(psigma) + '/'
    if long==False:
        infilename          = inlocation+'ths_h%i' % testh +filestext+'_analysis.txt'
    else:
        infilename          = inlocation+'ths_h%i' % testh +filestext+'_long_analysis.txt'
    
    ## Reading in values:
    infile = open(infilename,'r')
    lines  = infile.readlines()
    Nread  = int(lines[0].split()[2])
    Nth    = int(lines[1].split()[3])
    min_th = int(lines[3].split()[2])
    max_th = int(lines[4].split()[2])
    avg_th = int(lines[5].split()[4])
    rms_th = int(lines[5].split()[5])
    #
    min_zh = int(lines[7].split()[2])
    max_zh = int(lines[8].split()[2])
    avg_zh = int(lines[9].split()[4])
    rms_zh = int(lines[9].split()[5])
    infile.close()
    
    ## Appending values to list:
    Nread_list.append(Nread)
    Nth_list.append(Nth)
    # time
    th_avg_list.append(avg_th)
    th_rms_list.append(rms_th)
    th_max_list.append(max_th)
    th_min_list.append(min_th)
    # Distance (for verification):
    zh_avg_list.append(avg_zh)
    zh_rms_list.append(rms_zh)
    zh_max_list.append(max_zh)
    zh_min_list.append(min_zh)

plt.figure(figsize=[10,8])
plt.errorbar(spacings, th_avg_list, yerr=th_rms_list, fmt="none", capsize=2)
plt.plot(spacings,th_avg_list, '.', color='tab:blue')
plt.fill_between(spacings, th_min_list,th_max_list, alpha=0.3)
###plt.plot(spacings,th_min_list, '.', color='tab:blue')
###plt.plot(spacings,th_max_list, '.', color='tab:blue')
plt.xlabel('Spacing $d$',fontsize=15)
plt.ylabel(r'Arrival time [s]',fontsize=15)
plt.title('Arrival times at z=%i for a dynamic brush' %testh, fontsize=15)
plt.savefig(plot_th_av)

plt.figure(figsize=[10,8])
plt.errorbar(spacings, zh_avg_list, yerr=zh_rms_list, fmt="none", capsize=2)
plt.plot(spacings,zh_avg_list, '.', color='tab:blue')
plt.fill_between(spacings, zh_min_list,zh_max_list, alpha=0.3)
###plt.plot(spacings,zh_min_list, '.', color='tab:blue')
###plt.plot(spacings,zh_max_list, '.', color='tab:blue')
plt.xlabel('Spacing $d$',fontsize=15)
plt.ylabel(r'z-value at arrival [m]',fontsize=15)
plt.title('z-value at arrival for a dynamic brush' %testh, fontsize=15)
plt.savefig(plot_zh_av)

### For some sims, we shouldn't 
fracNth = np.zeros(Nsp-Nlong)
j = 0
for i in range(Nsp):
    long = longlst[i]
    if long==False:
        fracNth[j] = Nth_list[i]/Nread_list[i]
        j+=1


plt.figure(figsize=[10,8])
plt.plot(spacings, Nth_list_short, label=r'$N_{t,h}$')
plt.plot(spacings, Nread_list_short, label=r'$N_{all sims.}$')
plt.xlabel(r'Spacing $d$',fontsize=15)
plt.ylabel(r'',fontsize=15)
plt.title(r'Number of arrivals $N_{t,h}$', fontsize=15)
plt.legend(loc='upper right')
plt.savefig(plot_Nth)



plt.figure(figsize=[10,8])
plt.plot(spacings, fracNth)
plt.xlabel(r'Spacing $d$',fontsize=15)
plt.ylabel(r'$N_{t,h}$/$N_{all sims.}$',fontsize=15)
plt.title(r'Fraction of beads to arrive at z=%i, $t_{max}=%.1e$' % (testh, simlen), fontsize=15)
plt.savefig(plot_ftr)