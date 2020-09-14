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


long = False
damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
if long==True:
    spacings = [2.5,3.5,5]
else:
    spacings = [1,1.25,1.5,2,3,4,5,6,7,8,10,15,25,50,75,100]
Nsp = len(spacings)
psigma  = 1

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T         = 3
testh     = 20     # For transport measurements (default 50, might change for tests)
confignrs = np.arange(1,1001)
if long==True:
    Nsteps = 10001
else:
    Nsteps = 2001

N10 = spacings.index(10)

writeevery   = 100
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime*writeevery
simlen       = timestepsize*Nsteps

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
endlocation          = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/d_vs_tr/'
filestext            = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])

if long==True:
    plot_Nth   = endlocation + 'Ntr_vs_d_r%i_long.png' % testh
    plot_ftr   = endlocation + 'fractiontoreach_vs_d_r%i_long.png' % testh
    plot_th_av = endlocation + 'tr_vs_d_av_rms_maxmin_r%i_long.png' % testh
    plot_zh_av = endlocation + 'r_vs_d_av_rms_maxmin_r%i_long.png' % testh
    plot_th_av_short = endlocation + 'tr_vs_d_av_rms_maxmin_smallds_r%i_long.png' % testh
else:
    plot_Nth   = endlocation + 'Ntr_vs_d_r%i.png' % testh
    plot_ftr   = endlocation + 'fractiontoreach_vs_d_r%i.png' % testh
    plot_th_av = endlocation + 'tr_vs_d_av_rms_maxmin_r%i.png' % testh
    plot_zh_av = endlocation + 'r_vs_d_av_rms_maxmin_r%i.png' % testh
    plot_th_av_short = endlocation + 'tr_vs_d_av_rms_maxmin_smallds_r%i.png' % testh

# Looping over the spacings:
for spacing in spacings:
    inlocation          = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp + str(psigma) + '/'
    if long==False:
        infilename          = inlocation+'trs_r%i' % testh +filestext+'_analysis.txt'
    else:
        infilename          = inlocation+'trs_r%i' % testh +filestext+'_long_analysis.txt'
    
    ## Reading in values:
    infile = open(infilename,'r')
    lines  = infile.readlines()
    Nread  = int(lines[0].split()[2])
    Nth    = int(lines[1].split()[3])
    min_th = float(lines[3].split()[2])
    max_th = float(lines[4].split()[2])
    avg_th = float(lines[5].split()[4])
    rms_th = float(lines[5].split()[5])
    #
    min_zh = float(lines[7].split()[2])
    max_zh = float(lines[8].split()[2])
    avg_zh = float(lines[9].split()[4])
    rms_zh = float(lines[9].split()[5])
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
plt.title(r'Arrival times at $\sqrt{x^2+y^2}$=%i for a dynamic brush' %testh, fontsize=15)
plt.savefig(plot_th_av)

if N10>1:
    plt.figure(figsize=[10,8])
    plt.errorbar(spacings[0:N10], th_avg_list[0:N10], yerr=th_rms_list[0:N10], fmt="none", capsize=2)
    plt.plot(spacings[0:N10],th_avg_list[0:N10], '.', color='tab:blue')
    plt.fill_between(spacings[0:N10], th_min_list[0:N10],th_max_list[0:N10], alpha=0.3)
    plt.xlabel('Spacing $d$',fontsize=15)
    plt.ylabel(r'Arrival time [s]',fontsize=15)
    plt.title(r'Arrival times at $\sqrt{x^2+y^2}$=%i for a dynamic brush' %testh, fontsize=15)
    plt.savefig(plot_th_av_short)

plt.figure(figsize=[10,8])
plt.errorbar(spacings, zh_avg_list, yerr=zh_rms_list, fmt="none", capsize=2)
plt.plot(spacings,zh_avg_list, '.', color='tab:blue')
plt.fill_between(spacings, zh_min_list,zh_max_list, alpha=0.3)
###plt.plot(spacings,zh_min_list, '.', color='tab:blue')
###plt.plot(spacings,zh_max_list, '.', color='tab:blue')
plt.xlabel('Spacing $d$',fontsize=15)
plt.ylabel(r'$\sqrt{x^2+y^2}$ at arrival [m]',fontsize=15)
plt.title('$\sqrt{x^2+y^2}$ at arrival for a dynamic brush', fontsize=15)
plt.savefig(plot_zh_av)


plt.figure(figsize=[10,8])
plt.plot(spacings, Nth_list, label=r'$N_{t,r}$')
plt.plot(spacings, Nread_list, label=r'$N_{all sims.}$')
plt.xlabel(r'Spacing $d$',fontsize=15)
plt.ylabel(r'Number of arrivals $N_{t,r}$',fontsize=15)
plt.title(r'Number of arrivals $N_{t,r}$ at $\sqrt{x^2+y^2}$=%i' % testh, fontsize=15)
plt.legend(loc='lower right')
plt.savefig(plot_Nth)

fracNth = np.zeros(Nsp)
for i in range(Nsp):
    fracNth[i] = Nth_list[i]/Nread_list[i]

plt.figure(figsize=[10,8])
plt.plot(spacings, fracNth)
plt.xlabel(r'Spacing $d$',fontsize=15)
plt.ylabel(r'$N_{t,r}$/$N_{all sims.}$',fontsize=15)
plt.title(r'Fraction of beads to arrive at $\sqrt{x^2+y^2}$=%i, $t_{max}=%.1e$' % (testh, simlen), fontsize=15)
plt.savefig(plot_ftr)