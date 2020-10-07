import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
from scipy import ndimage                           # For Euclidean distance measurement
import data_treatment as datr
import numpy as np
import random
import math
import time
import os
import glob
import copy

damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
popup_plots = False
spacing    = 100
psigma     = 1
sigma_atom = psigma

spacings = [1,1.25,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,10,15,25,50,75,100]
Nd       = len(spacings)

basepath = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Singletrajectories/'
outfilen = basepath+'Ncoll_vs_d.txt'
plotname = basepath+'Ncoll_vs_d.png'
plotname_small = basepath+'Ncoll_vs_d_smallds.png'

outfile = open(outfilen,'w')

i = 0
ncoll_avg = np.zeros(Nd)
ncoll_rms = np.zeros(Nd)
for spacing in spacings:
    print('spacing:', spacing)
    outfolder_coll = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Singletrajectories/Spacing' + str(spacing) + '/'
    #ncollfilename = outfolder_coll + 'Ncoll_selectedconfigs.txt'
    ncollfilename = outfolder_coll + 'Ncoll_manyconfigs.txt'
    ncollfile     = open(ncollfilename,'r')
    lines = ncollfile.readlines()
    
    counter = 0
    ncoll_list = [] 
    for line in lines:
        ncoll_this    = float(line.split()[1])
        ncoll_avg[i] += ncoll_this
        ncoll_list.append(ncoll_this)
        counter+=1
    
    ncoll_avg[i]/=counter
    #Find rms
    for j in range(counter):
        ncoll_rms[i] += (ncoll_avg[i]-ncoll_list[j])**2
    ncoll_rms[i] = np.sqrt(ncoll_rms[i]/(counter-1))
    
    outfile.write('%.2f %.5e %.5e\n' % (spacing,ncoll_avg[i],ncoll_rms[i]))
    i+=1
outfile.close()

plt.figure(figsize=(6,5))
plt.errorbar(spacings,ncoll_avg,yerr=ncoll_rms, capsize=2)
plt.xlabel('$d$ [nm]')
plt.ylabel(r'$N_{coll}$')
plt.title(r'Average number of collisions in 200 000 time steps vs spacing $d$')
plt.tight_layout()
plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.errorbar(spacings[:13],ncoll_avg[:13],yerr=ncoll_rms[:13], capsize=2)
plt.xlabel('$d$ [nm]')
plt.ylabel(r'$N_{coll}$')
plt.title(r'Average number of collisions in 200 000 time steps vs spacing $d$')
plt.tight_layout()
plt.savefig(plotname_small)
