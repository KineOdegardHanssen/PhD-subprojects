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

Nconfigs    = 100#100
Nplacements = 10#10
#
damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
popup_plots = False
spacings = [1, 1.25, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 15, 25, 50, 75, 100]
psigma  = 1
density = 0.238732414637843 # Yields mass 1 for bead of radius 1 nm
print('spacing:', spacing)
print('psigma:', psigma)
N = len(spacings) # ... I don't use this (yet)

outpath = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/'
outfilename = outpath + 'exittimes_av_rms_range.txt'
plotname    = outpath + 'exittimes_av_rms_range.py'

outfile = open(outfilename,'w')
outfile.write('mean(exittimes)    rms(exittimes)    min(exittimes)    max(exittimes)\n')

av_exittimes  = []
rms_exittimes = []
min_exittimes = []
max_exittimes = []

for spacing in spacings:
    basepath        = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/Spacing'+str(spacing)+'/'
    basepath        = basepath + 'Radius' + str(psigma) + '/'
    location_config = basepath +'Initial_configs/Before_bead/'
    endlocation     = basepath + 'Results/'
    
    infilename = endlocation+'exittimes.txt'
    
    infile = open(infilename, 'r')
    header = infile.readline()
    lines  = infile.readlines()
    Nlines = len(lines) 
    
    exittimes = [] ##np.zeros(Nlines) # Or use a list?
    
    for line in lines:
        words = line.split()
        if len(words)>0:
            exittimes.append(float(words[0]))
    
    Nexits = len(exittimes)
    
    av_exittime  = 0
    rms_exittime = 0
    min_exittime = 0
    max_exittime = 0
    if Nexits!=0:
        av_exittime = np.mean(exittimes)
        for i in range(Nexits):
            rms_exittime += (av_exittime-exittimes[i])**2
        rms_exittime = np.sqrt(rms_exittime/(Nexits-1))
        min_exittime = min(exittimes)
        max_exittime = max(exittimes)
    
    outfile.write('%.6e %.6e %.6e %.6e\n' % (av_exittime,rms_exittime,min_exittime,max_exittime)) # Resolution on this thing?
    
    av_exittimes.append(av_exittime)
    rms_exittimes.append(rms_exittime)
    min_exittimes.append(min_exittime)
    max_exittimes.append(max_exittime)

outfile.close()

'''
for i in range(N):
    print('spacing:',spacings[i], '; av_exittimes:', av_exittimes[i], '; rms_exittimes:', rms_exittimes[i],  '; min_exittimes:', min_exittimes[i],  '; max_exittimes:', max_exittimes[i] )
'''

plt.figure(figsize=[10,8])
plt.errorbar(spacings, av_exittimes, yerr=rms_exittimes, fmt="none", capsize=2)
plt.plot(spacings,av_exittimes, '.', color='tab:blue')
plt.fill_between(spacings, min_exittimes,max_exittimes, alpha=0.3)
###plt.plot(spacings,min_exittimes, '.', color='tab:blue')
###plt.plot(spacings,max_exittimes, '.', color='tab:blue')
plt.xlabel('Spacing d',fontsize=15)
plt.ylabel(r'Exit time [s]',fontsize=15)
plt.title('Exit time for brush histogram',fontsize=15)
plt.show()

'''
print('spacings:',spacings)
print('av_exittimes:',av_exittimes)
print('rms_exittimes:',rms_exittimes)
print('min_exittimes:',min_exittimes)
print('max_exittimes:',max_exittimes)
'''