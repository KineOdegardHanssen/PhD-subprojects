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
spacing = 3
Kthetas = [0.1, 1, 14, 100, 1000]
psigma  = 1
density = 0.238732414637843 # Yields mass 1 for bead of radius 1 nm
print('spacing:', spacing)
print('psigma:', psigma)
N = len(Kthetas) # ... I don't use this (yet)

outpath = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Exittime_vs_Ktheta/Spacing'+str(spacing)+'/Sigma_bead_1/'
outfilename = outpath + 'exittimes_av_rms_range.txt'
nofilename  = outpath + 'Nexits_vs_Ktheta.txt'
plotname    = outpath + 'exittimes_av_rms_range_dynamic_varyKtheta.png'
logplotname = outpath + 'exittimes_av_rms_range_varyKtheta_log.png'
noplotname  = outpath + 'Nexits_vs_Ktheta.png'


outfile = open(outfilename,'w')
outfile.write('mean(exittimes)    rms(exittimes)    min(exittimes)    max(exittimes)\n')

nofile = open(nofilename,'w')

av_exittimes  = []
rms_exittimes = []
min_exittimes = []
max_exittimes = []
Nexits_list   = []

for Ktheta in Kthetas:
    basepath        = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/' % damp + 'Sigma_bead_' +str(psigma) + '/VaryKtheta/Ktheta'+str(Ktheta)+'/'
    endlocation     = basepath
    
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
    Nexits_list.append(Nexits)
    
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
    
    outfile.write('%.2f %.6e %.6e %.6e %.6e\n' % (Ktheta, av_exittime,rms_exittime,min_exittime,max_exittime)) # Resolution on this thing?
    nofile.write('%.2f %i\n' % (Ktheta,Nexits))
    
    av_exittimes.append(av_exittime)
    rms_exittimes.append(rms_exittime)
    min_exittimes.append(min_exittime)
    max_exittimes.append(max_exittime)

outfile.close()
nofile.close()

'''
for i in range(N):
    print('Ktheta:',Kthetas[i], '; av_exittimes:', av_exittimes[i], '; rms_exittimes:', rms_exittimes[i],  '; min_exittimes:', min_exittimes[i],  '; max_exittimes:', max_exittimes[i] )
'''

plt.figure(figsize=[10,8])
plt.errorbar(Kthetas, av_exittimes, yerr=rms_exittimes, fmt="none", capsize=2)
plt.plot(Kthetas,av_exittimes, '.', color='tab:blue')
plt.fill_between(Kthetas, min_exittimes,max_exittimes, alpha=0.3)
###plt.plot(Kthetas,min_exittimes, '.', color='tab:blue')
###plt.plot(Kthetas,max_exittimes, '.', color='tab:blue')
plt.xlabel(r'$K_\theta$',fontsize=15)
plt.ylabel(r'Exit time [s]',fontsize=15)
plt.title('Exit time for dynamic brush, histogram',fontsize=15)
plt.savefig(plotname)

fig = plt.figure(figsize=[10,8])
ax = fig.add_subplot(1,1,1)
plt.errorbar(Kthetas, av_exittimes, yerr=rms_exittimes, fmt="none", capsize=2)
plt.plot(Kthetas,av_exittimes, '.', color='tab:blue')
plt.fill_between(Kthetas, min_exittimes,max_exittimes, alpha=0.3)
ax.set_xscale('log')
#plt.set_xscale('log')
###plt.plot(Kthetas,min_exittimes, '.', color='tab:blue')
###plt.plot(Kthetas,max_exittimes, '.', color='tab:blue')
plt.xlabel(r'$K_\theta$',fontsize=15)
plt.ylabel(r'Exit time [s]',fontsize=15)
plt.title('Exit time for dynamic brush, histogram',fontsize=15)
plt.savefig(logplotname)

plt.figure(figsize=[10,8])
plt.plot(Kthetas,Nexits_list, '-o')
plt.xlabel(r'$K_\theta$',fontsize=15)
plt.ylabel(r'Number of exits',fontsize=15)
plt.title('Number of exits for dynamic brush, histogram',fontsize=15)
plt.savefig(noplotname)
plt.show()


'''
print('spacings:',spacings)
print('av_exittimes:',av_exittimes)
print('rms_exittimes:',rms_exittimes)
print('min_exittimes:',min_exittimes)
print('max_exittimes:',max_exittimes)
'''