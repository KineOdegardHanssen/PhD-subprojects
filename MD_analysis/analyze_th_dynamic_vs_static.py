import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
from scipy import ndimage                           # For Euclidean distance measurement
import numpy as np
import random
import math
import time
import os
import glob
import copy

# Big and ugly code. Not too much to gain from finesse in this case.

def avg_and_rms(x):
    N = len(x)
    avg = np.mean(x)
    rms = 0
    for i in range(N):
        rms += (x[i]-avg)*(x[i]-avg)
    rms = np.sqrt(rms/(N-1)) 
    return rms, avg

long = False
damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
spacing = 5
psigma  = 1
print('spacing:', spacing)
print('psigma:', psigma)

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T         = 3
testh     = 20     # For transport measurements (default 50, might change for tests)



#-----------------------------PLOTNAMES, ETC.------------------------------#
endlocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/th_stat_vs_dyn/d'+str(spacing)+'/'
if long==False:
    plotname             = endlocation+'ths_h%i_d' % testh + str(spacing) +'_stat_vs_dyn.png'
    plotname_dots_lines  = endlocation+'ths_h%i_d' % testh + str(spacing) +'_dotline_stat_vs_dyn.png'
    plotname_loglog      = endlocation+'ths_h%i_d' % testh + str(spacing) +'_loglog_stat_vs_dyn.png'
    plotname_ylog        = endlocation+'ths_h%i_d' % testh + str(spacing) +'_ylog_stat_vs_dyn.png'
else:
    plotname             = endlocation+'ths_h%i_d' % testh + str(spacing) +'_long_stat_vs_dyn.png'
    plotname_dots_lines  = endlocation+'ths_h%i_d' % testh + str(spacing) +'_dotline_long_stat_vs_dyn.png'
    plotname_loglog      = endlocation+'ths_h%i_d' % testh + str(spacing) +'_loglog_long_stat_vs_dyn.png'
    plotname_ylog        = endlocation+'ths_h%i_d' % testh + str(spacing) +'_ylog_long_stat_vs_dyn.png'

#--------------------------------DYNAMIC-----------------------------------#
confignrs = np.arange(1,1001)

inlocation          = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp + str(psigma) + '/'
filestext            = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])
if long==False:
    infilename           = inlocation+'ths_h%i' % testh +filestext+'.txt'
else:
    infilename           = inlocation+'ths_h%i' % testh +filestext+'_long.txt'

infile    = open(infilename,'r')
firstline = infile.readline()
Nreal     = float(firstline.split()[2]) # Nread in the original script
lines     = infile.readlines()

testh_times = []
exiths      = []

for line in lines:
    words = line.split()
    if len(words)>0:
        testh_times.append(float(words[0]))
        exiths.append(float(words[1]))
infile.close()

## Distribution of times when walker reaches height testh:
Nth = len(testh_times)
if Nth>0:
    maxth = max(testh_times)
    minth = min(testh_times)
    maxzh = max(exiths)
    minzh = min(exiths)
    
    avgth, rmsth = avg_and_rms(testh_times)
    avgzh, rmszh = avg_and_rms(exiths)
else: # Not sure this really helps with anything...
    maxth = 0
    minth = 0
    maxzh = 0
    minzh = 0
    
    avgth = 0; rmsth = 0
    avgzh = 0; rmszh = 0

hist_th, bin_edges_th = np.histogram(testh_times, bins=20)

#define width of each column
width_dynamic = bin_edges_th[1]-bin_edges_th[0]
#standardize each column by dividing with the maximum height
hist_th = hist_th/float(Nreal)


halfwidth = width_dynamic/2.
bin_edges_shifted = bin_edges_th[:-1]+halfwidth

hist_th_dynamic = hist_th
bin_edges_th_dynamic = bin_edges_shifted

#------------------------------STATIC------------------------------#
Nconfigs    = 100
Nplacements = 10
confignrs      = np.arange(1,Nconfigs+1)
beadplacements = np.arange(1,Nplacements+1)


basepath       = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/Spacing'+str(spacing)+'/'
inlocation     = basepath + 'Radius1/Results/'
filestext      = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])+'_placements'+str(beadplacements[0])+'to'+str(beadplacements[-1])

if long==False:
    infilename           = inlocation+'ths_h%i' % testh +filestext+'.txt'
else:
    infilename           = inlocation+'ths_h%i' % testh +filestext+'_long.txt'

infile    = open(infilename,'r')
firstline = infile.readline()
Nreal     = float(firstline.split()[2]) # Nread in the original script
lines     = infile.readlines()

testh_times = []
exiths      = []

for line in lines:
    words = line.split()
    if len(words)>0:
        testh_times.append(float(words[0]))
        exiths.append(float(words[1]))
infile.close()

## Distribution of times when walker reaches height testh:
Nth = len(testh_times)
if Nth>0:
    maxth = max(testh_times)
    minth = min(testh_times)
    maxzh = max(exiths)
    minzh = min(exiths)
    
    avgth, rmsth = avg_and_rms(testh_times)
    avgzh, rmszh = avg_and_rms(exiths)
else: # Not sure this really helps with anything... 
    maxth = 0
    minth = 0
    maxzh = 0
    minzh = 0
    
    avgth = 0; rmsth = 0
    avgzh = 0; rmszh = 0

hist_th, bin_edges_th = np.histogram(testh_times, bins=20)

#define width of each column
width_static = bin_edges_th[1]-bin_edges_th[0]
#standardize each column by dividing with the maximum height
hist_th = hist_th/float(Nreal)

halfwidth = width_static/2.
bin_edges_shifted = bin_edges_th[:-1]+halfwidth

hist_th_static = hist_th
bin_edges_th_static = bin_edges_shifted


#-----------------------------PLOTTING-----------------------------#
# Plotting dots, not histograms
figure()
plt.plot(bin_edges_th_dynamic,hist_th_dynamic, '-o', label='Dynamic brush')
plt.plot(bin_edges_th_static,hist_th_static, '-o', label='Static brush')
plt.xlabel(r'$t_h$ [s]')
plt.ylabel(r'No. of exits/$N_{sims}$') 
plt.title('$t_h$ in system size by d = %i nm, h = %.1f' % (spacing,testh))
plt.legend(loc='upper right')
plt.savefig(plotname_dots_lines)

figure()
plt.loglog(bin_edges_th_dynamic,hist_th_dynamic, '-o', label='Dynamic brush')
plt.loglog(bin_edges_th_static,hist_th_static, '-o', label='Static brush')
plt.xlabel(r'$t_h$ [s]')
plt.ylabel(r'No. of exits/$N_{sims}$') 
plt.title('$t_h$ in system size by d = %i nm, h = %.1f' % (spacing,testh)) 
plt.legend(loc='upper right')
plt.savefig(plotname_loglog)

fig = figure()
ax = fig.add_subplot(1,1,1)
plt.plot(bin_edges_th_dynamic,hist_th_dynamic, '-o', label='Dynamic brush')
plt.plot(bin_edges_th_static,hist_th_static, '-o', label='Static brush')
plt.xlabel(r'$t_h$ [s]')
plt.ylabel(r'No. of exits/$N_{sims}$')
ax.set_yscale('log')
plt.title('$t_h$ in system size by d = %i nm, h = %.1f' % (spacing,testh)) 
plt.legend(loc='upper right')
plt.savefig(plotname_ylog)


