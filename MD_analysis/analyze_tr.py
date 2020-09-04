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

def avg_and_rms(x):
    N = len(x)
    avg = np.mean(x)
    rms = 0
    for i in range(N):
        rms += (x[i]-avg)*(x[i]-avg)
    rms = np.sqrt(rms/(N-1)) 
    return rms, avg

long = True
damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
spacing = 5
psigma  = 1
print('spacing:', spacing)
print('psigma:', psigma)

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T         = 3
testr     = 20     # For transport measurements (default 50, might change for tests)
confignrs = np.arange(1,1001)


endlocation          = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp + str(psigma) + '/'
filestext            = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])
if long==False:
    infilename           = endlocation+'trs_r%i' % testr +filestext+'.txt'
    outfilename          = endlocation+'trs_r%i' % testr +filestext+'_analysis.txt'
    plotname             = endlocation+'trs_r%i' % testr +filestext+'_dynamic_correct.png'
    plotname_dots_lines  = endlocation+'trs_r%i' % testr +filestext+'_dotline_dynamic_correct.png'
    plotname_loglog      = endlocation+'trs_r%i' % testr +filestext+'_loglog_dynamic_correct.png'
    plotname_ylog        = endlocation+'trs_r%i' % testr +filestext+'_ylog_dynamic_correct.png'
else:
    infilename           = endlocation+'trs_r%i' % testr +filestext+'_long.txt'
    outfilename          = endlocation+'trs_r%i' % testr +filestext+'_long_analysis.txt'
    plotname             = endlocation+'trs_r%i' % testr +filestext+'_long_dynamic_correct.png'
    plotname_dots_lines  = endlocation+'trs_r%i' % testr +filestext+'_dotline_long_dynamic_correct.png'
    plotname_loglog  = endlocation+'trs_r%i' % testr +filestext+'_loglog_long_dynamic_correct.png'
    plotname_ylog    = endlocation+'trs_r%i' % testr +filestext+'_ylog_long_dynamic_correct.png'

infile    = open(infilename,'r')
firstline = infile.readline()
Nreal     = float(firstline.split()[2]) # Nread in the original script
lines     = infile.readlines()

testr_times = []
exitrs      = []

for line in lines:
    words = line.split()
    if len(words)>0:
        testr_times.append(float(words[0]))
        exitrs.append(float(words[1]))
infile.close()

## Distribution of times when walker reaches height testr:
Ntr = len(testr_times)
if Ntr>0:
    maxtr = max(testr_times)
    mintr = min(testr_times)
    maxr  = max(exitrs)
    minr  = min(exitrs)
    
    avgtr, rmstr = avg_and_rms(testr_times)
    avgr,  rmsr = avg_and_rms(exitrs)
else: # Not sure this really helps with anything...
    maxtr = 0
    mintr = 0
    maxzr = 0
    minzr = 0
    
    avgtr = 0; rmstr = 0
    avgr  = 0; rmsr  = 0
outfile = open(outfilename,'w')
outfile.write('Files read: %i\n' % Nreal)
outfile.write('Number of tr\'s: %i\n'% Ntr)
outfile.write('-----------tr------------\n')
outfile.write('Min tr: %.16f\n' % mintr )
outfile.write('Max tr: %.16f\n' % maxtr )
outfile.write('Avg tr, rms tr: %.16f %.16f\n' % (avgtr, rmstr))
outfile.write('-----------xh------------\n')
outfile.write('Min  r: %.16f\n' % minr )
outfile.write('Max  r: %.16f\n' % maxr )
outfile.write('Avg r, rms r: %.16f %.16f\n' % (avgr, rmsr))
outfile.close()

hist_tr, bin_edges_tr = np.histogram(testr_times, bins=20)

# Creating histogram II
# (https://stackoverflow.com/questions/22241240/how-to-normalize-a-histogram-in-python)
#define width of each column
width = bin_edges_tr[1]-bin_edges_tr[0]
#standardize each column by dividing with the maximum height
hist_tr = hist_tr/float(Nreal)
#plot
plt.bar(bin_edges_tr[:-1],hist_tr,width = width)
plt.xlabel(r'$t_r$ [s]')
plt.ylabel(r'No. of arrivals/$N_{sims}$') 
plt.title('$t_r$ in system size by d = %i nm, r = %.1f nm' % (spacing,testr)) 
plt.savefig(plotname)

# Plotting dots, not histograms
halfwidth = width/2.
bin_edges_shifted = bin_edges_tr[:-1]+halfwidth
figure()
plt.plot(bin_edges_shifted,hist_tr, '-o')
plt.xlabel(r'$t_r$ [s]')
plt.ylabel(r'No. of arrivals/$N_{sims}$') 
plt.title('$t_r$ in system size by d = %i nm, r = %.1f nm' % (spacing,testr)) 
plt.savefig(plotname_dots_lines)

figure()
plt.loglog(bin_edges_shifted,hist_tr, '-o')
plt.xlabel(r'$t_r$ [s]')
plt.ylabel(r'No. of arrivals/$N_{sims}$') 
plt.title('$t_h$ in system size by d = %i nm, r = %.1f nm' % (spacing,testr)) 
plt.savefig(plotname_loglog)

fig = figure()
ax = fig.add_subplot(1,1,1)
plt.plot(bin_edges_shifted,hist_tr, '-o')
plt.xlabel(r'$t_r$ [s]')
plt.ylabel(r'No. of arrivals/$N_{sims}$')
ax.set_yscale('log')
plt.title('$t_r$ in system size by d = %i nm, r = %.1f n,' % (spacing,testr)) 
plt.savefig(plotname_ylog)