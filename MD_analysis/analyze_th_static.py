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

long = False
damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
spacing = 8
psigma  = 1
print('spacing:', spacing)
print('psigma:', psigma)
Nconfigs    = 100
Nplacements = 10

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T         = 3
testh     = 50     # For transport measurements (default 50, might change for tests)
confignrs      = np.arange(1,Nconfigs+1)
beadplacements = np.arange(1,Nplacements+1)


basepath        = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/Spacing'+str(spacing)+'/'
endlocation     = basepath + 'Radius1/Results/'
filestext            = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])+'_placements'+str(beadplacements[0])+'to'+str(beadplacements[-1])
################################
infilename           = endlocation+'ths_h%i' % testh +filestext+'.txt'
outfilename          = endlocation+'ths_h%i' % testh +filestext+'_static_analysis.txt'
plotname             = endlocation+'ths_h%i' % testh +filestext+'_static_correct.png'
################################
if long==False:
    infilename           = endlocation+'ths_h%i' % testh +filestext+'.txt'
    outfilename          = endlocation+'ths_h%i' % testh +filestext+'_static_analysis.txt'
    plotname             = endlocation+'ths_h%i' % testh +filestext+'_static_correct.png'
    plotname_dots_lines  = endlocation+'ths_h%i' % testh +filestext+'_dotline_static_correct.png'
    plotname_loglog      = endlocation+'ths_h%i' % testh +filestext+'_loglog_static_correct.png'
    plotname_ylog        = endlocation+'ths_h%i' % testh +filestext+'_ylog_static_correct.png'
else:
    infilename           = endlocation+'ths_h%i' % testh +filestext+'_long.txt'
    outfilename          = endlocation+'ths_h%i' % testh +filestext+'_long_static_analysis.txt'
    plotname             = endlocation+'ths_h%i' % testh +filestext+'_long_static_correct.png'
    plotname_dots_lines  = endlocation+'ths_h%i' % testh +filestext+'_dotline_long_static_correct.png'
    plotname_loglog  = endlocation+'ths_h%i' % testh +filestext+'_loglog_long_static_correct.png'
    plotname_ylog        = endlocation+'ths_h%i' % testh +filestext+'_ylog_long_static_correct.png'

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

outfile = open(outfilename,'w')
outfile.write('Files read: %i\n' % Nreal)
outfile.write('Number of th\'s: %i\n'% Nth)
outfile.write('-----------th------------\n')
outfile.write('Min th: %.16f\n' % minth )
outfile.write('Max th: %.16f\n' % maxth )
outfile.write('Avg th, rms th: %.16f %.16f\n' % (avgth, rmsth))
outfile.write('-----------xh------------\n')
outfile.write('Min zh: %.16f\n' % minzh )
outfile.write('Max zh: %.16f\n' % maxzh )
outfile.write('Avg zh, rms zh: %.16f %.16f\n' % (avgzh, rmszh))
outfile.close()

hist_th, bin_edges_th = np.histogram(testh_times, bins=20)

# Creating histogram II
# (https://stackoverflow.com/questions/22241240/how-to-normalize-a-histogram-in-python)
#define width of each column
width = bin_edges_th[1]-bin_edges_th[0]
#standardize each column by dividing with the maximum height
hist_th = hist_th/float(Nreal)
#plot
plt.bar(bin_edges_th[:-1],hist_th,width = width)
plt.xlabel(r'$t_h$ [s]')
plt.ylabel(r'No. of exits/$N_{sims}$') 
plt.title('$t_h$ in system size by d = %i nm, h = %.1f' % (spacing,testh)) 
plt.savefig(plotname)

# Plotting dots, not histograms
halfwidth = width/2.
bin_edges_shifted = bin_edges_th[:-1]+halfwidth
figure()
plt.plot(bin_edges_shifted,hist_th, '-o')
plt.xlabel(r'$t_h$ [s]')
plt.ylabel(r'No. of exits/$N_{sims}$') 
plt.title('$t_h$ in system size by d = %i nm, h = %.1f' % (spacing,testh)) 
plt.savefig(plotname_dots_lines)

figure()
plt.loglog(bin_edges_shifted,hist_th, '-o')
plt.xlabel(r'$t_h$ [s]')
plt.ylabel(r'No. of exits/$N_{sims}$') 
plt.title('$t_h$ in system size by d = %i nm, h = %.1f' % (spacing,testh)) 
plt.savefig(plotname_loglog)

fig = figure()
ax = fig.add_subplot(1,1,1)
plt.plot(bin_edges_shifted,hist_th, '-o')
plt.xlabel(r'$t_h$ [s]')
plt.ylabel(r'No. of exits/$N_{sims}$')
ax.set_yscale('log')
plt.title('$t_h$ in system size by d = %i nm, h = %.1f' % (spacing,testh)) 
plt.savefig(plotname_ylog)