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

## NB NB NB!: In this script, I divide by the bin length as well as Nreal.

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
spacings = [3,5,7]
psigma  = 1

# Figuring out total simulation length
Nsteps       = 1000000
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime
totlen       = Nsteps*timestepsize

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T         = 3
testh     = 20     # For transport measurements (default 50, might change for tests)
hists     = []
edges     = []
confignrs = np.arange(1,1001)
Nd = len(spacings)

# Locations, text snippets, etc.
baselocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/'
outlocation  = baselocation + 'd_vs_tr/'
filestext    = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])

# Plot names:
hstring = 'r' + '%i_' % testh
for i in range(Nd):
    hstring = hstring + 'd'+str(spacings[i])+'_' 

if long==False:
    plotname_dots_lines  = outlocation+ 'trs_' + hstring + filestext+'_dotline_dynamic_correct.png'
    plotname_loglog      = outlocation+ 'trs_' + hstring + filestext+'_loglog_dynamic_correct.png'
    plotname_ylog        = outlocation+ 'trs_' + hstring + filestext+'_ylog_dynamic_correct.png'
else:
    plotname_dots_lines  = outlocation+ 'trs_' + hstring + filestext+'_dotline_long_dynamic_correct.png'
    plotname_loglog      = outlocation+'trs_'+ hstring + filestext +'_loglog_long_dynamic_correct.png'
    plotname_ylog        = outlocation+'trs_'+ hstring + filestext+'_ylog_long_dynamic_correct.png'


for i in range(Nd):
    spacing        = spacings[i]
    inlocation     = baselocation + 'Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp + str(psigma) + '/'
    if long==False:
        infilename = inlocation+'trs_r%i' % testh +filestext+'.txt'
    else:
        infilename = inlocation+'trs_r%i' % testh +filestext+'_long.txt'
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
    infile.close()
    
    hist_th, bin_edges_th = np.histogram(testh_times, bins=20)


    width = bin_edges_th[1]-bin_edges_th[0]
    hist_th = hist_th/float(Nreal) #*width)*totlen # Is this thing properly normalized?
    
    # Plotting dots, not histograms
    halfwidth = width/2.
    bin_edges_shifted = bin_edges_th[:-1]+halfwidth
    
    hists.append(hist_th)
    edges.append(bin_edges_shifted)

figure()
for i in range(Nd):
    plt.plot(edges[i],hists[i], '-o', label='d=%.2f' % spacings[i])
plt.xlabel(r'$t_r$ [s]')
plt.ylabel(r'No. of arrivals/$N_{sims}$') 
plt.title('$t_r$ in system size by h = %i nm' % testh) 
plt.legend(loc='upper right')
plt.savefig(plotname_dots_lines)

figure()
for i in range(Nd):
    plt.loglog(edges[i],hists[i], '-o', label='d=%.2f' % spacings[i])
plt.xlabel(r'$t_r$ [s]')
plt.ylabel(r'No. of arrivals/$N_{sims}$') 
plt.title('$t_r$ in system size by r = %i nm' % testh) 
plt.legend(loc='upper right')
plt.savefig(plotname_loglog)

fig = figure()
ax = fig.add_subplot(1,1,1)
for i in range(Nd):
    plt.plot(edges[i],hists[i], '-o', label='d=%.2f' % spacings[i])
plt.xlabel(r'$t_r$ [s]')
plt.ylabel(r'No. of arrivals/$N_{sims}$')
ax.set_yscale('log')
plt.title('$t_r$ in system size by r = %i nm' % testh) 
plt.legend(loc='upper right')
plt.savefig(plotname_ylog)