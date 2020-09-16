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

long = False
damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
if long==False:
    spacings = [1,1.25,1.5,2,3,4,5,6,7,8,10,15,25,50,75,100]
else:
    spacings = [2.5,3.5,5]
areas = []
psigma  = 1

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T         = 3
testr     = 20     # For transport measurements (default 50, might change for tests)
confignrs = np.arange(1,1001)

endlocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/d_vs_tr/'
filestext            = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])
outfilename   = endlocation+'trs_r%i' % testr +filestext+'_areas.txt'
plotname_all =  endlocation+'trs_r%i' % testr +filestext+'_areas.png'

outfile = open(outfilename, 'w')

for spacing in spacings:
    endlocation          = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp + str(psigma) + '/'  
    if long==False:
        infilename           = endlocation+'trs_r%i' % testr +filestext+'.txt'
        plotname             = endlocation+'trs_r%i' % testr +filestext+'_dynamic_withfit.png'
    else:
        infilename           = endlocation+'trs_r%i' % testr +filestext+'_long.txt'
        plotname             = endlocation+'trs_r%i' % testr +filestext+'_long_dynamic_withfit.png'
    
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
    
    hist_tr, bin_edges_tr = np.histogram(testr_times, bins=20)
    
    # Creating histogram and treating data:
    # (https://stackoverflow.com/questions/22241240/how-to-normalize-a-histogram-in-python)
    #define width of each column
    width = bin_edges_tr[1]-bin_edges_tr[0]
    midpoints = bin_edges_tr[:-1] + width/2.
    #standardize each column by dividing with the maximum height
    hist_tr = hist_tr/float(Nreal)
    Nhist   = len(hist_tr)
    
    # Base area:
    area = 0
    for i in range(Nhist):
        area += hist_tr[i]  # Will multiply by width at the end (currently constant bin width)
    
    # Find maximum for finding the tail: # Add a few points afterwards, maybe?
    maxind = np.argmax(hist_tr)
    
    # Find the tail through ylog:
    # Suitable indices after the tail:
    #histend  = hist_tr[maxind+1:]
    #log_ind  = np.where(histend>0) + maxind+1
    log_ind  = np.where(hist_tr[maxind+1:]>0)
    t_li     = midpoints[log_ind]-midpoints[maxind]
    hist_log = np.log(hist_tr[log_ind])
    
    print('log_ind:',log_ind)
    print('hist_log:',hist_log)
    print('hist_tr[log_ind]:',hist_tr[log_ind])
    print('hist_tr:',hist_tr)
    try:
        coeffs, covs = np.polyfit(t_li, hist_log, 1)
    except:
        print('Empty list. Exception raised. d = %.2f' % spacing)
        area *= width
        outfile.write('%.2f %.16e\n' % (spacing,area))
        areas.append(area)
        continue    
    print('coeffs:', coeffs)
    try:
        slope = coeffs[0]
    except:
        print('Negative infinity. Exception raised. d = %.2f' % spacing)
        area *= width
        outfile.write('%.2f %.16e\n' % (spacing,area))
        areas.append(area)
        continue
    const = coeffs[1]
    slope_rms = np.sqrt(covs[0,0])
    const_rms = np.sqrt(covs[1,1])
    fitfunc   = const*np.exp(slope*t_li)
    
    textr_array = []
    vextr_array = []
    
    nsteps = 10000 # Is this excessive?
    if slope>0:
        counter    = 0
        t_extr     = midpoints[-1]-midpoints[maxind] + width
        value_extr = const*np.exp(slope*t_extr) # Or should I just use a linear function?
        while value_extr>0 or counter>nsteps:
            area += value_extr # WARNING: This is for const. bin size
            # Appending value into array
            textr_array.append(t_extr+midpoints[maxind])
            vextr_array.append(value_extr)
            # Updating values for next check
            counter += 1
            t_extr  += width
            value_extr = const*np.exp(slope*t_extr)
        plt.figure()
        plt.plot(midpoints,hist_tr,'-o')
        plt.plot(textr_array,vextr_array,'--')
        plt.plot(t_li,fitfunc)
        plt.xlabel('t [s]')
        plt.ylabel(r'No. of arrivals/$N_{sims}$')
        plt.title('$t_r$ in system size by d = %i nm, r = %.1f n,' % (spacing,testr))
        plt.savefig(plotname)
        
        # Do more stuff:
        area *= width
        outfile.write('%.2f %.16e\n' % (spacing,area))
    else:
        area *= width
        outfile.write('%.2f %.16e WARNING: AREA NOT EXTRAPOLATED\n' % (spacing,area))
    area.append(areas)
outfile.close()

plt.figure(figsize=(6,5))
plt.plot(spacings, areas)
plt.xlabel(r'Spacing $d$ [nm]')
plt.ylabel(r'Area under No. of arrivals/$N_{sims}$ vs $t$')
plt.title(r'Area vs d under $t_r$ vs $N_{arrival}/N_{sims}$ graph')
plt.tight_layout()
plt.savefig(plotname_all)

print('done')
