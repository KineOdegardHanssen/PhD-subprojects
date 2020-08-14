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
psigma  = 1
density = 0.238732414637843 # Yields mass 1 for bead of radius 1 nm
print('spacing:', spacing)
print('psigma:', psigma)

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
Nbins  = int(Nexits/10)      # Don't know if this is the best solution...
hist, bin_edges = np.histogram(exittimes, bins=20)
dbin = bin_edges[1]-bin_edges[0]
dbinhalf = dbin/2

print('exittimes:',exittimes)
print('bin_edges:',bin_edges)

'''
plt.figure(figsize=[10,8])

plt.bar(bin_edges[:-1], hist, width = 0.01, color='#0504aa')
#plt.xlim(min(bin_edges), max(bin_edges))
plt.grid(axis='y')
plt.xlabel(r'Exit time [s]',fontsize=15)
plt.ylabel('Number of exits',fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.title('Exit time for brush, d = %i nm, histogram' % spacing,fontsize=15)
plt.show()
#plt.savefig(plotname_exittimes_binned)
'''

plt.figure(figsize=[10,8])
plt.plot(bin_edges[:-1],hist, 'o')
plt.xlabel(r'Exit time [s]',fontsize=15)
plt.ylabel('Number of exits',fontsize=15)
plt.title('Exit time for brush, d = %i nm, histogram' % spacing,fontsize=15)
plt.show()

'''
plt.figure(figsize=[10,8])

plt.bar(bin_edges[:-1], hist, width = 0.5, color='#0504aa',alpha=0.7)
plt.xlim(min(bin_edges), max(bin_edges))
#plt.grid(axis='y', alpha=0.75)
plt.xlabel('Value',fontsize=15)
plt.ylabel('Frequency',fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel('Frequency',fontsize=15)
plt.title('Normal Distribution Histogram',fontsize=15)
plt.show()
'''

# Creating histogram 
fig, axs = plt.subplots(1, 1, 
                        figsize =(10, 7),  
                        tight_layout = True) 
  
axs.hist(exittimes, bins = 20) 
plt.xlabel(r'Exit time [s]',fontsize=15)
plt.ylabel('Number of exits',fontsize=15)
plt.title('Exit time for brush, d = %i nm, histogram' % spacing,fontsize=15)
  
# Show plot 
plt.show() 
