import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob

# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
bulkdiffusion = False#False
substrate     = False
cutit         = False
maxh          = 55.940983199999756

spacings = [1,2,3,5,10,50,100] #[7,10,15,25]#[3,4,5,6]#[]#[3,4,5,7,10,15,25]
psigma   = 1
damp     = 10

Nsteps = 2001
bnrs      = np.arange(1,11)
confignrs = np.arange(1,101)

if bulkdiffusion==True:
    parentfolder = 'Pure_bulk/'
    filestext     = '_seed0to1000'
    if substrate==True:
        parentfolder = 'Bulk_substrate/'
else:
    parentfolder = 'Brush/'
    filestext     = '_config'+str(confignrs[0])+'to'+str(confignrs[-1])+'_placements'+str(bnrs[0])+'to'+str(bnrs[-1])

print('psigma:', psigma)
print('type(psigma):', type(psigma))
# Make file names
outlocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/D_vs_d/Nocut/'
plotname    = outlocation + 'dpardz_together.png'
outfilename = outlocation + 'dpardz_together.txt'

outfile     = open(outfilename, 'w')

counter = 0
allgraphs = []
for spacing in spacings:    
    counter += 1
    
    # Filenames:
    endlocation   = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/Spacing'+str(spacing)+'/Radius'+ str(psigma) + '/Nocut/'
    if cutit==True:
        endlocation = endlocation + 'maxh'+ str(maxh)+'/'
    infilename  = endlocation+'av_ds'+filestext+'.txt'

    # I need to set the file name in an easier way, but for now I just use this:  ## Might want to add a loop too, if I have more files...

    # Should divide into folders in a more thorough manner?
    # Extracting the correct names (all of them)
    plotseed = 0
    plotdirs = False
    test_sectioned = False
    unitlength   = 1e-9
    unittime     = 2.38e-11 # s
    timestepsize = 0.00045*unittime
    print('timestepsize:', timestepsize)
    Npartitions = 5 # For extracting more walks from one file (but is it really such a random walk here...?)

    # Set up arrays
    times  = np.zeros(Nsteps)
    dz2s   = np.zeros(Nsteps) 
    dpar2s = np.zeros(Nsteps)
    dpardz = np.zeros(Nsteps)

    # Read file
    infile = open(infilename,'r')
    header = infile.readline()
    lines  = infile.readlines()
    
    outfile.write('%.2f\n' % spacing)
    #    	0			1		2		3		4			5		6
    # (times_single[i], times_single_real[i], averageRs_SI[i], averagedxs_SI[i], averagedys_SI[i], averagedzs_SI[i], averagedparallel_SI[i]))
    #Loop over sections
    for i in range(1,Nsteps):
        line  = lines[i]
        words = line.split()
        times[i]  = float(words[1])
        dz2s[i]   = float(words[5])
        dpar2s[i] = float(words[6])
        dpardz[i] = dpar2s[i]/dz2s[i]
        outfile.write('%.16e ' % dpardz[i])
        
    outfile.write('\n')
    allgraphs.append(dpardz)
    infile.close()

    print('spacing:', spacing)
    print('psigma:', psigma)
    print('damp:', damp)

plt.figure(figsize=(6,5))
for i in range(counter):
    thisdpardz = allgraphs[i]
    plt.plot(times[1:], thisdpardz[1:], label='d=%.2f' % spacings[i])
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.xlabel(r'Time (s)')
plt.ylabel(r'$\frac{<dx^2+dy^2>}{<dz^2>}$')
plt.title(r'$\frac{<dx^2+dy^2>}{<dz^2>}$ vs time SI')
plt.tight_layout()
plt.legend(loc='upper right')
plt.savefig(plotname)

'''# If a lot of graphs:
plt.figure(figsize=(6,5))
ax = plt.subplot(111)
for i in range(counter):
    thisdpardz = allgraphs[i]
    ax.plot(times[1:], thisdpardz[1:], label='d=%.2f' % spacings[i])
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)#, prop={'size': 20})
plt.xlabel(r'Time (s)')
plt.ylabel(r'$\frac{<dx^2+dy^2>}{<dz^2>}$')
plt.title(r'$\frac{<dx^2+dy^2>}{<dz^2>}$ vs time SI')
plt.savefig(plotname)
'''