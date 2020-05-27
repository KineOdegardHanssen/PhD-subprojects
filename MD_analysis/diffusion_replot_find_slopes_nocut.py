from pylab import *
import matplotlib.pyplot as plt                     # To plot
import numpy as np
import math

damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
spacing = 1.5
psigma  = 1
print('spacing:', spacing)
print('psigma:', psigma)
# I need to set the file name in an easier way, but for now I just use this:  ## Might want to add a loop too, if I have more files...

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T            = 3
confignrs    = np.arange(1,1001)#300)#20)#101)#1001)#22) 
Nseeds       = len(confignrs)      # So that I don't have to change that much
Nsteps       = 2001 # 20001 before
unitlength   = 1e-9
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime#*writeevery
Npartitions  = 5 # For extracting more walks from one file (but is it really such a random walk here...?)
minlength    = int(floor(Nsteps/Npartitions)) # For sectioning the data
print('timestepsize:', timestepsize)

endlocation          = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp+str(psigma) + '/Nocut/'
filestext            = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])

# Text files
infilename_ds       = endlocation+'av_ds_'+filestext+'_nocut.txt'

# Plots
plotname             = endlocation+'par_ort_beginning_'+filestext+'_nocut.png'

infile_ds = open(infilename_ds, 'r')
firstline = infile_ds.readline()

#Time step; Time; <R^2>; <dx^2>; <dy^2>; <dz^2>; <dx^2+dy^2>
#  [0]       [1]   [2]     [3]    [4]      [5]       [6]
timesteps = []
times     = []
dR2       = []
dx2       = []
dy2       = []
dz2       = []
dpar2     = []

lines = infile_ds.readlines()

for line in lines:
    words = line.split()
    timesteps.append(int(words[0]))
    times.append(float(words[1]))
    dR2.append(float(words[2]))
    dx2.append(float(words[3]))
    dy2.append(float(words[4]))
    dz2.append(float(words[5]))
    dpar2.append(float(words[6]))
infile_ds.close()

timesteps = np.array(timesteps)
times     = np.array(times)
dR2       = np.array(dR2)
dx2       = np.array(dx2)
dy2       = np.array(dy2)
dz2       = np.array(dz2)
dpar2     = np.array(dpar2)

## To determine the range
# Interactive
plt.figure(figsize=(6,5))
plt.plot(timesteps, dR2, label=r'$<R^2>$')
plt.plot(timesteps, dx2, label=r'$<dx^2>$')
plt.plot(timesteps, dy2, label=r'$<dy^2>$')
plt.plot(timesteps, dz2, label=r'$<dz^2>$')
plt.plot(timesteps, dpar2, label=r'$<dx^2+dy^2>$')
plt.xlabel(r'Index (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %.2f nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.show()

# Saving
xmax_plot = 100
plt.figure(figsize=(6,5))
plt.plot(timesteps, dR2, label=r'$<R^2>$')
plt.plot(timesteps, dx2, label=r'$<dx^2>$')
plt.plot(timesteps, dy2, label=r'$<dy^2>$')
plt.plot(timesteps, dz2, label=r'$<dz^2>$')
plt.plot(timesteps, dpar2, label=r'$<dx^2+dy^2>$')
plt.xlabel(r'Index (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %.2f nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.axis([0, xmax_plot, 0, max(dR2[0:xmax_plot])])
plt.savefig(plotname)


