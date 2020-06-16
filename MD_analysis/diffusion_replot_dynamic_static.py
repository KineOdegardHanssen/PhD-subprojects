from pylab import *
import matplotlib.pyplot as plt                     # To plot
import numpy as np
import math

showplot = False

damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
spacing = 100
psigma  = 1
print('spacing:', spacing)
print('psigma:', psigma)
# I need to set the file name in an easier way, but for now I just use this:  ## Might want to add a loop too, if I have more files...

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T            = 3
Nsteps       = 2001 # 20001 before
unitlength   = 1e-9
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime#*writeevery
Npartitions  = 5 # For extracting more walks from one file (but is it really such a random walk here...?)
minlength    = int(floor(Nsteps/Npartitions)) # For sectioning the data
print('timestepsize:', timestepsize)
confignrs_dynamic = np.arange(1,1001)
confignrs_static  = np.arange(1,101)
placements        = np.arange(1,11)

endlocation_static  = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/Spacing' + str(spacing) + '/Radius' + str(psigma) + '/Nocut/'
endlocation_dynamic = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp+str(psigma) + '/Nocut/'
endlocation_out     = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/Dynamic_vs_static/'
filestext_static    = 'config'+str(confignrs_static[0])+'to'+str(confignrs_static[-1])+'_placements'+str(placements[0])+'to'+str(placements[-1])
filestext_dynamic   = 'config'+str(confignrs_dynamic[0])+'to'+str(confignrs_dynamic[-1])


# Text files
infilename_ds_static   = endlocation_static +'av_ds_'+filestext_static+'.txt'
infilename_ds_dynamic  = endlocation_dynamic +'av_ds_'+filestext_dynamic+'_nocut.txt'

# Plots
plotname               = endlocation_out+'d'+str(spacing)+'_par_ort_beginning_dynamic_and_static.png'

# Read static:
infile_ds_static = open(infilename_ds_static, 'r')
firstline = infile_ds_static.readline()

#Time step; Time; <R^2>; <dx^2>; <dy^2>; <dz^2>; <dx^2+dy^2>
#  [0]       [1]   [2]     [3]    [4]      [5]       [6]
timesteps = []
times     = []
dR2       = []
dx2       = []
dy2       = []
dz2       = []
dpar2     = []

lines = infile_ds_static.readlines()

for line in lines:
    words = line.split()
    timesteps.append(int(words[0]))
    times.append(float(words[1]))
    dR2.append(float(words[2]))
    dx2.append(float(words[3]))
    dy2.append(float(words[4]))
    dz2.append(float(words[5]))
    dpar2.append(float(words[6]))
infile_ds_static.close()

timesteps_static = np.array(timesteps)
times_static     = np.array(times)
dR2_static       = np.array(dR2)
dx2_static       = np.array(dx2)
dy2_static       = np.array(dy2)
dz2_static       = np.array(dz2)
dpar2_static     = np.array(dpar2)

# Read dynamic:
infile_ds_dynamic = open(infilename_ds_dynamic, 'r')
firstline = infile_ds_dynamic.readline()

#Time step; Time; <R^2>; <dx^2>; <dy^2>; <dz^2>; <dx^2+dy^2>
#  [0]       [1]   [2]     [3]    [4]      [5]       [6]
timesteps = []
times     = []
dR2       = []
dx2       = []
dy2       = []
dz2       = []
dpar2     = []

lines = infile_ds_dynamic.readlines()

for line in lines:
    words = line.split()
    timesteps.append(int(words[0]))
    times.append(float(words[1]))
    dR2.append(float(words[2]))
    dx2.append(float(words[3]))
    dy2.append(float(words[4]))
    dz2.append(float(words[5]))
    dpar2.append(float(words[6]))
infile_ds_dynamic.close()

timesteps_dynamic = np.array(timesteps)
times_dynamic     = np.array(times)
dR2_dynamic       = np.array(dR2)
dx2_dynamic       = np.array(dx2)
dy2_dynamic       = np.array(dy2)
dz2_dynamic       = np.array(dz2)
dpar2_dynamic     = np.array(dpar2)

## To determine the range
# Interactive
if showplot==True:
    plt.figure(figsize=(8,5))
    ax = plt.subplot(111)
    ax.plot(timesteps_dynamic, dR2_dynamic, color='r', label=r'$<R^2>$')
    ax.plot(timesteps_dynamic, dz2_dynamic, color='b', label=r'$<dz^2>$')
    ax.plot(timesteps_dynamic, dpar2_dynamic, color='g', label=r'$<dx^2+dy^2>$')
    ax.plot(timesteps_static, dR2_static, color='lightcoral', label=r'$<R^2>$')
    ax.plot(timesteps_static, dz2_static, color='c', label=r'$<dz^2>$')
    ax.plot(timesteps_static, dpar2_static, color='limegreen', label=r'$<dx^2+dy^2>$')
    plt.xlabel(r'Index (s)')
    plt.ylabel(r'Distance$^2$ [in unit length]')
    plt.title(r'RMSD in brush, d = %.2f nm, $\sigma_b=%.2f$' % (spacing,psigma))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.show()

# Saving
xmax_plot = 100
plt.figure(figsize=(8,5))
ax = plt.subplot(111)
ax.plot(timesteps_dynamic, dR2_dynamic, color='r', label=r'$<R^2>$')
ax.plot(timesteps_dynamic, dz2_dynamic, color='b', label=r'$<dz^2>$')
ax.plot(timesteps_dynamic, dpar2_dynamic, color='g', label=r'$<dx^2+dy^2>$')
ax.plot(timesteps_static, dR2_static, color='lightcoral', label=r'$<R^2>$')
ax.plot(timesteps_static, dz2_static, color='c', label=r'$<dz^2>$')
ax.plot(timesteps_static, dpar2_static, color='limegreen', label=r'$<dx^2+dy^2>$')
plt.xlabel(r'Index (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'MSD in brush, d = %.2f nm, $\sigma_b=%.2f$' % (spacing,psigma))
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname)

print('d:', spacing, '; mean(dR2_static)/mean(dR2_dynamic):', (np.mean(dR2_static)/np.mean(dR2_dynamic)))
