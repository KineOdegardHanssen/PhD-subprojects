from pylab import *
import matplotlib.pyplot as plt                     # To plot
import numpy as np
import math

damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
spacing = 1
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

endlocation_cut     = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp+str(psigma) + '/'
endlocation_nocut   = endlocation_cut + 'Nocut/'
filestext           = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])

# Text files
infilename_ds_nocut = endlocation_nocut+'av_ds_'+filestext+'_nocut.txt'
infilename_ds_cut   = endlocation_cut+'av_ds_'+filestext+'.txt'

# Plots
plotname_all        = endlocation_nocut+'par_ort_all_'+filestext+'_nocutandcuttogether.png'
plotname_beginning  = endlocation_nocut+'par_ort_beginning_'+filestext+'_nocutandcuttogether.png'


## Read in files
# No cutting
infile_ds_nocut = open(infilename_ds_nocut, 'r')
firstline       = infile_ds_nocut.readline()

#Time step; Time; <R^2>; <dx^2>; <dy^2>; <dz^2>; <dx^2+dy^2>
#  [0]       [1]   [2]     [3]    [4]      [5]       [6]
timesteps   = []
times       = []
dR2_nocut   = []
dx2_nocut   = []
dy2_nocut   = []
dz2_nocut   = []
dpar2_nocut = []

lines = infile_ds_nocut.readlines()

for line in lines:
    words = line.split()
    timesteps.append(int(words[0]))
    times.append(float(words[1]))
    dR2_nocut.append(float(words[2]))
    dx2_nocut.append(float(words[3]))
    dy2_nocut.append(float(words[4]))
    dz2_nocut.append(float(words[5]))
    dpar2_nocut.append(float(words[6]))
infile_ds_nocut.close()

timesteps   = np.array(timesteps)
times       = np.array(times)
dR2_nocut   = np.array(dR2_nocut)
dx2_nocut   = np.array(dx2_nocut)
dy2_nocut   = np.array(dy2_nocut)
dz2_nocut   = np.array(dz2_nocut)
dpar2_nocut = np.array(dpar2_nocut)


############
# Cutting
infile_ds_cut = open(infilename_ds_cut, 'r')
firstline     = infile_ds_cut.readline()

#Time step; Time; <R^2>; <dx^2>; <dy^2>; <dz^2>; <dx^2+dy^2>
#  [0]       [1]   [2]     [3]    [4]      [5]       [6]
# Assuming the times are the same.
dR2_cut   = []
dx2_cut   = []
dy2_cut   = []
dz2_cut   = []
dpar2_cut = []

lines = infile_ds_cut.readlines()

for line in lines:
    words = line.split()
    dR2_cut.append(float(words[2]))
    dx2_cut.append(float(words[3]))
    dy2_cut.append(float(words[4]))
    dz2_cut.append(float(words[5]))
    dpar2_cut.append(float(words[6]))
infile_ds_cut.close()

dR2_cut   = np.array(dR2_cut)
dx2_cut   = np.array(dx2_cut)
dy2_cut   = np.array(dy2_cut)
dz2_cut   = np.array(dz2_cut)
dpar2_cut = np.array(dpar2_cut)


## To determine the range
# Interactive
plt.figure(figsize=(8,5))
ax = plt.subplot(111)
# Cut
ax.plot(timesteps, dR2_cut, label=r'$<R^2>$, cut')
#plt.plot(timesteps, dx2_cut, label=r'$<dx^2>$, cut')
#plt.plot(timesteps, dy2_cut, label=r'$<dy^2>$, cut')
ax.plot(timesteps, dz2_cut, label=r'$<dz^2>$, cut')
ax.plot(timesteps, dpar2_cut, label=r'$<dx^2+dy^2>$, cut')
# No cut
ax.plot(timesteps, dR2_nocut, label=r'$<R^2>$, no cut')
#plt.plot(timesteps, dx2_nocut, label=r'$<dx^2>$, no cut')
#plt.plot(timesteps, dy2_nocut, label=r'$<dy^2>$, no cut')
ax.plot(timesteps, dz2_nocut, label=r'$<dz^2>$, no cut')
ax.plot(timesteps, dpar2_nocut, label=r'$<dx^2+dy^2>$, no cut')
plt.xlabel(r'Index (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %.2f nm, $\sigma_b=%.2f$' % (spacing,psigma))
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.show()

# Saving
plt.figure(figsize=(8,5))
ax = plt.subplot(111)
# Cut
ax.plot(timesteps, dR2_cut, label=r'$<R^2>$, cut')
#plt.plot(timesteps, dx2_cut, label=r'$<dx^2>$, cut')
#plt.plot(timesteps, dy2_cut, label=r'$<dy^2>$, cut')
ax.plot(timesteps, dz2_cut, label=r'$<dz^2>$, cut')
ax.plot(timesteps, dpar2_cut, label=r'$<dx^2+dy^2>$, cut')
# No cut
ax.plot(timesteps, dR2_nocut, label=r'$<R^2>$, no cut')
#plt.plot(timesteps, dx2_nocut, label=r'$<dx^2>$, no cut')
#plt.plot(timesteps, dy2_nocut, label=r'$<dy^2>$, no cut')
ax.plot(timesteps, dz2_nocut, label=r'$<dz^2>$, no cut')
ax.plot(timesteps, dpar2_nocut, label=r'$<dx^2+dy^2>$, no cut')
plt.xlabel(r'Index (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %.2f nm, $\sigma_b=%.2f$' % (spacing,psigma))
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname_all)


xmax_plot = 100
plt.figure(figsize=(8,5))
ax = plt.subplot(111)
# Cut
ax.plot(timesteps, dR2_cut, label=r'$<R^2>$, cut')
#plt.plot(timesteps, dx2_cut, label=r'$<dx^2>$, cut')
#plt.plot(timesteps, dy2_cut, label=r'$<dy^2>$, cut')
ax.plot(timesteps, dz2_cut, label=r'$<dz^2>$, cut')
ax.plot(timesteps, dpar2_cut, label=r'$<dx^2+dy^2>$, cut')
# No cut
ax.plot(timesteps, dR2_nocut, label=r'$<R^2>$, no cut')
#plt.plot(timesteps, dx2_nocut, label=r'$<dx^2>$, no cut')
#plt.plot(timesteps, dy2_nocut, label=r'$<dy^2>$, no cut')
ax.plot(timesteps, dz2_nocut, label=r'$<dz^2>$, no cut')
ax.plot(timesteps, dpar2_nocut, label=r'$<dx^2+dy^2>$, no cut')
plt.xlabel(r'Index (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %.2f nm, $\sigma_b=%.2f$' % (spacing,psigma))
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) # ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., prop={'size': 20})
ax.axis([0, xmax_plot, 0, max(dR2_nocut[0:xmax_plot])])
plt.savefig(plotname_beginning)

print('Done')
