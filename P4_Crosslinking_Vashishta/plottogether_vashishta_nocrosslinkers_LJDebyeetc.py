import matplotlib.pyplot as plt                     # To plot
import numpy as np
import random
import math
import time
import os
import glob
import copy

spacing   = 3
psigma    = 1
radius    = psigma
charge    = -1
confignrs = np.arange(1,1001)
damp = 10

vashishtalocation     = '/home/kine/Documents/P4_Crosslinking_Vashishta/Spacing'+str(spacing)+'/Sigma_free_'  + str(radius) + '/Charge%i/' % charge
filestext             = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])+'_nocrosslinkers'
# Text files
infilename_vashishta = vashishtalocation+'av_ds_'+filestext+'.txt'

LJDebyelocation     = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp+str(psigma) + '/Charge'+str(charge)+'/'
filestext           = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])
# Text files
infilename_LJDebye = LJDebyelocation+'av_ds_'+filestext+'.txt'

plotname = vashishtalocation+'compare_Vashishta_LJDebye.png'


### LJDebye
infile_LJDebye = open(infilename_LJDebye, 'r')
firstline     = infile_LJDebye.readline()

#Time step; Time; <R^2>; <dx^2>; <dy^2>; <dz^2>; <dx^2+dy^2>
#  [0]       [1]   [2]     [3]    [4]      [5]       [6]
timesteps_LJDebye = [] # The times should be the same, but being careful never hurts.
times_LJDebye     = []
dR2_LJDebye       = []
dx2_LJDebye       = []
dy2_LJDebye       = []
dz2_LJDebye       = []
dpar2_LJDebye     = []

lines = infile_LJDebye.readlines()

for line in lines:
    words = line.split()
    timesteps_LJDebye.append(int(words[0]))
    times_LJDebye.append(float(words[1]))
    dR2_LJDebye.append(float(words[2]))
    dx2_LJDebye.append(float(words[3]))
    dy2_LJDebye.append(float(words[4]))
    dz2_LJDebye.append(float(words[5]))
    dpar2_LJDebye.append(float(words[6]))
infile_LJDebye.close()

timesteps_LJDebye = np.array(timesteps_LJDebye)
times_LJDebye     = np.array(times_LJDebye)
dR2_LJDebye       = np.array(dR2_LJDebye)
dx2_LJDebye       = np.array(dx2_LJDebye)
dy2_LJDebye       = np.array(dy2_LJDebye)
dz2_LJDebye       = np.array(dz2_LJDebye)
dpar2_LJDebye     = np.array(dpar2_LJDebye)

### Vashishta
infile_vashishta = open(infilename_vashishta, 'r')
firstline     = infile_vashishta.readline()

#Time step; Time; <R^2>; <dx^2>; <dy^2>; <dz^2>; <dx^2+dy^2>
#  [0]       [1]   [2]     [3]    [4]      [5]       [6]
timesteps_vashishta = [] # The times should be the same, but being careful never hurts.
times_vashishta     = []
dR2_vashishta       = []
dx2_vashishta       = []
dy2_vashishta       = []
dz2_vashishta       = []
dpar2_vashishta     = []

lines = infile_vashishta.readlines()

for line in lines:
    words = line.split()
    timesteps_vashishta.append(int(words[0]))
    times_vashishta.append(float(words[1]))
    dR2_vashishta.append(float(words[2]))
    dx2_vashishta.append(float(words[3]))
    dy2_vashishta.append(float(words[4]))
    dz2_vashishta.append(float(words[5]))
    dpar2_vashishta.append(float(words[6]))
infile_vashishta.close()

timesteps_vashishta = np.array(timesteps_vashishta)
times_vashishta     = np.array(times_vashishta)
dR2_vashishta       = np.array(dR2_vashishta)
dx2_vashishta       = np.array(dx2_vashishta)
dy2_vashishta       = np.array(dy2_vashishta)
dz2_vashishta       = np.array(dz2_vashishta)
dpar2_vashishta     = np.array(dpar2_vashishta)


plt.figure(figsize=(6,5))
plt.plot(timesteps_LJDebye, dR2_LJDebye, label='dR2, LJDebye')
plt.plot(timesteps_LJDebye, dz2_LJDebye, label='dz2, LJDebye')
plt.plot(timesteps_LJDebye, dpar2_LJDebye, label='dpar2, LJDebye')
#
plt.plot(timesteps_vashishta, dR2_vashishta, '--', label='dR2, Vashishta')
plt.plot(timesteps_vashishta, dz2_vashishta, '--', label='dz2, Vashishta')
plt.plot(timesteps_vashishta, dpar2_vashishta, '--', label='dpar2, Vashishta')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname)
plt.show()
