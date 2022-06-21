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

vashishtalocation     = '/home/kine/Documents/P4_Crosslinking_Vashishta/Spacing'+str(spacing)+'/Sigma_free_'  + str(radius) + '/Charge%i' % charge
filestext             = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])+'_nocrosslinkers'
# Text files
outfilename_vashishta = vashishtalocation+'av_ds_'+filestext+'.txt'

LJDebyelocation     = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp+str(psigma) + '/Charge'+str(charge)+'/'
filestext           = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])
# Text files
outfilename_LJDebye = LJDebyelocation+'av_ds_'+filestext+'.txt'



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
    times_cut.append(float(words[1]))
    dR2_cut.append(float(words[2]))
    dx2_cut.append(float(words[3]))
    dy2_cut.append(float(words[4]))
    dz2_cut.append(float(words[5]))
    dpar2_cut.append(float(words[6]))
infile_ds_cut.close()

timesteps_cut = np.array(timesteps_cut)
times_cut     = np.array(times_cut)
dR2_cut       = np.array(dR2_cut)
dx2_cut       = np.array(dx2_cut)
dy2_cut       = np.array(dy2_cut)
dz2_cut       = np.array(dz2_cut)
dpar2_cut     = np.array(dpar2_cut)
