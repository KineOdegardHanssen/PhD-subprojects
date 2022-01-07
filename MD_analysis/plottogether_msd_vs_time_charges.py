import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob


dzs_all   = []
dpar_all  = []
times_all = []

spacing = 10

# Fixed parameters
psigma   = 1 # For instance 
damp     = 10
# Input booleans for file selection:
bulkdiffusion = False
substrate     = False
moresigmas    = False
big           = False
bulk_cut      = False
confignrs     = np.arange(1,1001)
filestext     = '_seed'+str(confignrs[0])+'to'+str(confignrs[-1])

charges  = [-2,-1,-0.25,0,1,2]
colors   = ['r','b','c','g','k','darkviolet']
markers  = ['-*','-X','-o',5,1] # Can I use this?
Ncharges = len(charges)

endlocation = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/Nocut/'

#Bulk:
bulklocation = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/'
bulkfilename  = bulklocation + 'diffusion_bulk'+filestext
if bulk_cut==True:
    bulkfilename = bulkfilename +'_cut.txt'
else:
    bulkfilename = bulkfilename +'_uncut.txt'

bulkfile  = open(bulkfilename, 'r')
# D_R2  sigmaD_R2 b_R2 sigmab_R2; D_z2  sigmaD_z2  b_z2  sigmaD_z2; D_par2 sigmaD_par2  b_par2  sigmab_par2
bulklines = bulkfile.readlines()
bulkline  = bulklines[1]
words     = bulkline.split()

# Ds
DRs_bulk = float(words[0])
Dzs_bulk = float(words[4])
Dparallel_bulk = float(words[8])

# Ds, stdv
DRs_stdv_bulk = float(words[1])
Dzs_stdv_bulk = float(words[5])
Dparallel_stdv_bulk = float(words[9])

bulkfile.close()

filestext           = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])

for charge in charges:
    if charge!=0:
        endlocation_cut     = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp+str(psigma) + '/Charge'+str(charge)+'/'
        endlocation_nocut   = endlocation_cut + 'Nocut/'
    else:
        endlocation_cut     = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp+str(psigma) + '/'
        endlocation_nocut   = endlocation_cut + 'Nocut/'
        
    # Text files
    try:
        infilename_ds_nocut = endlocation_nocut+'av_ds_'+filestext+'_nocut.txt'
        infile_ds_nocut = open(infilename_ds_nocut,'r')
    except:
        infilename_ds_nocut = endlocation_nocut+'av_ds_'+filestext+'_nocut_long.txt'
        infile_ds_nocut = open(infilename_ds_nocut,'r')
    
    
    firstline = infile_ds_nocut.readline()
    lines = infile_ds_nocut.readlines()
    N_dyn = len(lines)-1
    
    timesteps_nocut = []
    times_nocut     = []
    dR2_nocut       = []
    dx2_nocut       = []
    dy2_nocut       = []
    dz2_nocut       = []
    dpar2_nocut     = []
    
    
    for line in lines:
        words = line.split()
        timesteps_nocut.append(int(words[0]))
        times_nocut.append(float(words[1]))
        dR2_nocut.append(float(words[2]))
        dx2_nocut.append(float(words[3]))
        dy2_nocut.append(float(words[4]))
        dz2_nocut.append(float(words[5]))
        dpar2_nocut.append(float(words[6]))
    infile_ds_nocut.close()
    
    dR2_nocut = np.array(dR2_nocut)
    dx2_nocut = np.array(dx2_nocut)
    dy2_nocut = np.array(dy2_nocut)
    dz2_nocut = np.array(dz2_nocut)
    dpar2_nocut = np.array(dpar2_nocut)
    
    dzs_all.append(dz2_nocut)
    dpar_all.append(dpar2_nocut)
    times_all.append(times_nocut)

#print('dzs_all:',dzs_all)

plt.figure(figsize=(14,5))
ax = plt.subplot(111)
for i in range(Ncharges):
    ax.plot(times_all[i], dzs_all[i], color=colors[i], label=r'$q=$%.2f' % charges[i])
plt.xlabel(r'time step', fontsize=12)
plt.ylabel(r'<$dz^2$>', fontsize=12)
#plt.ylabel(r'<$dz^2$>, d=%.1f' % spacing, fontsize=12)
plt.title(r'$d$= %.1f nm' % spacing)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.show()

plt.figure(figsize=(14,5))
ax = plt.subplot(111)
for i in range(Ncharges):
    ax.plot(times_all[i], dpar_all[i], color=colors[i], label=r'$q=$%.2f' % charges[i])
plt.xlabel(r'time step', fontsize=12)
plt.ylabel(r'<$d\parallel^2$>', fontsize=12)
#plt.ylabel(r'<$d\parallel^2$>, d=%.1f' % spacing, fontsize=12)
plt.title(r'$d$ = %.1f nm' % spacing)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.show()
print('dzs_all[0]-dzs_all[1]:',dzs_all[0]-dzs_all[1])
print('dzs_all[1]-dzs_all[2]:',dzs_all[1]-dzs_all[2])
print('dzs_all[3]-dzs_all[2]:',dzs_all[3]-dzs_all[2])
print('dzs_all[0]-dzs_all[4]:',dzs_all[0]-dzs_all[4])
print('dzs_all[1]-dzs_all[4]:',dzs_all[1]-dzs_all[4])
print('dzs_all[2]-dzs_all[4]:',dzs_all[2]-dzs_all[4])
print('dzs_all[3]-dzs_all[4]:',dzs_all[3]-dzs_all[4])


print('np.sum(dzs_all[0]-dzs_all[1]):',np.sum(dzs_all[0]-dzs_all[1]))
print('np.sum(dzs_all[1]-dzs_all[2]):',np.sum(dzs_all[1]-dzs_all[2]))
print('np.sum(dzs_all[3]-dzs_all[2]):',np.sum(dzs_all[3]-dzs_all[2]))
print('np.sum(dzs_all[1]-dzs_all[4]):',np.sum(dzs_all[3]-dzs_all[4]))
print('np.sum(dzs_all[2]-dzs_all[4]):',np.sum(dzs_all[2]-dzs_all[4]))
#plt.savefig(plotname_d_dyn)
