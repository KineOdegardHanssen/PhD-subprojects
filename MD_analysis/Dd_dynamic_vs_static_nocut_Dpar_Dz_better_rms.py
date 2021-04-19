import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob

def rmsd(x,y):
    Nx = len(x)
    Ny = len(y)
    if Nx!=Ny:
        print('WARNING! Nx!=Ny. Could not calculate rmsd value')
        return 'WARNING! Nx!=Ny. Could not calculate rmsd value'
    delta = 0
    for i in range(Nx):
        delta += (x[i]-y[i])*(x[i]-y[i])
    delta = np.sqrt(delta/(Nx-1))
    return delta

Nintervals = 10

# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
psigma   = 1 # For instance 
pmass    = 1 # I don't use this anymore, but it got stuck in the file names. Have constant mass density now.
damp     = 10
# Input booleans for file selection:
bulkdiffusion = False
substrate     = False
moresigmas    = False
big           = False
bulk_cut      = False
confignrs     = np.arange(1,1001)
filestext     = '_seed'+str(confignrs[0])+'to'+str(confignrs[-1])

endlocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/Nocut/'
if moresigmas==True:
    endlocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Moresigmas/Nocut/'

basepath_base      = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/'
endlocation_static = basepath_base+'D_vs_d/Nocut/'

bulklocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/'

bulkfilename  = bulklocation + 'diffusion_bulk'+filestext
if bulk_cut==True:
    bulkfilename = bulkfilename +'_cut.txt'
else:
    bulkfilename = bulkfilename +'_uncut.txt'

## Files to read
brushfilename_dyn  = endlocation + 'D_vs_d_better_rms_Nintervals%i.txt' % Nintervals
brushfilename_stat = endlocation_static + 'D_vs_d_static_better_rms_Nintervals%i.txt' % Nintervals
## Files to write to
if big==False:
    plotname     = endlocation_static+'Dd_dyn_vs_stat_noDR_better_rms_Nintervals%i.png' % Nintervals
    plotname_cut = endlocation_static+'Dd_dyn_vs_stat_cut_noDR_better_rms_Nintervals%i.png' % Nintervals
    plotname_twoinone = endlocation_static+'Dd_dyn_vs_stat_noDR_twoinone_better_rms_Nintervals%i.png' % Nintervals
else:
    plotname     = endlocation_static+'Dd_dyn_vs_stat_big_noDR_better_rms_Nintervals%i.png' % Nintervals
    plotname_cut = endlocation_static+'Dd_dyn_vs_stat_cut_big_noDR_better_rms_Nintervals%i.png' % Nintervals
    plotname_twoinone = endlocation_static+'Dd_dyn_vs_stat_big_noDR_twoinone_better_rms_Nintervals%i.png' % Nintervals

# Dynamic sims:
brushfile_dyn = open(brushfilename_dyn, 'r')

lines = brushfile_dyn.readlines()
N_dyn = len(lines)-1

# ds
spacings_dyn = np.zeros(N_dyn)
# Ds
DRs_dyn = np.zeros(N_dyn)
Dxs_dyn = np.zeros(N_dyn)
Dys_dyn = np.zeros(N_dyn)
Dzs_dyn = np.zeros(N_dyn)
Dparallel_dyn = np.zeros(N_dyn)
# Ds, stdv
DRs_stdv_dyn = np.zeros(N_dyn)
Dxs_stdv_dyn = np.zeros(N_dyn)
Dys_stdv_dyn = np.zeros(N_dyn)
Dzs_stdv_dyn = np.zeros(N_dyn)
Dparallel_stdv_dyn = np.zeros(N_dyn)

for i in range(1,N_dyn+1):
    words = lines[i].split()
    j = i-1
    
    spacings_dyn[j] = float(words[0])
    DRs_dyn[j] = float(words[1])
    Dzs_dyn[j] = float(words[3])
    Dparallel_dyn[j] = float(words[5])
    # Ds, stdv
    DRs_stdv_dyn[j] = float(words[2])
    Dzs_stdv_dyn[j] = float(words[4])
    Dparallel_stdv_dyn[j] = float(words[6])
    
brushfile_dyn.close()

##########

# Static sims:
brushfile_stat = open(brushfilename_stat, 'r')

lines = brushfile_stat.readlines()
N_stat = len(lines)-1

# ds
spacings_stat = np.zeros(N_stat)
# Ds
DRs_stat = np.zeros(N_stat)
Dxs_stat = np.zeros(N_stat)
Dys_stat = np.zeros(N_stat)
Dzs_stat = np.zeros(N_stat)
Dparallel_stat = np.zeros(N_stat)
# Ds, stdv
DRs_stdv_stat = np.zeros(N_stat)
Dxs_stdv_stat = np.zeros(N_stat)
Dys_stdv_stat = np.zeros(N_stat)
Dzs_stdv_stat = np.zeros(N_stat)
Dparallel_stdv_stat = np.zeros(N_stat)

for i in range(1,N_stat+1):
    words = lines[i].split()
    j = i-1
    
    spacings_stat[j] = float(words[0])
    DRs_stat[j] = float(words[1])
    Dzs_stat[j] = float(words[3])
    Dparallel_stat[j] = float(words[5])
    # Ds, stdv
    DRs_stdv_stat[j] = float(words[2])
    Dzs_stdv_stat[j] = float(words[4])
    Dparallel_stdv_stat[j] = float(words[6])
    
brushfile_stat.close()

###
#Bulk:

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

# Divide by bulk:
for i in range(N_stat):
    DRnew = DRs_stat[i]/DRs_bulk
    Dznew = Dzs_stat[i]/DRs_bulk
    Dparnew = Dparallel_stat[i]/DRs_bulk
    # Ds, stdv
    DRs_stdv_stat[i] = abs(DRnew)*np.sqrt((DRs_stdv_stat[i]/DRs_stat[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dzs_stdv_stat[i] = abs(Dznew)*np.sqrt((Dzs_stdv_stat[i]/Dzs_stat[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dparallel_stdv_stat[i] = abs(Dparnew)*np.sqrt((Dparallel_stdv_stat[i]/Dparallel_stat[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    # Ds
    DRs_stat[i] = DRnew
    Dzs_stat[i] = Dznew
    Dparallel_stat[i] = Dparnew

for i in range(N_dyn):
    DRnew = DRs_dyn[i]/DRs_bulk
    Dznew = Dzs_dyn[i]/DRs_bulk
    Dparnew = Dparallel_dyn[i]/DRs_bulk
    # Ds, stdv
    DRs_stdv_dyn[i] = abs(DRnew)*np.sqrt((DRs_stdv_dyn[i]/DRs_dyn[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dzs_stdv_dyn[i] = abs(Dznew)*np.sqrt((Dzs_stdv_dyn[i]/Dzs_dyn[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dparallel_stdv_dyn[i] = abs(Dparnew)*np.sqrt((Dparallel_stdv_dyn[i]/Dparallel_dyn[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    # Ds
    DRs_dyn[i] = DRnew
    Dzs_dyn[i] = Dznew
    Dparallel_dyn[i] =  Dparnew


if big==False:
    plt.figure(figsize=(8,5))
    ax = plt.subplot(111)
    ax.errorbar(spacings_dyn, Dzs_dyn, yerr=Dzs_stdv_dyn, color='b', capsize=2, label=r'$D_\perp$, dyn.')
    ax.errorbar(spacings_dyn, Dparallel_dyn, yerr=Dparallel_stdv_dyn, color='g', capsize=2, label=r'$D_\parallel$, dyn.')
    ax.errorbar(spacings_stat, Dzs_stat, yerr=Dzs_stdv_stat, color='c', capsize=2, label=r'$D_\perp$, stat.')
    ax.errorbar(spacings_stat, Dparallel_stat, yerr=Dparallel_stdv_stat, color='limegreen', capsize=2, label=r'$D_\parallel$, stat.')
    if moresigmas==True:
        plt.xlabel(r'$d/\sigma_b$')
    else:
        plt.xlabel(r'$d$')
    plt.ylabel(r'Diffusion constant $D$')
    if moresigmas==True:
        plt.title('Diffusion constant $D$ vs $d/\sigma_b$ for dynamic and static brushes', y=1.03)
    else:
        plt.title('Diffusion constant $D$ vs $d$ for dynamic and static brushes', y=1.03)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
else:
    plt.figure(figsize=(16,10))
    plt.rc('xtick', labelsize=20) 
    plt.rc('ytick', labelsize=20) 
    ax = plt.subplot(111)
    ax.errorbar(spacings_dyn, Dzs_dyn, yerr=Dzs_stdv_dyn, color='b', linewidth=7.0, capsize=2, label=r'$D_\perp$, dyn.')
    ax.errorbar(spacings_dyn, Dparallel_dyn, yerr=Dparallel_stdv_dyn, color='g', linewidth=7.0, capsize=2, label=r'$D_\parallel$, dyn.')
    ax.errorbar(spacings_stat, Dzs_stat, yerr=Dzs_stdv_stat, color='c', linewidth=7.0, capsize=2, label=r'$D_\perp$, stat.')
    ax.errorbar(spacings_stat, Dparallel_stat, yerr=Dparallel_stdv_stat, color='limegreen', linewidth=7.0, capsize=2, label=r'$D_\parallel$, stat.')
    if moresigmas==True:
        plt.xlabel(r'$d/\sigma_b$', fontsize=20)
    else:
        plt.xlabel(r'$d$', fontsize=20)
    plt.ylabel(r'Diffusion constant $D$', fontsize=20)
    if moresigmas==True:
        plt.title('Diffusion constant $D$ vs $d/\sigma_b$ for dynamic and static brushes', fontsize=28, y=1.03)
    else:
        plt.title('Diffusion constant $D$ vs $d$ for dynamic and static brushes', fontsize=28, y=1.03)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., prop={'size': 20})
plt.savefig(plotname)

if big==False:
    plt.figure(figsize=(6.4,5))
    ax = plt.subplot(111)
    ax.errorbar(spacings_dyn, Dzs_dyn, yerr=Dzs_stdv_dyn, color='b', capsize=2, label=r'$D_\perp$, dyn.')
    ax.errorbar(spacings_dyn, Dparallel_dyn, yerr=Dparallel_stdv_dyn, color='g', capsize=2, label=r'$D_\parallel$, dyn.')
    ax.errorbar(spacings_stat, Dzs_stat, yerr=Dzs_stdv_stat, color='c', capsize=2, label=r'$D_\perp$, stat.')
    ax.errorbar(spacings_stat, Dparallel_stat, yerr=Dparallel_stdv_stat, color='limegreen', capsize=2, label=r'$D_\parallel$, stat.')
    if moresigmas==True:
        plt.xlabel(r'$d/\sigma_b$')
    else:
        plt.xlabel(r'$d$')
    plt.ylabel(r'Diffusion constant $D$')
    if moresigmas==True:
        plt.title('Diffusion constant $D$ vs $d/\sigma_b$ for dynamic and static brushes',  y=1.03)
    else:
        plt.title('Diffusion constant $D$ vs $d$ for dynamic and static brushes', y=1.03)
    ax.axis([0,10,0,0.7]) # 6e-7 before we divided by Dbulk
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
else:
    plt.figure(figsize=(12.8,10))
    ax = plt.subplot(111)
    ax.errorbar(spacings_dyn, Dzs_dyn, yerr=Dzs_stdv_dyn, color='b', linewidth=7.0, capsize=2, label=r'$D_\perp$, dyn.')
    ax.errorbar(spacings_dyn, Dparallel_dyn, yerr=Dparallel_stdv_dyn, color='g', linewidth=7.0, capsize=2, label=r'$D_\parallel$, dyn.')
    ax.errorbar(spacings_stat, Dzs_stat, yerr=Dzs_stdv_stat, color='c', linewidth=7.0, capsize=2, label=r'$D_\perp$, stat.')
    ax.errorbar(spacings_stat, Dparallel_stat, yerr=Dparallel_stdv_stat, color='limegreen', linewidth=7.0, capsize=2, label=r'$D_\parallel$, stat.')
    if moresigmas==True:
        plt.xlabel(r'$d/\sigma_b$', fontsize=20)
    else:
        plt.xlabel(r'$d$', fontsize=20)
    plt.ylabel(r'Diffusion constant $D$', fontsize=20)
    if moresigmas==True:
        plt.title('Diffusion constant $D$ vs $d/\sigma_b$ for dynamic and static brushes', fontsize=28,  y=1.03)
    else:
        plt.title('Diffusion constant $D$ vs $d$ for dynamic and static brushes', fontsize=28, y=1.03)
    ax.axis([0,10,0,0.7]) # 6e-7 before we divided by Dbulk
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.savefig(plotname_cut)

################## twoinone ####################################
if big==False:
    print("Plottin it")
    #plt.figure(figsize=(16,5))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
    if moresigmas==True:
        fig.suptitle('Diffusion constant $D/D_{bulk}$ vs $d/\sigma_b$ for dynamic and static brushes')#, y=1.03)
    else:
        fig.suptitle('Diffusion constant $D/D_{bulk}$ vs $d$ for dynamic and static brushes')#, y=1.03)
    ax1.errorbar(spacings_dyn, Dzs_dyn, yerr=Dzs_stdv_dyn, color='b', capsize=2, label=r'$D_\perp$, dyn.')
    ax1.errorbar(spacings_dyn, Dparallel_dyn, yerr=Dparallel_stdv_dyn, color='g', capsize=2, label=r'$D_\parallel$, dyn.')
    ax1.errorbar(spacings_stat, Dzs_stat, yerr=Dzs_stdv_stat, color='c', capsize=2, label=r'$D_\perp$, stat.')
    ax1.errorbar(spacings_stat, Dparallel_stat, yerr=Dparallel_stdv_stat, color='limegreen', capsize=2, label=r'$D_\parallel$, stat.')
    if moresigmas==True:
        ax1.set(xlabel=r'$d/\sigma_b$', ylabel='Diffusion constant $D/D_{bulk}$')
    else:
        ax1.set(xlabel=r'$d$', ylabel='Diffusion constant $D/D_{bulk}$')
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax1.legend(loc="lower right")
    
    ax2.errorbar(spacings_dyn, Dzs_dyn, yerr=Dzs_stdv_dyn, color='b', capsize=2, fmt='.', label=r'$D_\perp$, dyn.')
    ax2.errorbar(spacings_dyn, Dparallel_dyn, yerr=Dparallel_stdv_dyn, color='g', capsize=2, fmt='.', label=r'$D_\parallel$, dyn.')
    ax2.errorbar(spacings_stat, Dzs_stat, yerr=Dzs_stdv_stat, color='c', capsize=2, fmt='.', label=r'$D_\perp$, stat.')
    ax2.errorbar(spacings_stat, Dparallel_stat, yerr=Dparallel_stdv_stat, color='limegreen', capsize=2, fmt='.', label=r'$D_\parallel$, stat.')
    if moresigmas==True:
        ax2.set(xlabel=r'$d/\sigma_b$', ylabel=r'Diffusion constant $D/D_{bulk}$')
    else:
        ax2.set(xlabel=r'$d$', ylabel=r'Diffusion constant $D/D_{bulk}$')
    ax2.axis([0,10,0,0.7]) # 6e-7 before we divided by Dbulk
    #fig.tight_layout()
    plt.show()
    fig.savefig(plotname_twoinone)
else:
    print("Plottin it")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
    if moresigmas==True:
       fig.suptitle('Diffusion constant $D$ vs $d/\sigma_b$ for dynamic and static brushes', fontsize=28, y=1.03)
    else:
        fig.suptitle('Diffusion constant $D$ vs $d$ for dynamic and static brushes', fontsize=28, y=1.03)
    ax1.errorbar(spacings_dyn, Dzs_dyn, yerr=Dzs_stdv_dyn, color='b', linewidth=7.0, capsize=2, label=r'$D_\perp$, dyn.')
    ax1.errorbar(spacings_dyn, Dparallel_dyn, yerr=Dparallel_stdv_dyn, color='g', linewidth=7.0, capsize=2, label=r'$D_\parallel$, dyn.')
    ax1.errorbar(spacings_stat, Dzs_stat, yerr=Dzs_stdv_stat, color='c', linewidth=7.0, capsize=2, label=r'$D_\perp$, stat.')
    ax1.errorbar(spacings_stat, Dparallel_stat, yerr=Dparallel_stdv_stat, color='limegreen', linewidth=7.0, capsize=2, label=r'$D_\parallel$, stat.')
    if moresigmas==True:
        ax1.set(xlabel=r'$d/\sigma_b$', ylabel=r'Diffusion constant $D/D_{bulk}$', fontsize=20)
    else:
        ax1.set(xlabel=r'$d$', ylabel=r'Diffusion constant $D/D_{bulk}$', fontsize=20)
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    ax2.errorbar(spacings_dyn, Dzs_dyn, yerr=Dzs_stdv_dyn, color='b', linewidth=7.0, capsize=2, fmt='.', label=r'$D_\perp$, dyn.')
    ax2.errorbar(spacings_dyn, Dparallel_dyn, yerr=Dparallel_stdv_dyn, color='g', linewidth=7.0, capsize=2, fmt='.', label=r'$D_\parallel$, dyn.')
    ax2.errorbar(spacings_stat, Dzs_stat, yerr=Dzs_stdv_stat, color='c', linewidth=7.0, capsize=2, fmt='.', label=r'$D_\perp$, stat.')
    ax2.errorbar(spacings_stat, Dparallel_stat, yerr=Dparallel_stdv_stat, color='limegreen', linewidth=7.0, capsize=2, fmt='.', label=r'$D_\parallel$, stat.')
    if moresigmas==True:
        ax2.set(xlabel=r'$d/\sigma_b$', ylabel=r'Diffusion constant $D/D_{bulk}$', fontsize=20)
    else:
        ax2.set(xlabel=r'$d$', ylabel=r'Diffusion constant $D/D_{bulk}$', fontsize=20)
    ax2.axis([0,10,0,0.7]) # 6e-7 before we divided by Dbulk
    ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #fig.tight_layout()
    plt.show()
    fig.savefig(plotname_twoinone)

print("Done.")