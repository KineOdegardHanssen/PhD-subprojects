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
old_bulk      = False
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
if old_bulk==False:
    bulkfilename  = bulklocation + 'diffusion_bulk'+filestext+'_new.txt'

## Files to read
brushfilename_dyn  = endlocation + 'D_vs_d_better_rms_Nestimates%i.txt' % Nintervals
brushfilename_stat = endlocation_static + 'D_vs_d_static_better_rms_Nestimates%i.txt' % Nintervals
## Files to write to
if big==False:
    plotname     = endlocation_static+'Dd_dyn_vs_stat_noDR_better_rms_Nestimates%i.png' % Nintervals
    plotname_cut = endlocation_static+'Dd_dyn_vs_stat_cut_noDR_better_rms_Nestimates%i.png' % Nintervals
    plotname_twoinone = endlocation_static+'Dd_dyn_vs_stat_noDR_twoinone_better_rms_Nestimates%i.png' % Nintervals
else:
    plotname     = endlocation_static+'Dd_dyn_vs_stat_big_noDR_better_rms_Nestimates%i.png' % Nintervals
    plotname_cut = endlocation_static+'Dd_dyn_vs_stat_cut_big_noDR_better_rms_Nestimates%i.png' % Nintervals
    plotname_twoinone = endlocation_static+'Dd_dyn_vs_stat_big_noDR_twoinone_better_rms_Nestimates%i.png' % Nintervals

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
if old_bulk==True: # Worse rms
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
else:
    bulklines = bulkfile.readlines()
    bulkline  = bulklines[1]
    words     = bulkline.split()
    
    # Ds
    DRs_bulk = float(words[1])
    Dzs_bulk = float(words[3])
    Dparallel_bulk = float(words[5])
    
    # Ds, stdv
    DRs_stdv_bulk = float(words[2])
    Dzs_stdv_bulk = float(words[4])
    Dparallel_stdv_bulk = float(words[6])
    
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

minforsmall = 0
maxforsmall = 0
i = 0
d = 1
while d<11:
    # Max vals:
    Dthismax_zdyn = Dzs_dyn[i]+Dzs_stdv_dyn[i]
    Dthismax_pdyn = Dparallel_dyn[i]+Dparallel_stdv_dyn[i]
    Dthismax_zstat = Dzs_stat[i]+Dzs_stdv_stat[i]
    Dthismax_pstat = Dparallel_stat[i]+Dparallel_stdv_stat[i]
    # Min vals:
    Dthismin_zdyn = Dzs_dyn[i]-Dzs_stdv_dyn[i]
    Dthismin_pdyn = Dparallel_dyn[i]-Dparallel_stdv_dyn[i]
    Dthismin_zstat = Dzs_stat[i]-Dzs_stdv_stat[i]
    Dthismin_pstat = Dparallel_stat[i]-Dparallel_stdv_stat[i]
    # Testing if larger:
    if Dthismax_zdyn>maxforsmall:
        maxforsmall=Dthismax_zdyn
    if Dthismax_pdyn>maxforsmall:
        maxforsmall=Dthismax_pdyn
    if Dthismax_zstat>maxforsmall:
        maxforsmall=Dthismax_zstat
    if Dthismax_pstat>maxforsmall:
        maxforsmall=Dthismax_pstat
    # Testing if smaller:
    if Dthismin_zdyn<minforsmall:
        minforsmall=Dthismin_zdyn
    if Dthismin_pdyn<minforsmall:
        minforsmall=Dthismin_pdyn
    if Dthismin_zstat<minforsmall:
        minforsmall=Dthismin_zstat
    if Dthismin_pstat<minforsmall:
        minforsmall=Dthismin_pstat
    i+=1
    d = spacings_dyn[i] 

if minforsmall<0:
    minforsmall*=1.2
else:
    minforsmall*=0.8
maxforsmall*=1.05
print('maxforsmall:',maxforsmall)

## Put matplotlib setting here?: # No. I do want italics in some equations.
##params = {'mathtext.default': 'regular' }          
##plt.rcParams.update(params)

if big==False:
    plt.figure(figsize=(8,5))
    ax = plt.subplot(111)    
    ax.plot(spacings_dyn, Dzs_dyn, '-o' , color='b',label=r'$D_\perp$, dyn.')
    ax.plot(spacings_dyn, Dparallel_dyn, '-*', color='g', label=r'$D_\parallel$, dyn.')
    ax.plot(spacings_stat, Dzs_stat, '-s', color='c', label=r'$D_\perp$, stat.')
    ax.plot(spacings_stat, Dparallel_stat, '-d', color='limegreen',label=r'$D_\parallel$, stat.')
    ax.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='b',alpha=0.2)
    ax.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
    ax.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='c', alpha=0.2)
    ax.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen',alpha=0.2)
    if moresigmas==True:
        plt.xlabel(r'$d/\sigma_b$')
    else:
        plt.xlabel(r'$d$ (nm)')
    plt.ylabel(r'$D/D_{\mathregular{bulk}}$')
    #if moresigmas==True:
    #    plt.title('Diffusion constant $D$ vs $d/\sigma_b$ for dynamic and static brushes', y=1.03)
    #else:
    #    plt.title('Diffusion constant $D$ vs $d$ for dynamic and static brushes', y=1.03)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
else:
    plt.figure(figsize=(16,10))
    plt.rc('xtick', labelsize=20) 
    plt.rc('ytick', labelsize=20) 
    ax = plt.subplot(111)
    ax.plot(spacings_dyn, Dzs_dyn, '-o' , color='b',label=r'$D_\perp$, dyn.')
    ax.plot(spacings_dyn, Dparallel_dyn, '-*', color='g', label=r'$D_\parallel$, dyn.')
    ax.plot(spacings_stat, Dzs_stat, '-s', color='c', label=r'$D_\perp$, stat.')
    ax.plot(spacings_stat, Dparallel_stat, '-d', color='limegreen',label=r'$D_\parallel$, stat.')
    ax.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='b',alpha=0.2)
    ax.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
    ax.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='c', alpha=0.2)
    ax.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen',alpha=0.2)
    if moresigmas==True:
        plt.xlabel(r'$d/\sigma_b$', fontsize=20)
    else:
        plt.xlabel(r'$d$ (nm)', fontsize=20)
    plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=20)
    #if moresigmas==True:
    #    plt.title('Diffusion constant $D$ vs $d/\sigma_b$ for dynamic and static brushes', fontsize=28, y=1.03)
    #else:
    #    plt.title('Diffusion constant $D$ vs $d$ for dynamic and static brushes', fontsize=28, y=1.03)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., prop={'size': 20})
plt.savefig(plotname)

if big==False:
    plt.figure(figsize=(6.4,5))
    ax = plt.subplot(111)    
    ax.plot(spacings_dyn, Dzs_dyn, '-o' , color='b',label=r'$D_\perp$, dyn.')
    ax.plot(spacings_dyn, Dparallel_dyn, '-*', color='g', label=r'$D_\parallel$, dyn.')
    ax.plot(spacings_stat, Dzs_stat, '-s', color='c', label=r'$D_\perp$, stat.')
    ax.plot(spacings_stat, Dparallel_stat, '-d', color='limegreen',label=r'$D_\parallel$, stat.')
    ax.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='b',alpha=0.2)
    ax.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
    ax.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='c', alpha=0.2)
    ax.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen',alpha=0.2)
    if moresigmas==True:
        plt.xlabel(r'$d/\sigma_b$')
    else:
        plt.xlabel(r'$d$ (nm)')
    plt.ylabel(r'$D/D_{\mathregular{bulk}}$')
    #if moresigmas==True:
    #    plt.title('Diffusion constant $D$ vs $d/\sigma_b$ for dynamic and static brushes',  y=1.03)
    #else:
    #    plt.title('Diffusion constant $D$ vs $d$ for dynamic and static brushes', y=1.03)
    ax.axis([0,11,minforsmall,maxforsmall]) # 6e-7 before we divided by Dbulk
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
else:
    plt.figure(figsize=(12.8,10))
    ax = plt.subplot(111)
    ax.plot(spacings_dyn, Dzs_dyn, '-o' , color='b',label=r'$D_\perp$, dyn.')
    ax.plot(spacings_dyn, Dparallel_dyn, '-*', color='g', label=r'$D_\parallel$, dyn.')
    ax.plot(spacings_stat, Dzs_stat, '-s', color='c', label=r'$D_\perp$, stat.')
    ax.plot(spacings_stat, Dparallel_stat, '-d', color='limegreen',label=r'$D_\parallel$, stat.')
    ax.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='b',alpha=0.2)
    ax.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
    ax.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='c', alpha=0.2)
    ax.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen',alpha=0.2)
    if moresigmas==True:
        plt.xlabel(r'$d/\sigma_b$', fontsize=20)
    else:
        plt.xlabel(r'$d$ (nm)', fontsize=20)
    plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=20)
    #if moresigmas==True:
    #    plt.title('Diffusion constant $D$ vs $d/\sigma_b$ for dynamic and static brushes', fontsize=28,  y=1.03)
    #else:
    #    plt.title('Diffusion constant $D$ vs $d$ for dynamic and static brushes', fontsize=28, y=1.03)
    ax.axis([0,11,minforsmall,maxforsmall]) # 6e-7 before we divided by Dbulk
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.savefig(plotname_cut)

################## twoinone ####################################
if big==False:
    print("Plottin it")
    #plt.figure(figsize=(16,5))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
    #if moresigmas==True:
    #    fig.suptitle('Diffusion constant $D/D_{\mathregular{bulk}}$ vs $d/\sigma_b$ for dynamic and static brushes')#, y=1.03)
    #else:
    #    fig.suptitle('Diffusion constant $D/D_{\mathregular{bulk}}$ vs $d$ for dynamic and static brushes')#, y=1.03)
    ax1.plot(spacings_dyn, Dzs_dyn, '-o' , color='b',label=r'$D_\perp$, dyn.')
    ax1.plot(spacings_dyn, Dparallel_dyn, '-*', color='g', label=r'$D_\parallel$, dyn.')
    ax1.plot(spacings_stat, Dzs_stat, '-s', color='c', label=r'$D_\perp$, stat.')
    ax1.plot(spacings_stat, Dparallel_stat, '-d', color='limegreen',label=r'$D_\parallel$, stat.')
    ax1.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='b',alpha=0.2)
    ax1.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
    ax1.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='c', alpha=0.2)
    ax1.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen',alpha=0.2)
    if moresigmas==True:
        ax1.set(xlabel=r'$d/\sigma_b$', ylabel='$D/D_{\mathregular{bulk}}$')
    else:
        ax1.set(xlabel=r'$d$ (nm)', ylabel='$D/D_{\mathregular{bulk}}$')
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax1.legend(loc="lower right")
    ax2.plot(spacings_dyn, Dzs_dyn, '-o' , color='b',label=r'$D_\perp$, dyn.')
    ax2.plot(spacings_dyn, Dparallel_dyn, '-*', color='g', label=r'$D_\parallel$, dyn.')
    ax2.plot(spacings_stat, Dzs_stat, '-s', color='c', label=r'$D_\perp$, stat.')
    ax2.plot(spacings_stat, Dparallel_stat, '-d', color='limegreen',label=r'$D_\parallel$, stat.')
    ax2.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='b',alpha=0.2)
    ax2.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
    ax2.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='c', alpha=0.2)
    ax2.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen',alpha=0.2)
    if moresigmas==True:
        ax2.set(xlabel=r'$d/\sigma_b$', ylabel=r'$D/D_{\mathregular{bulk}}$')
    else:
        ax2.set(xlabel=r'$d$ (nm)', ylabel=r'$D/D_{\mathregular{bulk}}$')
    ax2.axis([0,11,minforsmall,maxforsmall]) # 6e-7 before we divided by Dbulk
    #fig.tight_layout()
    plt.show()
    fig.savefig(plotname_twoinone)
else:
    print("Plottin it")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
    #if moresigmas==True:
    #   fig.suptitle('Diffusion constant $D$ vs $d/\sigma_b$ for dynamic and static brushes', fontsize=28, y=1.03)
    #else:
    #    fig.suptitle('Diffusion constant $D$ vs $d$ for dynamic and static brushes', fontsize=28, y=1.03)
    ax1.plot(spacings_dyn, Dzs_dyn, '-o' , color='b',label=r'$D_\perp$, dyn.')
    ax1.plot(spacings_dyn, Dparallel_dyn, '-*', color='g', label=r'$D_\parallel$, dyn.')
    ax1.plot(spacings_stat, Dzs_stat, '-s', color='c', label=r'$D_\perp$, stat.')
    ax1.plot(spacings_stat, Dparallel_stat, '-d', color='limegreen',label=r'$D_\parallel$, stat.')
    ax1.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='b',alpha=0.2)
    ax1.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
    ax1.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='c', alpha=0.2)
    ax1.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen',alpha=0.2)
    if moresigmas==True:
        ax1.set(xlabel=r'$d/\sigma_b$', ylabel=r'$D/D_{\mathregular{bulk}}$', fontsize=20)
    else:
        ax1.set(xlabel=r'$d$ (nm)', ylabel=r'$D/D_{\mathregular{bulk}}$', fontsize=20)
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax2.plot(spacings_dyn, Dzs_dyn, '-o' , color='b',label=r'$D_\perp$, dyn.')
    ax2.plot(spacings_dyn, Dparallel_dyn, '-*', color='g', label=r'$D_\parallel$, dyn.')
    ax2.plot(spacings_stat, Dzs_stat, '-s', color='c', label=r'$D_\perp$, stat.')
    ax2.plot(spacings_stat, Dparallel_stat, '-d', color='limegreen',label=r'$D_\parallel$, stat.')
    ax2.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='b',alpha=0.2)
    ax2.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
    ax2.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='c', alpha=0.2)
    ax2.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen',alpha=0.2)
    if moresigmas==True:
        ax2.set(xlabel=r'$d/\sigma_b$', ylabel=r'$D/D_{\mathregular{bulk}}$', fontsize=20)
    else:
        ax2.set(xlabel=r'$d$ (nm)', ylabel=r'$D/D_{\mathregular{bulk}}$', fontsize=20)
    ax2.axis([0,11,minforsmall,maxforsmall]) # 6e-7 before we divided by Dbulk
    ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #fig.tight_layout()
    plt.show()
    fig.savefig(plotname_twoinone)

print("Done.")
plt.show()