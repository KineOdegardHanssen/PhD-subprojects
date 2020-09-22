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

def avg_and_rms(x):
    N = len(x)
    avg = np.mean(x)
    rms = 0
    for i in range(N):
        rms += (x[i]-avg)*(x[i]-avg)
    rms = np.sqrt(rms/(N-1)) 
    return avg, rms

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
brushfilename_dyn  = endlocation + 'D_vs_d.txt'
brushfilename_stat = endlocation_static + 'D_vs_d_static.txt'

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
    Dzs_dyn[j] = float(words[5])
    Dparallel_dyn[j] = float(words[9])
    # Ds, stdv
    DRs_stdv_dyn[j] = float(words[2])
    Dzs_stdv_dyn[j] = float(words[6])
    Dparallel_stdv_dyn[j] = float(words[10])
    
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

testr     = 20     # For transport measurements (default 50, might change for tests)
confignrs = np.arange(1,1001)


endlocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/tD_vs_th/'
plotname1   = endlocation+'tD_vs_th.png'
plotname2   = endlocation+'tD_div_th.png'
plotname1_smalld = endlocation+'tD_vs_th_smalld.png'
plotname2_smalld = endlocation+'tD_div_th_smalld.png'


spacings = [1,1.25,1.5,2,2.5,3,3.5,4,5,6,7,8,10,15,25,50,75,100]
avtrs    = []
avths    = []
rmstrs    = []
rmsths    = []
for spacing in spacings:
    # For tr (movement in xy-plane): 
    endlocation          = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp + str(psigma) + '/'
    filestext            = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])
    try:
        infilename = endlocation+'trs_r%i' % testr +filestext+'.txt'
        infile     = open(infilename,'r')
    except:
        infilename  = endlocation+'trs_r%i' % testr +filestext+'_long.txt'
        infile      = open(infilename,'r')
    
    firstline = infile.readline()
    Nreal     = float(firstline.split()[2]) # Nread in the original script
    lines     = infile.readlines()
    
    testr_times = []
    exitrs      = []
    
    for line in lines:
        words = line.split()
        if len(words)>0:
            testr_times.append(float(words[0]))
    infile.close()
    avtr, rmstr = avg_and_rms(testr_times) 
    avtrs.append(avtr)
    rmstrs.append(rmstr)
    
    # For th (movement in z-dir):
    try:
        infilename = endlocation+'ths_h%i' % testr +filestext+'.txt'
        infile     = open(infilename,'r')
    except:
        infilename  = endlocation+'ths_h%i' % testr +filestext+'_long.txt'
        infile      = open(infilename,'r')
    
    firstline = infile.readline()
    Nreal     = float(firstline.split()[2]) # Nread in the original script
    lines     = infile.readlines()
    
    for line in lines:
        words = line.split()
        if len(words)>0:
            testr_times.append(float(words[0]))
    infile.close()
    avth, rmsth = avg_and_rms(testr_times)
    avths.append(avth)
    rmsths.append(rmsth)


testr_nm = testr*1e-9 # Length in nm
tDlh = np.zeros(N_dyn)
tDlr = np.zeros(N_dyn)
th_div_tDlh = np.zeros(N_dyn)
tr_div_tDlr = np.zeros(N_dyn)
for i in range(N_dyn):
    tDlh[i] = testr_nm**2/Dzs_dyn[i]
    tDlr[i] = testr_nm**2/Dparallel_dyn[i]
    th_div_tDlh[i] = avths[i]/tDlh[i]
    tr_div_tDlr[i] = avtrs[i]/tDlr[i]

plt.figure()
#plt.errorbar(spacings, avtrs, yerr=rmstrs, label=r'$<t_r>$', capsize=2)
#plt.errorbar(spacings, avths, yerr=rmsths, label=r'$<t_h>$', capsize=2)
plt.plot(spacings, avtrs, label=r'$<t_r>$')
plt.plot(spacings, avths, label=r'$<t_h>$')
plt.plot(spacings_dyn, tDlh, label=r'$t_{D,h}$')
plt.plot(spacings_dyn, tDlr, label=r'$t_{D,r}$')
plt.xlabel('d [nm]')
plt.ylabel('t [s]')
plt.title('$t_{D,h}$ and $<t_h>$ vs d')
plt.legend(loc='upper right')
plt.savefig(plotname1)

plt.figure()
plt.plot(spacings_dyn, th_div_tDlh, label=r'$<t_h>/t_{D,h}$')
plt.plot(spacings_dyn, tr_div_tDlr, label=r'$<t_r>/t_{D,r}$')
plt.xlabel('d [nm]')
plt.ylabel('Fraction')
plt.title('$<t_h>/t_{D,h}$ vs d')
plt.legend(loc='lower right')
plt.savefig(plotname2)


# Could have picked out small d in a neater way...
i = 0
smalld_spacings    = []
smalld_avtrs       = []
smalld_avths       = []
smalld_rmstr       = []
smalld_rmsth       = []
smalld_tDlh        = []
smalld_tDlr        = []
smalld_th_div_tDlh = []
smalld_tr_div_tDlr = []

while spacings[i]<15:
    smalld_spacings.append(spacings[i])
    smalld_avtrs.append(avtrs[i])
    smalld_avths.append(avths[i])
    smalld_rmstr.append(rmstrs[i])
    smalld_rmsth.append(rmsths[i])
    smalld_tDlr.append(tDlr[i])
    smalld_tDlh.append(tDlh[i])
    smalld_th_div_tDlh.append(th_div_tDlh[i])
    smalld_tr_div_tDlr.append(tr_div_tDlr[i])
    i+=1

plt.figure()
plt.errorbar(smalld_spacings, smalld_avtrs, yerr=smalld_rmstr, label=r'$<t_r>$', capsize=2)
plt.errorbar(smalld_spacings, smalld_avths, yerr=smalld_rmsth, label=r'$<t_h>$', capsize=2)
plt.plot(smalld_spacings, smalld_tDlh, label=r'$t_{D,h}$')
plt.plot(smalld_spacings, smalld_tDlr, label=r'$t_{D,r}$')
plt.xlabel('d [nm]')
plt.ylabel('t [s]')
plt.title('$t_{D,h}$ and $<t_h>$ vs d')
plt.legend(loc='upper right')
plt.savefig(plotname1_smalld)

plt.figure()
plt.plot(smalld_spacings, smalld_th_div_tDlh, label=r'$<t_h>/t_{D,h}$')
plt.plot(smalld_spacings, smalld_tr_div_tDlr, label=r'$<t_r>/t_{D,r}$')
plt.xlabel('d [nm]')
plt.ylabel('Fraction')
plt.title('$<t_h>/t_{D,h}$ vs d')
plt.legend(loc='lower right')
plt.savefig(plotname2_smalld)


print('------------------------------------------------------------------------------------')
print('avths:',avths)
print('---')
print('tDlh:',tDlh)
print('------------')
print('th_div_tDlh:',th_div_tDlh)
print('------------------------------------------------------------------------------------')
print('spacings:', spacings)
print('---')
print('avtrs:',avtrs)
print('---')
print('tDlr:',tDlr)
print('------------')
print('tr_div_tDlr:',tr_div_tDlr)
print('------------------------------------------------------------------------------------')