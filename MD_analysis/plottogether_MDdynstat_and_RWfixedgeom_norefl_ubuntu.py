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

showplots = False

# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
rsphere  = 1.0 #0.8
sigma    = 13 # Change this for the modelling part

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

endlocation = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/Nocut/'
if moresigmas==True:
    endlocation = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Moresigmas/Nocut/'

basepath_base      = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/'
endlocation_static = basepath_base+'D_vs_d/Nocut/'

bulklocation = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/'

bulkfilename  = bulklocation + 'diffusion_bulk'+filestext
if bulk_cut==True:
    bulkfilename = bulkfilename +'_cut.txt'
else:
    bulkfilename = bulkfilename +'_uncut.txt'

## Files to read
brushfilename_dyn  = endlocation + 'D_vs_d.txt'
brushfilename_stat = endlocation_static + 'D_vs_d_static.txt'

## Files to write to
rmsdfilename_RW = endlocation_static+'D_rmsd_close_packing_vs_RWgeom_norefl_rsphere'+str(rsphere)+'_maxlen1.0_modelr'+str(sigma)+'.txt'
rmsdfilename_forest = endlocation_static+'D_rmsd_close_packing_vs_forest_modelr'+str(sigma)+'.txt'
rmsdfilename_dyn = endlocation_static+'D_rmsd_close_packing_vs_dyn_modelr'+str(sigma)+'.txt'
rmsdfilename_stat = endlocation_static+'D_rmsd_close_packing_vs_stat_modelr'+str(sigma)+'.txt'
if big==False:
    plotname     = endlocation_static+'D_graftdens_dyn_vs_stat_vs_RWgeom_norefl_rsphere'+str(rsphere)+'.png'
    plotname_cut = endlocation_static+'D_graftdens_dyn_vs_stat_vs_RWgeom_cut_norefl_rsphere'+str(rsphere)+'.png' ## Use this?
    plotname_loglog = endlocation_static+'D_graftdens_dyn_vs_stat_vs_RWgeom_loglog_norefl_rsphere'+str(rsphere)+'.png'
    plotname_ylog = endlocation_static+'D_graftdens_dyn_vs_stat_vs_RWgeom_ylog_norefl_rsphere'+str(rsphere)+'.png'
    ### d
    plotname_d     = endlocation_static+'D_vs_d_dyn_vs_stat_vs_RWgeom_norefl_rsphere'+str(rsphere)+'_modelr'+str(sigma)+'.png'
    plotname_cut_d = endlocation_static+'D_vs_d_dyn_vs_stat_vs_RWgeom_cut_norefl_rsphere'+str(rsphere)+'_modelr'+str(sigma)+'.png' ## Use this?
    plotname_loglog_d = endlocation_static+'D_vs_d_dyn_vs_stat_vs_RWgeom_loglog_norefl_rsphere'+str(rsphere)+'.png'
    plotname_ylog_d = endlocation_static+'D_vs_d_dyn_vs_stat_vs_RWgeom_ylog_norefl_rsphere'+str(rsphere)+'.png'
else:
    plotname     = endlocation_static+'D_graftdens_dyn_vs_stat_vs_RWgeom_big_norefl_rsphere'+str(rsphere)+'.png'
    plotname_cut = endlocation_static+'D_graftdens_dyn_vs_stat_cut_vs_RWgeom_big_norefl_rsphere'+str(rsphere)+'.png' ## Use this?
    plotname_loglog = endlocation_static+'D_graftdens_dyn_vs_stat_vs_RWgeom_big_loglog_norefl_rsphere'+str(rsphere)+'.png'
    plotname_ylog = endlocation_static+'D_graftdens_dyn_vs_stat_vs_RWgeom_big_ylog_norefl_rsphere'+str(rsphere)+'.png'
    ### d 
    plotname_d     = endlocation_static+'D_vs_d_dyn_vs_stat_vs_RWgeom_big_norefl_rsphere'+str(rsphere)+'_modelr'+str(sigma)+'.png'
    plotname_cut_d = endlocation_static+'D_vs_d_dyn_vs_stat_cut_vs_RWgeom_big_norefl_rsphere'+str(rsphere)+'_modelr'+str(sigma)+'.png' ## Use this?
    plotname_loglog_d = endlocation_static+'D_vs_d_dyn_vs_stat_vs_RWgeom_big_loglog_norefl_rsphere'+str(rsphere)+'.png'
    plotname_ylog_d = endlocation_static+'D_vs_d_dyn_vs_stat_vs_RWgeom_big_ylog_norefl_rsphere'+str(rsphere)+'.png'

# Dynamic sims:
brushfile_dyn = open(brushfilename_dyn, 'r')

lines = brushfile_dyn.readlines()
N_dyn = len(lines)-1

# ds
dens_dyn = np.zeros(N_dyn)
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
    
    d = float(words[0])
    spacings_dyn[j] = d
    dens_dyn[j] = np.pi/float(d**2)
    DRs_dyn[j]  = float(words[1])
    Dzs_dyn[j]  = float(words[5])
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
dens_stat     = np.zeros(N_stat)
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
    
    d = float(words[0])
    spacings_stat[j] = d
    dens_stat[j] = np.pi/float(d**2)
    DRs_stat[j]  = float(words[1])
    Dzs_stat[j]  = float(words[3])
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
    DRs_stat[i] = DRs_stat[i]/DRs_bulk
    Dzs_stat[i] = Dzs_stat[i]/DRs_bulk
    Dparallel_stat[i] =  Dparallel_stat[i]/DRs_bulk
    # Ds, stdv
    DRs_stdv_stat[i] = abs(DRs_bulk)*np.sqrt((DRs_stdv_stat[i]/DRs_stat[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dzs_stdv_stat[i] = abs(DRs_bulk)*np.sqrt((Dzs_stdv_stat[i]/Dzs_stat[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dparallel_stdv_stat[i] = abs(DRs_bulk)*np.sqrt((Dparallel_stdv_stat[i]/Dparallel_stat[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)

for i in range(N_dyn):
    DRs_dyn[i] = DRs_dyn[i]/DRs_bulk
    Dzs_dyn[i] = Dzs_dyn[i]/DRs_bulk
    Dparallel_dyn[i] =  Dparallel_dyn[i]/DRs_bulk
    # Ds, stdv
    DRs_stdv_dyn[i] = abs(DRs_bulk)*np.sqrt((DRs_stdv_dyn[i]/DRs_dyn[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dzs_stdv_dyn[i] = abs(DRs_bulk)*np.sqrt((Dzs_stdv_dyn[i]/Dzs_dyn[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dparallel_stdv_dyn[i] = abs(DRs_bulk)*np.sqrt((Dparallel_stdv_dyn[i]/Dparallel_dyn[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)

## Forest (a new addition):
endlocation_forest = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_forest/D_vs_d/Nocut/'
infilename_f  = endlocation_forest+'Dd_div_Dbulk_vs_d_2_forest_uncut.txt'
infile_f = open(infilename_f,'r')
lines = infile_f.readlines()
N_f   = len(lines)-1

# ds
dens_f     = np.zeros(N_f)
spacings_f = np.zeros(N_f)
# Ds
DRs_f = np.zeros(N_f)
Dxs_f = np.zeros(N_f)
Dys_f = np.zeros(N_f)
Dzs_f = np.zeros(N_f)
Dparallel_f = np.zeros(N_f)
# Ds, stdv
DRs_stdv_f = np.zeros(N_f)
Dxs_stdv_f = np.zeros(N_f)
Dys_stdv_f = np.zeros(N_f)
Dzs_stdv_f = np.zeros(N_f)
Dparallel_stdv_f = np.zeros(N_f)

for i in range(1,N_f+1):
    words = lines[i].split()
    j = i-1
    
    d = float(words[0])
    spacings_f[j] = d
    dens_f[j] = np.pi/float(d**2)
    DRs_f[j]  = float(words[1])
    Dzs_f[j]  = float(words[3])
    Dparallel_f[j] = float(words[5])
    # Ds, stdv
    DRs_stdv_f[j] = float(words[2])
    Dzs_stdv_f[j] = float(words[4])
    Dparallel_stdv_f[j] = float(words[6])
infile_f.close()

## RW part:
Nsteps   = 1000   # Increase #Maybe for other values of hitprob
Nreal    = 10000
folder   = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Randomwalk_fixedgeometry/'

infilename_Dbulk = folder+'2D_bulk_Nsteps%i_Nreal%i_D.txt' %(Nsteps,Nreal)
infile_Ddbulk = open(infilename_Dbulk,'r')
line = infile_Ddbulk.readline()
Dbulk = float(line.split()[0])
infile_Ddbulk.close()

folder   = folder+'No_reflection/'
infilename_Ddbulk = folder+'D_vs_d_Nsteps%i_Nreal%i_rsphere'%(Nsteps,Nreal)+str(rsphere)+'_norefl.txt' 
infile_Ddbulk = open(infilename_Ddbulk,'r')
lines = infile_Ddbulk.readlines()
print('lines:',lines)

print('infilename_Ddbulk:',infilename_Ddbulk)
phit = []
DRW  = []
dg   = []
for line in lines:
    words = line.split()
    print('words:',words)
    #if len(words)!=0:
    #    phit.append(float(words[0]))
    #    DRW.append(float(words[1]))
    d = float(words[0])
    dg.append(d)
    phit.append(np.pi*(rsphere**2)/d**2)
    DRW.append(float(words[1])/Dbulk)
infile_Ddbulk.close()
print('phit:',phit)


ds = np.linspace(1,100,51)
porosity = 1-(np.pi*(sigma/ds)**2)
Ds_op = 2.0/(3.0-porosity) # Diffusion constant for ordered packings
print('porosity:',porosity)
print('Ds_op:',Ds_op)

# Density:
plt.figure(figsize=(6,5))
plt.errorbar(dens_dyn, Dparallel_dyn, yerr=Dparallel_stdv_dyn, color='g', capsize=2, label=r'$D_\parallel$, dyn.')
plt.errorbar(dens_stat, Dparallel_stat, yerr=Dparallel_stdv_stat, color='limegreen', capsize=2, label=r'$D_\parallel$, stat.')
plt.errorbar(dens_f, Dparallel_f, yerr=Dparallel_stdv_f, color='springgreen', capsize=2, label=r'$D_\parallel$, straight')
plt.plot(phit,DRW, color='b', label='Random walk, fixed geom., no refl.')
plt.xlabel(r'Density $\phi$')
plt.ylabel(r'Normalized diffusion coefficient $D/D_{bulk}$')
plt.title(r'$D/D_{bulk}$ vs $\phi$')
plt.legend(loc='upper right', prop={'size': 12})
plt.tight_layout()
plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.loglog(dens_dyn[2:], Dparallel_dyn[2:], '-o', color='g', label=r'$D_\parallel$, dyn.')
plt.loglog(dens_stat[2:], Dparallel_stat[2:], '-o', color='limegreen', label=r'$D_\parallel$, stat.')
plt.loglog(dens_f[2:], Dparallel_f[2:], '-o', color='springgreen', label=r'$D_\parallel$, straight')
plt.loglog(phit,DRW, '-o', color='b', label='Random walk, fixed geometry')
plt.xlabel(r'Density $\phi$')
plt.ylabel(r'Normalized diffusion coefficient $D/D_{bulk}$')
plt.title(r'$D/D_{bulk}$ vs $\phi$')
plt.legend(loc='lower left', prop={'size': 12})
plt.tight_layout()
plt.savefig(plotname_loglog)

plt.figure(figsize=(6,5))
plt.plot(dens_dyn[2:], Dparallel_dyn[2:], '-o', color='g', label=r'$D_\parallel$, dyn.')
plt.plot(dens_stat[2:], Dparallel_stat[2:], '-o', color='limegreen', label=r'$D_\parallel$, stat.')
plt.plot(dens_f, Dparallel_f, '-o', color='springgreen',label=r'$D_\parallel$, straight')
plt.plot(phit,DRW, '-o', color='b', label='Random walk, fixed geom., no refl.')
plt.xlabel(r'Density $\phi$')
plt.ylabel(r'Normalized diffusion coefficient $D/D_{bulk}$')
plt.title(r'$D/D_{bulk}$ vs $\phi$')
plt.legend(loc='lower left', prop={'size': 12})
plt.yscale('log')
plt.tight_layout()
plt.savefig(plotname_ylog)

### Spacing
plt.figure(figsize=(6,5))
plt.errorbar(spacings_dyn, Dparallel_dyn, yerr=Dparallel_stdv_dyn, color='g', capsize=2, label=r'$D_\parallel$, dyn.')
plt.errorbar(spacings_stat, Dparallel_stat, yerr=Dparallel_stdv_stat, color='limegreen', capsize=2, label=r'$D_\parallel$, stat.')
plt.errorbar(spacings_f, Dparallel_f, yerr=Dparallel_stdv_f, color='springgreen', capsize=2, label=r'$D_\parallel$, straight')
plt.plot(dg,DRW, color='b', label='Random walk, fixed geom., no refl.')
plt.plot(ds,Ds_op, '--', label='Ordered packings')
plt.xlabel(r'Spacing $d$')
plt.ylabel(r'Normalized diffusion coefficient $D/D_{bulk}$')
plt.title(r'$D/D_{bulk}$ vs $d$')
plt.legend(loc='lower right', prop={'size': 12})
plt.tight_layout()
plt.savefig(plotname_d)

plt.figure(figsize=(6,5))
plt.loglog(spacings_dyn[2:], Dparallel_dyn[2:], '-o', color='g', label=r'$D_\parallel$, dyn.')
plt.loglog(spacings_stat[2:], Dparallel_stat[2:], '-o', color='limegreen', label=r'$D_\parallel$, stat.')
plt.loglog(spacings_f, Dparallel_f, '-o', color='springgreen', label=r'$D_\parallel$, straight')
plt.loglog(dg,DRW, '-o', color='b', label='Random walk, fixed geom., no refl.')
plt.xlabel(r'Spacing $d$')
plt.ylabel(r'Normalized diffusion coefficient $D/D_{bulk}$')
plt.title(r'$D/D_{bulk}$ vs $d$')
plt.legend(loc='lower right', prop={'size': 12})
plt.tight_layout()
plt.savefig(plotname_loglog_d)

plt.figure(figsize=(6,5))
plt.plot(spacings_dyn[2:], Dparallel_dyn[2:], '-o', color='g', label=r'$D_\parallel$, dyn.')
plt.plot(spacings_stat[2:], Dparallel_stat[2:], '-o', color='limegreen', label=r'$D_\parallel$, stat.')
plt.plot(spacings_f[2:], Dparallel_f[2:], '-o',color='springgreen', label=r'$D_\parallel$, straight')
plt.plot(dg,DRW, '-o', color='b', label='Random walk, geom., no refl.')
plt.xlabel('Density (or probability of hitting obstacle $p_{hit}$)')
plt.xlabel(r'Spacing $d$')
plt.title(r'$D/D_{bulk}$ vs $d$')
plt.legend(loc='lower right', prop={'size': 12})
plt.yscale('log')
plt.tight_layout()
plt.savefig(plotname_ylog_d)


dg = np.array(dg)
spacings_stat = np.array(spacings_stat)
spacings_dyn  = np.array(spacings_dyn)
spacings_f    = np.array(spacings_f)
porosity_rw   = 1-(np.pi*(sigma/dg)**2)
porosity_stat = 1-(np.pi*(sigma/spacings_stat)**2)
porosity_dyn  = 1-(np.pi*(sigma/spacings_dyn)**2)
porosity_f    = 1-(np.pi*(sigma/spacings_f)**2)
Ds_op_rw      = 2.0/(3.0-porosity_rw) # Diffusion constant for ordered packings
Ds_op_stat    = 2.0/(3.0-porosity_stat)
Ds_op_dyn     = 2.0/(3.0-porosity_dyn)
Ds_op_f       = 2.0/(3.0-porosity_f)

D_rmsd_rw   = rmsd(DRW,Ds_op_rw)
D_rmsd_stat = rmsd(Dparallel_stat,Ds_op_stat)
D_rmsd_dyn  = rmsd(Dparallel_dyn,Ds_op_dyn)
D_rmsd_f    = rmsd(Dparallel_f,Ds_op_f)

# Write rmsd to file
'''
rmsdfilename_RW
rmsdfilename_forest
rmsdfilename_dyn
rmsdfilename_stat
'''

rmsdfile = open(rmsdfilename_RW,'w')
rmsdfile.write('%.16f' % D_rmsd_rw)
rmsdfile.close()

rmsdfile = open(rmsdfilename_forest,'w')
rmsdfile.write('%.16f' % D_rmsd_f)
rmsdfile.close()

rmsdfile = open(rmsdfilename_dyn,'w')
rmsdfile.write('%.16f' % D_rmsd_dyn)
rmsdfile.close()

rmsdfile = open(rmsdfilename_stat,'w')
rmsdfile.write('%.16f' % D_rmsd_stat)
rmsdfile.close()

if showplots==True:
    plt.show()
