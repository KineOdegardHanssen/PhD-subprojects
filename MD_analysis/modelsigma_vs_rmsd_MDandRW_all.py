import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob

showfig = True

rsphere = 1.0 #0.8
sigmas  = np.arange(1,16) # Change this for the modelling part

# Fixed parameters
psigma   = 1 # For instance
damp     = 10

## Paths:
basepath_base      = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/'
endlocation_static = basepath_base+'D_vs_d/Nocut/'
## Figure name:
plotname_rw   = endlocation_static+'D_rmsd_vs_sigma_allmodels_RW_rsphere'+str(rsphere)+'.png'
plotname_f    = endlocation_static+'D_rmsd_vs_sigma_allmodels_forest.png'
plotname_dyn  = endlocation_static+'D_rmsd_vs_sigma_allmodels_dynamic.png'
plotname_stat = endlocation_static+'D_rmsd_vs_sigma_allmodels_static.png'

# Ordered packing
rmsds_op_RW   = []
rmsds_op_f    = []
rmsds_op_dyn  = []
rmsds_op_stat = []
# HYperbola of revolution
rmsds_hr_RW   = []
rmsds_hr_f    = []
rmsds_hr_dyn  = []
rmsds_hr_stat = []
# Not monosized spheres
rmsds_nm_RW   = []
rmsds_nm_f    = []
rmsds_nm_dyn  = []
rmsds_nm_stat = []
# Partly saturated homogeneous isotropic monodisperse sphere packings
rmsds_ps_RW   = []
rmsds_ps_f    = []
rmsds_ps_dyn  = []
rmsds_ps_stat = []
# Overlapping spheres
rmsds_os_RW   = []
rmsds_os_f    = []
rmsds_os_dyn  = []
rmsds_os_stat = []
# Overlapping cylinders
rmsds_oc_RW   = []
rmsds_oc_f    = []
rmsds_oc_dyn  = []
rmsds_oc_stat = []
# Heterogeneous catalyst
rmsds_hc_RW   = []
rmsds_hc_f    = []
rmsds_hc_dyn  = []
rmsds_hc_stat = []
# Cation-exchange resin membrane
rmsds_rm_RW   = []
rmsds_rm_f    = []
rmsds_rm_dyn  = []
rmsds_rm_stat = []


for sigma in sigmas:
    ## Files to read from: # I messed up the naming a bit...
    rmsdfilename_RW = endlocation_static+'D_rmsd_close_packing_vs_RWgeom_norefl_rsphere'+str(rsphere)+'_maxlen1.0_modelr'+str(sigma)+'.txt'
    rmsdfilename_f = endlocation_static+'D_rmsd_close_packing_vs_forest_modelr'+str(sigma)+'.txt'
    rmsdfilename_dyn = endlocation_static+'D_rmsd_close_packing_vs_dyn_modelr'+str(sigma)+'.txt'
    rmsdfilename_stat = endlocation_static+'D_rmsd_close_packing_vs_stat_modelr'+str(sigma)+'.txt'
    ## Opening files:
    rmsdfile_RW   = open(rmsdfilename_RW,'r')
    rmsdfile_f    = open(rmsdfilename_f,'r')
    rmsdfile_dyn  = open(rmsdfilename_dyn,'r')
    rmsdfile_stat = open(rmsdfilename_stat,'r')
    ## Reading content:
    line_RW   = rmsdfile_RW.readline()
    line_f    = rmsdfile_f.readline()
    line_dyn  = rmsdfile_dyn.readline()
    line_stat = rmsdfile_stat.readline()
    #
    words_RW   = line_RW.split()
    words_f    = line_f.split()
    words_dyn  = line_dyn.split()
    words_stat = line_stat.split()
    ## Append to list:
    # Ordered packing (op)
    rmsds_op_RW.append(float(words_RW[0]))
    rmsds_op_f.append(float(words_f[0]))
    rmsds_op_dyn.append(float(words_dyn[0]))
    rmsds_op_stat.append(float(words_stat[0]))
    # Hyperbola of revolution (hr)
    rmsds_hr_RW.append(float(words_RW[1]))
    rmsds_hr_f.append(float(words_f[1]))
    rmsds_hr_dyn.append(float(words_dyn[1]))
    rmsds_hr_stat.append(float(words_stat[1]))
    # Not for monosized spheres (nm)
    rmsds_nm_RW.append(float(words_RW[2]))
    rmsds_nm_f.append(float(words_f[2]))
    rmsds_nm_dyn.append(float(words_dyn[2]))
    rmsds_nm_stat.append(float(words_stat[2]))
    # Partly saturated homogeneous isotropic monodisperse sphere packings (ps)
    rmsds_ps_RW.append(float(words_RW[3]))
    rmsds_ps_f.append(float(words_f[3]))
    rmsds_ps_dyn.append(float(words_dyn[3]))
    rmsds_ps_stat.append(float(words_stat[3]))
    # Overlapping spheres (os)
    rmsds_os_RW.append(float(words_RW[4]))
    rmsds_os_f.append(float(words_f[4]))
    rmsds_os_dyn.append(float(words_dyn[4]))
    rmsds_os_stat.append(float(words_stat[4]))
    # Overlapping cylinders (oc)
    rmsds_oc_RW.append(float(words_RW[5]))
    rmsds_oc_f.append(float(words_f[5]))
    rmsds_oc_dyn.append(float(words_dyn[5]))
    rmsds_oc_stat.append(float(words_stat[5]))
    # Heterogeneous catalyst (hc)
    rmsds_hc_RW.append(float(words_RW[6]))
    rmsds_hc_f.append(float(words_f[6]))
    rmsds_hc_dyn.append(float(words_dyn[6]))
    rmsds_hc_stat.append(float(words_stat[6]))
    # Cation-exchange resin membrane (rm)
    rmsds_rm_RW.append(float(words_RW[7]))
    rmsds_rm_f.append(float(words_f[7]))
    rmsds_rm_dyn.append(float(words_dyn[7]))
    rmsds_rm_stat.append(float(words_stat[7]))
    ## Close files:
    rmsdfile_f.close()
    rmsdfile_RW.close()
    rmsdfile_dyn.close()
    rmsdfile_stat.close()

plt.figure(figsize=(14,5))
ax = plt.subplot(111)
ax.plot(sigmas,rmsds_op_RW,label='Ordered packing')
ax.plot(sigmas,rmsds_hr_RW,label='A hyperbola of revolution')
ax.plot(sigmas,rmsds_nm_RW,label='Not for monosized spheres')
ax.plot(sigmas,rmsds_ps_RW,label='Monodisperse sphere packings')
ax.plot(sigmas,rmsds_os_RW,label='Overlapping spheres')
ax.plot(sigmas,rmsds_oc_RW,label='Overlapping cylinders')
ax.plot(sigmas,rmsds_hc_RW,label='Heterogeneous catalyst')
ax.plot(sigmas,rmsds_rm_RW,label='Cation-exchange resin membrane')
plt.xlabel(r'Obstacle radius $\sigma_m$, model')
plt.ylabel('RMSD D, models vs random walk')
plt.title('RMSD D, models vs random walk')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname_rw)


plt.figure(figsize=(14,5))
ax = plt.subplot(111)
ax.plot(sigmas,rmsds_op_f,label='Ordered packing')
ax.plot(sigmas,rmsds_hr_f,label='A hyperbola of revolution')
ax.plot(sigmas,rmsds_nm_f,label='Not for monosized spheres')
ax.plot(sigmas,rmsds_ps_f,label='Monodisperse sphere packings')
ax.plot(sigmas,rmsds_os_f,label='Overlapping spheres')
ax.plot(sigmas,rmsds_oc_f,label='Overlapping cylinders')
ax.plot(sigmas,rmsds_hc_f,label='Heterogeneous catalyst')
ax.plot(sigmas,rmsds_rm_f,label='Cation-exchange resin membrane')
plt.xlabel(r'Obstacle radius $\sigma_m$, model')
plt.ylabel('RMSD D, models vs straight brush MD')
plt.title('RMSD D, models vs straight brush MD')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname_f)



plt.figure(figsize=(14,5))
ax = plt.subplot(111)
ax.plot(sigmas,rmsds_op_dyn,label='Ordered packing')
ax.plot(sigmas,rmsds_hr_dyn,label='A hyperbola of revolution')
ax.plot(sigmas,rmsds_nm_dyn,label='Not for monosized spheres')
ax.plot(sigmas,rmsds_ps_dyn,label='Monodisperse sphere packings')
ax.plot(sigmas,rmsds_os_dyn,label='Overlapping spheres')
ax.plot(sigmas,rmsds_oc_dyn,label='Overlapping cylinders')
ax.plot(sigmas,rmsds_hc_dyn,label='Heterogeneous catalyst')
ax.plot(sigmas,rmsds_rm_dyn,label='Cation-exchange resin membrane')
plt.xlabel(r'Obstacle radius $\sigma_m$, model')
plt.ylabel('RMSD D, models vs dynamic brush MD')
plt.title('RMSD D, models vs dynamic brush MD')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname_dyn)

plt.figure(figsize=(14,5))
ax = plt.subplot(111)
ax.plot(sigmas,rmsds_op_stat,label='Ordered packing')
ax.plot(sigmas,rmsds_hr_stat,label='A hyperbola of revolution')
ax.plot(sigmas,rmsds_nm_stat,label='Not for monosized spheres')
ax.plot(sigmas,rmsds_ps_stat,label='Monodisperse sphere packings')
ax.plot(sigmas,rmsds_os_stat,label='Overlapping spheres')
ax.plot(sigmas,rmsds_oc_stat,label='Overlapping cylinders')
ax.plot(sigmas,rmsds_hc_stat,label='Heterogeneous catalyst')
ax.plot(sigmas,rmsds_rm_stat,label='Cation-exchange resin membrane')
plt.xlabel(r'Obstacle radius $\sigma_m$, model')
plt.ylabel('RMSD D, models vs static brush MD')
plt.title('RMSD D, models vs static brush MD')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname_stat)

if showfig==True:
    plt.show()

