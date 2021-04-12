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
plotname = endlocation_static+'D_rmsd_vs_sigma_close_packing_models_rsphere'+str(rsphere)+'.png'

rmsds_RW   = []
rmsds_f    = []
rmsds_dyn  = []
rmsds_stat = []

for sigma in sigmas:
    ## Files to read from:
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
    ## Append to list:
    rmsds_RW.append(float(line_RW.split()[0]))
    rmsds_f.append(float(line_f.split()[0]))
    rmsds_dyn.append(float(line_dyn.split()[0]))
    rmsds_stat.append(float(line_stat.split()[0]))
    ## Close files:
    rmsdfile_f.close()
    rmsdfile_RW.close()
    rmsdfile_dyn.close()
    rmsdfile_stat.close()

plt.figure(figsize=(6,5))
plt.plot(sigmas,rmsds_RW,label='RW')
plt.plot(sigmas,rmsds_dyn,label='MD dyn.')
plt.plot(sigmas,rmsds_stat,label='MD static')
plt.plot(sigmas,rmsds_f,label='MD straight')
plt.xlabel(r'Obstacle radius $\sigma_m$, model')
plt.ylabel('RMSD D, close packing vs data')
plt.title('RMSD D, close packing vs data')
plt.legend(loc='center right')
plt.tight_layout()
plt.savefig(plotname)

if showfig==True:
    plt.show()

