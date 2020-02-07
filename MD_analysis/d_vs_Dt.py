from skimage import morphology, measure, io, util   # For performing morphology operations (should maybe import more from skimage...)
from mpl_toolkits.mplot3d import Axes3D             # Plotting in 3D
import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
from scipy.ndimage import measurements, convolve    # Should have used the command below, but changing now will lead to confusion
from scipy import ndimage                           # For Euclidean distance measurement
import maintools_percolation as perctools
import numpy as np
import random
import math
import time
import os
import glob


damp = 15
spacing = 10
sigma_atom = 1


# Conversion factors
writeevery   = 10
unitlength   = 1e-9
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime

configfirst  = 1
configlast   = 1000
filestext    = '_config'+str(configfirst)+'to'+str(configlast)
endlocation  = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing%i/damp%i_diffseedLgv/Brush/Sigma_bead_' % (spacing,damp)+str(sigma_atom) + '/'

infilename_D  = endlocation+'diffusion'+filestext+'.txt'
infilename_t  = endlocation+'average_collisioninterval.txt'
infilename_dt = endlocation+'dt.txt'
outfilename   = endlocation+'Dt.txt'
#plotname     = endlocation+'d_vs_sqrtDtcoll.png' # I need more if I'm gonna plot

infile_dt = open(infilename_dt, 'r')
line      = infile_dt.readline()
words     = line.split()
dt        = float(words[0])
infile_dt.close()

infile_t     = open(infilename_t, 'r')
line         = infile_t.readline()
words        = line.split()
delta_t_av   = float(words[0])*dt
delta_t_rms  = float(words[1])*dt
infile_t.close()

'''
D_R2  sigmaD_R2 b_R2 sigmab_R2; D_z2  sigmaD_z2  b_z2  sigmaD_z2; D_par2 sigmaD_par2  b_par2  sigmab_par2
7.90903e-08 5.39804e-11 -1.67672e-15 5.35587e-18 2.73310e-08 8.89148e-11 -5.89414e-16 5.35587e-18 6.83171e-08 7.28447e-11 -1.96395e-15 5.35587e-18
'''

infile_D  = open(infilename_D, 'r')
lines     = infile_D.readlines()
words     = lines[1].split()
DR2_av    = float(words[0])
DR2_rms   = float(words[1])
bR2_av    = float(words[2])
bR2_rms   = float(words[3])
Dz2_av    = float(words[4])
Dz2_rms   = float(words[5])
bz2_av    = float(words[6])
bz2_rms   = float(words[7])
Dpar2_av  = float(words[8])
Dpar2_rms = float(words[9])
bpar2_av  = float(words[10])
bpar2_rms = float(words[11])
infile_D.close()

alpha = np.sqrt(delta_t_av*DR2_av) # Should I use this D? Probably, since the chain beads can be anywhere.
alpha_rms = 1/(2*alpha)*np.sqrt((DR2_av*delta_t_rms)**2+(delta_t_av*DR2_rms)**2) # because alpha=sqrt(D*t)

outfile = open(outfilename,'w')
outfile.write('%.16e %.16e' % (alpha, alpha_rms))
outfile.close()

print('Alpha, alpharms: %.16e %.16e' % (alpha, alpha_rms))
