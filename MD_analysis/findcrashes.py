from skimage import morphology, measure, io, util   # For performing morphology operations (should maybe import more from skimage...)
from mpl_toolkits.mplot3d import Axes3D             # Plotting in 3D
import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
from scipy.ndimage import measurements, convolve    # Should have used the command below, but changing now will lead to confusion
from scipy import ndimage                           # For Euclidean distance measurement
import maintools_percolation as perctools
import data_treatment as datr
import numpy as np
import random
import math
import time
import os
import glob
import copy

spacing = 10
sigma_atom = 1
configs = np.arange(1,5)

basefolder = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Brush/Quadr_M9N101_ljunits_Langevin_scaled_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_ljcut1p122/Spacing%i/Sigma_bead_'% (spacing)+ str(sigma_atom)+'/'
folder = basefolder + 'log-files/'
outfolder = basefolder
outfilename = outfolder + 'collisiontimes.txt'

all_crashtimes = []

outfile = open(outfilename, 'w')
outfile.write('Results written by findcrashes.py. Order of elements: time of collision, potential energy of bead\n')

for config in configs: # Write to file inside loop?
    print('config:', config)
    outfile.write('Config nr: %i\n' % config)
    
    filename_infolder = 'log.confignr%i' % config
    infilename = folder + filename_infolder
    
    infile = open(infilename, 'r')
    lines = infile.readlines()
    
    startit = 0
    times   = []
    crashtimes   = []
    pot_energies = []
    for line in lines:
        words = line.split()
        if len(words)>1:
            if words[0]=='Loop':
                break
            if startit==1:
                times.append(int(words[0]))
                pot_energies.append(float(words[1]))
            if words[0]=='Step':
                startit = 1
    infile.close()
    
    N = len(times)
    crashbefore = 0
    for i in range(N):
        pot_en   = pot_energies[i]
        thistime = times[i]
        if (pot_en!=0 and crashbefore==0):
            crashtimes.append(times[i])
            outfile.write('%i %.16f\n' % (thistime, pot_en))
            print('crashtime: %i, pot.en.:, %.5f' % (thistime, pot_en))
            crashbefore = 1
        elif pot_en==0:
            crashbefore = 0
    
    all_crashtimes.append(crashtimes)

outfile.close()
